#define _XOPEN_SOURCE 700 // For zlib/POSIX functions if needed, though we strictly follow C99 for logic
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* Dependencies: kseq.h and khashl.h must be available */
#include "kseq.h"
#include "khashl.h"

// Initialize kseq with zlib
KSEQ_INIT(gzFile, gzread)

// Define Hash Map: Key=uint64_t, Value=char* (read name)
KHASHL_MAP_INIT(KH_LOCAL, kmer_map_t, kmap, uint64_t, char*, kh_hash_uint64, kh_eq_generic)

// C99 compliant implementation of strdup
char *c99_strdup(const char *s) {
    if (!s) return NULL;
    size_t len = strlen(s) + 1;
    char *new_s = malloc(len);
    if (new_s) memcpy(new_s, s, len);
    return new_s;
}

// DNA encoding table: A=0, C=1, G=2, T=3, Others=4
unsigned char seq_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, // A, C, G
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, // T, U
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, // a, c, g
    4, 4, 4, 4,  3, 3, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, // t, u
    // ... rest default to 4 (initialized implicitly if static, but explicit here for clarity)
};

// Helper to encode a sequence of specific length into uint64
int encode_kmer(const char *seq, int len, uint64_t *val) {
    uint64_t x = 0;
    for (int i = 0; i < len; ++i) {
        uint8_t c = seq_nt4_table[(uint8_t)seq[i]];
        if (c > 3) return 0; // Invalid base (N)
        x = (x << 2) | c;
    }
    *val = x;
    return 1;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: ./find_kmers <queries.fasta> <genome.fasta>\n");
        return 1;
    }

    gzFile fp_query, fp_genome;
    kseq_t *ks;
    kmer_map_t *h = kmap_init();
    khint_t itr;
    int absent;
    int k_len = 0; 

    // ----------------------------------------------------------------
    // 1. Load Queries and Detect K
    // ----------------------------------------------------------------
    fprintf(stderr, "[INFO] Reading queries from %s...\n", argv[1]);
    fp_query = gzopen(argv[1], "r");
    if (!fp_query) { perror("Failed to open query file"); return 1; }
    
    ks = kseq_init(fp_query);
    uint64_t loaded_count = 0;

    // Detect K from first record
    if (kseq_read(ks) >= 0) {
        k_len = ks->seq.l;
        if (k_len > 31) {
            fprintf(stderr, "[ERROR] K-mer length %d exceeds limit of 31.\n", k_len);
            return 1;
        }
        if (k_len == 0) {
            fprintf(stderr, "[ERROR] Empty sequence detected.\n");
            return 1;
        }
        fprintf(stderr, "[INFO] Detected k-mer length: %d\n", k_len);

        uint64_t packed;
        if (encode_kmer(ks->seq.s, k_len, &packed)) {
            itr = kmap_put(h, packed, &absent);
            if (absent) kh_val(h, itr) = c99_strdup(ks->name.s);
            loaded_count++;
        }
    } else {
        fprintf(stderr, "[ERROR] Query file is empty.\n");
        return 1;
    }

    // Load remaining records
    while (kseq_read(ks) >= 0) {
        if (ks->seq.l != k_len) {
            continue; // Skip queries of different length
        }
        uint64_t packed;
        if (encode_kmer(ks->seq.s, k_len, &packed)) {
            itr = kmap_put(h, packed, &absent);
            if (absent) {
                kh_val(h, itr) = c99_strdup(ks->name.s); 
                loaded_count++;
            }
        }
    }
    kseq_destroy(ks);
    gzclose(fp_query);
    
    fprintf(stderr, "[INFO] Loaded %lu unique valid %d-mers.\n", (unsigned long)loaded_count, k_len);

    // ----------------------------------------------------------------
    // 2. Scan Genome and Output BED
    // ----------------------------------------------------------------
    fprintf(stderr, "[INFO] Scanning genome %s...\n", argv[2]);
    fp_genome = gzopen(argv[2], "r");
    if (!fp_genome) { perror("Failed to open genome file"); return 1; }
    
    ks = kseq_init(fp_genome);
    
    const uint64_t mask = (1ULL << (k_len * 2)) - 1;
    const int rc_shift = (k_len - 1) * 2;

    while (kseq_read(ks) >= 0) {
        uint64_t fwd_k = 0; 
        uint64_t rev_k = 0; 
        int len = 0;        
        
        for (int i = 0; i < ks->seq.l; ++i) {
            uint8_t c = seq_nt4_table[(uint8_t)ks->seq.s[i]];
            
            if (c < 4) {
                fwd_k = ((fwd_k << 2) | c) & mask;
                rev_k = (rev_k >> 2) | (((uint64_t)(c ^ 3)) << rc_shift);
                len++;
                
                if (len >= k_len) {
                    // BED Coordinates:
                    // chromStart (0-based inclusive) = i - k_len + 1
                    // chromEnd   (0-based exclusive) = i + 1
                    int start_pos = i - k_len + 1;
                    int end_pos   = i + 1;

                    // Check Forward
                    itr = kmap_get(h, fwd_k);
                    if (itr != kh_end(h)) {
                        // chrom, start, end, name, score, strand
                        printf("%s\t%d\t%d\t%s\t0\t+\n", 
                               ks->name.s, start_pos, end_pos, kh_val(h, itr));
                    }
                    
                    // Check Reverse
                    itr = kmap_get(h, rev_k);
                    if (itr != kh_end(h)) {
                        printf("%s\t%d\t%d\t%s\t0\t-\n", 
                               ks->name.s, start_pos, end_pos, kh_val(h, itr));
                    }
                }
            } else {
                len = 0;
                fwd_k = 0;
                rev_k = 0;
            }
        }
    }

    kseq_destroy(ks);
    gzclose(fp_genome);
    
    // Free map
    for (itr = 0; itr < kh_end(h); ++itr) {
        if (kh_exist(h, itr)) free(kh_val(h, itr));
    }
    kmap_destroy(h);

    return 0;
}

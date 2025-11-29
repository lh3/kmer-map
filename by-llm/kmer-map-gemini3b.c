#define _XOPEN_SOURCE 700 // Expose POSIX symbols like strdup
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
#include "khashl.h"

// Initialize kseq with zlib's gzFile
KSEQ_INIT(gzFile, gzread)

// Define the hash map using khashl macros
// Key: uint64_t (encoded k-mer sequence)
// Value: char* (k-mer name)
KHASHL_MAP_INIT(KH_LOCAL, kmer_map_t, kmer_map, uint64_t, char*, kh_hash_uint64, kh_eq_generic)

// Helper for C99 compliance if strdup is not available
static char *my_strdup(const char *s) {
    size_t len = strlen(s);
    char *d = malloc(len + 1);
    if (d) memcpy(d, s, len + 1);
    return d;
}

// Global table for nucleotide encoding (A=0, C=1, G=2, T=3)
static unsigned char seq_nt4_table[256];

void init_seq_nt4_table() {
    memset(seq_nt4_table, 4, 256); // Initialize all to 4 (invalid)
    seq_nt4_table['A'] = seq_nt4_table['a'] = 0;
    seq_nt4_table['C'] = seq_nt4_table['c'] = 1;
    seq_nt4_table['G'] = seq_nt4_table['g'] = 2;
    seq_nt4_table['T'] = seq_nt4_table['t'] = 3;
}

// Helper to free map values and destroy map
void cleanup_map(kmer_map_t *h) {
    if (!h) return;
    for (khint_t k = 0; k < kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            free(kh_val(h, k));
        }
    }
    kmer_map_destroy(h);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <kmers.fa/fq> <genome.fa/fq>\n", argv[0]);
        return 1;
    }

    init_seq_nt4_table();

    // -------------------------------------------------------------------------
    // Step 1: Load and Index K-mers
    // -------------------------------------------------------------------------
    gzFile fp_kmer = gzopen(argv[1], "r");
    if (!fp_kmer) {
        perror("Error opening k-mer file");
        return 1;
    }

    kseq_t *ks_kmer = kseq_init(fp_kmer);
    kmer_map_t *h = kmer_map_init();
    
    int K = -1; // To store the expected k-mer length
    int ret;
    
    while ((ret = kseq_read(ks_kmer)) >= 0) {
        int len = ks_kmer->seq.l;

        // Constraint: k-mers up to 31bp
        if (len > 31) {
            fprintf(stderr, "Warning: k-mer '%s' length %d > 31. Skipped.\n", ks_kmer->name.s, len);
            continue;
        }

        // Constraint: Equal length
        if (K == -1) {
            K = len;
        } else if (len != K) {
            fprintf(stderr, "Warning: k-mer '%s' length %d differs from expected %d. Skipped.\n", ks_kmer->name.s, len, K);
            continue;
        }

        // Encode Sequence
        uint64_t key = 0;
        int valid = 1;
        for (int i = 0; i < len; ++i) {
            uint8_t c = seq_nt4_table[(uint8_t)ks_kmer->seq.s[i]];
            if (c > 3) {
                valid = 0;
                break;
            }
            key = (key << 2) | c;
        }

        if (!valid) {
            fprintf(stderr, "Warning: k-mer '%s' contains non-ACGT characters. Skipped.\n", ks_kmer->name.s);
            continue;
        }

        // Constraint: No duplicates
        int absent;
        khint_t k_idx = kmer_map_put(h, key, &absent);
        if (!absent) {
            fprintf(stderr, "Warning: Duplicate k-mer sequence '%s' (name: %s). First instance kept.\n", ks_kmer->seq.s, ks_kmer->name.s);
        } else {
            kh_val(h, k_idx) = my_strdup(ks_kmer->name.s);
        }
    }
    
    kseq_destroy(ks_kmer);
    gzclose(fp_kmer);

    if (K == -1) {
        fprintf(stderr, "Error: No valid k-mers found.\n");
        cleanup_map(h);
        return 1;
    }

    // -------------------------------------------------------------------------
    // Step 2: Scan Genome
    // -------------------------------------------------------------------------
    gzFile fp_genome = gzopen(argv[2], "r");
    if (!fp_genome) {
        perror("Error opening genome file");
        cleanup_map(h);
        return 1;
    }

    kseq_t *ks_genome = kseq_init(fp_genome);

    // Prepare rolling hash constants
    // Mask: Keep only the lower 2*K bits (K <= 31, so 2*K <= 62)
    uint64_t mask = (1ULL << (2 * K)) - 1;
    // Shift for RC: New base complement goes to the most significant position
    int shift_rc = 2 * (K - 1);

    while ((ret = kseq_read(ks_genome)) >= 0) {
        int len = ks_genome->seq.l;
        if (len < K) continue; // Sequence too short

        uint64_t fwd = 0; // Forward strand rolling hash
        uint64_t rev = 0; // Reverse complement rolling hash
        int valid_bases = 0; // Count of consecutive valid ACGT bases

        for (int i = 0; i < len; ++i) {
            uint8_t code = seq_nt4_table[(uint8_t)ks_genome->seq.s[i]];

            if (code > 3) {
                // Non-ACGT base resets the rolling window
                valid_bases = 0;
                fwd = 0;
                rev = 0;
            } else {
                // Update Forward Hash: Shift left, add new base, mask
                fwd = ((fwd << 2) | code) & mask;

                // Update Reverse Hash: Shift right, add complement of new base at MSB
                // Complement: A(0)<->T(3), C(1)<->G(2) => code ^ 3
                rev = (rev >> 2) | ((uint64_t)(code ^ 3) << shift_rc);
                
                valid_bases++;
            }

            if (valid_bases >= K) {
                // 1. Check if Forward window matches a k-mer
                khint_t k_idx = kmer_map_get(h, fwd);
                if (k_idx != kh_end(h)) {
                    // Hit on Positive Strand
                    printf("%s\t%d\t%d\t%s\t0\t+\n", 
                           ks_genome->name.s, i - K + 1, i + 1, kh_val(h, k_idx));
                }

                // 2. Check if Reverse Complement window matches a k-mer
                // If RC(Window) is in our list, it means the list k-mer was found on the negative strand
                k_idx = kmer_map_get(h, rev);
                if (k_idx != kh_end(h)) {
                    // Hit on Negative Strand
                    printf("%s\t%d\t%d\t%s\t0\t-\n", 
                           ks_genome->name.s, i - K + 1, i + 1, kh_val(h, k_idx));
                }
            }
        }
    }

    // -------------------------------------------------------------------------
    // Cleanup
    // -------------------------------------------------------------------------
    kseq_destroy(ks_genome);
    gzclose(fp_genome);
    cleanup_map(h);

    return 0;
}

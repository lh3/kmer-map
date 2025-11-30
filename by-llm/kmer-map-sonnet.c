/*
 * K-mer Position Finder with 2-bit encoding
 * Finds positions of k-mers in a genome sequence
 * Both input files can be FASTA/FASTQ format, optionally gzip compressed
 * 
 * Compile: gcc -O3 -Wall -std=c99 -o kmer_find kmer_find.c -lz
 * Usage: ./kmer_find <kmers.fa[.gz]> <genome.fa[.gz]>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

// Include khashl.h for hash table
#define KHASHL_IMPLEMENTATION
#include "khashl.h"

// Hash table: key=encoded k-mer (64-bit), value=k-mer index
KHASHL_SET_INIT(KH_LOCAL, kh_kmer_t, kh_kmer, uint64_t, kh_hash_uint64, kh_eq_generic)

typedef struct {
    char **kmer_names;      // Original k-mer strings for output
    uint64_t *kmer_codes;   // 2-bit encoded k-mers
    uint32_t *kmer_ids;     // K-mer indices
    int n_kmers;
    int kmer_len;
    kh_kmer_t *h;
} kmer_set_t;

// 2-bit encoding: A=0, C=1, G=2, T=3
static inline int nt2bit(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;  // Invalid nucleotide
    }
}

// Encode a k-mer string to 64-bit integer
// Returns 1 on success, 0 if invalid nucleotides found
static inline int encode_kmer(const char *seq, int k, uint64_t *code) {
    uint64_t val = 0;
    for (int i = 0; i < k; i++) {
        int b = nt2bit(seq[i]);
        if (b < 0) return 0;  // Invalid nucleotide
        val = (val << 2) | b;
    }
    *code = val;
    return 1;
}

// Free k-mer set
void kmer_set_destroy(kmer_set_t *ks) {
    if (!ks) return;
    if (ks->kmer_names) {
        for (int i = 0; i < ks->n_kmers; i++) {
            free(ks->kmer_names[i]);
        }
        free(ks->kmer_names);
    }
    free(ks->kmer_codes);
    free(ks->kmer_ids);
    kh_kmer_destroy(ks->h);
    free(ks);
}

// Load k-mers from FASTA/FASTQ file
kmer_set_t *load_kmers(const char *fn) {
    gzFile fp = gzopen(fn, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Failed to open k-mer file: %s\n", fn);
        return NULL;
    }
    
    kseq_t *seq = kseq_init(fp);
    kmer_set_t *ks = calloc(1, sizeof(kmer_set_t));
    ks->h = kh_kmer_init();
    
    int capacity = 1024;
    ks->kmer_names = malloc(capacity * sizeof(char*));
    ks->kmer_codes = malloc(capacity * sizeof(uint64_t));
    ks->kmer_ids = malloc(capacity * sizeof(uint32_t));
    ks->n_kmers = 0;
    ks->kmer_len = -1;
    
    int has_length_warning = 0;
    int has_dup_warning = 0;
    int has_invalid_warning = 0;
    
    while (kseq_read(seq) >= 0) {
        int len = seq->seq.l;
        
        // Check k-mer length
        if (ks->kmer_len < 0) {
            ks->kmer_len = len;
            if (len > 31) {
                fprintf(stderr, "WARNING: K-mer length %d exceeds 31bp\n", len);
            }
        } else if (len != ks->kmer_len && !has_length_warning) {
            fprintf(stderr, "WARNING: K-mers have different lengths (%d vs %d)\n", 
                    ks->kmer_len, len);
            has_length_warning = 1;
        }
        
        // Encode k-mer
        uint64_t code;
        if (!encode_kmer(seq->seq.s, len, &code)) {
            if (!has_invalid_warning) {
                fprintf(stderr, "WARNING: K-mer with invalid nucleotides: %s\n", seq->seq.s);
                has_invalid_warning = 1;
            }
            continue;
        }
        
        // Check for duplicates
        khint_t k = kh_kmer_get(ks->h, code);
        if (k != kh_end(ks->h)) {
            if (!has_dup_warning) {
                fprintf(stderr, "WARNING: Duplicate k-mer found: %s\n", seq->seq.s);
                has_dup_warning = 1;
            }
            continue;
        }
        
        // Expand arrays if needed
        if (ks->n_kmers >= capacity) {
            capacity *= 2;
            ks->kmer_names = realloc(ks->kmer_names, capacity * sizeof(char*));
            ks->kmer_codes = realloc(ks->kmer_codes, capacity * sizeof(uint64_t));
            ks->kmer_ids = realloc(ks->kmer_ids, capacity * sizeof(uint32_t));
        }
        
        // Add to hash table
        int absent;
        k = kh_kmer_put(ks->h, code, &absent);
        
        // Store k-mer
        ks->kmer_names[ks->n_kmers] = strdup(seq->seq.s);
        ks->kmer_codes[ks->n_kmers] = code;
        ks->kmer_ids[ks->n_kmers] = ks->n_kmers;
        ks->n_kmers++;
    }
    
    kseq_destroy(seq);
    gzclose(fp);
    
    if (ks->n_kmers == 0) {
        fprintf(stderr, "ERROR: No k-mers loaded\n");
        kmer_set_destroy(ks);
        return NULL;
    }
    
    fprintf(stderr, "Loaded %d k-mers of length %d\n", ks->n_kmers, ks->kmer_len);
    return ks;
}

// Search k-mers in genome and output BED format
void search_kmers(const char *genome_fn, kmer_set_t *ks) {
    gzFile fp = gzopen(genome_fn, "r");
    if (!fp) {
        fprintf(stderr, "ERROR: Failed to open genome file: %s\n", genome_fn);
        return;
    }
    
    kseq_t *seq = kseq_init(fp);
    int k = ks->kmer_len;
    uint64_t mask = (1ULL << (2 * k)) - 1;  // Mask for k-mer bits
    
    while (kseq_read(seq) >= 0) {
        const char *chrom = seq->name.s;
        const char *s = seq->seq.s;
        int seq_len = seq->seq.l;
        
        if (seq_len < k) continue;
        
        // Encode first k-mer
        uint64_t fwd = 0, rev = 0;
        int valid_len = 0;
        
        for (int i = 0; i < seq_len; i++) {
            int b = nt2bit(s[i]);
            
            if (b < 0) {
                // Invalid nucleotide, reset
                valid_len = 0;
                fwd = rev = 0;
                continue;
            }
            
            // Update forward k-mer (shift left, add new base)
            fwd = ((fwd << 2) | b) & mask;
            
            // Update reverse k-mer (shift right, add complement to left)
            rev = (rev >> 2) | ((uint64_t)(b ^ 3) << (2 * (k - 1)));
            
            valid_len++;
            
            if (valid_len >= k) {
                // Check forward strand
                khint_t h = kh_kmer_get(ks->h, fwd);
                if (h != kh_end(ks->h)) {
                    int pos = i - k + 1;
                    // Find which k-mer this is
                    for (int j = 0; j < ks->n_kmers; j++) {
                        if (ks->kmer_codes[j] == fwd) {
                            printf("%s\t%d\t%d\t%s\t0\t+\n", 
                                   chrom, pos, pos + k, ks->kmer_names[j]);
                            break;
                        }
                    }
                }
                
                // Check reverse strand
                h = kh_kmer_get(ks->h, rev);
                if (h != kh_end(ks->h)) {
                    int pos = i - k + 1;
                    // Find which k-mer this is
                    for (int j = 0; j < ks->n_kmers; j++) {
                        if (ks->kmer_codes[j] == rev) {
                            printf("%s\t%d\t%d\t%s\t0\t-\n", 
                                   chrom, pos, pos + k, ks->kmer_names[j]);
                            break;
                        }
                    }
                }
            }
        }
    }
    
    kseq_destroy(seq);
    gzclose(fp);
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <kmers.fa[.gz]> <genome.fa[.gz]>\n", argv[0]);
        fprintf(stderr, "Output: BED format (chrom, start, end, name, score, strand)\n");
        return 1;
    }
    
    // Load k-mers
    kmer_set_t *ks = load_kmers(argv[1]);
    if (!ks) return 1;
    
    // Search k-mers in genome
    search_kmers(argv[2], ks);
    
    // Cleanup
    kmer_set_destroy(ks);
    
    return 0;
}

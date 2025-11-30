#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"

KSEQ_INIT(gzFile, gzread)

KHASH_INIT(kmer, uint64_t, char *, 1, kh_int64_hash_func, kh_int64_hash_equal)

static uint64_t encode_kmer(const char *s, int k, int *valid) {
    uint64_t key = 0;
    *valid = 1;
    for (int i = 0; i < k; ++i) {
        char c = toupper(s[i]);
        int b;
        if (c == 'A') b = 0;
        else if (c == 'C') b = 1;
        else if (c == 'G') b = 2;
        else if (c == 'T') b = 3;
        else { *valid = 0; return 0; }
        key = (key << 2) | b;
    }
    return key;
}

static uint64_t revcomp_key(uint64_t key, int k) {
    uint64_t rc = 0;
    for (int i = 0; i < k; ++i) {
        int b = key & 3;
        int cb = b ^ 3;
        rc = (rc << 2) | cb;
        key >>= 2;
    }
    return rc;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <genome.fa/q[.gz]> <kmers.fa/q[.gz]>\n", argv[0]);
        return 1;
    }

    const char *genome_file = argv[1];
    const char *kmers_file = argv[2];

    // Read k-mers
    gzFile fp_k = gzopen(kmers_file, "r");
    if (!fp_k) {
        fprintf(stderr, "Error: cannot open %s\n", kmers_file);
        return 1;
    }
    kseq_t *sk = kseq_init(fp_k);

    int k = -1;
    khash_t(kmer) *h = kh_init(kmer);
    while (kseq_read(sk) >= 0) {
        int l = (int)sk->seq.l;
        if (k == -1) k = l;
        else if (l != k) {
            fprintf(stderr, "Warning: unequal k-mer length %d != %d for %s\n", l, k, sk->name.s);
            continue;  // Skip unequal lengths
        }

        int valid;
        uint64_t key = encode_kmer(sk->seq.s, k, &valid);
        if (!valid) {
            fprintf(stderr, "Warning: invalid bases in k-mer %s\n", sk->name.s);
            continue;
        }

        int absent;
        khint_t it = kh_put(kmer, h, key, &absent);
        if (!absent) {
            fprintf(stderr, "Warning: duplicate k-mer sequence for %s (previous: %s)\n", sk->name.s, kh_val(h, it));
            free(kh_val(h, it));  // Overwrite with new name
        }
        kh_val(h, it) = strdup(sk->name.s);
    }

    kseq_destroy(sk);
    gzclose(fp_k);

    if (k == -1) {
        fprintf(stderr, "Error: no k-mers found\n");
        kh_destroy(kmer, h);
        return 1;
    }
    if (k > 31) {
        fprintf(stderr, "Error: k-mer length %d > 31\n", k);
        kh_destroy(kmer, h);
        return 1;
    }

    // Read genome
    gzFile fp_g = gzopen(genome_file, "r");
    if (!fp_g) {
        fprintf(stderr, "Error: cannot open %s\n", genome_file);
        kh_destroy(kmer, h);
        return 1;
    }
    kseq_t *sg = kseq_init(fp_g);

    while (kseq_read(sg) >= 0) {
        char *seq = sg->seq.s;
        int len = (int)sg->seq.l;
        char *chrom = sg->name.s;

        for (int i = 0; i <= len - k; ++i) {
            int valid;
            uint64_t fwd_key = encode_kmer(seq + i, k, &valid);
            if (!valid) continue;

            // Check forward strand
            khint_t it = kh_get(kmer, h, fwd_key);
            if (it != kh_end(h)) {
                printf("%s\t%d\t%d\t%s\t0\t+\n", chrom, i, i + k, kh_val(h, it));
            }

            // Check reverse complement strand
            uint64_t rev_key = revcomp_key(fwd_key, k);
            it = kh_get(kmer, h, rev_key);
            if (it != kh_end(h)) {
                printf("%s\t%d\t%d\t%s\t0\t-\n", chrom, i, i + k, kh_val(h, it));
            }
        }
    }

    kseq_destroy(sg);
    gzclose(fp_g);

    // Cleanup
    for (khint_t i = 0; i < kh_end(h); ++i) {
        if (kh_exist(h, i)) free(kh_val(h, i));
    }
    kh_destroy(kmer, h);

    return 0;
}

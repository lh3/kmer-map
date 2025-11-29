/*
 * find_kmers.c
 *
 * Find genomic positions of k-mers (<=31 bp) from FASTA/FASTQ input.
 *
 * Features:
 *  - FASTA/FASTQ or gzip-compressed inputs via kseq.h + zlib
 *  - k-mers: must be same length (warn and skip if inconsistent)
 *  - k <= 31 (2-bit packed in uint64_t)
 *  - canonical matching: forward and reverse-complement are both reported
 *  - duplicates: warned & ignored
 *  - output: BED6 (0-based)
 *
 * Compile:
 *   gcc -std=c99 -O3 -Wall -o find_kmers find_kmers.c -lz
 *
 * Requires:
 *   kseq.h   (from klib)
 *   khashl.h (from klib)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <zlib.h>

#include "kseq.h"
#include "khashl.h"

KSEQ_INIT(gzFile, gzread)

/* khashl: instantiate uint64_t → int map
 *   Table type: km64_t
 *   API prefix: k64_
 */
KHASHL_MAP_INIT(KH_LOCAL, km64_t, k64, uint64_t, int,
                kh_hash_uint64, kh_eq_generic)

char *kom_strdup(const char *src) // strdup() doesn't conform to C99
{
	size_t len;
	char *dst;
	len = strlen(src);
	dst = malloc(len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

/* Map base to 2-bit */
static inline int base2bit(char c) {
    switch (c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return 4;
    }
}

static inline uint64_t complement2bit(uint64_t b) {
    return (3u - b);  /* 0↔3, 1↔2 */
}

/* K-mer storage */
typedef struct {
    char **seqs;      /* original sequences */
    uint64_t *enc;    /* encoded 2-bit */
    int n, m;
    int k;
} KmerSet;

static void kset_init(KmerSet *s) {
    s->seqs = NULL;
    s->enc = NULL;
    s->n = s->m = 0;
    s->k = 0;
}

static void ks_push(KmerSet *s, const char *seq, uint64_t code) {
    if (s->n == s->m) {
        s->m = s->m ? s->m * 2 : 256;
        s->seqs = (char**)realloc(s->seqs, s->m * sizeof(char*));
        s->enc  = (uint64_t*)realloc(s->enc,  s->m * sizeof(uint64_t));
    }
    s->seqs[s->n] = kom_strdup(seq);
    s->enc[s->n]  = code;
    s->n++;
}

static void ks_free(KmerSet *s) {
    for (int i = 0; i < s->n; i++) free(s->seqs[i]);
    free(s->seqs);
    free(s->enc);
}

/* Encode k-mer into uint64_t (2 bits/base) */
static int encode_kmer(const char *s, int k, uint64_t *out) {
    uint64_t v = 0;
    for (int i = 0; i < k; i++) {
        int b = base2bit(s[i]);
        if (b > 3) return 0;
        v = (v << 2) | (uint64_t)b;
    }
    *out = v;
    return 1;
}

static void usage(const char *p) {
    fprintf(stderr, "Usage: %s <genome.fa[.gz]> <kmers.fa[.gz]>\n", p);
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc != 3)
        usage(argv[0]);

    const char *genome_path = argv[1];
    const char *kmer_path   = argv[2];

    /* ---------------------------
     * Load k-mers
     * --------------------------- */
    gzFile fk = gzopen(kmer_path, "r");
    if (!fk) {
        fprintf(stderr, "ERROR: cannot open %s\n", kmer_path);
        return 1;
    }
    kseq_t *ks = kseq_init(fk);

    KmerSet set;
    kset_init(&set);
    km64_t *map = k64_init();

    int expected_k = -1;
    int rec = 0;
    int skip_nonacgt = 0, skip_badlen = 0, dups = 0;

    while (kseq_read(ks) >= 0) {
        rec++;
        int len = ks->seq.l;
        if (len <= 0) {
            fprintf(stderr, "WARNING: empty k-mer sequence at record %d\n", rec);
            continue;
        }
        if (expected_k < 0)
            expected_k = len;

        if (len != expected_k) {
            fprintf(stderr,
                "WARNING: k-mer length %d inconsistent with expected %d at record %d — skipping\n",
                len, expected_k, rec);
            skip_badlen++;
            continue;
        }
        if (expected_k > 31) {
            fprintf(stderr, "ERROR: k=%d > 31 not supported\n", expected_k);
            return 1;
        }

        uint64_t code;
        if (!encode_kmer(ks->seq.s, expected_k, &code)) {
            fprintf(stderr, "WARNING: non-ACGT k-mer '%s' — skipping\n", ks->seq.s);
            skip_nonacgt++;
            continue;
        }

        int absent;
        khint_t it = k64_put(map, code, &absent);
        if (!absent) {
            dups++;
            continue;
        }
        ks_push(&set, ks->seq.s, code);
        kh_val(map, it) = set.n - 1;
    }

    kseq_destroy(ks);
    gzclose(fk);

    if (set.n == 0) {
        fprintf(stderr, "ERROR: no valid k-mers loaded\n");
        return 1;
    }

    if (dups) fprintf(stderr, "WARNING: %d duplicate k-mers ignored\n", dups);
    if (skip_badlen) fprintf(stderr, "WARNING: %d k-mers had inconsistent length\n", skip_badlen);
    if (skip_nonacgt) fprintf(stderr, "WARNING: %d k-mers contained non-ACGT bases\n", skip_nonacgt);

    int k = expected_k;
    uint64_t mask = (k == 32 ? ~(uint64_t)0 : ((1ULL << (2*k)) - 1ULL));

    /* ---------------------------
     * Scan genome
     * --------------------------- */
    gzFile fg = gzopen(genome_path, "r");
    if (!fg) {
        fprintf(stderr, "ERROR: cannot open %s\n", genome_path);
        return 1;
    }
    kseq_t *g = kseq_init(fg);

    while (kseq_read(g) >= 0) {
        const char *chrom = g->name.s;
        const char *seq = g->seq.s;
        int L = g->seq.l;
        if (L < k) continue;

        uint64_t f = 0, r = 0;
        int filled = 0;

        for (int i = 0; i < L; i++) {
            int b = base2bit(seq[i]);

            if (b > 3) {
                /* Reset on N or other base */
                f = r = 0;
                filled = 0;
                continue;
            }

            /* update forward */
            f = ((f << 2) | (uint64_t)b) & mask;

            /* update reverse-complement */
            uint64_t cb = complement2bit(b);
            r = (r >> 2) | (cb << (2*(k-1)));

            if (filled < k) {
                filled++;
                if (filled < k) continue;
            }

            int start = i - k + 1;

            /* forward match */
            khint_t itf = k64_get(map, f);
            if (itf < kh_end(map)) {
                int idx = kh_val(map, itf);
                printf("%s\t%d\t%d\t%s\t0\t+\n",
                    chrom, start, start + k, set.seqs[idx]);
            }

            /* reverse-complement match */
            khint_t itr = k64_get(map, r);
            if (itr < kh_end(map)) {
                int idx = kh_val(map, itr);
                printf("%s\t%d\t%d\t%s\t0\t-\n",
                    chrom, start, start + k, set.seqs[idx]);
            }
        }
    }

    kseq_destroy(g);
    gzclose(fg);
    k64_destroy(map);
    ks_free(&set);

    return 0;
}

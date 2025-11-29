#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, kmermap_t, kmap, uint64_t, char*, kh_hash_uint64, kh_eq_generic)

#define kom_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define kom_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))

uint8_t kom_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

typedef struct {
	int32_t ksize;
	kmermap_t *h;
} kmerdb_t;

char *kom_strdup(const char *src) // strdup() doesn't conform to C99
{
	size_t len;
	char *dst;
	len = strlen(src);
	dst = kom_malloc(char, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

kmerdb_t *km_read_kmer(const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	kmerdb_t *db;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(0, "rb");
	if (fp == 0) return 0;
	ks = kseq_init(fp);
	db = kom_calloc(kmerdb_t, 1);
	db->h = kmap_init();
	while (kseq_read(ks) >= 0) {
		int32_t i, absent;
		uint64_t x = 0;
		khint_t k;
		if (ks->seq.l > 31) {
			fprintf(stderr, "Warning: ignore %s as it is longer than 31bp\n", ks->name.s);
			continue;
		}
		if (db->ksize <= 0) {
			db->ksize = ks->seq.l;
		} else if (db->ksize != ks->seq.l) {
			fprintf(stderr, "Warning: ignore %s as it is of different lengths\n", ks->name.s);
			continue;
		}
		for (i = 0; i < ks->seq.l; ++i) {
			int c = kom_nt4_table[(uint8_t)ks->seq.s[i]];
			if (c > 3) break;
			x = x << 2 | c;
		}
		if (i < ks->seq.l) {
			fprintf(stderr, "Warning: ignore %s as it contains ambiguous bases\n", ks->name.s);
			continue;
		}
		k = kmap_put(db->h, x, &absent);
		if (!absent) {
			fprintf(stderr, "Warning: ignore %s as it is a duplicate\n", ks->name.s);
			continue;
		}
		kh_val(db->h, k) = kom_strdup(ks->name.s);
	}
	kseq_destroy(ks);
	gzclose(fp);
	return db;
}

void km_find_match(const char *fn, const kmerdb_t *db)
{
	gzFile fp;
	kseq_t *ks;
	uint64_t mask = (1ULL << db->ksize*2) - 1;
	int32_t shift = (db->ksize - 1) * 2;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(0, "rb");
	if (fp == 0) return;
	ks = kseq_init(fp);
	while (kseq_read(ks) >= 0) {
		int64_t i, l;
		uint64_t xf = 0, xr = 0;
		for (i = l = 0; i < ks->seq.l; ++i) {
			int c = kom_nt4_table[(uint8_t)ks->seq.s[i]];
			if (c < 4) {
				xf = (xf << 2 | c) & mask;
				xr = xr >> 2 | (uint64_t)(3 - c) << shift;
				if (++l >= db->ksize) {
					khint_t kf, kr;
					kf = kmap_get(db->h, xf);
					kr = kmap_get(db->h, xr);
					if (kf != kh_end(db->h))
						printf("%s\t%ld\t%ld\t%s\t%d\t+\n", ks->name.s, (long)(i + 1 - db->ksize), (long)(i + 1), kh_val(db->h, kf), db->ksize);
					if (kr != kh_end(db->h))
						printf("%s\t%ld\t%ld\t%s\t%d\t-\n", ks->name.s, (long)(i + 1 - db->ksize), (long)(i + 1), kh_val(db->h, kr), db->ksize);
				}
			} else l = 0, xf = xr = 0;
		}
	}
	kseq_destroy(ks);
	gzclose(fp);
}

int main(int argc, char *argv[])
{
	kmerdb_t *db;
	if (argc < 3) {
		fprintf(stderr, "Usage: map-kmer <kmer.fa> <genome.fa>\n");
		return 1;
	}
	db = km_read_kmer(argv[1]);
	if (db == 0) {
		fprintf(stderr, "ERROR: failed to load the kmer file\n");
		return 1;
	}
	km_find_match(argv[2], db);
	// free db here, but the OS will take care of that anyway
	return 0;
}

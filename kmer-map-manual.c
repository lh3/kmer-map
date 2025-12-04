#include <zlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "khashl.h"
KHASHL_MAP_INIT(KH_LOCAL, kmermap_t, kmap, uint64_t, char*, kh_hash_uint64, kh_eq_generic)

#define kom_malloc(type, cnt)       ((type*)malloc((cnt) * sizeof(type)))
#define kom_calloc(type, cnt)       ((type*)calloc((cnt), sizeof(type)))
#define kom_realloc(type, ptr, cnt) ((type*)realloc((ptr), (cnt) * sizeof(type)))

#define kom_grow(type, ptr, __i, __m) do { \
		if ((__i) >= (__m)) { \
			(__m) = (__i) + 1; \
			(__m) += ((__m)>>1) + 16; \
			(ptr) = kom_realloc(type, (ptr), (__m)); \
		} \
	} while (0)

#define kom_roundup64(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, (x)|=(x)>>32, ++(x))

static inline void str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kom_roundup64(s->m);
		s->s = kom_realloc(char, s->s, s->m);
	}
}

static inline void str_copy(kstring_t *s, const char *st, const char *en)
{
	str_enlarge(s, en - st);
	memcpy(&s->s[s->l], st, en - st);
	s->l += en - st;
}

int64_t kom_sprintf_lite(kstring_t *s, const char *fmt, ...)
{
	char buf[32]; // for integer to string conversion
	const char *p, *q;
	int64_t len = 0;
	va_list ap;
	va_start(ap, fmt);
	for (q = p = fmt; *p; ++p) {
		if (*p == '%') {
			if (p > q) {
				len += p - q;
				if (s) str_copy(s, q, p);
			}
			++p;
			if (*p == 'd') {
				int c, i, l = 0;
				unsigned int x;
				c = va_arg(ap, int);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
			} else if (*p == 'l' && *(p+1) == 'd') {
				int i, l = 0;
				long c;
				unsigned long x;
				c = va_arg(ap, long);
				x = c >= 0? c : -c;
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				if (c < 0) buf[l++] = '-';
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
				++p;
			} else if (*p == 'u') {
				int i, l = 0;
				uint32_t x;
				x = va_arg(ap, uint32_t);
				do { buf[l++] = x%10 + '0'; x /= 10; } while (x > 0);
				len += l;
				if (s) {
					str_enlarge(s, l);
					for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
				}
			} else if (*p == 's') {
				char *r = va_arg(ap, char*);
				int l;
				l = strlen(r);
				len += l;
				if (s) str_copy(s, r, r + l);
			} else if (*p == 'c') {
				++len;
				if (s) {
					str_enlarge(s, 1);
					s->s[s->l++] = va_arg(ap, int);
				}
			} else {
				fprintf(stderr, "ERROR: unrecognized type '%%%c'\n", *p);
				abort();
			}
			q = p + 1;
		}
	}
	if (p > q) {
		len += p - q;
		if (s) str_copy(s, q, p);
	}
	va_end(ap);
	if (s) s->s[s->l] = 0;
	return len;
}

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
	int32_t m, len;
	int32_t bid, boff;
	uint8_t **a;
} block_arena_t;

block_arena_t *ba_init(int32_t block_len)
{
	block_arena_t *ba;
	ba = kom_calloc(block_arena_t, 1);
	ba->len = block_len;
	ba->m = 4;
	ba->a = kom_malloc(uint8_t*, ba->m);
	ba->a[0] = kom_calloc(uint8_t, ba->len);
	return ba;
}

void *ba_alloc(block_arena_t *ba, int32_t size)
{
	void *p;
	assert(size < ba->len);
	if (ba->boff + size > ba->len) {
		ba->boff = 0;
		++ba->bid;
		kom_grow(uint8_t*, ba->a, ba->bid, ba->m);
		ba->a[ba->bid] = kom_calloc(uint8_t, ba->len);
	}
	p = (void*)&ba->a[ba->bid][ba->boff];
	ba->boff += size;
	return p;
}

char *ba_strdup(block_arena_t *ba, const char *src) // strdup() doesn't conform to C99
{
	size_t len;
	char *dst;
	len = strlen(src);
	dst = (char*)ba_alloc(ba, len + 1);
	memcpy(dst, src, len + 1);
	return dst;
}

typedef struct {
	int32_t ksize;
	kmermap_t *h;
	block_arena_t *ba;
} kmerdb_t;

kmerdb_t *km_read_kmer(const char *fn)
{
	gzFile fp;
	kseq_t *ks;
	kmerdb_t *db;
	fp = fn && strcmp(fn, "-")? gzopen(fn, "rb") : gzdopen(0, "rb");
	if (fp == 0) return 0;
	ks = kseq_init(fp);
	db = kom_calloc(kmerdb_t, 1);
	db->ba = ba_init(0x10000); // 64kb per block
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
		kh_val(db->h, k) = ba_strdup(db->ba, ks->name.s);
	}
	kseq_destroy(ks);
	gzclose(fp);
	return db;
}

void km_find_match(const char *fn, const kmerdb_t *db)
{
	gzFile fp;
	kseq_t *ks;
	kstring_t out = {0,0,0};
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
					out.l = 0;
					if (kf != kh_end(db->h))
						kom_sprintf_lite(&out, "%s\t%ld\t%ld\t%s\t%d\t+\n", ks->name.s, (long)(i + 1 - db->ksize), (long)(i + 1), kh_val(db->h, kf), db->ksize);
					if (kr != kh_end(db->h))
						kom_sprintf_lite(&out, "%s\t%ld\t%ld\t%s\t%d\t-\n", ks->name.s, (long)(i + 1 - db->ksize), (long)(i + 1), kh_val(db->h, kr), db->ksize);
					if (out.l > 0)
						fwrite(out.s, 1, out.l, stdout);
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

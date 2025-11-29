CC=			gcc
CFLAGS=		-g -Wall -O3 -std=c99
CXXFLAGS=	$(CFLAGS)
PROG=		kmer-map-gemini3a kmer-map-gemini3b kmer-map-gpt5a

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -I.. -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

kmer-map-gemini3a:kmer-map-gemini3a.o
		$(CC) -o $@ $< -lz

kmer-map-gemini3b:kmer-map-gemini3b.o
		$(CC) -o $@ $< -lz

kmer-map-gpt5a:kmer-map-gpt5a.o
		$(CC) -o $@ $< -lz

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

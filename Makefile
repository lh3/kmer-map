CC=			gcc
CFLAGS=		-g -Wall -O3 -std=c99
CXXFLAGS=	$(CFLAGS)
PROG=		kmer-map

ifneq ($(asan),)
	CFLAGS+=-fsanitize=address
	LIBS+=-fsanitize=address -ldl
endif

.SUFFIXES:.c .cpp .o
.PHONY:all clean depend

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

kmer-map:kmer-map.o
		$(CC) -o $@ $< -lz

clean:
		rm -fr *.o a.out $(PROG) *~ *.a *.dSYM

depend:
		(LC_ALL=C; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c *.cpp)

# DO NOT DELETE

kmer-map.o: kseq.h khashl.h

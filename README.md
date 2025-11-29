Prompt:

> Write a C program to find the positions of a list of k-mers in a genome. Both
> k-mers and the genome are provided in the FASTA/FASTQ format, optionally
> gzip'd. The k-mers are up to 31bp in length. They are supposed to be of equal
> length and have no duplicates. Give warnings if the input file doesn't meet
> these conditions. Output in the 6-column BED format.
>
> Use kseq.h and zlib for parsing FASTA/FASTQ with seamless gzip support. Use
> khashl.h (not khash.h) for k-mer hash table. Try not to use other
> non-standard libraries. Pay attention to code clarity and code efficiency.
> Conform to C99.

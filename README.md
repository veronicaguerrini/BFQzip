# BFQzip

We propose the first lossy reference-free and assembly-free compression approach for FASTQ files, which combines both DNA bases and quality score information in the reads to smooth the quality scores and to apply a noise reduction of the bases, while keeping variant calling performance comparable to that with original data.

The strategy is based on the Extended Burrows-Wheeler Transform (EBWT) and positional clustering, and we present implementations in both internal memory and external memory.

## Install

```sh
git clone --recursive https://github.com/veronicaguerrini/BFQzip
cd BFQzip 
make
```

## Run

Given a string collection in FASTQ format (e.g. example/reads.fastq), to run BFQzip in internal memory

```sh
python3 BFQzip.py example/reads.fastq -o output_reads
```
while to run BFQzip in external memory

```sh
python3 BFQzip_ext.py example/reads.fastq -o output_reads
```

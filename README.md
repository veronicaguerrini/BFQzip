# BFQzip

We propose the first lossy reference-free and assembly-free compression approach for FASTQ file that combines both DNA bases and quality information in the reads to smooth quality scores and to apply a noise reduction of bases, while keeping variant calling performance comparable to the original data.

The strategy is based on the Extended Burrows-Wheeler Transform (EBWT) and positional clustering.
We present implementations in both internal and external memory.

## Install

```sh
git clone --recursive https://github.com/bfqzip/BFQzip
cd BFQzip 
make
```

## Run

Given a string collection dataset.fastq, to run BFQzip in internal memory

```sh
python3 BFQzip.py dataset.fastq -o output_reads
```
while to run BFQzip in external memory

```sh
python3 BFQzip_ext.py dataset.fastq -o output_reads
```

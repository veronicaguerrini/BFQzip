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

#### References:

    *** BFQzip
    
    Veronica Guerrini, Felipe A. Louza, Giovanna Rosone
    Lossy Compressor preserving variant calling through Extended BWT
    BIOINFORMATICS 2022 (accepted) 
    13th International Conference on Bioinformatics Models, Methods and Algorithms

---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html

<small> Supported by the University of Pisa under the ``PRA – Progetti di Ricerca di Ateneo'' (Institutional Research Grants) - Project no. PRA\_2020\-2021\_26 ``Metodi Informatici Integrati per la Biomedica''.

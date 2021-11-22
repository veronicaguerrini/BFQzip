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

    *** EBWT
    
    Sabrina Mantaci, Antonio Restivo, Giovanna Rosone, Marinella Sciortino,
    An extension of the Burrows-Wheeler Transform.
    Theoretical Computer Science (2007) 387(3): 298-312,
    doi: 10.1016/j.tcs.2007.07.014
    
    Markus J. Bauer, Anthony J. Cox, Giovanna Rosone,
    Lightweight algorithms for constructing and inverting the BWT of string collections. 
    Theoretical Computer Science (2013) 483: 134-148,
    doi: 10.1016/j.tcs.2012.02.002
    
    *** Positional Clustering
    
    Nicola Prezza, Nadia Pisanti, Marinella Sciortino, Giovanna Rosone,
    SNPs detection by eBWT positional clustering,
    Algorithms for Molecular Biology (2019) 14(1):3,
    doi: 10.1186/s13015-019-0137-8
    
    Nicola Prezza, Nadia Pisanti, Marinella Sciortino, Giovanna Rosone,
    Variable-order reference-free variant discovery with the Burrows-Wheeler Transform.
    BMC Bioinformatics (2020) 21,
    doi: 10.1186/s12859-020-03586-3
    
    *** BFQzip
    
    Veronica Guerrini, Felipe A. Louza, Giovanna Rosone,
    Lossy Compressor preserving variant calling through Extended BWT
    BIOINFORMATICS 2022 (accepted) 
    13th International Conference on Bioinformatics Models, Methods and Algorithms

---
<small> Supported by the project Italian MIUR-SIR [CMACBioSeq][240fb5f5] ("_Combinatorial methods for analysis and compression of biological sequences_") grant n.~RBSI146R5L. P.I. Giovanna Rosone</small>

[240fb5f5]: http://pages.di.unipi.it/rosone/CMACBioSeq.html

<small> Supported by the University of Pisa under the "_PRA – Progetti di Ricerca di Ateneo_" (Institutional Research Grants) - Project no. PRA 2020-2021 "_Metodi Informatici Integrati per la Biomedica_".</small>

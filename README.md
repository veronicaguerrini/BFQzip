# BFQzip

We propose the first **lossy** **reference-free** and **assembly-free** compression approach for FASTQ files, which combines both DNA bases and quality score information in the reads to smooth the quality scores and to apply a noise reduction of the bases, while keeping variant calling performance comparable to that with original data.

The strategy is based on the Extended Burrows-Wheeler Transform (**EBWT**) and the **positional clustering** framework, and it can be summarized in four main steps:

1. *Data structures building*, for which one could use any state-of-the art tool that computes both the EBWT and its associated permutation of quality scores,

2. *Positional cluster detecting*, according to the positions of local minima in the LCP array,

3. *Noise reduction and quality score smoothing*, 

4. *FASTQ reconstruction* and *Compression*.

Given as input a FASTQ file containing a collection *S* of reads, in step 1 we build the ebwt(*S*) (EBWT output string) and the string qs(*S*), which is the concatenation of quality scores associated with the symbols appearing in the ebwt(*S*). In step 2, we exploit the positional clustering framework to detect blocks in the ebwt(*S*), and thus in qs(*S*).
In step 3, these blocks allow us not only to smooth their quality scores, but also the corresponding bases, replacing those that are believed to be noise introduced during the sequencing process. In step 4, we use the LF mapping on the ebwt(*S*) to output a new (modified) FASTQ file and compress it by using state-of-the-art compressors. 

We present implementations in both **internal memory** and **external memory**.

Note that Step 1 can be performed by any tool according to the resources available. For example, [gsufsort](https://github.com/felipelouza/gsufsort) runs in internal memory, while [egap](https://github.com/felipelouza/egap) and [BCR](https://github.com/giovannarosone/BCR_LCP_GSA) run in external memory.
Both implementations of step 2 need the ebwt(*S*) and qs(*S*), while the external memory version needs in addition the LCP array lcp(*S*). 
Indeed, the internal memory approach represents ebwt(*S*) via the compressed suffix tree described in [Prezza and Rosone, 2021](https://doi.org/10.1016/j.tcs.2020.11.024), where the lcp(*S*) is deduced from the ebwt(*S*). 
During the FASTQ reconstruction of step 4, the LF-mapping is implemented either in internal memory (via suffix-tree navigation) or in external memory.

## Install

```sh
git clone --recursive https://github.com/veronicaguerrini/BFQzip
cd BFQzip 
make
```

## Run

We propose three different modes to compress FASTQ files by using two well-known compressors: [PPMd](https://www.7-zip.org/7z.html) and [BSC](http://libbsc.com/).

Mode 1 (option *-1* or *--m1*) consists in compressing the whole FASTQ file reconstructed (bases and quality scores components are interleaved as in the FASTQ format). Note that we are not interested in compressing the headers component, for which one can use any state-of-the-art strategy that tokenizes it. Nevertheless, to report the original headers component in the FASTQ file reconstructed, please use option *--headers*.

Mode 2 (option *-2* or *--m2*) consists in compressing the bases and the quality scores components separately. 

Mode 3 (option *-3* or *--m3*) consists in compressing the bases, the quality scores, and the headers components separately.

Then, given a string collection in FASTQ format (e.g. example/reads.fastq), to compress all its components separately, run BFQzip

- in internal memory 

```sh
python3 BFQzip.py example/reads.fastq -o output_reads --m3
```
-or, in external memory

```sh
python3 BFQzip_ext.py example/reads.fastq -o output_reads --m3
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
    
    *** QS permutation
    
    Lilian Janin, Giovanna Rosone, and Anthony J. Cox: 
    Adaptive reference-free compression of sequence quality scores. 
    Bioinformatics (2014) 30 (1): 24-30, 
    doi: 10.1093/bioinformatics/btt257
    
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

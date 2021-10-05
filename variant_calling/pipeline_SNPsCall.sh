#!/bin/bash
####### Script for performing GATK best practices Variant Calling ########


##Please insert the correct file paths
### Software tools
bwa_program="./bwa-0.7.17/bwa"
gatk_program="./gatk-4.2.0.0/gatk"

## Reference FASTA
ref="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

#######

## input args
fastq_1=$1
fastq_2=$2
data=$(echo $fastq_1| cut -d'_' -f 1)


#run bwa index if files do not exist
if test ! -f $ref.amb;
then
	$bwa_program index $ref > "bwa_index_"$ref".stdout" 2> "bwa_index_"$ref".stdout"
fi


##Step 1 Alignment – Map to Reference 
$bwa_program mem -Y -R '@RG\tID:sample1\tLB:lib1\tPL:ILLUMINA\tPM:HISEQ\tSM:sample1' $ref $fastq_1 $fastq_2 > "aligned_"$data".sam" 2> "bwa_mem_"$data".stderr"

##Step 2 Mark Duplicates + Sort
$gatk_program MarkDuplicatesSpark -I "aligned_"$data".sam" -M "dedup_"$data".txt" -O "sorted_dedup_"$data".bam" 2> "gatk_mark_"$data".stderr"


#only run if .fai file does not exist
if test ! -f $ref.fai;
then
	samtools faidx $ref
fi

$gatk_program CreateSequenceDictionary -R $ref

##Step 3 Call Variants
$gatk_program HaplotypeCaller -R $ref -I "sorted_dedup_"$data".bam" -O "raw_variants_"$data".vcf" 2> "gatk_hapl_"$data".stderr"

##Step 5 Extract SNPs
$gatk_program SelectVariants -R $ref -V "raw_variants_"$data".vcf" -select-type SNP -O "raw_snps_"$data".vcf" 2> "gatk_select_"$data".stderr"

##Step 6 Filter SNPs
$gatk_program VariantFiltration -R $ref -V "raw_snps_"$data".vcf" -O "filtered_snps_"$data".vcf" -filter-name "QD_filter" -filter "QD < 2.0" -filter-name "FS_filter" -filter "FS > 60.0" -filter-name "MQ_filter" -filter "MQ < 40.0" -filter-name "SOR_filter" -filter "SOR > 4.0" -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" 2> "gatk_filtr_"$data".stderr"







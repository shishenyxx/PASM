# Snakemake pipeline for annotations and the calculations of exact binomial confidence intervals

## The annotation pipeline takes a list of candidate variants as input, calculates and annotates the information needed for downstream analysis and calculate the number of reference counts, alternative counts, mutant allelic fractions (MAFs), and the exact binomial confidence intervals of the MAFs from the "tumor" bam. Details about the input and output are listed below:

The pipeline is initially written by Renee D. George, and re-implemented and re-worte by Xin Xu, with help form Xiaoxu Yang and Martin Breuss.
All rights reserved.

## Input:
#### sample
User defiend name for the specific sample.
#### tumor
The name of the "tumor" sample from the tags in the bam header.
#### normal
The name of the "normal" sample from the tags in the bam header.
#### tumor_path
Path to the "tumor" bam file.
#### normal_path
Path to the "normal" bam file.
#### vcf_path
Path to the list of variants you want to annotate and calculate for the "tumor" bam. If you want to compare the same variant in tumor and normal, add another input row and switch the "tumor" and "normal".
#### gvcf_path
Path to the variant quality score recalibrated gvcfs from haplotype caller. This will help the calculation such as "near indel".

## Output:
  ID
  CHROM
  POS
  REF
  ALT
  ANNO
  GENE
  GNOMAD_FREQ
  REPEAT_MASKER
  SEGDUP
  HOMOPOLYMER
  REF_SEQ
  DINUCLEOTIDE
  NEAR_INDEL
  UCSC_RPMSK
  REF_COUNT
  ALT_COUNT
  MAF
  LOWER_CI
  UPPER_CI
  CI_IS_GREATER
  NORMAL_REF_COUNT
  NORMAL_ALT_COUNT
  NORMAL_MAF
  NORMAL_LOWER_CI
  NORMAL_UPPER_CI
  NORMAL_CI_IS_GREATER
  TUMOR_IS_BLOOD
  TUMOR_IS_SPERM
  
  

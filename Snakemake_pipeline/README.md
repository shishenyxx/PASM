# Snakemake pipeline for annotations and the calculations of exact binomial confidence intervals

## The annotation pipeline takes a list of candidate variants as input, calculates and annotates the information needed for downstream analysis and calculate the number of reference counts, alternative counts, mutant allelic fractions (MAFs), and the exact binomial confidence intervals of the MAFs from the "tumor" bam. Details about the input and output are listed below:

The pipeline is initially written by Renee D. George, and re-implemented and re-worte by Xin Xu, with help form Xiaoxu Yang and Martin Breuss.
All rights reserved.

## Before start, you have to install:
#### [Annovar](https://doc-openbio.readthedocs.io/projects/annovar/en/latest/)

## Input:
These are headers of the input file list.
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

## Config files:
These are files you need to prepare for the annotation scripts, saved in the file snake_conf.yaml
#### input_files
Path to the input file list.
#### out_dir
Path to the output directory.
#### scratch_dir
Path to the scratch files. Note that the number of temporary files will be euqals to two- or three-times the number of total listed vriants.
#### annovar
Path to your annovar annotate_variation.pl.
#### annovar_db
Path to your annovar databases.
#### homopolymer_script, ci_script, and summerize_script
Path to the helper scripts, by default the scripts are provided in the helper_scripts folder.
#### bed_file
The total length of the chromosome corresponding to your reference genome file.
#### ref_fasta
Your reference genome.
#### gnomad_af
vcf.gz file with only AF information from gnomAD (corresponding to your reference genome file).
#### ucsc_rpmsk
Repeatmask track from UCSC genome browser (corresponding to your reference genome file). the default format is #bin	swScore	milliDiv	milliDel	milliIns	genoName	genoStart	genoEnd	genoLeft	strand	repName	repClass	repFamily	repStart	repEndrepLeft	id
#### repeat_masker
Additional genomic repeats annotated to the output table in bed format (corresponding to your reference genome file).
#### segdup
Segmental duplications annotated to the output table in bed format (corresponding to your reference genome file).

## Output:
#### ID
The name of the "tumor" bam tag where the information of this line is calculated from.
#### CHROM
Chromosome number of each variant.
#### POS
Genomic position of each variant.
#### REF
Reference allele from according to the reference genome file given as "ref_fasta".
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
  
  

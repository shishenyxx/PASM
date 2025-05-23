# PASM, TAS/TASeq, and MPAS/snMPAS

## 1. Overview:

Here are scripts initially written for Postzygotic Amplicon Sequencing for Mosaicism (PASM), Targeted Amplicon Sequencing (TAS/TASeq), and some codes for the method we now define as Massive Parallel Amplicon Sequencing (MPAS). We provided a Perl+R version, two standalone Python versions, and a Snakemake pipeline. The scripts and pipelines are useful for the calculation of variant allelic fraction (AF) based on not only amplicon-based deep sequencing data but also the AF estimation as well as annotations for postzygotic mosaic variant studies from all kinds of Next Generation Sequencing (NGS) data.

--------------
## 2. Versions and updates:

For the calculation of confidence intervals, you can choose [exact binomial confidence interval in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/binom.test) ([Clopper-Pearson interval](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval) by default), or the iterative methods which consider the base qualities from each base, described in <i>[Xu, Yang, and Wu et al. Wei and Zhang. 2015](https://doi.org/10.1002/humu.22819)</i> and <i>[Yang and Liu et al. Wu, Wei, and Zhang. 2017](https://doi.org/10.1038/s41598-017-15814-7)</i>, different versions of scripts are available.

### Snakemake pipelines:

A [Snakemake pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) with exact binomial CIs and detailed annotations was implemented by Xin Xu and Xiaoxu Yang, with great input from Martin Breuss, the pipeline is based on the Python scripts written by Xin Xu and a previous Snakemake pipeline written by Martin Breuss and Renee D. George. (2019-08-12)



### Python versions:

A [new python version](https://github.com/shishenyxx/PASM/tree/master/2019-09-25-new-python-MAF-binom-calculator) with exact binomial CIs was implemented by Xin Xu supervised by Xiaoxu Yang (2019-07-24) and fixed by Jiawei Shen (2022-04-22).


A [python version](https://github.com/shishenyxx/PASM/blob/master/CI_calculator.py) was implemented by Xianing Zheng and supervised by Xiaoxu Yang. (2016-04-17)


### Perl and R versions:
 #### Before you start the Perl + R version:
Note that the Perl package [Statistics::R](https://metacpan.org/pod/Statistics::R) is used to call the [yyxMosaicHunter](https://github.com/Yyx2626/yyxMosaicHunter) package in R written by Adam Yongxin Ye.
<br/>Dependencies of yyxMosaicHunter 0.1.4 are: `Rcpp`
`pryr`. (2014-11-11)

 #### Instructions for the Perl + R version:
The [first part of the Perl version](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_calculate_only_2016_12_03.pl) is a pileup filter, it takes in SAMTools mpileup results and calculates different characters to count the bases, written by Jiarui Li, modified by Xiaoxu Yang and Xianing Zheng. (2015-03-24)

You can also only [output the base quality](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03_output_basequality.pl) and deal with the base qualities in R. 


If you want to calculate the CIs with PASM Bayesian model, you can use [this Perl script](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03.pl), or a [older version perl script](https://github.com/shishenyxx/PASM/blob/master/old_get_ref_alt_baseQ_corrected_2016_07_14.pl). 

--------------

## 3. Example usage:

For the Perl version: `samtools mpileup -r ${chr}:${pos}-${pos} -f <reference_file> -Q 0 -q 0 -AB -d 5000000 <input_bam> | ./get_ref_alt_baseQ_corrected_calculate_only_2016_12_03.pl`

--------------

## 4. Experimental design:


![image](https://github.com/user-attachments/assets/86e1da03-23a3-4537-975e-6b7f831e3ea0)

1. Primers for the amplicons are designed based on the [Primer3 Command Line version](https://github.com/shishenyxx/primer3), or [CREPE](https://www.biorxiv.org/content/10.1101/2025.03.28.646040v1). Alternative homozygous SNPs with 0 reference bases from preliminary sequencing data are recommended as negative controls. Heterozygous SNPs with bulk AF around 50% or known heterozygous in the proband can be used as positive controls.
2. After amplicon design, libraries are made for Illumina or Ion Torrent platforms.
3. After sequencing, data is mapped, indel realigned and the base quality scores recalibrated.
4. Decision boundaries are based on a] the 1-(95% binomial lower CI of alt homo) (95% percentile for an FDR of 0.05 for homozygous); b] assuming no reverse mutation, the 95% binomial lower CI of hets (95% percentile for an FDR of 0.05 for heterozygous).

Percentiles for decision boundaries could be adjusted to get different accuracy.
   
--------------

## 5. Related publications:

* [Identification of Novel Mosaic Variants in Focal Epilepsy-Associated Patients’ Brain Lesions](https://www.mdpi.com/2073-4425/16/4/421) 2025. <i>(Genes)</i>
* [Cell-type-resolved mosaicism reveals clonal dynamics of the human forebrain](https://www.nature.com/articles/s41586-024-07292-5) 2024. <i>(Nature)</i>
* [Control-independent mosaic single nucleotide variant detection with DeepMosaic.](https://www.nature.com/articles/s41587-022-01559-w) 2023. <i>(Nature Biotechnology)</i>
* [Comprehensive multi-omic profiling of somatic mutations in malformations of cortical development.](https://doi.org/10.1038/s41558-022-01276-9) 2023. <i>(Nature Genetics)</i>
* [Somatic mosaicism reveals clonal distributions of neocortical development.](https://www.nature.com/articles/s41586-022-04602-7) 2022. <i>(Nature)</i>
* [Unbiased mosaic variant assessment in sperm: a cohort study to test predictability of transmission.](https://elifesciences.org/articles/78459) 2022. <i>(eLife)</i>
* [Developmental and temporal characteristics of clonal sperm mosaicism.](http://www.doi.org/10.1016/j.cell.2021.07.024) 2021. <i>(Cell)</i>
* [Comprehensive identification of somatic nucleotide variants in human brain tissue.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02285-3) 2021. <i>(Genome Biology)</i>
* [Autism risk in offspring can be assessed through quantification of male sperm mosaicism.](https://doi.org/10.1038/s41591-019-0711-0) 2020. <i>(Nature Medicine)</i>
* [mTOR pathway somatic variants and the molecular pathogenesis of hemimegalencephaly.](https://doi.org/10.1002/epi4.12377) 2020. <i>(Epilepsia Open)</i>
* [Genomic mosaicism in the pathogenesis and inheritance of a Rett syndrome cohort.](https://doi.org/10.1038/s41436-018-0348-2) 2019. <i>(Genetics in Medicine)</i>
* [Mosaicism and incomplete penetrance of <i>PCDH19</i> mutations.](http://dx.doi.org/10.1136/jmedgenet-2017-105235) 2019. <i>(Journal of Medical Genetics)</i>
* [Somatic double-hit in <i>MTOR</i> and <i>RPS6</i> in hemimegalencephaly with intractable epilepsy.](https://doi.org/10.1093/hmg/ddz194) 2019. <i>(Human Molecular Genetics)</i>
* [<i>ATP1A3</i> mosaicism in families with alternating hemiplegia of childhood.](https://doi.org/10.1111/cge.13539) 2019. <i>(Clinical Genetics)</i>
* [Distinctive types of postzygotic single-nucleotide mosaicisms in healthy individuals revealed by genome-wide profiling of multiple organs.](https://doi.org/10.1371/journal.pgen.1007395) 2018. <i>(PLoS Genetics)</i>
* [MosaicHunter: accurate detection of postzygotic single-nucleotide mosaicism through next-generation sequencing of unpaired, trio, and paired samples.](https://doi.org/10.1093/nar/gkx024) 2017. <i>(Nucleic Acids Research)</i>
* [Ultrasensitive and high-efficiency screen of de novo low-frequency mutations by o2n-seq.](https://doi.org/10.1038/ncomms15335) 2017. <i>(Nature Communications)</i>
* [Postzygotic single‐nucleotide mosaicisms contribute to the etiology of autism spectrum disorder and autistic traits and the origin of mutations.](https://doi.org/10.1002/humu.23255) 2017. <i>(Human Mutation)</i>
* [Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort.](https://doi.org/10.1038/s41598-017-15814-7) 2017. <i>(Scientific Reports)</i>
* [Amplicon resequencing identified parental mosaicism for approximately 10% of <i>“de novo” SCN1A</i> mutations in children with Dravet syndrome.](https://doi.org/10.1002/humu.22819) 2015. <i>(Human Mutation)</i>

-----------------------------------

## 6. Contributors:

Xiaoxu Yang, Jiarui Li, Adam Yongxin Ye, Xianing Zheng, Sheng Wang, Martin W Breuss, Xin Xu, Weizhen Zhou.

-----------------------------------

## 7. Contact:

:email: Xiaoxu Yang: [xiaoxu.yang@genetics.utah.edu](mailto:xiaoxu.yang@genetics.utah.edu), [xiaoxuyanglab@gmail.com](mailto:xiaoxuyanglab@gmail.com)

-----------------------------------

## 8. Cite the code:
* Cite the Python version or the snakemake wrapper and the binomial model:

    Yang X*,#, Xu X*, <i>et al.</i> [Control-independent mosaic single nucleotide variant detection with DeepMosaic.](https://www.nature.com/articles/s41587-022-01559-w) 2023. <i>Nature Biotechnology</i> 

    Yang X* and Breuss MW*, <i>et al.</i> [Developmental and temporal characteristics of clonal sperm mosaicism.](http://www.doi.org/10.1016/j.cell.2021.07.024) 2021. <i>Cell</i> 
           
* Cite the Perl+R version and the Bayesian model: 

    Yang X* and Liu A*, <i>et al.</i> [Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort.](https://doi.org/10.1038/s41598-017-15814-7) 2017. <i>Scientific Reports</i>

    Xu X*, Yang X*, and Wu Q*, <i>et al.</i> [Amplicon resequencing identified parental mosaicism for approximately 10% of <i>“de novo” SCN1A</i> mutations in children with Dravet syndrome.](https://doi.org/10.1002/humu.22819) 2015. <i>Human Mutation</i>


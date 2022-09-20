# PASM and MPAS

## 1. Overview:

Here are scripts initially written for Postzygotic Amplicon Sequencing for Mosaicism (PASM) and some codes for the method we now define as Massive Parallel Amplicon Sequencing (MPAS). We provided a perl+R version, two standalone python versions, and a Snakemake pipeline. The scripts and pipelines are useful for the calculation of variant allelic fraction (AF) based on not only amplicon based deep sequencing data, but also the AF estimation as well as annotations for postzygotic mosaic variant studies from all kinds of Next Generation Sequencing (NGS) data.

--------------
## 2. Versions and updates:

For the calculation of confidence intervals, you can choose [exact binomial confidence interval in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/binom.test) ([Clopper-Pearson interval](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval) by default), or the iterative methods which considers the base qualities from each base, described in <i>[Xu , Yang, and Wu et al. Wei and Zhang. 2015](https://doi.org/10.1002/humu.22819)</i> and <i>[Yang and Liu et al. Wu, Wei, and Zhang. 2017](https://doi.org/10.1038/s41598-017-15814-7)</i>, different versions of scripts are available.

### Snakemake pipelines:

A [Snakemake pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) with exact binomial CIs and detailed annotations was implemented by Xin Xu and Xiaoxu Yang, with great input form Martin Breuss, the pipeline is based on the Python scripts written by Xin Xu and a previous Snakemake pipeline written by Martin Breuss and Renee D. George. (2019-08-12)



### Python versions:

A [new python version](https://github.com/shishenyxx/PASM/tree/master/2019-09-25-new-python-MAF-binom-calculator) with exact binomial CIs was implemented by Xin Xu supervised by Xiaoxu Yang (2019-07-24) and fixed by Jiawei Shen (2022-04-22).


A [python version](https://github.com/shishenyxx/PASM/blob/master/CI_calculator.py) was implemented by Xianing Zheng supervised by Xiaoxu Yang. (2016-04-17)


### Perl and R versions:
 #### Before you start the Perl + R version:
Note that the Perl package [Statistics::R](https://metacpan.org/pod/Statistics::R) is used to call the [yyxMosaicHunter](https://github.com/Yyx2626/yyxMosaicHunter) package in R written by Adam Yongxin Ye.
<br/>Dependencies of yyxMosaicHunter 0.1.4 are: `Rcpp`
`pryr`. (2014-11-11)

 #### Instructions for the Perl + R version:
The [first part of the perl version](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_calculate_only_2016_12_03.pl) is a pileup filter, it takes in SAMTools mpileup results and calculate different characters to count the bases, written by Jiarui Li, modified by Xiaoxu Yang and Xianing Zheng. (2015-03-24)

You can also only [output the base qualtiy](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03_output_basequality.pl) and deal with the base qualities in R. 


If you want to calculate the CIs with PASM Bayesian model, you can use [this perl script](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03.pl), or a [older version perl script](https://github.com/shishenyxx/PASM/blob/master/old_get_ref_alt_baseQ_corrected_2016_07_14.pl). 

--------------



## 3. Example usage:

For the Perl version: `samtools mpileup -r ${chr}:${pos}-${pos} -f <reference_file> -Q0 -q0 -AB -d3000 <input_bam> | ./get_ref_alt_baseQ_corrected_calculate_only_2016_12_03.pl`

--------------

## 4. Experimental design:

Primers for the amplicons are designed based on the [Primer3 Command Line version](https://github.com/shishenyxx/primer3).

--------------

## 5. Related publications:

* [Somatic mosaicism reveals clonal distributions of neocortical development.](https://www.nature.com/articles/s41586-022-04602-7) 2022. <i>(Nature)</i>
* [Comprehensive multiomic profiling of somatic mutations in malformations of cortical development.](https://www.biorxiv.org/content/10.1101/2022.04.07.487401v2.full) 2022. <i>(bioRxiv)</i>
* [Unbiased mosaic variant assessment in sperm: a cohort study to test predictability of transmission.](https://elifesciences.org/articles/78459) 2022. <i>(eLife)</i>
* [Developmental and temporal characteristics of clonal sperm mosaicism.](http://www.doi.org/10.1016/j.cell.2021.07.024) 2021. <i>(Cell)</i>
* [Comprehensive identification of somatic nucleotide variants in human brain tissue.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02285-3) 2021. <i>(Genome Bio.)</i>
* [DeepMosaic: Control-independent mosaic single nucleotide variant detection using deep convolutional neural networks.](https://www.biorxiv.org/content/10.1101/2020.11.14.382473v2.full) 2021. <i>(bioRxiv)</i>
* [Autism risk in offspring can be assessed through quantification of male sperm mosaicism.](https://doi.org/10.1038/s41591-019-0711-0) 2020. <i>(Nat. Med.)</i>
* [mTOR pathway somatic variants and the molecular pathogenesis of hemimegalencephaly.](https://doi.org/10.1002/epi4.12377) 2020. <i>(Epi. Open)</i>
* [Genomic mosaicism in the pathogenesis and inheritance of a Rett syndrome cohort.](https://doi.org/10.1038/s41436-018-0348-2) 2019. <i>(Genet. Med.)</i>
* [Mosaicism and incomplete penetrance of <i>PCDH19</i> mutations.](http://dx.doi.org/10.1136/jmedgenet-2017-105235) 2019. <i>(J. Med. Genet.)</i>
* [Somatic double-hit in <i>MTOR</i> and <i>RPS6</i> in hemimegalencephaly with intractable epilepsy.](https://doi.org/10.1093/hmg/ddz194) 2019. <i>(Hum. Mol. Genet.)</i>
* [<i>ATP1A3</i> mosaicism in families with alternating hemiplegia of childhood.](https://doi.org/10.1111/cge.13539) 2019. <i>(Clin. Genet.)</i>
* [Distinctive types of postzygotic single-nucleotide mosaicisms in healthy individuals revealed by genome-wide profiling of multiple organs.](https://doi.org/10.1371/journal.pgen.1007395) 2018. <i>(PLoS Genet.)</i>
* [MosaicHunter: accurate detection of postzygotic single-nucleotide mosaicism through next-generation sequencing of unpaired, trio, and paired samples.](https://doi.org/10.1093/nar/gkx024) 2017. <i>(Nucleic acids res.)</i>
* [Ultrasensitive and high-efficiency screen of de novo low-frequency mutations by o2n-seq.](https://doi.org/10.1038/ncomms15335) 2017. <i>(Nat. Comms.)</i>
* [Postzygotic single‐nucleotide mosaicisms contribute to the etiology of autism spectrum disorder and autistic traits and the origin of mutations.](https://doi.org/10.1002/humu.23255) 2017. <i>(Hum. Mutat.)</i>
* [Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort.](https://doi.org/10.1038/s41598-017-15814-7) 2017. <i>(Sci. Rep.)</i>
* [Amplicon resequencing identified parental mosaicism for approximately 10% of <i>“de novo” SCN1A</i> mutations in children with Dravet syndrome.](https://doi.org/10.1002/humu.22819) 2015. <i>(Hum. Mutat.)</i>

-----------------------------------

## 6. Contact:

:email: Dr. Xiaoxu Yang: [xiy010@health.ucsd.edu](mailto:xiy010@health.ucsd.edu), [yangxiaoxu-shishen@hotmail.com](mailto:yangxiaoxu-shishen@hotmail.com)



-----------------------------------

## 7. Cite the code:
* Cite the new Python version or the snakemake wrapper:

    Yang X and Breuss MW, <i>et al.</i> Gleeson JG. [Developmental and temporal characteristics of clonal sperm mosaicism.](http://www.doi.org/10.1016/j.cell.2021.07.024) 2021. <i>Cell</i> 

* Cite the Perl+R version: 

    Yang X and Liu A, <i>et al.</i> Wei L and Zhang Y. [Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort.](https://doi.org/10.1038/s41598-017-15814-7) 2017. <i>Sci. Rep.</i>

    Xu X, Yang X, and Wu Q, <i>et al.</i> Wei L and Zhang Y. [Amplicon resequencing identified parental mosaicism for approximately 10% of <i>“de novo” SCN1A</i> mutations in children with Dravet syndrome.](https://doi.org/10.1002/humu.22819) 2015. <i>Hum. Mutat.</i>


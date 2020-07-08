# PASM
--------------
## Overview
Here are scripts initially written for Postzygotic Amplicon Sequencing for Mosaicism (PASM). We provided a perl+R version, two standalone python versions, and a Snakemake pipeline. The scripts and pipelines are useful for the calculation of variant allelic fraction (AF) based on not only amplicon based deep sequencing data, but also the AF estimation as well as annotations for postzygotic mosaic variant studies from all kinds of different NGS data.

--------------
## Versions and updates
For the calculation of confidence intervals, you can choose [exact binomial confidence interval in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/binom.test) ([Clopper-Pearson interval](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval) by default), or the methods described in [Yang and Liu et al. 2017](https://doi.org/10.1038/s41598-017-15814-7).

The [first part](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_calculate_only_2016_12_03.pl) is a pileup filter, it takes in SAMTools mpileup results and calculate different characters to count the bases, written by Jiarui Li, modified by Xiaoxu Yang and Xianing Zheng. (2015-03-24)

--------------

You can also only [output the base qualtiy](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03_output_basequality.pl) and deal with the base qualities in R. 

--------------

If you also want to calculate the CIs with PASM Bayesian model, you can use [this script](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03.pl), or a [older version](https://github.com/shishenyxx/PASM/blob/master/old_get_ref_alt_baseQ_corrected_2016_07_14.pl). Note that the Perl package Statistics::R is used to call the [yyxMosaicHunter](https://github.com/Yyx2626/yyxMosaicHunter) package in R written by Adam Yongxin Ye.
Dependencies of yyxMosaicHunter 0.1.4 are: `Rcpp`
`pryr`. (2014-11-11)

--------------

A [python version](https://github.com/shishenyxx/PASM/blob/master/CI_calculator.py) was implemented by Xianing Zheng. (2016-04-17)

--------------

A [new python version](https://github.com/shishenyxx/PASM/blob/master/2019-09-25-new-python-MAF-binom-calculator/compute_maf_binom.py) with exact binomial CIs was implemented by Xin Xu. (2019-07-24)

--------------

A [Snakemake pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) with exact binomial CIs was implemented by Xin Xu. (2019-08-12)

--------------
## Related publications:
* [Autism risk in offspring can be assessed through quantification of male sperm mosaicism](https://doi.org/10.1038/s41591-019-0711-0) <i>(Nat. Med. 2020)</i>
* [mTOR pathway somatic variants and the molecular pathogenesis of hemimegalencephaly](https://doi.org/10.1002/epi4.12377) <i>(Epi. Open 2020)</i>
* [Genomic mosaicism in the pathogenesis and inheritance of a Rett syndrome cohort](https://doi.org/10.1038/s41436-018-0348-2) <i>(Genet. Med. 2019)</i>
* [Mosaicism and incomplete penetrance of PCDH19 mutations](http://dx.doi.org/10.1136/jmedgenet-2017-105235) <i>(J. Med. Genet. 2019)</i>
* [Somatic double-hit in MTOR and RPS6 in hemimegalencephaly with intractable epilepsy](https://doi.org/10.1093/hmg/ddz194) <i>(Hum. Mol. Genet. 2019)</i>
* [ATP1A3 mosaicism in families with alternating hemiplegia of childhood](https://doi.org/10.1111/cge.13539) <i>(Clin. Genet. 2019)</i>
* [Distinctive types of postzygotic single-nucleotide mosaicisms in healthy individuals revealed by genome-wide profiling of multiple organs](https://doi.org/10.1371/journal.pgen.1007395) <i>(PLoS Genet. 2018)</i>
* [Ultrasensitive and high-efficiency screen of de novo low-frequency mutations by o2n-seq](https://doi.org/10.1038/ncomms15335) <i>(Nat. Comms. 2017)</i>
* [Postzygotic single‐nucleotide mosaicisms contribute to the etiology of autism spectrum disorder and autistic traits and the origin of mutations.](https://doi.org/10.1002/humu.23255) <i>(Hum. Mutat. 2017)</i>
* [Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort](https://doi.org/10.1038/s41598-017-15814-7) <i>(Sci. Rep. 2017)</i>
* [Amplicon resequencing identified parental mosaicism for approximately 10% of “de novo” SCN1A mutations in children with Dravet syndrome.](https://doi.org/10.1002/humu.22819) <i>(Hum. Mutat. 2015)</i>

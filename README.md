# PASM
Here are scripts initially written for Postzygotic Amplicon Sequencing for Mosaicism (PASM). We provided a perl+R version, two standalone python versions, and a Snakemake pipeline. The scripts and pipelines are useful for the calculation of variant allelic fraction (AF) based on not only amplicon based deep sequencing data, but also the AF estimation as well as annotations for postzygotic mosaic variant studies from all kinds of different NGS data.

--------------
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
[Amplicon resequencing identified parental mosaicism for approximately 10% of “de novo” SCN1A mutations in children with Dravet syndrome (2015)](https://doi.org/10.1002/humu.22819)

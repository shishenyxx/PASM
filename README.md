# PASM
## Here are scripts for PGM amplicon sequencing for Mosaicism (PASM)
--------------

The [first part](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_calculate_only_2016_12_03.pl) is a pileup filter, it takes in SAMTools mpileup results and calculate different characters to count the bases, written by Jiarui Li, modified by Xiaoxu Yang and Xianing Zheng.

--------------

You can also only [output the base qualtiy](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03_output_basequality.pl) and deal with the base qalities in R.

--------------

If you also want to calculate the CIs with PASM Bayesian model, you can use [this script](https://github.com/shishenyxx/PASM/blob/master/get_ref_alt_baseQ_corrected_2016_12_03.pl), or a [older version](https://github.com/shishenyxx/PASM/blob/master/old_get_ref_alt_baseQ_corrected_2016_07_14.pl). Note that the Perl package Statistics::R is used to call the [yyxMosaicHunter](https://github.com/Yyx2626/yyxMosaicHunter) package in R written by Adam Yongxin Ye.
Dependencisof yyxMosaicHunter 0.1.4 are: Rcpp
pryr

--------------

A [python version](https://github.com/shishenyxx/PASM/blob/master/CI_calculator.py) was implemented by Xianing Zheng.

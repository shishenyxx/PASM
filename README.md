# PASM and MPAS

## 1. Overview
Here are scripts initially written for Postzygotic Amplicon Sequencing for Mosaicism (PASM) and some codes for the method we now define as Massive Parallel Amplicon Sequencing (MPAS). We provided a perl+R version, two standalone python versions, and a Snakemake pipeline. The scripts and pipelines are useful for the calculation of variant allelic fraction (AF) based on not only amplicon based deep sequencing data, but also the AF estimation as well as annotations for postzygotic mosaic variant studies from all kinds of Next Generation Sequencing (NGS) data.

--------------
## 2. Versions and updates
For the calculation of confidence intervals, you can choose [exact binomial confidence interval in R](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/binom.test) ([Clopper-Pearson interval](https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval) by default), or the iterative methods which considers the base qualities from each base, described in <i>[Xu , Yang, and Wu et al. Wei and Zhang. 2015](https://doi.org/10.1002/humu.22819)</i> and <i>[Yang and Liu et al. Wu, Wei, and Zhang. 2017](https://doi.org/10.1038/s41598-017-15814-7)</i>, different versions of scripts are available.

### Snakemake pipelines:

A [Snakemake pipeline](https://github.com/shishenyxx/PASM/tree/master/Snakemake_pipeline) with exact binomial CIs and detailed annotations was implemented by Xin Xu and Xiaoxu Yang, with great input form Martin Breuss, the pipeline is based on the Python scripts written by Xin Xu and a previous Snakemake pipeline written by Martin Breuss and Renee D. George. (2019-08-12)



### Python versions:

A [new python version](https://github.com/shishenyxx/PASM/tree/master/2019-09-25-new-python-MAF-binom-calculator) with exact binomial CIs was implemented by Xin Xu and Xiaoxu Yang. (2019-07-24)


A [python version](https://github.com/shishenyxx/PASM/blob/master/CI_calculator.py) was implemented by Xianing Zheng with help from Xiaoxu Yang. (2016-04-17)


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



## 3. Related publications:
* Developmental and temporal characteristics of clonal sperm mosaicism. <i>([Cell 2021](http://www.doi.org/10.1016/j.cell.2021.07.024 ))</i>
* Comprehensive identification of somatic nucleotide variants in human brain tissue. <i>([Genome Bio. 2021](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02285-3))</i>
* DeepMosaic: Control-independent mosaic single nucleotide variant detection using deep convolutional neural networks. <i>([bioRxiv 2021](https://www.biorxiv.org/content/10.1101/2020.11.14.382473v2.full))</i>
* Somatic mosaicism in the mature brain reveals clonal cellular distributions during cortical development. <i>([bioRxiv 2020](https://www.biorxiv.org/content/10.1101/2020.08.10.244814v1.full))</i>
* Autism risk in offspring can be assessed through quantification of male sperm mosaicism. <i>([Nat. Med. 2020](https://doi.org/10.1038/s41591-019-0711-0))</i>
* mTOR pathway somatic variants and the molecular pathogenesis of hemimegalencephaly. <i>([Epi. Open 2020](https://doi.org/10.1002/epi4.12377))</i>
* Genomic mosaicism in the pathogenesis and inheritance of a Rett syndrome cohort. <i>([Genet. Med. 2019](https://doi.org/10.1038/s41436-018-0348-2))</i>
* Mosaicism and incomplete penetrance of <i>PCDH19</i> mutations. <i>([J. Med. Genet. 2019](http://dx.doi.org/10.1136/jmedgenet-2017-105235))</i>
* Somatic double-hit in <i>MTOR</i> and <i>RPS6</i> in hemimegalencephaly with intractable epilepsy. <i>([Hum. Mol. Genet. 2019](https://doi.org/10.1093/hmg/ddz194))</i>
* <i>ATP1A3</i> mosaicism in families with alternating hemiplegia of childhood. <i>([Clin. Genet. 2019](https://doi.org/10.1111/cge.13539))</i>
* Distinctive types of postzygotic single-nucleotide mosaicisms in healthy individuals revealed by genome-wide profiling of multiple organs. <i>([PLoS Genet. 2018](https://doi.org/10.1371/journal.pgen.1007395))</i>
* MosaicHunter: accurate detection of postzygotic single-nucleotide mosaicism through next-generation sequencing of unpaired, trio, and paired samples. <i>([Nucleic acids res. 2017](https://doi.org/10.1093/nar/gkx024))</i>
* Ultrasensitive and high-efficiency screen of de novo low-frequency mutations by o2n-seq. <i>([Nat. Comms. 2017](https://doi.org/10.1038/ncomms15335))</i>
* Postzygotic single‐nucleotide mosaicisms contribute to the etiology of autism spectrum disorder and autistic traits and the origin of mutations. <i>([Hum. Mutat. 2017](https://doi.org/10.1002/humu.23255))</i>
* Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort. <i>([Sci. Rep. 2017](https://doi.org/10.1038/s41598-017-15814-7))</i>
* Amplicon resequencing identified parental mosaicism for approximately 10% of <i>“de novo” SCN1A</i> mutations in children with Dravet syndrome. <i>([Hum. Mutat. 2015](https://doi.org/10.1002/humu.22819))</i>

-----------------------------------
### 4. Contact:

:email: Dr. Xiaoxu Yang: [xiy010@health.ucsd.edu](mailto:xiy010@health.ucsd.edu), [yangxiaoxu-shishen@hotmail.com](mailto:yangxiaoxu-shishen@hotmail.com)



-----------------------------------
### 5. Cite the code:
Cite the new Python version or the snakemake wrapper:

Yang X and Breuss MW, <i>et al.</i> Gleeson JG. 2021. Developmental and temporal characteristics of clonal sperm mosaicism. <i>[Cell](http://www.doi.org/10.1016/j.cell.2021.07.024)</i> 

Cite the R version: 

Yang X and Liu A, <i>et al.</i> Wei L and Zhang Y. 2017. Genomic mosaicism in paternal sperm and multiple parental tissues in a Dravet syndrome cohort. <i>[Sci. Rep.](https://doi.org/10.1038/s41598-017-15814-7)</i>

Xu X, Yang X, and Wu Q, <i>et al.</i> Wei L and Zhang Y. 2015. Amplicon resequencing identified parental mosaicism for approximately 10% of <i>“de novo” SCN1A</i> mutations in children with Dravet syndrome. <i>[Hum. Mutat.](https://doi.org/10.1002/humu.22819)</i>


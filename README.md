Code for all sequence processing and data analysis associated with the manuscript: "Plant Secondary Compound- and Antibiotic-Induced Community Disturbances Improve the Establishment of Foreign Gut Microbiota"

Contents of the files are as follows:

"Alpha diversity covariate testing.R": R script containing linear mixed effects models examining the effect of previous treatment and water intake on the alpha divserity metrics of animals in the antibiotic experiment. 

"deseq2 analysis.R": R script containing all DESeq2 analysis for both the PSC and antibiotic experiment.

"diversity analysis qiime2.txt": QIIME2 script containing all alpha and beta-diversity analysis for both the PSC and antibiotic experiment. Sequence processing which generated these final files can be found in "sequence-processing-establishment.txt". 

"indicator species analysis.R": R script for Indicator Species analysis, final indicator numbers were exported and analyzed using the "sourcetracker and indicspeices analysis.R" code. 

"sequence-processing-establishment.txt": contains QIIME2 script from raw sequences to rarefied feature table which is used as input for "diversity analysis qiime2.txt". 

"sourcetracker and indicspecies analysis.R": R script containing comparisons of microbes sourced from donors and indicator species abundance for both the PSC and antibiotic experiment. Inputs for the indicspecies analysis can be found in "indicator species analysis.R". 



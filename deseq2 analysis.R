#DESeq2 Analysis for Establishment paper 
#data from qiime2 version 2019.4 
#R version 4.1.1 (2021-08-10)

#load packages
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)
library(ggplot2)
setwd("pathway/to/data")


#PSC/creosote experiment
#read in qiime2 data as a phyloseq object
physeq <-qza_to_phyloseq(
  features="final-table-pseudo.qza",
  tree="rooted-tree.qza",
  taxonomy="taxonomy.qza",
  metadata = "novel-resource-metadata.txt")

#load phyloseq to run deseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(phyloseq)
library(DESeq2)
library(ashr)
 
#look at differentially abundant microbes between donors and recipients 
day0 <- subset_samples(physeq, Sampletype=="day0")

deseq.data.0 <- phyloseq_to_deseq2(day0, ~Population)
deseq.data.0$Population <- relevel(deseq.data.0$Population, "whiterock")
as.data.frame( colData(deseq.data.0) )
diagdds.0 = DESeq(deseq.data.0, test="Wald", fitType="parametric")
resultsNames(diagdds.0)
res.corrected.0 <- lfcShrink(diagdds.0, coef="Population_lytleranch_vs_whiterock", type="ashr")
alpha = 0.01
sigtab.0 = res.corrected.0[which(res.corrected.0$padj < alpha), ]
sigtab.0 = cbind(as(sigtab.0, "data.frame"), as(tax_table(physeq)[rownames(sigtab.0), ], "matrix"))
head(sigtab.0)
dim(sigtab.0)

#order results by smallest p-value 
res.0.ordered <- res.corrected.0[order(res.corrected.0$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.0 <- subset(res.corrected.0, padj < 0.01)
resSig.0


#+resin vs -resin animals after fecal transplant
#subset data
day5 <- subset_samples(physeq, Sampletype=="day5")
day5.recipients <- subset_samples(physeq, Population=="whiterock")

deseq.data.5.recipients <- phyloseq_to_deseq2(day5.recipients, ~Treatment)
diagdds.5.r = DESeq(deseq.data.5.recipients, test="Wald", fitType="parametric")
res.5.r = results(diagdds.5.r, cooksCutoff = FALSE)
#treatment vs control 
resultsNames(diagdds.5.r)
res.corrected.5.r <- lfcShrink(diagdds.5.r, coef="Treatment_treatment_vs_control", type="ashr")
alpha = 0.01
sigtab.5.r = res.corrected.5.r[which(res.corrected.5.r$padj < alpha), ]
sigtab.5.r = cbind(as(sigtab.5.r, "data.frame"), as(tax_table(physeq)[rownames(sigtab.5.r), ], "matrix"))
head(sigtab.5.r)
dim(sigtab.5.r)

#order results by smallest p-value 
res.5.ordered <- res.corrected.5.r[order(res.corrected.5.r$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.5 <- subset(res.5.ordered, padj < 0.01)
resSig.5


#+resin animals vs +resin animals before and after fecal transplants
#subset for treatment group
treatment <- subset_samples(physeq, Treatment=="treatment")
treatment2 <- subset_samples(treatment, Thirty=="nay")
deseq.data.treatment <- phyloseq_to_deseq2(treatment2, ~Sampletype)
diagdds.treatment = DESeq(deseq.data.treatment, test="Wald", fitType="parametric")
res.treatment = results(diagdds.treatment, cooksCutoff = FALSE)
resultsNames(diagdds.treatment)
res.corrected.treatment <- lfcShrink(diagdds.treatment, coef="Sampletype_day5_vs_day0", type="ashr")
alpha = 0.01
sigtab.treatment = res.corrected.treatment[which(res.corrected.treatment$padj < alpha), ]
sigtab.treatment = cbind(as(sigtab.treatment, "data.frame"), as(tax_table(physeq)[rownames(sigtab.treatment), ], "matrix"))
head(sigtab.treatment)
dim(sigtab.treatment)

res.treatment.ordered <- res.corrected.treatment[order(res.corrected.treatment$pvalue),]

#order results by smallest p-value 
res.treatment.ordered <- res.corrected.treatment[order(res.corrected.treatment$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.treatment <- subset(res.treatment.ordered, padj < 0.01)
resSig.treatment

#-resin animals compared to -resin animals before and after transplant
#subset data
control <- subset_samples(physeq, Treatment=="control")
control2 <- subset_samples(control, Thirty=="nay")
deseq.data.control <- phyloseq_to_deseq2(control2, ~Sampletype)
deseq.data.control$Sampletype <- relevel(deseq.data.control$Sampletype, "day0")
diagdds.control = DESeq(deseq.data.control, test="Wald", fitType="parametric")
res.control = results(diagdds.control, cooksCutoff = FALSE)
resultsNames(diagdds.control)
res.corrected.control <- lfcShrink(diagdds.control, coef="Sampletype_day5_vs_day0", type="ashr")
alpha = 0.01
sigtab.control = res.corrected.control[which(res.corrected.control$padj < alpha), ]
sigtab.control = cbind(as(sigtab.control, "data.frame"), as(tax_table(physeq)[rownames(sigtab.control), ], "matrix"))
head(sigtab.control)
dim(sigtab.control)

#order results by smallest p-value 
res.control.ordered <- res.corrected.control[order(res.corrected.control$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.control <- subset(res.control.ordered, padj < 0.01)
resSig.control

#antibiotic experiment
library(qiime2R)
setwd("pathway/to/data")
#read in qiime2 data as a phyloseq object
physeq.anti <-qza_to_phyloseq(
  features="single-final-table.qza",
  tree="rooted-tree-single.qza",
  taxonomy="taxonomy-single.qza",
  metadata = "antibiotic-metadata.tsv")
#load phyloseq to run deseq2
#look at differentially abundant microbes between donors and recipients 
day0.abx <- subset_samples(physeq.anti, deseq.0=="yes")

deseq.data.0.abx <- phyloseq_to_deseq2(day0.abx, ~Population)
diagdds.0.abx = DESeq(deseq.data.0.abx, test="Wald", fitType="parametric")
res.0.abx = results(diagdds.0.abx, cooksCutoff = FALSE)
resultsNames(diagdds.0.abx)
res.corrected.0.abx <- lfcShrink(diagdds.0.abx, coef="Population_whiterock_vs_lytleranch", type="ashr")
alpha = 0.01
sigtab.0.abx = res.0.abx[which(res.0.abx$padj < alpha), ]
sigtab.0.abx = cbind(as(sigtab.0.abx, "data.frame"), as(tax_table(physeq.anti)[rownames(sigtab.0.abx), ], "matrix"))
head(sigtab.0.abx)
dim(sigtab.0.abx)

#order results by smallest p-value 
res.0.ordered.abx <- res.0.abx[order(res.0.abx$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.0.abx <- subset(res.0.ordered.abx, padj < 0.01)
resSig.0.abx


#+antibiotic vs -antibiotic after transplant
#subset data 
day5.abx <- subset_samples(physeq.anti, Sampletype=="day5")
day5.abx2 <- subset_samples(day5.abx, Population=="whiterock")

deseq.data.5.abx <- phyloseq_to_deseq2(day5.abx2, ~Treatment)
diagdds.5.abx = DESeq(deseq.data.5.abx, test="Wald", fitType="parametric")
res.5.abx = results(diagdds.5.abx, cooksCutoff = FALSE)
resultsNames(diagdds.5.abx)
res.corrected.5.abx <- lfcShrink(diagdds.5.abx, coef="Treatment_treatment_vs_control", type="ashr")
alpha = 0.01
sigtab.5.abx = res.corrected.5.abx[which(res.corrected.5.abx$padj < alpha), ]
sigtab.5.abx = cbind(as(sigtab.5.abx, "data.frame"), as(tax_table(physeq.anti)[rownames(sigtab.5.abx), ], "matrix"))
head(sigtab.5.abx)
dim(sigtab.5.abx)

#order results by smallest p-value
res.5.ordered.abx <- res.corrected.5.abx[order(res.corrected.5.abx$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.5.abx <- subset(res.5.ordered.abx, padj < 0.01)
resSig.5.abx

#+antibiotic animals vs +antibiotic animals before and after transplant
#subset data
treatment.abx <- subset_samples(physeq.anti, Treatment=="treatment")
treatment.abxpart2 <- subset_samples(treatment.abx , abx=="no")
treatment2.abx <- subset_samples(treatment.abxpart2, Thirty=="no")
deseq.data.treatment.abx <- phyloseq_to_deseq2(treatment2.abx, ~Sampletype)
diagdds.treatment.abx = DESeq(deseq.data.treatment.abx, test="Wald", fitType="parametric")
res.treatment.abx = results(diagdds.treatment.abx, cooksCutoff = FALSE)
resultsNames(diagdds.treatment.abx)
res.corrected.treatment.abx <- lfcShrink(diagdds.treatment.abx, coef="Sampletype_day5_vs_day0", type="ashr")
alpha = 0.01
sigtab.treatment.abx = res.corrected.treatment.abx[which(res.corrected.treatment.abx$padj < alpha), ]
sigtab.treatment.abx = cbind(as(sigtab.treatment.abx, "data.frame"), as(tax_table(physeq.anti)[rownames(sigtab.treatment.abx), ], "matrix"))
head(sigtab.treatment.abx)
dim(sigtab.treatment.abx)

#order results by smallest p-value 
res.treatment.ordered.abx <- res.corrected.treatment.abx[order(res.corrected.treatment.abx$pvalue),]

#subset data for only p-values with adjusted p that is significant 
resSig.treatment.abx <- subset(res.treatment.ordered.abx, padj < 0.01)
resSig.treatment.abx


#-antibiotic animals compared to - antibiotic animals before and after translpant
#subset data
control <- subset_samples(physeq.anti, Treatment=="control")
control2 <- subset_samples(control, Thirty=="nay")
deseq.data.control.abx <- phyloseq_to_deseq2(control2, ~Sampletype)
deseq.data.control.abx$Sampletype <- relevel(deseq.data.control.abx$Sampletype, "day0")
diagdds.control.abx = DESeq(deseq.data.control.abx, test="Wald", fitType="parametric")
res.control.abx = results(diagdds.control.abx, cooksCutoff = FALSE)
resultsNames(diagdds.control.abx)
res.corrected.control.abx <- lfcShrink(diagdds.control.abx, coef="Sampletype_day5_vs_day0", type="ashr")
alpha = 0.01
sigtab.control.abx = res.corrected.control.abx[which(res.corrected.control.abx$padj < alpha), ]
sigtab.control.abx = cbind(as(sigtab.control.abx, "data.frame"), as(tax_table(physeq.anti)[rownames(sigtab.control.abx), ], "matrix"))
head(sigtab.control.abx)
dim(sigtab.control.abx)

#order results by smallest p-value 
res.control.ordered.abx <- res.corrected.control.abx[order(res.corrected.control.abx$pvalue),]

#then subset data for only p-values with adjusted p that is significant 
resSig.control.abx <- subset(res.control.ordered.abx, padj < 0.01)
resSig.control.abx


#+antibiotic animals vs +antibiotic animals immediately after administration of abx
abx.compare <- subset_samples(physeq.anti, abx.compare=="yes")

deseq.data.abx.compare <- phyloseq_to_deseq2(abx.compare, ~Sampletype)
diagdds.abx.compare = DESeq(deseq.data.abx.compare, test="Wald", fitType="parametric")
res.abx.compare = results(diagdds.abx.compare, cooksCutoff = FALSE)
resultsNames(diagdds.abx.compare)
res.corrected.abx.compare <- lfcShrink(diagdds.abx.compare, coef="Sampletype_day0_vs_antibiotic", type="ashr")
alpha = 0.01
sigtab.abx.compare = res.corrected.abx.compare[which(res.corrected.abx.compare$padj < alpha), ]
sigtab.abx.compare = cbind(as(sigtab.abx.compare, "data.frame"), as(tax_table(physeq.anti)[rownames(sigtab.abx.compare), ], "matrix"))
head(sigtab.abx.compare)
dim(sigtab.abx.compare)

#order results by smallest p-value to make life easier
res.0.ordered.abx.compare <- res.corrected.abx.compare[order(res.corrected.abx.compare$pvalue),]

#then subset data for only p-values with adjusted p that is significant 
resSig.abx.compare<- subset(res.abx.compare, padj < 0.01)
resSig.abx.compare
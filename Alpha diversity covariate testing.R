#Running  Linear Mixed Effects Models on the alpha diversity data to determine the effect of total water intake on diversity measurements, as well as the effect of previous treatment 
#We did not see a significant effect of water intake 
#all LME models showed that previous treatment did not account for a large portion of the variance in our data, nor did any of our significance values change.
library(rstatix)
library(lmerTest)
library(lme4)
setwd("pathway-to-file-containing-all-alpha-metrics")
data <- read.csv("table-containing-all-alpha-metrics")

#subsetting for the different time points
day.0 <- data[which(data$Day=='0'),]
day.5 <- data[which(data$Day=='5'),]
day.30 <- data[which(data$Day=='30'),]
abx <- data[which(data$Day=='3'),]

#subsetting for different treatment groups
treatment <- data[which(data$Treatment=='treatment'),]
#retaining +antibiotic animals day 0 and day 3, post abx 
treatment.abx <- data[which(data$abx.compare.1=='yes'),]

#LMM for each alpha-diversity metric day 0, pre transplant
shannon.lm <- lmer(Shannon ~ Treatment + (1|PreviousTreatment), data=day.0)
summary(shannon.lm)
#singular fitting model, not significant, previous treatment explains 1.851e-15 of the variance

#observed otus
observed.lm <- lmer(ObservedOtus ~ Treatment + (1|PreviousTreatment), data=day.0)
summary(observed.lm)
#singular fitting model, not significant, previous treatment explains 0 of the variance. 0 is not actually possible, but the number is low enough to say previous treatment does not explain significant variance.

#faith's pd
faith.lm <- lmer(Faiths ~ Treatment + (1|PreviousTreatment), data=day.0)
summary(faith.lm)
#singular fitting model, not significant, previous treatment explains 0 of the variance. 

#day 5, post transplant
shannon.lm.5 <- lmer(Shannon ~ Treatment + (1|PreviousTreatment), data=day.5)
summary(shannon.lm.5)
#significance matches wilcoxon test, previous treatment variance is >0.01

#observed otus
observed.lm.5 <- lmer(ObservedOtus ~ Treatment + (1|PreviousTreatment), data=day.5)
summary(observed.lm.5)
#singular fitting model, previous treatment accounts for 2.014e-10 of the variance. Significance is not different from previous wilcoxon test. 

#faith's pd
faith.lm.5 <- lmer(Faiths ~ Treatment + (1|PreviousTreatment), data=day.5)
summary(faith.lm.5)
#not significant, previous treatment explains little variance, significance matches wilcoxon tests.  

#30 days post transplant
shannon.lm.30 <- lmer(Shannon ~ Treatment + (1|PreviousTreatment), data=day.30)
summary(shannon.lm.30)
#singular fit, significance matches wilcoxon test, previous treatment variance is 0. Can't actually be 0 but we can assume it is not contributing significant variance. 

#observed otus
observed.lm.30 <- lmer(ObservedOtus ~ Treatment + (1|PreviousTreatment), data=day.30)
summary(observed.lm.30)
#singular fitting model, previous treatment accounts for 0 of the variance. Significance is not different from wilcoxon test. 

#faith's pd
faith.lm.30 <- lmer(Faiths ~ Treatment + (1|PreviousTreatment), data=day.30)
summary(faith.lm.30)
#singular fitting model, previous treatment explains 0 of the variance, significane matches reported tests. 

#did water intake affect the shannon/otus/faiths? 
#right after receiving neomycin and after transplants

#Faith's PD post abx
cor.test(abx$water_intake, abx$Faiths, method = c("spearman"))
#not significant

#Faith's PD day 5, post transplant
cor.test(day.5$water_intake, day.5$Faiths, method = c("spearman"))
#not significant

#Shannon diversity index post abx
cor.test(abx$water_intake, abx$Shannon, method = c("spearman"))
#not significant

#Shannon diversity index day 5, post transplant
cor.test(day.5$water_intake, day.5$Shannon, method = c("spearman"))

#Observed ASVs post abx
cor.test(abx$water_intake, abx$ObservedOtus, method = c("spearman"))
#not significant

#Observed ASVs day 5, post transplant
cor.test(day.5$water_intake, day.5$ObservedOtus, method = c("spearman"))

#wilcoxon tests comparing the +antibiotic animals to themselves before and after administration of antibiotics to see if antibiotic effected alpha-diverstiy metrics
shannon.cox <- wilcox.test(Shannon ~ as.factor(Day), data = treatment.abx, paired = TRUE)
shannon.cox

observed.cox <- wilcox.test(ObservedOtus ~ as.factor(Day), data = treatment.abx, paired = TRUE)
observed.cox

faiths.cox <- wilcox.test(Faiths ~ as.factor(Day), data = treatment.abx, paired = TRUE)
faiths.cox

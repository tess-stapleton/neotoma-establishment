#Sourcetracker and indicspecies analysis
#R version 4.1.1 (2021-08-10)

setwd("pathway/to/sourcetracker/data")
source.data <- read.csv("sourcetracker.data")

library(ggplot2)


#PSC experiment
#-resin vs +resin % microbes sourced from donors after transplant
#mann whitney u-test 
wilcox.test(P.Donor~Treatment, data=source.data)
# w = 4, p-value = 0.00373

#-resin vs +resin % microbes sourced from donors 30 days after transplant
thirty <- read.csv("sourcetracker-data-here")
thirty$ID <- factor(thirty$ID)

#mann-whitney u-test
wilcox.test(sourced~treatment, data=thirty)
#W = 10, p-value = 0.04009

#antibiotic experiment
setwd("pathway/to/data")
anti.data <- read.csv("antibiotic.sourcetracker.data")

#was there significant difference in the % donor microbes between groups?
#wilcoxin' 
wilcox.test(P.Donor~Treatment, data=anti.data)
#w=11, p=0.09732

#adding previoustreatment as a covariate
lme.source <- lmer(P.Donor ~ Treatment + (1|PreviousTreatment), data=anti.data)
summary(lme.source)
#previous treatment accounts for extremely low variance, 9.849e-06, should not be included as covariate

means.plot.anti <- ggplot(anti.data, aes(x=Treatment, y=Abs.P.Change)) + 
  geom_boxplot()

#do the groups differ in % microbes sourced from donors 30 days after transplant?
thirty.abx <- read.csv("source-data-here")
wilcox.test(sourced~treatment, data=thirty.abx)
#W = 25, p-value = 1
lme.source.30 <- lmer(sourced ~ treatment + (1|PreviousTreatment), data=thirty.abx)
summary(lme.source.30)
#singular fit, previous treatment 0 variance, not significant. 


#Change in Indicator species from Indicspecies, multipatt analysis
setwd("pathway/to/file/here")
relative <- read.csv("abundance.of.indicator.species.data")
#subset for the different time points
relative.0 <- relative[which(relative$day=='0'), ]
relative.5 <- relative[which(relative$day=='5'), ]
relative.recipients.0 <- relative.0[which(relative.0$group=='recipients'), ]
relative.recipients.5 <- relative.5[which(relative.5$group=='recipients'), ]
library(ggplot2)
abundance.plot <- ggplot(relative.0, aes(x=as.factor(day), y=abundance, fill=treatment)) +
  geom_boxplot()

#testing for differences between donors and recipients before transplant
wilcox.test(abundance~group, data=relative.0)
#mann whitney u-test instead
wilcox.test(abundance~treatment, data=relative.recipients.0)
#W = 48, p-value = 0.02051
kruskal.test(abundance~treatment, data=relative.0)
library(FSA)
post.hoc = dunnTest(abundance ~ treatment,
              data=relative.0,
              method="bh")

res.aov <- aov(abundance ~ treatment, data = relative.0)
summary(res.aov)
TukeyHSD(res.aov)

pairwise.wilcox.test(relative.0$abundance, relative.0$treatment,
                     p.adjust.method = "BH")

wilcox.test(abundance~group, data=relative.0)
#donors vs recipients
#W = 75, p-value = 0.000129

#wilcox test abundance of indicator features b/w +resin/-resin animals after transplant
wilcox.test(abundance~treatment, data=relative.recipients.5)
#W = 9, p-value = 0.0289

#paired test to see if treatment/control significantly increased it's abundance of indicator features
#subset for -resin 
control.indic.abundance <- relative[ which(relative$treatment=='control'), ]
wilcox.test(abundance ~ day, data = control.indic.abundance, paired = TRUE)
#W = 4, p-value = 0.1094

#subset for +resin animals
treatment.indic.abundance <- relative[ which(relative$treatment=='treatment'), ]
#test for significant change over time
wilcox.test(abundance ~ day, data = treatment.indic.abundance, paired = TRUE)
#V = 0, p-value = 0.007813

#is there a difference b/w +resin and -resin animals 30 days after transplant?
setwd("pathway/to/data")
indic.30 <- read.csv("indicator species 30 data")
#subset for day 30, recipients
indic.30.recipients <- indic.30[which(indic.30$Day=='30'), ]

#do +resin animals have higher abundance than -resin?
wilcox.test(p.indicator ~ Treatment, data = indic.30.recipients, paired = FALSE)
#W = 19, p-value = 0.3357

#antibiotic experiment
relative.anti <- read.csv("antibiotic indicator species data")
#subset for before and after transplant
relative.0.anti <- relative.anti[which(relative.anti$day=='0'), ]
relative.5.anti <- relative.anti[which(relative.anti$day=='5'), ]
relative.recipients.0.anti <- relative.0.anti[which(relative.0.anti$group=='recipients'), ]
library(ggplot2)
abundance.plot2 <- ggplot(relative.anti, aes(x=as.factor(day), y=abundance, fill=treatment)) +
  geom_boxplot()

#do donors harbor more indicators than recipients prior to fecal transplant
wilcox.test(abundance~group, data=relative.0.anti)
#W = 68, p-value = 0.000688

#do +antibiotic animals harbor more indicators than -antibiotic prior to transplant
wilcox.test(abundance~treatment, data=relative.recipients.0.anti)
#W = 32, p-value = 0.3829

#previous treatment as covariate
lme.indic <- lmer(abundance ~ treatment + (1|previoustreatment), data=relative.recipients.0.anti)
summary(lme.indic)
#not significant, previous treatment variance 8.314e-07

#no sig difference in relative abundance of donor indicators in treatment/control but significant in donor vs the two recipient groups

#did -antibiotic animals show increase in indicators?
#subset for -antibiotic
control.indic.abundance.anti <- relative.anti[ which(relative.anti$treatment=='control'), ]
control.indic.abundance.anti <- relative.anti[ which(control.indic.abundance.anti$day=='control'), ]

#mann whitney u
wilcox.test(abundance ~ day, data = control.indic.abundance.anti, paired = TRUE)
#V=14, p-value = 1

#effect of previous treatment?
lme.indic.control <- lmer(abundance ~ as.factor(day) + (1|previoustreatment) + (1|SampleID), data=control.indic.abundance.anti)
summary(lme.indic.control)
#singular fit, variance from previous treatment very low, 6.081e-15,no significant difference

#now treatment
treatment.indic.abundance.anti <- relative.anti[ which(relative.anti$treatment=='treatment'), ]

lme.indic.treatment <- lmer(abundance ~ as.factor(day) + (1|previoustreatment) + (1|SampleID), data=treatment.indic.abundance.anti)
summary(lme.indic.treatment)
#singular fit, variance from previous treatment low, 2.72e-05, no significant difference b/w any dates

#trying a man whitney
wilcox.test(abundance ~ day, data = treatment.indic.abundance.anti, paired = TRUE)
#V = 6, p-value = 0.2188

#30 days 
indic.abx.30 <- read.csv("30 day indicator species data")
indic.30.recipients.abx <- indic.abx.30[which(indic.abx.30$day=='30'), ]

#difference in abundance of indicator species for +antibiotic vs -antibiotic?
wilcox.test(p.indicators ~ treatment, data = indic.30.recipients.abx, paired = FALSE)
#W = 24, p-value = 1
#no, no difference b/w control and treatment group after 30 days 

#effect of previous treatment?
lme.indic.compare <- lmer(p.indicators ~ treatment + (1|previoustreatment), data=indic.30.recipients.abx)
summary(lme.indic.compare)
#singular fit, not significant, variance from previous treatment 0 

#Indic species analysis 
#data from qiime2 v 2019.4 
#data are ASVs clustered at the 100% sequence identity level

library(indicspecies)
setwd("pathway/to/data")
indic.data <- read.csv("metadata")
indic.wild <- indic.data[ which(indic.data$Day=='0'), ]

abund = indic.wild[,6:ncol(indic.wild)]
abund.5 = indic.5[,6:ncol(indic.5)]

group = indic.wild$group
group.5 = indic.5$group

inv = multipatt(abund, group, func = "IndVal.g", control = how(nperm=9999))

summary(inv, indvalcomp=TRUE)
inv$sign
dat.multipatt.summary<-capture.output(summary(inv, indvalcomp=TRUE))
write.csv(dat.multipatt.summary, "indic.species.csv")


#correctto p-value for multiple comparisons
#first make into data table
library(data.table)
novel.indic <- read.csv("indic.species.formatted.csv")
novel.table<-as.data.table(novel.indic, keep.rownames=TRUE)
#add adjusted p-value
novel.table[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
novel.table[p.value.bh<=0.05, ]
write.csv(novel.table, "indic.species.corrected.csv")

#format in excel
#you can extract parts of a cell using 
#=LEFT(text,LEN(text)-n) where text is the name of the cell you want to extract characters from and n is the number of characters to remove. You can also you =RIGHT to go from the other direction. It takes some work, but you can do it. 
#now run for white rock animals and then remove features that exist in both 
sc2 = indicators(X=abund, cluster=group, group="WhiteRock",
                 max.order = TRUE, verbose=TRUE,
                 At=0.8, Bt=0.35)
print(sc2, sqrtIVt = 0.6)
white.rock.indic <-capture.output(print(sc2, indvalcomp=TRUE))
write.csv(white.rock.indic, "white.rock.indic.csv")

options(max.print=5.5E5)


#antibiotic experiment
setwd("pathway/to/data")
indic.data.a <- read.csv("metadata")
indic.wild.a <- indic.data.a[ which(indic.data.a$day=='0'), ]

abund.a = indic.wild.a[,6:ncol(indic.wild.a)]
group.a = indic.wild.a$group
inv.a = multipatt(abund.a, group.a, func = "IndVal.g", control = how(nperm=9999))
summary(inv.a, indvalcomp=TRUE)
inv.a$sign
dat.multipatt.summary.a<-capture.output(summary(inv.a, indvalcomp=TRUE))
write.csv(dat.multipatt.summary.a, "antibiotic.indic.species.csv")

#correct p-value for multiple comparisons
#first make into data table
library(data.table)
antibiotic.indic <- read.csv("antibiotic.indic.species.formatted.csv")
antibiotic.table<-as.data.table(antibiotic.indic, keep.rownames=TRUE)
#add adjusted p-value
antibiotic.table[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
antibiotic.table[p.value.bh<=0.05, ]
write.csv(antibiotic.table, "antibiotict.indic.species.corrected.csv")
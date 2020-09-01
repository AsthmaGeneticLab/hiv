#################################################################
########### Differential Expression Analysis (DGE) for HIV ######
#################################################################

## This file contains computational R code for DGE for HIV with confoudning variables
## with all confounding variables

## in this version mode of delivery, Celect & Cemerg have been merged and delivery_merged_is_natural is generated Yes/No field. Yes=Natural and No= Celect or Cemerg   

# Change the path to the directory
# On Windows use a forward slash / instead of the usual /.
setwd("D:/R-scripts/hiv_chord_gene")

#####DIFFERENTIAL EXPRESSION ANALYSIS WITH LIMMA ######
library(limma)
# 'ExprsMat_GitHub.txt' has been taken from https://github.com/BreenMS/Umbilical-cord-blood-transcriptome-analysis
exprs <- read.delim("Data/DrakensteinUCBExprs_updated/ExprsMat_GitHub.txt", check.names=FALSE, stringsAsFactors=FALSE, header=T, row.names=1)

targets <- read.delim("Data/DrakensteinUCBExprs_updated/MetadataHiv_Phenotype_Only_Ordered_Delvr_Merged.csv", sep = ',', check.names=FALSE, stringsAsFactors=FALSE)
Group <- factor(targets$child_hiv_birth, levels=c("HIV unexposed","HIV exposed uninfected")) # On which we want to measure, that column should be in the right side, e.g. Yes = HIVexposed 

#Add covariates
Nicotine <- factor(targets$smokecat)
Alcohol <- factor(targets$alcohol)
RIN <- (targets$rin)
Batch <- factor(targets$batch)
Gender <- factor(targets$gender)
GestationalAge <- (targets$gestation_delivery)
Delivery <- factor(targets$delivery_merged_is_natural)
Ethnicity <- factor(targets$ethnicity)
PTSD <-factor(targets$ptsd)

design <- model.matrix(~Group+Nicotine+Alcohol+RIN+Batch+Gender+GestationalAge+Delivery+Ethnicity+PTSD)
dim(design)
colnames(design)
#design

fit <- lmFit(exprs, design)
fit2 <- eBayes(fit)
topTable(fit2, adjust="BH")

Diff<- topTable(fit2, coef=2, n=25000) # Here coef=2 means Yes = HIVexposed 
write.table(Diff, file="DGE/v1.1-on-updated-data/DGE_Hiv_result4.2.txt", quote=FALSE, row.names=TRUE)
write.csv(Diff, file = 'DGE/v1.1-on-updated-data/DGE_Hiv_result4.2.csv', quote=FALSE, row.names=TRUE)



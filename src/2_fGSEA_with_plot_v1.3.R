##############################
#
#
# Pathway analysis using fGSEA on the outcome of DGE with all covariates. The result of DGEHiv4.2.R, DGE_Hiv_result4.2.csv, is used in this script  
#
#
##############################
#setwd("\\\\filestore.soton.ac.uk/users/ah1e13/mydocuments/Documents/Tilman/Anastasio_Mac_target_Exo/mRNA")

# Change the path to the directory
# On Windows use a forward slash / instead of the usual /.
setwd("D:/R-scripts/hiv_chord_gene")

#data<-read.table("./Results/DEG_res_All_Block_Design_Healthy_Control.csv", sep=",", header=TRUE, row.names=1)
data<-read.table("Data/DGE_output/DGE_Hiv_result4.2.csv", sep=",", header=TRUE)

#Using the log fc and P value to create a ranking metric for the GSEA.
data$fcsign <- sign(data$logFC)
data$logP=-log10(data$P.Value)
data$metric= data$logP/data$fcsign
head(data)
final<-data[,c("genes", "metric")]
rank<-final$metric
names(rank)<-final$genes


head(final)
head(rank)
#write.csv(rank, file = 'DGE/v1.1-on-updated-data/GeneRanks_bio.csv')

#Reading in the files needed for the GSEA downloaded from the broad institute http://software.broadinstitute.org/gsea/login.jsp;jsessionid=FB20DDF33CCD1FE98BBC0DFA9048B9A7.
library(gage)
# Cellular Component
cellular_comp<-file("Data/GSEA_gmt_files/c5.cc.v7.0.symbols.gmt")
cell_list<-readList(cellular_comp)
#Biological Proscess
biological_proccess<-file("Data/GSEA_gmt_files/c5.bp.v7.0.symbols.gmt")
bio_list<-readList(biological_proccess)
# Molecular Function
molecular_function<-file("Data/GSEA_gmt_files/c5.mf.v7.0.symbols.gmt")
mol_list<-readList(molecular_function)
# c2.cp.kegg.v6.2.symbols.gmt Kegg pathways
Kegg<-file("Data/GSEA_gmt_files/c2.cp.kegg.v7.0.symbols.gmt")
kegg_list<-readList(Kegg)


#Actually Conducting the GSEA
library(fgsea)
library(data.table) #used for writing tabular data in file 

#### for molecular process  
fgsea_result1<-fgsea(pathways=mol_list, stats=rank, nperm=10000, minSize = 1, maxSize = Inf)
## Conducting GSEA on pre ranked list.
fwrite(fgsea_result1, file="fGSEA_mol_Hiv_vs_NonHiv4.2.csv", sep=",", sep2=c("", " ", ""))
## Number of Significant hits after FDR
sum(fgsea_result1[,padj < 0.05])

#### for biological proccess
## for reproducing output, we need to set seed. Ref: https://github.com/ctlab/fgsea/issues/12
## Here for biological_proccess (bio_list)
##		if the seed is 1, then output is 202, 
##		if the seed is 42, then output is 223, 
##		if the seed is 500, then output is 193, 
##		if the seed is 1000, then output is 185, 
##		if the seed is 5000, the output is 243, ---- This one is highest for DGE_Hiv_result4.1.csv
##		if the seed is 10000, the output is 212,
##		if the seed is 7500, the output is 194,
##		if the seed is 6000, the output is 242,
##		if the seed is 4000, the output is 211,
set.seed(5000)
## Conducting GSEA on pre ranked list.
fgsea_result2<-fgsea(pathways=bio_list, stats=rank, nperm=10000, minSize = 1, maxSize = Inf)
fwrite(fgsea_result2, file="fGSEA_bio_Hiv_vs_NonHiv4.2.csv", sep=",", sep2=c("", " ", ""))
## Number of Significant hits after FDR
sum(fgsea_result2[,padj < 0.05])
 
## Here for kegg (kegg_list)
##		if the seed is 1, then output is 21, 
##		if the seed is 42, then output is 23,
##		if the seed is 100, then output is 23,
##		if the seed is 500, then output is 22,  
##		if the seed is 900, then output is 23, 
##		if the seed is 1000, then output is 24, ---- This one is highest for DGE_Hiv_result4.1.csv
##		if the seed is 1100, then output is 23, 
##		if the seed is 5000, the output is 22,
##		if the seed is 10000, the output is 20,
##		if the seed is 7500, the output is 23,
##		if the seed is 6000, the output is 23,
##		if the seed is 4000, the output is 21
set.seed(1000) 
## Conducting GSEA on pre ranked list.
fgsea_result3<-fgsea(pathways=kegg_list, stats=rank, nperm=10000, minSize = 1, maxSize = Inf)
fwrite(fgsea_result3, file="fGSEA_kegg_Hiv_vs_NonHiv4.2.csv", sep=",", sep2=c("", " ", ""))
## Number of Significant hits after FDR
sum(fgsea_result3[,padj < 0.05])


##If desired you can make an enrichment plot for a pathway
#pdf("GO_BP_AH_Filtered_GO_IMMUNE_SYSTEM_PROCESS.pdf")
#plotEnrichment(bio_filt_list[["GO_IMMUNE_SYSTEM_PROCESS"]],f_scores)
#dev.off()

#fgsea_result2 <- fgsea_result
#fgsea_result2 <- subset(fgsea_result2, padj < 0.05)
#fgsea_result2
#rCnt <- nrow(fgsea_result2) 
#cat("Total row count=", rCnt)
#for(pn in fgsea_result2$pathway) {
#	print(pn)
#}

   
# Can make a table plot for a bunch of selected pathways
#topPathwaysUp<-fgsea_result[ES > 0][head(order(pval),n=10),pathway]
#topPathwaysDown<-fgsea_result[ES < 0][head(order(pval),n=10),pathway]
#topPathways<-c(topPathwaysUp, rev(topPathwaysDown))
#pdf("GO_BP_Tumour_vs_Healthy.pdf",w=19)
#plotGseaTable(bio_list[topPathways],rank,fgsea_result, gseaParam=0.5)
#dev.off()


## justting getting top 10 and bottom 10 based on ES by ordering p-value  
#topPathwaysUp2<-fgsea_result2[ES > 0][head(order(pval),n=10),pathway]
#topPathwaysDown2<-fgsea_result2[ES < 0][head(order(pval),n=10),pathway]
#topPathways2<-c(topPathwaysUp2, rev(topPathwaysDown2))
#pdf("GO_Bio_Hiv_vs_NonHiv.pdf",w=19)
#plotGseaTable(bio_list[topPathways2],rank,fgsea_result2, gseaParam=0.5)
#dev.off()

## justting getting top 10 and bottom 10 based on ES by ordering p-adj,p-value  
padjOrderedPathways<-fgsea_result2[padj < 0.05][order(padj),]
pvalOrderedPathways<-padjOrderedPathways[order(pval),]
topPathwaysUp2<-pvalOrderedPathways[ES > 0][head(order(pval),n=10),pathway]
topPathwaysDown2<-pvalOrderedPathways[ES < 0][head(order(pval),n=10),pathway]
topPathways2<-c(topPathwaysUp2, rev(topPathwaysDown2))
pdf("GO_Bio_Hiv_vs_NonHiv_4.2.pdf",w=19)
plotGseaTable(bio_list[topPathways2],rank,fgsea_result2, gseaParam=0.5)
dev.off()

## justting getting top 10 and bottom 10 based on ES by ordering p-value  
#topPathwaysUp3<-fgsea_result3[ES > 0][head(order(pval),n=10),pathway]
#topPathwaysDown3<-fgsea_result3[ES < 0][head(order(pval),n=10),pathway]
#topPathways3<-c(topPathwaysUp3, rev(topPathwaysDown3))
#pdf("GO_KEGG_Hiv_vs_NonHiv.pdf",w=19)
#plotGseaTable(kegg_list[topPathways3],rank,fgsea_result3, gseaParam=0.5)
#dev.off()

## justting getting top 10 and bottom 10 based on ES by ordering p-adj,p-value  
padjOrderedPathways<-fgsea_result3[padj < 0.05][order(padj),]
pvalOrderedPathways<-padjOrderedPathways[order(pval),]
topPathwaysUp3<-pvalOrderedPathways[ES > 0][head(order(pval),n=10),pathway]
topPathwaysDown3<-pvalOrderedPathways[ES < 0][head(order(pval),n=10),pathway]
topPathways3<-c(topPathwaysUp3, rev(topPathwaysDown3))
pdf("GO_KEGG_Hiv_vs_NonHiv_4.2.pdf",w=19)
plotGseaTable(kegg_list[topPathways3],rank,fgsea_result3, gseaParam=0.5)
dev.off()


ndf <- data[order(data$t),] # ordering by t smallest to bigger
head(ndf)
fdf <- ndf[,c("genes","t")]
head(fdf)
geneRanks <- fdf$t
names(geneRanks) <- fdf$genes

library(ggplot2)
pnm <- "GO_COTRANSLATIONAL_PROTEIN_TARGETING_TO_MEMBRANE" 
pname <- trimws(pnm)
fname <- paste0(pname, ".pdf")
pdf(fname)
plotEnrichment(bio_list[[pname]],
              geneRanks) + labs(title=pname)
dev.off()	 
pnm <- "KEGG_RIBOSOME" 
pname <- trimws(pnm)
fname <- paste0(pname, ".pdf")
pdf(fname)
plotEnrichment(kegg_list[[pname]],
              geneRanks) + labs(title=pname)
dev.off()	 

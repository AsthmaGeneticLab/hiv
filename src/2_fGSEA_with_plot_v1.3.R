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
set.seed(5000)
## Conducting GSEA on pre ranked list.
fgsea_result2<-fgsea(pathways=bio_list, stats=rank, nperm=10000, minSize = 1, maxSize = Inf)
fwrite(fgsea_result2, file="fGSEA_bio_Hiv_vs_NonHiv4.2.csv", sep=",", sep2=c("", " ", ""))
## Number of Significant hits after FDR
sum(fgsea_result2[,padj < 0.05])
 
## Here for kegg (kegg_list)
set.seed(1000) 
## Conducting GSEA on pre ranked list.
fgsea_result3<-fgsea(pathways=kegg_list, stats=rank, nperm=10000, minSize = 1, maxSize = Inf)
fwrite(fgsea_result3, file="fGSEA_kegg_Hiv_vs_NonHiv4.2.csv", sep=",", sep2=c("", " ", ""))
## Number of Significant hits after FDR
sum(fgsea_result3[,padj < 0.05])
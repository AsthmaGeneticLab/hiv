#rm(list=ls())

# Change the path to the directory
# On Windows use a forward slash / instead of the usual /.
setwd("D:/R-scripts/hiv_chord_gene")

load('WGCNA/Drakenstein_WGCNA_pow12_all_on_DGE_with_all_covar_using_DGE_rslt_v4.2.RData')
options(stringsAsFactors = FALSE);

colorh1 = modulesPRE;
 
modules<-cbind(colnames(dataExpr0),colorh1)
modules<-as.data.frame(modules)
colnames(modules)<-c("gene","color")

unique_color<-unique(colorh1)
X<-nrow(modules)

HivFu = "HivFu_v1.2"
#result<-read.csv('i_frc_litres__6wks_cov.csv',row.names=1,check.names=F)
#gene sorting code are available in dge_gene_sort.R 
result <- read.delim("Data/DGE_output/DGE_Hiv_result4.2_Sorted.txt", check.names=TRUE, stringsAsFactors=FALSE, header=TRUE) #THIS FILE CONTAINS A VECTOR OF DGE P-VALUES FROM DIFFERENTIAL GENE EXPRESSION ANALYSES, GENES NEED TO BE SORTED IN THE EXACT SAME ORDER AS MATRIX INPUT
colnames(result) = names(result)

res<-result[result$P.Value<=0.01,]
Y<-nrow(res)
mat<-matrix(,nrow=length(unique_color),ncol=2)
for(i in 1:length(unique_color)){
	module_color<-modules[modules$color%in%unique_color[i],]
	m<-nrow(module_color)
	res_color<-module_color[module_color$gene%in%res$genes,]
	n<-nrow(res_color)
	
	#Hypergeometric test (phyper) 
	#URL: https://seqqc.wordpress.com/2019/07/25/how-to-use-phyper-in-r/
	#Test for under-representation (depletion)
	#phyper(Overlap, group2, Total-group2, group1, lower.tail= TRUE) 
	#Test for over-representation (enrichment)
	#phyper(Overlap-1, group2, Total-group2, group1, lower.tail= TRUE) 
	
	hyp<-phyper(n-1,m,X-m,Y,lower.tail=F) 
	
	mat[i,1]<-signif(hyp,3)
	mat[i,2]<-n
	
}
rownames(mat)<-unique_color
mat<-cbind(mat,p.adjust(mat[,1],method="BH"))
colnames(mat)<-c("p-value","N","FDR")

mat<-as.data.frame(mat)
FDRmat<-mat[mat$FDR<=0.05,]

overep_filename<-paste(HivFu,"OverRep.csv",sep=".")
write.csv(mat,file=overep_filename)
overepFDR_filename<-paste(HivFu,"OverRepFDR.csv",sep=".")
write.csv(FDRmat,file=overepFDR_filename)

for(j in 1:nrow(FDRmat)){
	rowname<-rownames(FDRmat)[j]
	var<-paste("kME",rowname,sep="")
	FilterGenes= abs(datKME[[var]])>.7
	hubgenes<-dimnames(data.frame(dataExpr0))[[2]][FilterGenes]
	
	gene_list<-unique(c(hubgenes,res$genes))
	filename<-paste(HivFu,rowname,"gene_list.csv",sep=".")
	write.csv(gene_list,file=filename)
}



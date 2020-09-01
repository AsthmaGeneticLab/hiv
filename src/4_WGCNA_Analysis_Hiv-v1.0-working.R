#This file contains computational R code for weighted gene co-expression network analysis (WGCNA).
# WGCNA - functions - https://www.rdocumentation.org/packages/WGCNA/versions/1.25-1

# Change the path to the directory
# On Windows use a forward slash / instead of the usual /.
setwd("D:/R-scripts/wgcna-analysis/Drakenstein")

library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)
allowWGCNAThreads(n=34) # in practice based on core size of the machine, it should be set, e.g., Breen's paper used 34

#Input normalized, quality controlled gene expression matrix
M1 <- read.delim("data/DrakensteinUCBExprs_updated/ExprsMat_GitHub.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE)
colnames(M1)=names(M1)
dataExpr0<-as.data.frame(t(M1))
#head(dataExpr0)
dim(dataExpr0)
gsg=goodSamplesGenes(dataExpr0,verbose=3)
gsg$allOK

nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

#Input MetaData File
traitData <- read.csv("data/DrakensteinUCBExprs_updated/MetadataHiv_Phenotype_Only_Ordered_WGCNA.csv")
dim(traitData)

#MATCH TRAITS TO Gene Expression Data
rowsExpr <- rownames(dataExpr0)
traitRows <- match(rowsExpr,traitData$id)
datTraits = traitData[traitRows, -1];
#datTraits = traitData[traitRows,];
rownames (datTraits) = traitData[traitRows, 1];
table(rownames(datTraits)==rownames(dataExpr0)) # EVERYTHING OK?
names(datTraits)

powers=c(1:20) # in practice this should include powers up to 20.
sft0=pickSoftThreshold(dataExpr0,powerVector=powers, networkType="signed")

pdf("Soft_Threshold.pdf")
par(mfrow=c(1,2))
plot(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2], xlab="Soft Threshold (power)",ylab="SFT, signed R^2",type="n",main=paste("Scale independence"))
text(sft0$fitIndices[,1],-sign(sft0$fitIndices[,3])*sft0$fitIndices[,2],labels=powers,col="red")
abline(h=0.90,col="red")    #CHOOSE A  R^2 CUT-OFF OF H
plot(sft0$fitIndices[,1],sft0$fitIndices[,5],type="n",
xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft0$fitIndices[,1],sft0$fitIndices[,5],labels=powers,col="red")
dev.off()

adjacencyPre = adjacency((dataExpr0),power=12 ,type="signed") # for signed network the defualt value is 12.
diag(adjacencyPre)=0
dissTOMPre   = 1-TOMsimilarity(adjacencyPre, TOMType="signed")
geneTreePre  = hclust(as.dist(dissTOMPre), method="average")

#Plot gene tree
pdf("GeneTree.pdf",height=6,width=12)
par(mfrow=c(1,1))
plot(geneTreePre,xlab="",sub="",main="Gene clustering on TOM-based dissimilarity", labels=FALSE,hang=0.04);
dev.off()

#MODULE ASSIGNMENTS with cutreeHybrid algorithm
mColorh=NULL
#a medium sensitivity is (deepSplit=2) - from hovarth book
#If you are concerned about apparently distinct branches being lumped into a single module, 
#you can (1) decrease mergeCutHeight (this should also get rid of modules with genes in different blocks) and (2) increase deepSplit. - ref: https://support.bioconductor.org/p/94677/ 
for (ds in 0:4){ 
 tree = cutreeHybrid(dendro = geneTreePre, pamStage=FALSE,
   minClusterSize = (50), cutHeight = 0.99, 
   deepSplit = ds, distM = dissTOMPre)
 mColorh=cbind(mColorh,labels2colors(tree$labels));
}

#Plot dendrogram and module assignments with deep-split options
pdf("DeepSplit_Choices.pdf", height=10,width=25); 
plotDendroAndColors(geneTreePre, mColorh, paste("dpSplt =",0:4), main = "Co-Expression Network",dendroLabels=FALSE);
dev.off()

#SET DEEP SPLIT CHOICE AND NAME OUR COLORS
modulesPRE =  mColorh[,3]
table(modulesPRE)

#Check to see if network modules should be cut and merged...
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")

pdf("Module_Relationshipsq.pdf")
plot(METree, main ="Clustering of Module Eigengenes",xlab ="",sub="")
abline(h=0.1,col="red")
abline(h=0.15,col="red")
abline(h=0.2,col="red")
abline(h=0.25,col="red")
dev.off()
MEDissThres = 0.25

####################################### SHOULD WE MERGE MODULES BASED ON THE ABOVE? #######################################
merge = mergeCloseModules(dataExpr0, modulesPRE, cutHeight = MEDissThres, verbose=3)
mergedColors = merge$colors
table(mergedColors)
table(modulesPRE)
mergedMEs = merge$newME

#COMPARE UNMERGED MODULES TO MERGED MODULES
pdf("PreDeploy_Network_Unmerged_Merged.pdf", w=9)
#plotDendroAndColors(geneTreeControlPre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(geneTreePre, cbind(modulesPRE, mergedColors),c("Dynamic Tree Cut","Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()

modulesPRE = merge$colors
table(modulesPRE)
MEList = moduleEigengenes(dataExpr0, colors=modulesPRE)
MEs=MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method ="average")
##########################################################################################################################
##Our final set of parameters
# Beta value = 12
# deep split = 4
# Min Module Size = 50
# cutHeight = 0.99
# Merge modules past height = 0.25

#CALCULATE PC FOR VISUALIZATION
PCsPD    = moduleEigengenes((dataExpr0),  colors=modulesPRE) 
ME_PD    = PCsPD$eigengenes
distPCPD = 1-abs(cor(ME_PD,use="p"))
distPCPD = ifelse(is.na(distPCPD), 0, distPCPD)
pcTreePD = hclust(as.dist(distPCPD),method="average") 
MDS_PD   = cmdscale(as.dist(distPCPD),2)
colorsPD = names(table(modulesPRE))
names = row.names((dataExpr0))

#Plot diagnostic module plots
pdf("Module_Visualizationq.pdf",height=8,width=8)
par(mfrow=c(1,1), mar=c(5.1, 4.1, 4.1, 3) + 0.1, cex=1)
plot(pcTreePD, xlab="",ylab="",main="",sub="")
plot(MDS_PD, col= colorsPD,  main="MDS plot", cex=2, pch=19)

for (which.module in names(table(modulesPRE)))
{
 par(mfrow=c(2,1), mar=c(4, 4.1, 4.1, 2))
 plotMat(t(scale(dataExpr0[,modulesPRE==which.module])),
,cex.axis=2,nrgcols=100,rlabels=F,tck=0, rcols=which.module,main=paste("Heatmap",which.module,"Module"))

  ME = ME_PD[, paste("ME",which.module, sep="")] 
  n<- barplot(ME, col=which.module, cex.main=1, ylab="Eigengene Expression",xlab="")
  axis(1,at=n, labels=row.names(dataExpr0), las=2, cex.axis=0.4, font=2)
};
dev.off();

#DETERMINE GENE-TRAIT RELATIONSHIPS ACROSS ALL GENES AND MODULES
names(datTraits)
HIV = as.data.frame(datTraits$HIVn)
names(HIV)="HIV"
GS.HIV=as.numeric(cor(dataExpr0,HIV,use="p"))
GS.HIVColor=numbers2colors(GS.HIV,signed=T)

datColors0=data.frame(modulesPRE, GS.HIVColor)

#Plot gene tree and color bars
pdf("Drakenstein_ColoredModules.pdf",height=8,width=14)
plotDendroAndColors(geneTreePre, colors=datColors0, main="Gene Dendrogram and Module Colors (Beta=12)", groupLabels=c("Module colors", "HIV"), dendroLabels=FALSE, hang=0.03, addGuide=FALSE, guideHang=0.05) 
dev.off()

# Write out ME table
MEs0  =  moduleEigengenes(dataExpr0,  modulesPRE)$eigengenes
MEs = orderMEs(MEs0)
write.table(MEs, "Drakenstein_MEs_q1.txt", sep="\t")
write.table(MEs, "Drakenstein_MEs_q1.csv", sep=",")


#Quantify ME-trait asssociations 
moduleTraitCor = cor(MEs, datTraits, use = "p") #uses pearson correlation -may consider method='spearman' for particular experiments, variables)
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
data = cbind(moduleTraitCor,moduleTraitPvalue)
pdf("Modules_vs_Traits.pdf")
textMatrix = paste(signif(moduleTraitCor, 3), "\n(",
signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 5, 5));
labeledHeatmap(Matrix = moduleTraitCor,
xLabels = names(datTraits),
yLabels = names(MEs),
ySymbols = names(MEs),
colorLabels = FALSE,
colors = blueWhiteRed(50),
textMatrix = textMatrix,
setStdMargins = FALSE,
cex.text = 0.4,
zlim = c(-.3,.3),
main = paste("ME-Trait Relationships"))
dev.off()

#DETERMINE MODULE SIG. USING -LOG PVALUE FOR EACH INDIVIDUAL GENE
#gene sorting code are available in dge_gene_sort.R 
datSummary <- read.delim("Data/DGE_output/DGE_Hiv_result4.2_Sorted.txt", check.names=TRUE, stringsAsFactors=FALSE, row.names=1, header=TRUE) #THIS FILE CONTAINS A VECTOR OF DGE P-VALUES FROM DIFFERENTIAL GENE EXPRESSION ANALYSES, GENES NEED TO BE SORTED IN THE EXACT SAME ORDER AS MATRIX INPUT
#head(datSummary)
GS_HIV=-log10(datSummary$P.Value) #HIV DGE P-VALUES

#PLOT MODULE SIGNIFICANCE
pdf("ModuleSignificance.pdf")
par(mfrow=c(1,1))
verboseBarplot(GS_HIV,modulesPRE,color=colorsPD, border=colorsPD,main="HIV" ,xlab="",ylab="MS (-log10 P-value)" ,las=2, cex.axis=0.7,cex.lab=0.7, cex=0.7,cex.main=0.9, ylim=c(0,1.2))
abline(h=0.6,col="black",lty = 2)   
dev.off()

##CELL-TYPE AND GENE-SET ENRICHMENT
Gene  = colnames(dataExpr0)
enrichments = userListEnrichment(Gene, modulesPRE,
fnIn = NULL,
catNmIn = NULL,
#fnIn = c("GeneList","ModuleColors"),
#catNmIn =  c("Genes","Modules"),
nameOut = "ModuleEnrichment2.csv", 
useBrainLists = TRUE,
useBloodAtlases = TRUE,
omitCategories = "grey", 
outputCorrectedPvalues = TRUE,
useStemCellLists = TRUE, 
outputGenes = TRUE, 
minGenesInCategory = 2, 
useBrainRegionMarkers = TRUE,
useImmunePathwayLists = TRUE,
usePalazzoloWang = TRUE)

#PLOT CONNECTIVITY VS GS.TRAIT MEASURES
CM_Pre=intramodularConnectivity(adjacencyPre,colors=modulesPRE)
names(CM_Pre) 
whichmodule="black";
restrict1=modulesPRE==whichmodule
pdf("Black_Connectivity.pdf")
verboseScatterplot (CM_Pre$kWithin[restrict1], GS.HIV[restrict1],col=modulesPRE[restrict1],
xlab="Connectivity (k) ",ylab="G.S. Eventual HIV(p)") 
dev.off()

whichmodule="purple";
restrict1=modulesPRE==whichmodule
pdf("Purple_Connectivity.pdf")
verboseScatterplot (CM_Pre$kWithin[restrict1], GS.HIV[restrict1],col=modulesPRE[restrict1],
xlab="Connectivity (k) ",ylab="G.S. Eventual HIV(p)") 
dev.off()

# CALCULATE MODULE MEMBERSHIP VALUES (aka. module eigengene based connectivity kME)
datKME=signedKME(dataExpr0, MEs)
colorOfColumn=substring(names(datKME),4)

pdf("Regress_Modules.pdf",w=13,h=5)
par(mfrow = c(2,5))
 for (module in names(table(modulesPRE))) 
  {
  column = match(module,colorOfColumn) 
  restModule=modulesPRE==module
  
  verboseScatterplot(datKME[restModule,column],GS.HIV[restModule],
  xlab=paste("MM ",module,""),ylab="GS.HIV",
  main=paste("kME vs GS"),col=module,abline=T,las=1, cex.axis=0.9, pch=16)  
  }
dev.off()

#MODULE MEMBERSHIP (kME) KME is defined as correlation between expression and modules
geneModuleMembership1 = signedKME((dataExpr0), MEs)
colnames(geneModuleMembership1)=paste("PC",colorsPD,".cor",sep=""); 
MMPvalue1=corPvalueStudent(as.matrix(geneModuleMembership1),dim(dataExpr0)[[2]]); 
colnames(MMPvalue1)=paste("PC",colorsPD,".pval",sep="");

Gene       = rownames(t(dataExpr0))
kMEtable1  = cbind(Gene,Gene,modulesPRE)
for (i in 1:length(colorsPD))
kMEtable1 = cbind(kMEtable1, geneModuleMembership1[,i], MMPvalue1[,i])
colnames(kMEtable1)=c("PSID","Gene","Module",sort(c(colnames(geneModuleMembership1), colnames(MMPvalue1))))
write.csv(kMEtable1,"Drakenstein_kMEtable.csv",row.names=FALSE)

topGenesKME = NULL
for (c in 1:length(colorsPD)){
 kMErank1    = rank(-geneModuleMembership1[,c])
 maxKMErank  = rank(apply(cbind(kMErank1+.00001),1,max))
 topGenesKME = cbind(topGenesKME,Gene[maxKMErank<=20])
}; 
colnames(topGenesKME) = colorsPD
topGenesKME #OUTPUT TOP HUB GENES

#WRITE OUTPUT
#datTraits = datTraits[, -1] #if the id column in not deleted earlier then we need to delete here.
Pre<-as.data.frame(datTraits)
names(Pre)<-"Pre"
modNames = substring(names(MEs), 3)
nGenes = ncol(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
nSamples = nrow(dataExpr0) #  MODULE MEMBERSHIP FROM CORRELATION BETWEEN EIGENGENE AND GENE
geneModuleMembership<-as.data.frame(cor(dataExpr0,MEs,use = "p")) #PEARSON CORRELATION MEMBERSHIP OF EACH GENE TO A MODULE 
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples)) #P VALUE SIGNIFICANCE FOR EACH GENE IN EACH MODULE

MM_names<-names(geneModuleMembership)
MM_names<-substring(MM_names,3,length(MM_names))
names(geneModuleMembership)<-paste("MM.",MM_names,sep="")
names(MMPvalue)<-paste("Mp.",MM_names,sep="")

geneTraitSignificance = as.data.frame(cor(dataExpr0, Pre, use = "p")) #  Generate correlations and p-value for each gene against the trait.  
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Pre), sep="")
names(GSPvalue) = paste("p.GS.", names(Pre), sep="")

Colors=modulesPRE
tempout<-cbind(geneModuleMembership,MMPvalue,Colors)
sortedout<-tempout[,sort(names(tempout))]
geneinfo<-sortedout[order(sortedout[,1]),]
write.csv(geneinfo,file="Network_WGCNA_OUTPUT_q.csv")

#HERE WE EXPORT 20 TOP INTRAMODULAR HUB GENES PER MODULE OF INTEREST

modules = c("black") #"purple", "yellow", "lightgreen" #,"yellow", "lightgreen", "black"
probes=names(dataExpr0)
inModule=is.finite(match(modulesPRE, modules));
modProbes=probes[inModule] ##
modTom=dissTOMPre[inModule, inModule]
dimnames(modTom)=list(modProbes, modProbes)
nTopHubs = 50
kIN = softConnectivity(dataExpr0[, modProbes])
selectHubs = (rank (-kIN) <= nTopHubs)

# saving the image of the data so that it can be used later
save.image(file = "WGCNA/Drakenstein_WGCNA_pow12_all_on_DGE_with_all_covar_using_DGE_rslt_v4.2.RData")
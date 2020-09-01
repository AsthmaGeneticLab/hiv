#This scripts find out Linear Regresssion on Tidal Volume (TV) for Hub Genes of Black module
#This work is based on DGE results v4.2 with all covariates/confounding where mode of delivery, Celect & Cemerg have been merged and delivery_merged_is_natural is generated Yes/No field. Yes=Natural and No=Celect or Cemerg
#1) It collects Hub Genes intramodular connectivity > 0.7 which returns 177 Hub genes
#2) Filters ExprsMat with these 177 genes
#3) Read LungMetaData2yr4Hiv for TV (Tidal volume at 24 months, i.e., 2 year)
#4) Read LungMetaData4HivCovar for covariates/confounding variables
#5) Combine 177 genes ExprsMat data with covariates/confounding variables for TV
#6) Prepare TV linear regression for each gene with all covariates  

# Change the path to the directory
# On Windows use a forward slash / instead of the usual /.
setwd("D:/R-scripts/hiv_chord_gene")
options(stringsAsFactors  =  FALSE)

library(WGCNA)
library(cluster)
options(stringsAsFactors  =  FALSE)

load('WGCNA/Drakenstein_WGCNA_pow12_all_on_DGE_with_all_covar_using_DGE_rslt_v4.2.RData')
head(GS.HIV) 
head(datKME) 

var<-paste("kME","black",sep="")
FilterGenes=abs(datKME[[var]])>.7 #returns 177 genes
hubgenes<-dimnames(data.frame(dataExpr0))[[2]][FilterGenes] 
#hubgenes

ExprsMat <- read.delim("Data/DrakensteinUCBExprs_updated/ExprsMat_GitHub.txt", check.names=TRUE, stringsAsFactors=FALSE, header=TRUE)
colnames(ExprsMat)=names(ExprsMat)
rownames(ExprsMat)=ExprsMat$SYMBOL;

matchedRows <- match(hubgenes,ExprsMat$SYMBOL)
head(matchedRows)
blackHubGeneExprsMat <- ExprsMat[matchedRows,]
head(colnames(blackHubGeneExprsMat))
head(rownames(blackHubGeneExprsMat))

sampleWiseBlackHubGeneExpr<-as.data.frame(t(blackHubGeneExprsMat))
head(colnames(sampleWiseBlackHubGeneExpr))
head(rownames(sampleWiseBlackHubGeneExpr))
dim(sampleWiseBlackHubGeneExpr) #it return 145 177, so need to remove the first row it make it 144 
#rownames(sampleWiseBlackHubGeneExpr)
sampleWiseBlackHubGeneExpr = sampleWiseBlackHubGeneExpr[-1,]
#head(sampleWiseBlackHubGeneExpr)
sampleWiseBlackHubGeneExpr2 <- as.data.frame(sapply(sampleWiseBlackHubGeneExpr, as.numeric)) #<- sapply is here to convert string to number - but loss the rownames 
rownames(sampleWiseBlackHubGeneExpr2) <- rownames(sampleWiseBlackHubGeneExpr)
head(sampleWiseBlackHubGeneExpr2$C11orf10)
head(sampleWiseBlackHubGeneExpr2[,1])

lungMetaData <- read.csv("Data/DrakensteinUCBExprs_updated/LungMetaData2yr4Hiv.csv")
colnames(lungMetaData)=names(lungMetaData)
dim(lungMetaData)

lungMetaDataCovar <- read.csv("Data/DrakensteinUCBExprs_updated/LungMetaData4HivCovar.csv")
colnames(lungMetaDataCovar)=names(lungMetaDataCovar)
dim(lungMetaDataCovar)

totalGenes <- ncol(sampleWiseBlackHubGeneExpr2)
geneCoefMat <- NULL
for (geneIndex in c(1:totalGenes)){
	geneName <- colnames(sampleWiseBlackHubGeneExpr2)[geneIndex]
	print(geneName)
	singleBlackHubGeneExpr = cbind(sampleWiseBlackHubGeneExpr2[,geneIndex], lungMetaDataCovar)
	#colnames(singleBlackHubGeneExpr)[1] <- geneName
	colnames(singleBlackHubGeneExpr)[1] <- "Gene"
	rownames(singleBlackHubGeneExpr) <- rownames(sampleWiseBlackHubGeneExpr2)
	singleBlackHubGeneExpr <- singleBlackHubGeneExpr[ , -which(names(singleBlackHubGeneExpr) %in% c("Barcode"))]

	ctable <- cbind(singleBlackHubGeneExpr, lungMetaData)
	ctable <- ctable[ , -which(names(ctable) %in% c("Barcode"))]

	ZWeight <- (ctable$zwei_6M)
	Gender <- factor(ctable$gender_n)
	GestationalAge <- (ctable$gestation_delivery)
	SmokeCat <- factor(ctable$smokecat_n)
	
	tvBlackHubGeneLm <- lm(tidal_vol_24M ~ Gene + ZWeight + Gender + GestationalAge + SmokeCat, ctable)	
	geneNameWithExt <- paste(geneName,"_with_covar.txt",sep="")
	outputFileName <- paste("WGCNA/HubGenes177/LuFu_TV_2yr_lm_output_4_", geneNameWithExt, sep="")

	geneNameWithExt2 <- paste(geneName,"_with_covar.csv",sep="")
	outputFileName2 <- paste("WGCNA/HubGenes177/LuFu_TV_2yr_lm_CoefMat_4_", geneNameWithExt2, sep="")
	CoefMat <- summary(tvBlackHubGeneLm)$coefficients
	#write.csv (CoefMat, outputFileName2)  
	
	if (geneIndex == 1) {
		geneCoefMat <- CoefMat[2,]
	} else {
		geneCoefMat <- rbind(geneCoefMat, CoefMat[2,])
	}
}
rownames(geneCoefMat) <- colnames(sampleWiseBlackHubGeneExpr2)
outputFileNameAll <- paste("WGCNA/HubGenes177/LuFu_TV_2yr_lm_CoefMat_4_", "All_genes_with_covar.csv", sep="")
write.csv (geneCoefMat, outputFileNameAll) 


geneCoefMatWithFDR<-cbind(geneCoefMat,p.adjust(geneCoefMat[,4],method="BH"))
colnames(geneCoefMatWithFDR)[5] <- "FDR"
write.csv (geneCoefMatWithFDR, "WGCNA/HubGenes177/LuFu_TV_2yr_lm_CoefMat_4_All_genes_with_covar_WithFDR.csv");

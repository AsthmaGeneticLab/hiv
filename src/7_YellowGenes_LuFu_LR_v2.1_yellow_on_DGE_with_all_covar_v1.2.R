#This scripts find out Linear Regresssion on Tidal Volume (TV) for Hub Genes of Yellow module
#This work is based on DGE results with all covariates/confounding
#1) It collects Hub Genes intramodular connectivity > 0.7 which returns 191 Hub genes
#2) Filters ExprsMat with these 191 genes
#3) Read LungMetaData4Hiv for TV
#4) Read LungMetaData4HivCovar for covariates/confounding variables
#5) Combine 191 genes ExprsMat data with covariates/confounding variables for TV
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

var<-paste("kME","yellow",sep="")
FilterGenes=abs(datKME[[var]])>.7 #returns 191 genes
hubgenes<-dimnames(data.frame(dataExpr0))[[2]][FilterGenes] 
#hubgenes

ExprsMat <- read.delim("Data/DrakensteinUCBExprs_updated/ExprsMat_GitHub.txt", check.names=TRUE, stringsAsFactors=FALSE, header=TRUE)
colnames(ExprsMat)=names(ExprsMat)
rownames(ExprsMat)=ExprsMat$SYMBOL;
ExprsMat$SYMBOL[9658]<-'TMEM189.UBE2V1' #9658 positioned Gene name is 'TMEM189-UBE2V1' which is normally converted TMEM189.UBE2V1 in WGCNA
rownames(ExprsMat)[9658]<-'TMEM189.UBE2V1'

matchedRows <- match(hubgenes,ExprsMat$SYMBOL)
head(matchedRows)
yellowHubGeneExprsMat <- ExprsMat[matchedRows,]
head(colnames(yellowHubGeneExprsMat))
head(rownames(yellowHubGeneExprsMat))

sampleWiseYellowHubGeneExpr<-as.data.frame(t(yellowHubGeneExprsMat))
head(colnames(sampleWiseYellowHubGeneExpr))
head(rownames(sampleWiseYellowHubGeneExpr))
dim(sampleWiseYellowHubGeneExpr) #it return 145 191, so need to remove the first row it make it 144 
#rownames(sampleWiseYellowHubGeneExpr)
sampleWiseYellowHubGeneExpr = sampleWiseYellowHubGeneExpr[-1,]
#head(sampleWiseYellowHubGeneExpr)
sampleWiseYellowHubGeneExpr2 <- as.data.frame(sapply(sampleWiseYellowHubGeneExpr, as.numeric)) #<- sapply is here to convert string to number - but loss the rownames 
rownames(sampleWiseYellowHubGeneExpr2) <- rownames(sampleWiseYellowHubGeneExpr)
head(sampleWiseYellowHubGeneExpr2$MARCH8)
head(sampleWiseYellowHubGeneExpr2[,1])

lungMetaData <- read.csv("Data/DrakensteinUCBExprs_updated/LungMetaData4Hiv.csv")
colnames(lungMetaData)=names(lungMetaData)
dim(lungMetaData)

lungMetaDataCovar <- read.csv("Data/DrakensteinUCBExprs_updated/LungMetaData4HivCovar.csv")
colnames(lungMetaDataCovar)=names(lungMetaDataCovar)
dim(lungMetaDataCovar)

totalGenes <- ncol(sampleWiseYellowHubGeneExpr2)
geneCoefMat <- NULL
for (geneIndex in c(1:totalGenes)){
	geneName <- colnames(sampleWiseYellowHubGeneExpr2)[geneIndex]
	print(geneName)
	singleBlackHubGeneExpr = cbind(sampleWiseYellowHubGeneExpr2[,geneIndex], lungMetaDataCovar)
	#colnames(singleBlackHubGeneExpr)[1] <- geneName
	colnames(singleBlackHubGeneExpr)[1] <- "Gene"
	rownames(singleBlackHubGeneExpr) <- rownames(sampleWiseYellowHubGeneExpr2)
	singleBlackHubGeneExpr <- singleBlackHubGeneExpr[ , -which(names(singleBlackHubGeneExpr) %in% c("Barcode"))]

	ctable <- cbind(singleBlackHubGeneExpr, lungMetaData)
	ctable <- ctable[ , -which(names(ctable) %in% c("Barcode"))]

	ZWeight <- (ctable$zwei_6M)
	Gender <- factor(ctable$gender_n)
	GestationalAge <- (ctable$gestation_delivery)
	SmokeCat <- factor(ctable$smokecat_n)
	
	tvBlackHubGeneLm <- lm(tidal_vol__6wks ~ Gene + ZWeight + Gender + GestationalAge + SmokeCat, ctable)	
	geneNameWithExt <- paste(geneName,"_with_covar.txt",sep="")
	outputFileName <- paste("WGCNA/YellowHubGenes191/LuFu_TV_lm_output_4_", geneNameWithExt, sep="")

	geneNameWithExt2 <- paste(geneName,"_with_covar.csv",sep="")
	outputFileName2 <- paste("WGCNA/YellowHubGenes191/LuFu_TV_lm_CoefMat_4_", geneNameWithExt2, sep="")
	CoefMat <- summary(tvBlackHubGeneLm)$coefficients
	#write.csv (CoefMat, outputFileName2)  
	
	if (geneIndex == 1) {
		geneCoefMat <- CoefMat[2,]
	} else {
		geneCoefMat <- rbind(geneCoefMat, CoefMat[2,])
	}
}
rownames(geneCoefMat) <- colnames(sampleWiseYellowHubGeneExpr2)
outputFileNameAll <- paste("WGCNA/YellowHubGenes191/LuFu_TV_lm_CoefMat_4_", "All_genes_with_covar.csv", sep="")
write.csv (geneCoefMat, outputFileNameAll) 

#############################################################################################
########### Sorting genes on the output of DGE analysis with all confounding variables ######
#############################################################################################

# Change the path to the directory
# On Windows use a forward slash / instead of the usual /.
setwd("D:/R-scripts/hiv_chord_gene")

exprsTable <- read.delim("Data/DrakensteinUCBExprs_updated/ExprsMat_GitHub.txt", check.names=FALSE, stringsAsFactors=FALSE, header=T)
dim(exprsTable)

dgeTable <- read.csv(file = 'Data/DGE_output/DGE_Hiv_result4.2.csv')
colnames(dgeTable) = names(dgeTable)
dim(exprsTable)

matchedRows <- match(exprsTable$SYMBOL,dgeTable$genes)
head(matchedRows)

sortedTable = dgeTable[matchedRows,]
head(sortedTable)

write.table(sortedTable, file="Data/DGE_output/DGE_Hiv_result4.2_Sorted.txt", quote=FALSE, row.names=FALSE, sep="\t")


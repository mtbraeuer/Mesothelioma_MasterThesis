## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 27.TCGA-MESO_stratification_hubgenes.R
##
## Purpose of script: Stratification of hub genes
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 22. September 2020
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
## http://r-addict.com/2016/11/21/Optimal-Cutpoint-maxstat.html
## https://rpkgs.datanovia.com/survminer/index.html
##
## ------------------------------------------------------
source("Initialize_TCGA-MESO.R")
#--------------------------------------------------
scriptName <- "27.TCGA-MESO_stratification_hubgenes"

#----------------------------------------------
# DATASET loading 
#----------------------------------------------
# Expressionset
Eset <- readRDS(file = file.path(resultsdir,
                                 "10.TCGA-MESO_Eset.ExpressionSet.RDS"))

EsetAnalysis <- Eset
pData <- pData(EsetAnalysis)
expressionMatrix <- exprs(EsetAnalysis)

#----------------------------------------------
# DATASET loading 
#----------------------------------------------
# Expressionset
Eset <- readRDS(file = file.path(resultsdir,
                                 "10.TCGA-MESO_Eset.ExpressionSet.RDS"))

EsetAnalysis <- Eset
pData <- pData(EsetAnalysis)
expressionMatrix <- exprs(EsetAnalysis)

GOI <- read.table(file.path("../hub_genes.txt"),
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  header = TRUE)

# patient ID, overall survival, vital status
stratificationDF <- as.data.frame(pData[c(2,3,4)])

# save pval of survival curves
pvalDF <- as.data.frame(matrix(data = NA,
                               nrow = nrow(GOI),
                               ncol = 2))
colnames(pvalDF) <- c("GeneSymbol", "pval")
pvalDF$GeneSymbol <- GOI$GeneSymbol
#i=1
for (i in 1:nrow(GOI)) {
    
    goiEntrezID <- as.character(GOI$EntrezID[i])
    goiName <- as.character(GOI$GeneSymbol[i])
    
    x <- StratificationLowHighExpressors_SurvivialCurves(expressionMatrix = expressionMatrix,
                                                         stratificationDF = stratificationDF,
                                                         goiEntrezID = goiEntrezID,
                                                         goiName = goiName)
    pData <- cbind(pData, x[[1]])
    colnames(pData)[ncol(pData)] <- paste("strat_", goiName, sep = "")
    
    pvalDF$pval[i] <- x[[2]]
    
}

#----------------------------------------------
# write overview table for survival data

# data frame
overviewDF <- as.data.frame(matrix(data = NA,
                                   nrow = nrow(GOI),
                                   ncol = 5))
colnames(overviewDF) <- c("GeneSymbol",
                          "Number_High",
                          "Number_Low",
                          "Median_OS_High_days",
                          "Median_OS_Low_days")
overviewDF$GeneSymbol <- GOI$GeneSymbol

for (i in 1:nrow(overviewDF)) { 
    log <- read.table(file.path(resultsdir,
                                paste("27.TCGA-MESO_stratification_hubgenes.",
                                      overviewDF[i,1],
                                      ".log.txt",
                                      sep="")),
                      sep = "\t")
    overviewDF$Number_High[i] <- log[1,1]
    overviewDF$Number_Low[i] <- log[2,1]
    overviewDF$Median_OS_High_days[i] <- log[3,1]
    overviewDF$Median_OS_Low_days[i] <- log[4,1]
    
}
overviewDF["Percent_High"] <- round((overviewDF$Number_High*100)/nrow(pData))
overviewDF["PErcent_Low"] <- round((overviewDF$Number_Low*100)/nrow(pData))

write.table(overviewDF,
            file.path(resultsdir,
                      paste(scriptName,
                            "overview_OS",
                            "txt",
                            sep = ".")),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

# save pvalues of survival curves
pvalDF$pval <- as.numeric(substr(pvalDF$pval, start = 4, stop = nchar(pvalDF$pval)))

# get siginificance levels
pvalDF$GeneSymbol[pvalDF$pval<=0.05]
pvalDF$GeneSymbol[pvalDF$pval<=0.01]
pvalDF$GeneSymbol[pvalDF$pval<=0.001]

write.table(pvalDF,
            file.path(resultsdir,
                      paste(scriptName,
                            "pvalues",
                            "txt",
                            sep = ".")),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE)

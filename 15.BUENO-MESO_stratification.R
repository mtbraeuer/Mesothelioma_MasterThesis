## ------------------------------------------------------
## Project: BUENO-MESO
## Script name: 15.BUENO-MESO_stratification.R
##
## Purpose of script: stratification into high and low expressors of GOI
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 25. January 2019
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
##
## ------------------------------------------------------
source("Initialize_BUENO-MESO.R")
scriptName <- "15.BUENO-MESO_stratification"

#----------------------------------------------
# DATASET loading 
#----------------------------------------------
# Expressionset
Eset <- readRDS(file = file.path(resultsdir,
                                 "10.BUENO-MESO_readKallisto.ExpressionSet.RDS"))

EsetAnalysis <- Eset
pData <- pData(EsetAnalysis)
expressionMatrix <- exprs(EsetAnalysis)

GOI <- read.table(file.path("../GOI.txt"),
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  header = TRUE)



# patient ID, overall survival, vital status
pData$VitalStatus <- as.factor(pData$VitalStatus)
colnames(pData)
stratificationDF <- as.data.frame(pData[c(1, 16, 3)])
colnames(stratificationDF)
str(stratificationDF)

# save pval of survival curves
pvalDF <- as.data.frame(matrix(data = NA,
                               nrow = nrow(GOI),
                               ncol = 2))
colnames(pvalDF) <- c("GeneSymbol", "pval")
pvalDF$GeneSymbol <- GOI$GeneSymbol




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
# resultsdir/15.TCGA-MESO_Stratification.TERT.log.txt

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
                                paste(scriptName,
                                      overviewDF[i,1],
                                      "log.txt",
                                      sep=".")),
                      sep = "\t")
    overviewDF$Number_High[i] <- log[1,1]
    overviewDF$Number_Low[i] <- log[2,1]
    overviewDF$Mean_OS_High_days[i] <- round(log[3,1],1)
    overviewDF$Mean_OS_Low_days[i] <- round(log[4,1],1)
    
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
            row.names = FALSE,
            sep = "\t")
#----------------------------------------------
# save Eset with stratification infos
pData(EsetAnalysis) <- pData
EsetAnalysis
colnames(pData(EsetAnalysis))

saveRDS(EsetAnalysis, 
        file.path(resultsdir,
                  paste(scriptName,
                        class(Eset),
                        "RDS", 
                        sep=".")))

## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 10.TCGA-MESO_expressionset_construction.R
##
## Purpose of script: Expressionset Construction
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 17. October 2019
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
##
## ------------------------------------------------------
source("Initialize_TCGA-MESO.R")
#--------------------------------------------------
scriptName <- "10.TCGA-MESO_Eset"
minMutEvents <- 4

#----------------------------------------------
# DATASET loading
#----------------------------------------------
# A1: normalized by voom
#voom_exprs <- readRDS(file.path(resultsdir, "10.TCGA-MESO_normalization.EList.Voom.RDS"))
#Exprs <- as.matrix(ExtractCounts(voom_exprs, return_class = NULL))

# A2: normalized by limma-trend
limma_trend <- readRDS(file.path(resultsdir,
                                 "05.TCGA-MESO_normalization.matrix.LogCPM.RDS"))
Exprs <- as.matrix(limma_trend)

# traitdata & sample & clinical annotation information
ClinicalData <- read.delim(file.path(datadir,
                                     "TCGA-MESO.ClinicalInfo.txt"),
                           sep="\t",
                           header = TRUE,
                           stringsAsFactors = FALSE)
pData <- ClinicalData

all(colnames(Exprs)==make.names(pData$barcode))

# feature data information
load(file.path(datadir, 
               "HumanGenes_GRCh38_p12.rda"))

# mutation data
MutationData <- read.table(file.path(datadir,
                                     "TCGA-MESO.Mutationdata_Mutect2.maf"))

#----------------------------------------------
# ExpressionSet construction
#----------------------------------------------

#----------------------------------------------
# 1. pData 
#----------------------------------------------
# change rownames and create column patient_ID
rownames(pData) <- sub("TCGA.", "", make.names(pData$patient))

pData["PatientID"] <- rownames(pData)

pData$barcode <- make.names(pData$barcode)
pData$morphology <- sub("X", "", make.names(pData$morphology))
# get relevant clinical information
pData <- pData[c(58, 41, 14, 12, 13, 3)]
colnames(pData)
colnames(pData) <- c("PatientID",
                     "Gender",
                     "OS_days",
                     "VitalStatus",
                     "Morphology",
                     "Barcode")
colnames(pData)

# morphology
pData$Morphology[pData$Morphology=="9050.3"] <- "meso"
pData$Morphology[pData$Morphology=="9051.3"] <- "sarc"
pData$Morphology[pData$Morphology=="9052.3"] <- "epi"
pData$Morphology[pData$Morphology=="9053.3"] <- "bi"
table(pData$Morphology)

# vital status
pData$VitalStatus <- as.factor(pData$VitalStatus)

str(pData)
#----------------------------------------------
# Mutation Data
#----------------------------------------------
mData <- MutationData[c(2, 16, 1, 5:10)]
colnames(mData)
colnames(mData) <- c("EntrezID",
                     "Barcode",
                     "Symbol",
                     "Chromosome",
                     "StartPos",
                     "EndPos",
                     "Strand",
                     "VariantClassification",
                     "VariantType")
colnames(mData)
str(mData)
# make barcode to names
mData$Barcode <- make.names(mData$Barcode)
mData["PatientID"] <- substr(mData$Barcode, 6,12)

# get all EntrezIDs with >= minMutEvents and their frequencies
mutTable <- as.data.frame(table(mData$EntrezID))
colnames(mutTable) <- c("EntrezID", "Frequencies")
mutTable <- setorder(x = mutTable, Frequencies)
mutTable <- subset(mutTable, Frequencies>=minMutEvents)

table(mutTable$Frequencies)
dim(mutTable)

keep <- match(mutTable$EntrezID, mData$EntrezID)
mutTable["Symbol"] <- mData$Symbol[keep]

#create subset for pData with mut info
subpData <- data.frame(matrix(data = 0,
                              nrow = nrow(pData),
                              ncol = nrow(mutTable)))

rownames(subpData) <- rownames(pData)
colnames(subpData) <- mutTable$EntrezID
dim(subpData)

# row-wise through mData
for (i in 1:nrow(mData)) {
    
    # row-wise throgh subpData - search for patient id
    for (x in 1:nrow(subpData)) {
        
        #check if patient ids matches
        if (rownames(subpData)[x]==mData$PatientID[i] ) {
            
            # col-wise - search for entrezid
            for (y in 1:ncol(subpData)) {
                
                # check if EntrezIDs matches
                if (colnames(subpData)[y]==mData$EntrezID[i]) {
                    
                    # correct mutation entry
                    subpData[x,y] <- 1
                    
                }
            }
        }
    }
}

write.table(mutTable, 
            quote = F,
            sep = "\t",
            row.names = F,
            file.path(resultsdir,
                      paste(scriptName,
                            "mData",
                            "txt",
                            sep=".")))

# subpData colnames
colnames(subpData) <- paste("Mut_", colnames(subpData), sep = "")

# if TRUE then bind to each other
all(rownames(pData)==rownames(subpData))
pData <- cbind(pData, subpData)



#----------------------------------------------
# Exprs - check if everything matches
#----------------------------------------------
#need to return true
all(colnames(Exprs)==pData$Barcode)

# change col names in Exprs to patient ID 
for (i in 1:ncol(Exprs)) {
    if(colnames(Exprs)[i] == pData$Barcode[i])
    {
        colnames(Exprs)[i] <- pData$PatientID[i]
    }else{
        print(paste(i, ":", colnames(Exprs)[i], " - ", pData$Barcode[i]), sep="")
    }
}

# check for rownames pData and colnames Exprs
for (i in 1:ncol(Exprs)) {
    if (colnames(Exprs)[i]==rownames(pData)[i]) {
    }else{
        print(paste("NO match ", i, ":", colnames(Exprs)[i], " - ", rownames(pData)[i]), sep="")
    }
}

# align pData to eset 
keep <- match(colnames(Exprs), rownames(pData))
pData <- pData[keep, ]

# check if number of rows of pData == numbers of expression set
# AND names have to be identical
# both need to return TRUE
all(rownames(pData)==colnames(Exprs))
all(nrow(pData)==ncol(Exprs))

#----------------------------------------------
# fData 
#----------------------------------------------
# mapping expression matrix and fData
fData <- gene.location

gene_entrez <- as.data.frame(mapIds(org.Hs.eg.db, 
                                    keys = fData$ensembl_gene_id,
                                    column = "ENTREZID",
                                    keytype = "ENSEMBL")) #get ENTREZIDs

# vector of entrez ids
colnames(gene_entrez)[1] <- "EntrezID"
# add new col to fdata with entrez id
fData$EntrezID <- as.character(gene_entrez$EntrezID)

# align fData to eset
keep <- match(rownames(Exprs), fData$EntrezID)

fData <- fData[keep,]
colnames(fData)
colnames(fData) <- c("Chromosome",
                     "StartPos",
                     "EndPos",
                     "Strand",
                     "EnsemblGeneID",
                     "Symbol",
                     "EntrezID")
rownames(fData) <- c(1:nrow(fData))

# cehck again for NAs
dim(fData)

table(is.na(fData$EntrezID))
table(is.na(rownames(Exprs)))

rownames(fData) <- fData$EntrezID

# check if align was correct
# both need to return TRUE
all(nrow(fData)==nrow(Exprs))
all(rownames(fData)==rownames(Exprs))

#----------------------------------------------
# Expressionset creation
#----------------------------------------------
Eset <- ExpressionSet(assayData=Exprs,
                      phenoData= AnnotatedDataFrame(pData),
                      featureData = AnnotatedDataFrame(fData),
                      annotation = "org.Hs.eg.db")
Eset


saveRDS(Eset, 
        file.path(resultsdir,
                  paste(scriptName,
                        class(Eset), 
                        "RDS", 
                        sep=".")))

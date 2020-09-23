## ------------------------------------------------------
## Project: BUENO-MESO
## Script name: 05.BUENO-MESO_checkAlignment.R
##
## Purpose of script: check completness of pseudo-alignment with kallisto
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 24.January 2020
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
## additional info
## check RNA-Seq data completness after alignment (bcbionextgen pipeline)
## files in data/rawData/EGAR00001406XXX/abundance.tsv 
## XXX is 109 to 319
## abundance.tsv: target_id(Ensembl transcript ID), length, eff_length, est_counts, tpm
## ckeck columns est_counts & tpm, they should not contain only 0 and NANs
## 
## ------------------------------------------------------
# libibary loading
source("Initialize_BUENO-MESO.R")

#--------------------------------------------------
# scriptName
scriptName <- "05.BUENO-MESO_checkAlignment"

#--------------------------------------------------
# data set loading

pData <- read.table(file.path(datadir,
                              "BUENO-MESO_pData.txt"),
                    sep = "\t",
                    stringsAsFactors = FALSE,
                    header = T)

################################################################################################
# check alignment
################################################################################################

# create checkerDF
checker <- as.data.frame(matrix(data = NA,
                                nrow = nrow(pData),
                                ncol = 6))
colnames(checker) <- c("RNATumorID",
                       "EGAR_ID",
                       "length",
                       "eff_length",
                       "est_counts",
                       "tpm")
checker$RNATumorID <- pData$RNATumor_ID
checker$EGAR_ID <- pData$EGAR_ID


# read all files into checkerDF
for (i in 1:nrow(checker)) {
    
    file <- read.delim(file = paste("data/rawData/",
                                    checker$EGAR_ID[i],
                                    "/",
                                    "abundance.tsv",
                                    sep = ""))
    checker$length <- sum(as.numeric(file[,2]))
    checker$eff_length <- sum(as.numeric(file[,3]))
    checker$est_counts <- sum(as.numeric(file[,4]))
    checker$tpm <- sum(as.numeric(file[,5]))
    print(checker$EGAR_ID[i])
}

table(checker$tpm)

# save checkerDF
write.table(checker,
            file.path(resultsdir,
                      paste(scriptName,
                            "_checkedAlignment.txt",
                            sep = "")),
            quote = FALSE,
            row.names = FALSE)

# copy files to data dir RNATumor_ID_abundance.tsv
for( i in 1:nrow(pData)){
    file.copy(from = paste("data/rawData/",
                           pData$EGAR_ID[i],
                           "/",
                           "abundance.tsv",
                           sep = ""),
              to = paste("data/",
                         pData$RNATumor_ID[i],
                         "_abundance.tsv",
                         sep = ""),
              recursive = FALSE)
    print(pData$RNATumor_ID[i])
}

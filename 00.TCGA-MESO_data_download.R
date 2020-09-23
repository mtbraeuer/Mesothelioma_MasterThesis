## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 00.TCGA-MESO_data_download.R
##
## Purpose of script: data download from TCGA (https://portal.gdc.cancer.gov/)
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
scriptName <- "00.TCGA-MESO_data_download"


#--------------------------------------------------
# download raw counts HTSeq for TCGA-COAD - aligned to hg38
query <- GDCquery(project = project,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query,
            method = "api",
            directory = datadir)
data <- GDCprepare(query,
                   directory = datadir)

# extract the count matrix
count.matrix <- assays(data)$`HTSeq - Counts`

# save the count matrix
write.table(count.matrix,
            file.path(datadir,
                      paste(project, "RawCountMatrix", "txt", sep = ".")),
            row.names = TRUE,
            quote = FALSE,
            sep = "\t")

# save the clinical information
write.table(colData(data),
            file.path(datadir,
                      paste(project, "ClinicalInfo", "txt", sep = ".")),
            row.names = TRUE,
            quote = FALSE,
            sep = "\t")

# save the gene data
write.table(rowData(data),
            file.path(datadir,
                      paste(project, "RawData", "txt", sep = ".")),
            row.names = TRUE,
            quote = FALSE,
            sep = "\t")
# save the summarized experiment
saveRDS(data, file.path(datadir, paste(project, class(data), "RDS", sep = ".")))

#--------------------------------------------------
# Mutation data (hg38)
# This exmaple will download MAF (mutation annotation files) f
# or variant calling pipeline muse. Pipelines options are: muse,  
# varscan2, somaticsniper, mutect. For more information please access GDC docs.

mafMuse <- GDCquery_Maf("MESO", 
                        pipelines = "muse",
                        save.csv = FALSE,
                        directory = datadir)

write.table(mafMuse,
            file.path(datadir,
                      paste(project,
                            "Mutationdata_Muse",
                            "maf",
                            sep=".")))

#--------------------------------------------------
mafVarscan2 <- GDCquery_Maf("MESO", 
                            pipelines = "varscan2",
                            save.csv = FALSE,
                            directory = datadir)

write.table(mafVarscan2,
            file.path(datadir,
                      paste(project,
                            "Mutationdata_Varscan2",
                            "maf",
                            sep=".")))

#--------------------------------------------------
mafSomaticsniper <- GDCquery_Maf("MESO", 
                                 pipelines = "somaticsniper",
                                 save.csv = FALSE,
                                 directory = datadir)

write.table(mafSomaticsniper,
            file.path(datadir,
                      paste(project,
                            "Mutationdata_Somaticsniper",
                            "maf",
                            sep=".")))

#--------------------------------------------------
mafMutect2 <- GDCquery_Maf("MESO", 
                           pipelines = "mutect2",
                           save.csv = FALSE,
                           directory = datadir)

write.table(mafMutect2,
            file.path(datadir,
                      paste(project,
                            "Mutationdata_Mutect2",
                            "maf",
                            sep=".")))

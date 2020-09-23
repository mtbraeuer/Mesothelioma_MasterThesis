## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 05.TCGA-MESO_normalization.R
##
## Purpose of script: NORMALIZATION of RNA HTcounts
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
scriptName <- "05.TCGA-MESO_normalization"

#----------------------------------------------
# DATASET loading
#----------------------------------------------
# expression data
Expressionset <- as.matrix(read.table(file.path(datadir,
                                                "TCGA-MESO.RawCountMatrix.txt"),
                                      sep="\t",
                                      header = TRUE))
GOI <- read.table(file.path("../GOI.txt"),
                  sep = "\t",
                  stringsAsFactors = FALSE,
                  header = TRUE)
#----------------------------------------------
# Mapping IDs
#----------------------------------------------
# attention check of last update of annotation org.Hs.eg!!!!!
Exprs <- MappingIDs(inmatrix = Expressionset, 
                    annotation = "org.Hs.eg", 
                    map_origin = "ENSEMBL", 
                    map_target = "ENTREZID")
dim(Exprs)

# check if GOI are still in
all(GOI$EntrezID %in% rownames(Exprs))
#----------------------------------------------
# Filtering
#----------------------------------------------
ExprsFiltered <- FilterByCPM(inmatrix = Exprs,
                             cutoff = 0.1, 
                             threshold = 0.2)
dim(ExprsFiltered)

# check if GOI are still in
all(GOI$EntrezID %in% rownames(ExprsFiltered))


str(ExprsFiltered)
write.table(ExprsFiltered,
            file = file.path(resultsdir, 
                             paste(scriptName,
                                   "_countmatrix_only_filteredmapped.txt",
                                   sep = "")),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


#----------------------------------------------
# NORMALIZATION
#----------------------------------------------
# prepare data input
Expressionset <- ExprsFiltered

#load table into DGEList object (rows = genes, col=samples, libsize, feature annotation)
dge <- DGEList(counts=Expressionset)
dim(dge)

# normalization after filtering with TMM normalization method (scale normalization)
dge <- edgeR::calcNormFactors(dge)

#save normalized count matrix into RDS file
saveRDS(dge, 
        file.path(resultsdir,
                  paste(scriptName,
                        class(dge), 
                        "RDS", 
                        sep=".")))

# choose which normalization method
minlib <- min(dge$samples$lib.size)
maxlib <- max(dge$samples$lib.size) 
checklibsize <- maxlib/minlib
print(paste("Min Lib Size:", minlib, sep = ""))
print(paste("Max Lib Size:", maxlib, sep = ""))
print(paste("max/min = ", checklibsize, sep = ""))

#----------------------------------------------
# 1. limma-trend
#----------------------------------------------
# criteria: max log fold change is <= 3 (libsize)
# sequencing depth is consistent across the RNA samples

# convert counts into logCPM (log2)
logCPM <- cpm(dge, 
              log=TRUE, 
              prior.count=3)

dim(logCPM)

# check if GOI are still in
all(GOI$EntrezID %in% rownames(logCPM))

#save log-transformed (limma-trend) into matrix file
saveRDS(logCPM, 
        file=file.path(resultsdir,
                       paste(scriptName,
                             class(logCPM), 
                             "LogCPM", 
                             "RDS", 
                             sep=".")))

#----------------------------------------------
# 2. voom
#----------------------------------------------
# critera: max log fold change is > 3 (libsize)
# voom transformation is applied to normalized & filtered DGEList object
v <- voom(dge,
          design=NULL,
          plot=TRUE)

dim(v)

#extract count matrix from v
voom <- ExtractCounts(v)
dim(voom)

# check if GOI are still in
all(GOI$EntrezID %in% rownames(voom))

# save log-transformed (voom) into RDS file
saveRDS(v,
        file=file.path(resultsdir, 
                       paste(scriptName,
                             class(v),
                             "Voom",
                             "RDS",
                             sep=".")))

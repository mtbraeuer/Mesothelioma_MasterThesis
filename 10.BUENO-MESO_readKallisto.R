## ------------------------------------------------------
## Project: BUENO-MESO
## Script name: 10.BUENO-MESO_readKallisto.R
##
## Purpose of script: read kallisto abundance.tsv files into expression matrix and then Eset
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
##
## ------------------------------------------------------
# libibary loading
source("Initialize_BUENO-MESO.R")
library(vsn)
library(tximport)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(AnnotationDbi)
library(DESeq2)
library(Biobase)

#--------------------------------------------------
# scriptName
scriptName <- "10.BUENO-MESO_readKallisto"

#--------------------------------------------------
# data set loading
pData <- read.delim(file.path(datadir,
                              "BUENO-MESO_pData.txt"),
                    sep = "\t",
                    stringsAsFactors = FALSE)
rownames(pData) <- pData$RNATumor_ID

################################################################################################
# read kallisto abundance files into expr matrix
################################################################################################
# get vector of files with "abundance.tsv" in data dir
files <- dir(datadir,
             pattern = "abundance.tsv",
             recursive = FALSE)
# get rnatumor_id of abundance.tsv files
samples <- sub(pattern = "_abundance.tsv", replacement = "", x = files)
# set file.path for each file
files <- file.path(datadir, files)
# set names
names(files) <- samples

# load mapping data base into txdb
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# get ENST IDs of txdb
k <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
# match ENST ID to Entrez ID
tx2gene <- select(txdb, k, "GENEID", "TXNAME")
# remove version of ENST gene id
# e.g. ENST00000456328.2 -> ENST00000456328
tx2gene$TXNAME <- substr(tx2gene$TXNAME, 1, 15)
# remove all entries which do not have an EntrezID
tx2gene <- tx2gene[!is.na(tx2gene$GENEID),]

# read all abundance.tsv into a further usable transcript-level matrix
txi.kallisto.tsv <- tximport::tximport(files, 
                                       type = "kallisto", 
                                       tx2gene = tx2gene, 
                                       ignoreAfterBar = FALSE,
                                       existenceOptional = FALSE)

# create annotation/fData with customized function
fdata <- CreateAnnotation("org.Hs.eg", 
                          keytype = "ENTREZID")
class(fdata)
colnames(fdata)
dim(fdata)
# gnerate the DESeqObject
ddsTxi <- DESeq2::DESeqDataSetFromTximport(txi.kallisto.tsv,
                                           colData = pData,
                                           design = ~1)


# filter for low reads
keep <- rowSums(counts(ddsTxi)) >= 10
dds <- ddsTxi[keep,]
dds <- dds[intersect(rownames(dds), rownames(fdata)),]
rowData(dds) <- fdata[rownames(dds),]

countMatrix <- counts(dds)
dim(countMatrix)


#----------------------------------------------
# save DESeq2 data object
saveRDS(object = dds,
        file.path(resultsdir,
                  paste(scriptName,
                        "_",
                        class(dds),
                        ".RDS",
                        sep = "")))


# check if count matrix has numeric entries (NOT integers!!!)
class(countMatrix[1,1])
countMatrix <- as.matrix(sapply(as.data.frame(countMatrix), as.numeric))
# add rownames
rownames(countMatrix) <- rownames(counts(dds))
str(countMatrix)
dim(countMatrix)


#----------------------------------------------
# save count matrix
write.table(x = countMatrix,
            file.path(resultsdir,
                      paste(scriptName,
                            "_rawCountMatrix",
                            ".txt",
                            sep = "")),
            sep = "\t",
            quote = FALSE,
            col.names = TRUE,
            row.names = TRUE)


################################################################################################
# create expression set
################################################################################################
# Variance stabilizing transformation
# https://rdrr.io/bioc/DESeq2/man/varianceStabilizingTransformation.html
vst_countMatrix <- DESeq2::vst(object = countMatrix, blind = FALSE)

# prior vst
pdf(file=file.path(plotsdir, 
                   paste(scriptName,
                         "_VarianceStabilizingTransformation_meanSDPlot_prior",
                         ".pdf",
                         sep="")),
    onefile = FALSE)
vsn::meanSdPlot(countMatrix)
dev.off()

# post vst
pdf(file=file.path(plotsdir, 
                   paste(scriptName,
                         "_VarianceStabilizingTransformation_meanSDPlot_post",
                         ".pdf",
                         sep="")),
    onefile = FALSE)
vsn::meanSdPlot(vst_countMatrix) 
dev.off()

class(vst_countMatrix[1,1])
exprs <- vst_countMatrix
#------------------------------------------------
# edit fData
fData <- fdata
colnames(fData)
fData <- fData[,-c(1, 2, 5, 6, 9)]

# align fData to eset
keep <- match(rownames(exprs), fData$ENTREZID)
fData <- fData[keep,]

# check again for NAs
dim(fData)

table(is.na(fData$ENTREZID))
table(is.na(rownames(exprs)))

rownames(fData) <- fData$ENTREZID


# check if align was correct
# both need to return TRUE
all(nrow(fData)==nrow(exprs))
all(rownames(fData)==rownames(exprs))


# pData
keep <- match(colnames(exprs), rownames(pData))
pData <- pData[keep,]

# check if number of rows of pData == numbers of expression set
# AND names have to be identical
# both need to return TRUE
all(rownames(pData)==colnames(exprs))
all(nrow(pData)==ncol(exprs))


#----------------------------------------------
# Expressionset creation
#----------------------------------------------
Eset <- ExpressionSet(assayData=exprs,
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

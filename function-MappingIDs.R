## ------------------------------------------------------
## Project: Mesothelioma Masterthesis
## Script name: function-MappingIDs.R
##
## Purpose of script: MappingIDsToEntrez
## 
##
## Author: Marie-Theres Bräuer
##
## Version: Jan 2019
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
## input1: raw count matrix (rows = genes, col = samples)
## input2: annotation file eg. "org.Hs.eg.db"
## if in doubt check the db out!
## input3: map_origin "from ID" eg ENSEMBL (uppercase!!!) 
## input4: map_target "to ID" eg ENTREZID (uppercase!!!)
## output: filtered matrix for requested ID
##
## ------------------------------------------------------
# LIBRARY loading
library(limma)
library(Biobase)
library(AnnotationDbi)

#----------------------------------------------
# function 
#----------------------------------------------
MappingIDs <- function(inmatrix, annotation, map_origin, map_target)
{
    # prepare input for analysis
    countmatrix <- as.data.frame(inmatrix)
    dim(countmatrix)
    
    # get the annotation file
    package <- paste(annotation, "db", sep = ".")
    library(package, character.only = TRUE)
    DBObject <- get(package)
    
    # mapping step
    # return data frame 
    ProbeID <- as.data.frame(mapIds(DBObject, 
                                    keys = rownames(countmatrix),
                                    column = map_target,
                                    keytype = map_origin))
    # rename 1 col
    colnames(ProbeID)[1] <- "ID"
    
    # add column with ID at end of inmatrix
    countmatrix$NEWID <- ProbeID$ID
    
    # remove all rows with NAs
    print(paste("Matrix with NAs ", nrow(countmatrix), sep = ":"))
    countmatrix <- na.omit(countmatrix)
    print(paste("Matrix without NAs ", nrow(countmatrix), sep = ":"))
    
    #create vector with new ids for later
    newids <- countmatrix$NEWID
    
    # remove the last col "NEWID"
    countmatrix <- as.matrix(countmatrix[, 1:length(countmatrix)-1])
    
    # rename rows with new ID
    rownames(countmatrix) <- newids
    
    # summarize multiple IDs which point to 1 NEWID with average --> avereps
    countmatrix <- avereps(countmatrix,
                           ID=rownames(countmatrix))
    print(paste("Matrix after avereps ", nrow(countmatrix), sep = ":"))
    
    # return mapped matrix
    countmatrix <- as.matrix(countmatrix)
    return(countmatrix)
}
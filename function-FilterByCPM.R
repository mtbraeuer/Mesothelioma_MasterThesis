## ------------------------------------------------------
## Project: Mesothelioma Masterthesis
## Script name: function-FilterByCPM.R
##
## Purpose of script: FilterByCPM
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
## input2: cutoff for filtering eg 0.10, 0.15, ...
## input3: threshold --> min count per CPM eg  0.10, 0.20, ...
## output: filtered matrix
## calculation steps
## 1 create cpm matrix 
## 2 create boolean matrix with threshold
## 3 create boolean matrix with cutoff
## 4 filter matrix with boolean vector and return result

## function call: x <- FilterByCPM(matrix, cutoff, threshold)
##
## ------------------------------------------------------
# LIBRARY loading
library(limma)
library(edgeR)
library(DESeq2)
library(Biobase)

#----------------------------------------------
# function 
#----------------------------------------------

FilterByCPM <- function(inmatrix, cutoff, threshold)
{
    # create cpm matrix
    CPM <- cpm(inmatrix,
               normalized.lib.sizes = TRUE,
               log = FALSE,
               prior.count=2)
    
    # create boolean matrix
    bm <- matrix(data = NA, nrow = nrow(CPM), ncol = ncol(CPM))
    # i = row
    for (i in 1:nrow(CPM)) {
        
        # j = col
        for (j in 1:ncol(CPM)) {
            if (CPM[i,j]>=threshold) {
                bm[i,j] <- 1
                
            }
            else{
                bm[i,j] <- 0
            }
        }
    }
   
    # calculate boolean vector 0/1 --> "(rowsum/n)>cutoff"
    bool <- vector(mode = "logical", length = nrow(CPM))
    n <- ncol(CPM)
    for (i in 1:nrow(bm)) {
        
        entry <- (sum(bm[i,])/n)
        if (entry>=cutoff) {
            bool[i] <- TRUE
            
        }
        else{
            bool[i] <- FALSE
        }
    }
    
    # filter matrix with boolean vector
    outmatrix <- subset(inmatrix, bool)
    
    return(outmatrix)

}

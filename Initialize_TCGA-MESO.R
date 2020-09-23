## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: Initialize_TCGA-MESO.R
##
## Purpose of script: Initialize directories for input, output, scripts
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 16. October 2019
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
##
## ------------------------------------------------------
print("Initializing default directories...")
directories <- c("data", 
                 "results", 
                 "plots", 
                 "scripts",
                 "systeminfo")
names(directories) <- directories
for (i in directories){
  assign(paste(i,"dir", sep = ""), i)
  if (!file.exists(i)){
    dir.create(i) 
  }
}

#--------------------------------------------------
# Load needed libaries for project with pacman
#--------------------------------------------------
# pacman checks if the needed package is installed
# if not it will download and install it from CRAN and Bioconductor
print("Loading libraries...")
if (!require("pacman", character.only = TRUE)){
  install.packages("pacman")
} 

pacman::p_load("WGCNA",
               "data.table",
               "dplyr",
               "edgeR",
               "genefilter",
               "ggfortify",
               "genefilter",
               "igraph",
               "limma",
               "magrittr",
               "maxstat",
               "org.Hs.eg.db",
               "rms",
               "SummarizedExperiment",
               "survminer",
               "survival",
               "stringr")
# some packages can not be downloaded and installed with pacman
# use the old way
bioc <- c("DESeq2",
          "TCGAbiolinks",
          "AnnotationDbi",
          "Biobase",
          "car",
          "ComplexHeatmap")
for(p in bioc){
  if(!require(p, character.only = TRUE)){
    BiocManager::install(p)
  }
  library(p, character.only = TRUE)

}

# check for loaded libs
print("List of loaded packages...")
print((.packages()))

#--------------------------------------------------
# Load customized functions for project in functionsdir
#--------------------------------------------------
print("Customized functions were loaded...")
functionsdir = "../CustomizedFunctions/"
files <- dir(functionsdir, pattern = "function", full.names = TRUE)
for (file in files){
  source(file)
}

#--------------------------------------------------
# Initialize project-wide parameter
#--------------------------------------------------
project <- "TCGA-MESO"

#--------------------------------------------------
# Clean-up unnecessary suff
#--------------------------------------------------
rm(directories, i, bioc, p, files, file)

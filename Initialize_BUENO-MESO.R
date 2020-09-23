## ------------------------------------------------------
## Project: BUENO-MESO
## Script name: Initialize_BUENO-MESO.R
##
## Purpose of script: set and define global parameter for project
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
## data downloaded from european genome archive via thomas server
## 211 MPM patients RNA-Seq(EGA: dataset 15, ~1.9TB) 
## Rafael Bueno et al. 2016, Comprehensive genomic analysis of malignant pleural mesothelioma 
## identifies recurrent mutations, gene fusions and splicing alterations
##
## ------------------------------------------------------

#--------------------------------------------------
# global settings
options(stringsAsFactors = FALSE)

#--------------------------------------------------
# Initialize directories for input, output, scripts
#--------------------------------------------------
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
project <- "BUENO-MESO"

GOI <- read.delim(file.path("../GOI.txt"),
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

#--------------------------------------------------
# Initialize project-wide plotting parameters
#--------------------------------------------------
DEG_volcanoplot_x_axis <- c(-10,10)
DEG_volcanoplot_y_axis <- c(0,45)

#--------------------------------------------------
# Clean-up unnecessary suff
#--------------------------------------------------
rm(directories, i, files, file)
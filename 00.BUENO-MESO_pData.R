## ------------------------------------------------------
## Project: BUENO-MESO
## Script name: 00.BUENO-MESO_pData.R
##
## Purpose of script: create pData of BUENO data set
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 24. January 2020
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
##
## ------------------------------------------------------
## libibary loading
source("Initialize_BUENO-MESO.R")

scriptName <- "00.BUENO-MESO_pData"

## ------------------------------------------------------
# data set loading

# pData
clinicalInfo <- read.csv(file.path(datadir, 
                                   "BUENO-MESO.ClinicalInfo.csv"),
                         stringsAsFactors = FALSE,
                         header = TRUE)
pData <- clinicalInfo

# egar id to identify the RNASeq data for each patient
egar <- read.table(file.path(datadir,
                             "BUENO-MESO-EGA_Dataset15_RNASeq.txt"),
                   header = TRUE,
                   stringsAsFactors = FALSE)


################################################################################################
# create pData
################################################################################################
# set coolumn names
# set these columns as written here: PATIENTID = RNATumor_ID, OS_days, VitalStatus, Sex !!!!!
colnames(pData)
# remove first row
pData <- pData[-c(1),]
colnames(pData) <- c("Number",                                                                              
                     "Tumor_ID",
                     "MatchedNormal_ID",
                     "DNATumor_ID",
                     "DNAMatchedNormal_ID",
                     "RNATumor_ID",
                     "PrimaryTissue",
                     "TissueDiagnosis",
                     "Histology",
                     "Histology_Bueno_expr_based",
                     "Exome",
                     "Targeted_exome_gene_panel_SPET",
                     "RNA_seq",
                     "X2_5M_Array",
                     "WGS_tumor_and_normal_sequenced",
                     "Percent_Epithelioid",
                     "Percent_Sarcomatoid",
                     "Sex",
                     "Smoking_History",
                     "AsbestosExposure_Study",
                     "Fibers_per_g_lung",
                     "VitalStatus",
                     "OS_years",
                     "Age_at_surgery_years",
                     "Stage",
                     "Type_of_Surgery",
                     "preop_Treatment_1_treated_0_untreated",
                     "Type_of_pre_op_Treatment",
                     "FISH_chrom3",
                     "FISH_chrom22",
                     "NF2_FISH",
                     "NF2_status_2_5M_array",
                     "NF2_fusion", 
                     "PD_L1_expression_RPKM")

colnames(pData)
dim(pData)

# reduce pData to
#RNATumor_ID, Sex, VitalStatus, OS_years,
#Histology, Histology_Bueno_expr_based, Age_at_surgery_years, AsbestosExposure_Study,  
#Percent_Epithelioid, Percent_Sarcomatoid, Tumor_ID, MatchedNormal_ID, DNATumor_ID, DNAMatchedNormal_ID 
pData <- pData[c(6, 18, 22, 23, 9, 10, 24, 20, 16, 17, 2, 3, 4, 5)]
colnames(pData)
# set rownames for deleting empty entries
rownames(pData) <- c(1:nrow(pData))

# there are some fake entries with no RNATumor_ID --> delete these rows
# 211 patients in total
rownames(pData[pData$RNATumor_ID=="",])
pData <- pData[-as.vector(as.numeric(rownames(pData[pData$RNATumor_ID=="",]))),]
dim(pData)
rownames(pData) <- c(1:nrow(pData))



# get EGAR00001406XXX ID for each patient
pData["EGAR_ID"] <- NA
egar["RNATumor_ID"] <- NA

# get RNATumor_ID for each corresponding egar ID
for (i in 1:nrow(egar) ) {
    x <- stringr::str_match(string = egar$File_name[i],
                            pattern = "(.*?)/.*/(.*?)_") # EGAR00001406\d*)\/.*\/(\d*)_.*
    egar$EGAR_ID[i] <- x[2]
    egar$RNATumor_ID[i] <- x[3]
}

# match EGAR ID and RNATumor ID in pData
keep <- match(pData$RNATumor_ID, egar$RNATumor_ID)
pData$EGAR_ID <- egar$EGAR_ID[keep]

# correct pData entries in some columns
# check which colmns need adjustment
colnames(pData)

# ok
table(pData$Sex)

# not ok
table(pData$VitalStatus)
pData$VitalStatus[pData$VitalStatus=="a"] <- 0
pData$VitalStatus[pData$VitalStatus=="d"] <- 1
table(pData$VitalStatus)

# make, numeric & add another col with days, multiply by 365.25
table(pData$OS_years)
pData$OS_years <- as.numeric(pData$OS_years)
pData["OS_days"] <- round(pData$OS_years*365.25, digits = 1)

# ok
table(pData$Histology)
table(pData$Histology_Bueno_expr_based)

# 
table(pData$Age_at_surgery_years)
hist(as.numeric(pData$Age_at_surgery_years), labels = T)

# convertable YES, NO, NA
table(pData$AsbestosExposure_Study)
pData["AsbestosExposure"] <- NA

pData$AsbestosExposure[pData$AsbestosExposure_Study=="39 years known"] <- "Y"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="none known"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="none known, but was in Navy for 2 years"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="possible"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="questionable"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="yes"] <- "Y"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="none"] <- "N"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="none known, but husband used to install asbestos pipe"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="none mentioned"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="probable"] <- "NA"
pData$AsbestosExposure[pData$AsbestosExposure_Study=="unknown"] <- "NA"

table(pData$AsbestosExposure)

# save pData
write.table(pData,
            file.path(datadir,
                      "BUENO-MESO_pData.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

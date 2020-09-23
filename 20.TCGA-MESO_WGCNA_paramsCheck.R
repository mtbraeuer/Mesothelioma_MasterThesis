## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 20.TCGA-MESO_WGCNA_paramsCheck.R
##
## Purpose of script: check for network parameters 
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 13. November 2019
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
## Weighted Gene Co-expression Network Analysis
## https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
## Langfelder & Horvath et al
##
## ------------------------------------------------------
source("Initialize_TCGA-MESO.R")
#--------------------------------------------------
scriptName <- "20.TCGA-MESO_WGCNA_paramsCheck"

library(WGCNA)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()

#----------------------------------------------
# DATASET loading 
#----------------------------------------------
# Expressionset
Eset <- readRDS(file = file.path(resultsdir,
                                 "15.TCGA-MESO_Stratification.ExpressionSet.RDS"))

EsetAnalysis <- Eset

EsetAnalysis

#----------------------------------------------
# PARAMETER SETTINGS
#----------------------------------------------
# GeneralParameters 
GeneralParameters <- list(
    # number of genes
    NumberOfGenes <- nrow(exprs(EsetAnalysis)), 
    # network construction
    Simple <- TRUE,
    # col names of trait data file
    # cluster tree method
    ClusterTreeMethod <- "average",
    CutHeightClusterTree <- 15
)

# Network Parameter Setting
NetworkParameters <- list(
    verbose = 3,
    # options for splitting data into blocks
    maxBlockSize = 25000,
    # adjacency function options
    networkType = "signed",
    # TOM
    label = "singleWGCNA",
    TOMType = "signed",
    saveTOMs = FALSE,
    saveTOMFileBase = file.path(resultsdir,paste(scriptName, "blockwiseTOM", sep=".")),
    # tree cut options
    minModuleSize = 30,
    # gene reassignment, module trimming, and module "significance" criteria
    reassignThreshold = 1e-6 ,
    # module merging
    numericLabels = TRUE,
    pamStage = FALSE,
    pamRespectsDendro = TRUE,
    #correlation type
    corType = "pearson"
)

#----------------------------------------------
# WGCNA 
#----------------------------------------------
# additional parameters
#-----------------------
GeneNames <- as.character(EsetAnalysis@featureData[["EntrezID"]])
DatExpr <- t(exprs(EsetAnalysis))
NumberOfSamples <- nrow(pData(EsetAnalysis))
NumberOfGenes <- ncol(DatExpr)

gsg = goodSamplesGenes(datExpr = DatExpr, verbose = 3);
gsg$allOK

#----------------------
# Network Construction
#----------------------
# calculation beta
powers <- c(seq(1,20,by=1))

sft <- pickSoftThreshold(DatExpr,
                         # RSquaredCut: desired minimum scale free topology 
                         # fitting index R^2 #(*** un/comment)
                         RsquaredCut = 0.80,
                         powerVector=powers,
                         networkType = NetworkParameters$networkType,
                         verbose = 6)

beta <- sft$powerEstimate
if (is.na(beta)){
    if (NetworkParameters$networkType == "unsigned"){
        beta <- 6
    } else {
        beta <- 12
    }
}
# append beta and esetAnalysis to networkParameter list
NetworkParameters$power <- beta
NetworkParameters$datExpr <- DatExpr

names(NetworkParameters)

collectGarbage()


#-------------
# for loop over multiple parameters for network optimization
#-------------
counter <- 1

# deppsplit 
for (deepsplit in c(2, 3)) {
    
    # merge cut height 
    for (mergecutheight in c(0.15, 0.20, 0.25, 0.30)) {
        
        # detectcutheight ,
        for (detectcutheight in c(0.990, 0.995, 0.998)) {
            
            label <- paste(paste("d", deepsplit, sep = ":"),
                           paste("mch", mergecutheight, sep = ":"),
                           paste("dch", detectcutheight, sep = ":"),
                           sep = "_")
            
            NetworkParameters$deepSplit <- deepsplit
            NetworkParameters$mergeCutHeight <- mergecutheight
            NetworkParameters$detectCutHeight <- detectcutheight
            
            # Network Construction
            Net <- do.call("blockwiseModules", NetworkParameters)
            
            Net$NetworkParameters <- as.list(args(blockwiseModules))
            Net$NetworkParameters[names(NetworkParameters)] <- NetworkParameters
            Net$eset <- EsetAnalysis
            collectGarbage()
            
            # define some global parameters for further analysis
            ModuleColors <- labels2colors(Net$colors)
            ModuleLabels <- Net$colors
            MEs <- Net$MEs
            GeneTree <- Net$dendrograms[[1]]
            ModNames <- substring(names(MEs), 3)
            ModuleNames <- list(substring(names(MEs),3))
            
            # Number of modules with number of genes
            # module "0" is for genes outside of all modules
            table(labels2colors(Net$colors))
            
            #----------------------------------------------
            # Cluster Dendrogram 
            #----------------------------------------------
            pdf(file=file.path(plotsdir,
                               paste(scriptName,
                                     counter,
                                     "Dendrogram", 
                                     "pdf", 
                                     sep = ".")),
                paper="a4",
                onefile = FALSE)
            
            plotDendroAndColors(GeneTree,
                                ModuleColors[Net$blockGenes[[1]]],
                                main = paste(counter,
                                             "-",
                                             scriptName,
                                             "\n",
                                             "Beta:",
                                             NetworkParameters$power,
                                             "- MergeCutHeight:",
                                             NetworkParameters$mergeCutHeight,
                                             "- DeepSplit:",
                                             NetworkParameters$deepSplit,
                                             "\n",
                                             "- DetectCutHeight",
                                             NetworkParameters$detectCutHeight),
                                dendroLabels = FALSE,
                                hang = 0.02,
                                addGuide = TRUE, 
                                guideHang = 0.05)
            dev.off()
            collectGarbage()

            counter = counter + 1
        }
    }
}

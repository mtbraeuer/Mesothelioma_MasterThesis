## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 25.TCGA-MESO_WGCNA.R
##
## Purpose of script: WGCNA
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
## Weighted Gene Co-expression Network Analysis
## https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
## Langfelder & Horvath et al
##
## ------------------------------------------------------
source("Initialize_TCGA-MESO.R")
#--------------------------------------------------
scriptName <- "25.TCGA-MESO_WGCNA"

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

GOI <- read.delim(file.path("../GOI.txt"),
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

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
    corType = "pearson",
    # module merging
    mergeCutHeight = 0.25,
    #tree cut option
    deepSplit = 2,
    detectCutHeight = 0.995
)

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

# Network Construction
# set cor function to WGCNA::cor function
cor <- WGCNA::cor
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
write.table(table(labels2colors(Net$colors)), 
            sep = "\t",
            file.path(resultsdir, 
                      paste(scriptName,
                            "NumberofGenesPerModule",
                            ".txt",
                            sep = "")))

saveRDS(Net,
        file = file.path(resultsdir,
                         paste(scriptName,
                               "Network",
                               ".RDS",
                               sep = "")))

save.image(file.path(resultsdir,
                     paste(scriptName,
                           "_workingImage",
                           ".RData",
                           sep = "")))

########################################################################################
# Plots
#----------------------------------------------
# Cluster Dendrogram 
#----------------------------------------------
pdf(file=file.path(plotsdir,
                   paste(scriptName,
                         "_Dendrogram", 
                         ".pdf", 
                         sep = "")),
    #paper="a4",
    onefile = FALSE,
    height = 8,
    width = 6)

plotDendroAndColors(GeneTree,
                    ModuleColors[Net$blockGenes[[1]]],
                    main = paste("Single WGCNA MPM-TCGA\n",
                                 "Cluster Dendrogram",
                                 # "NWTop:",
                                 # NetworkParameters$networkType,
                                 # "#Genes:",
                                 # NumberOfGenes,
                                 # "\n",
                                 # "Beta:",
                                 # NetworkParameters$power,
                                 # "MergeCutHeight:",
                                 # NetworkParameters$mergeCutHeight,
                                 # "DeepSplit:",
                                 # NetworkParameters$deepSplit,
                                 # "\n",
                                 # "DetectCutHeight",
                                 # NetworkParameters$detectCutHeight,
                                 # "Correlation:",
                                 # NetworkParameters$corType,
                                 sep = " "),
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = TRUE, 
                    guideHang = 0.05)
dev.off()
collectGarbage()


# traits of interest
pData <- pData(EsetAnalysis)
colnames(pData)

#----------------------------------------------
# Sample Dendogram & Trait Heatmap of general parameter
#----------------------------------------------
traits <- data.matrix(subset(pData(EsetAnalysis), select = c(4)))

pdf(file=file.path(plotsdir,
                   paste(scriptName,
                         "_SampleDendogram_TraitHeatmap_",
                         ".pdf", 
                         sep = "")),
    paper = "a4r",
    onefile = FALSE,
    width = 12)

traitcolors <- numbers2colors(traits,
                              signed = FALSE)

sampletree <- hclust(dist(DatExpr), 
                     method = ClusterTreeMethod)

plotDendroAndColors(sampletree, 
                    traitcolors, 
                    groupLabels = colnames(traits), 
                    main = paste(scriptName, 
                                 ":", 
                                 "Sample Dendrogram and Trait Heatmap"))
dev.off()

collectGarbage()


#----------------------------------------------
# Scale Independence Plot
#----------------------------------------------
pdf(file=file.path(plotsdir,
                   paste(scriptName, 
                         "_ScaleIndependence", 
                         ".pdf", 
                         sep = "")),
    paper = "a4",
    onefile = FALSE)

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n",
     main = paste(scriptName, "Scale independence", sep = " "))

text(sft$fitIndices[,1], 
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels = powers,
     cex = 0.9,
     col = "red") 

dev.off()

#----------------------------------------------
# Mean Connectivity Plot
#----------------------------------------------
pdf(file=file.path(plotsdir,
                   paste(scriptName, 
                         "_MeanConnectivity", 
                         ".pdf", 
                         sep = "")),
    paper = "a4",
    onefile = FALSE)


plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste(scriptName, "Mean connenctivity", sep = " "))

text(sft$fitIndices[,1], 
     sft$fitIndices[,5],
     labels = powers,
     cex = 0.9,
     col = "red")

dev.off()


########################################################################################
# Plots MOI
#----------------------------------------------
# Mutations Of Interest - mutation events
#----------------------------------------------
# define threshold for module-trait relationship
posthresholdModuleTrait <- 0.4
negthresholdModuleTrait <- -0.4

MEs0 = moduleEigengenes(DatExpr, ModuleColors)$eigengenes
MEs = orderMEs(MEs0)

# Mutsubset of pData
colnames(pData)
MOITraits <- data.matrix(subset(pData, select = c(7:31)))
colnames(MOITraits)
colnames(MOITraits) <- gsub("(.*?\\_)(.*)", "\\2", colnames(MOITraits))


# create correlation matrix for all mutTraits 
ModuleTraitCor <- cor(MEs, as.matrix(MOITraits), use = "p")
ModuleTraitPvalue = corPvalueStudent(ModuleTraitCor, NumberOfSamples)

# plot correlation of values within heatmap plot
# function
pdf(file=file.path(paste(plotsdir,
                         paste(scriptName,
                               "_CorrelationHeatmap_MOI",
                               ".pdf",
                               sep = ""),
                         sep = "/")),
    #paper="a4",
    onefile = FALSE,
    height = 10)


textMatrix <- paste(signif(ModuleTraitCor, 2), "\n(",
                    signif(ModuleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(ModuleTraitCor)

labeledHeatmap(Matrix = ModuleTraitCor,
               xLabels = colnames(MOITraits),
               yColorLabels = T,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               cex.lab.y = 0.4,
               zlim = c(-1,1),
               main = paste(scriptName, 
                            "CorrelationHeatmap - MOI"))
dev.off()


########################################################################################
# Plots GOI
#----------------------------------------------
# Genes Of Interest - expression
#----------------------------------------------
dim(DatExpr)

# not in DatExprs inculuded
table(GOI$EntrezID %in% colnames(DatExpr))

# data frame with GOI Traits
keep <- match(as.character(GOI$EntrezID),colnames(DatExpr))

length(keep)
keep
# subset of only GOI expression matrix
GOITraits <- DatExpr[,keep]
dim(GOITraits)
class(GOITraits)

# change col names for biologists
for (i in 1:ncol(GOITraits)) {
    for (j in 1:nrow(GOI)) {
        if (colnames(GOITraits)[i]==as.character(GOI$EntrezID)[j]) {
            
            colnames(GOITraits)[i] <- paste(colnames(GOITraits)[i], 
                                            GOI$GeneSymbol[j],
                                            sep = "_")
        }
    }
}
colnames(GOITraits)

# recalculate MEs with color labels
ModuleTraitCor = cor(MEs, as.matrix(GOITraits), use = "p")
ModuleTraitPvalue = corPvalueStudent(ModuleTraitCor, NumberOfSamples)


# plot correlation of values within heatmap plot
# function
pdf(file=file.path(plotsdir,
                   paste(scriptName,
                         "_CorrelationHeatmap_GOIexpr",
                         ".pdf", 
                         sep = "")),
    #paper="a4",
    onefile = FALSE,
    height = 10)


textMatrix <- paste(signif(ModuleTraitCor, 2), "\n(",
                    signif(ModuleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(ModuleTraitCor)

labeledHeatmap(Matrix = ModuleTraitCor,
               xLabels = colnames(GOITraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               cex.lab.y = 0.4,
               zlim = c(-1,1),
               main = paste(scriptName,
                            "Correlation Heatmap GOI expr", 
                            sep = " "))
dev.off()
collectGarbage()

########################################################################################
# Plots GOI
#----------------------------------------------
# stratification and overall survival (Traits of interest)
#----------------------------------------------
colnames(pData)
TOI <- pData[,c(3, 32:54)]
colnames(TOI)
class(TOI)
str(TOI)

# factor into numeric
for (i in 2:ncol(TOI)) {
    TOI[,i] <- ifelse(TOI[,i]=="high", 1, 0)
    # high =1, low = 0
    
    print(table(TOI[,i]))
    
}
# recalculate MEs with color labels
ModuleTraitCor = cor(MEs, as.matrix(TOI), use = "p")
ModuleTraitPvalue = corPvalueStudent(ModuleTraitCor, NumberOfSamples)

# plot correlation of values within heatmap plot
# function
pdf(file=file.path(plotsdir,
                   paste(scriptName,
                         "_CorrelationHeatmap_stratification",
                         ".pdf", 
                         sep = "")),
    #paper="a4",
    onefile = FALSE,
    height = 10)


textMatrix <- paste(signif(ModuleTraitCor, 2), "\n(",
                    signif(ModuleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(ModuleTraitCor)

labeledHeatmap(Matrix = ModuleTraitCor,
               xLabels = colnames(TOI),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               cex.lab.y = 0.4,
               zlim = c(-1,1),
               main = paste(scriptName,
                            "Correlation Heatmap stratification", 
                            sep = " "))
dev.off()
collectGarbage()

########################################################################################
# Plots GOI and strat
#----------------------------------------------
# stratification and overall survival (Traits of interest)
#----------------------------------------------
colnames(pData)
colnames(GOITraits)
# TOI <- pData[,c(3, 35, 42, 48, 50, 53, 54)]
# colnames(TOI)
# class(TOI)
# str(TOI)
# 
# # factor into numeric
# for (i in 2:ncol(TOI)) {
#     TOI[,i] <- ifelse(TOI[,i]=="high", 1, 0) 
#     # high =1, low = 0
#     print(table(TOI[,i]))
#     
# }

TOI <- cbind(pData$OS_days, GOITraits[,c(4,11,17,19,22,23)])
str(TOI)
colnames(TOI)

colnames(TOI) <- c("OS_days",
                   "CDKN2A_expr",
                   "IL6R_expr",
                   "NF2_expr",
                   "TERT_expr",
                   "TPH1_expr",
                   "YBX1_expr")


# recalculate MEs with color labels
ModuleTraitCor = cor(MEs, as.matrix(TOI), use = "p")
ModuleTraitPvalue = corPvalueStudent(ModuleTraitCor, NumberOfSamples)

# plot correlation of values within heatmap plot
# function
pdf(file=file.path(plotsdir,
                   paste(scriptName,
                         "_CorrelationHeatmap_GOI_start_expr",
                         ".pdf", 
                         sep = "")),
    #paper="a4",
    onefile = FALSE,
    height = 8,
    width = 6)


textMatrix <- paste(signif(ModuleTraitCor, 2), "\n(",
                    signif(ModuleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) <- dim(ModuleTraitCor)

labeledHeatmap(Matrix = ModuleTraitCor,
               xLabels = colnames(TOI),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.4,
               cex.lab.y = 0.4,
               zlim = c(-1,1),
               main = "Single WGCNA MPM-TCGA Correlation Heatmap")
dev.off()
collectGarbage()


#----------------------------------------------
# Get Module Membership of GOI
#----------------------------------------------
# module colors
colors <- unique(ModuleColors)

# list of variable names of each module where the entrezids are saved in
ListModules <- c()

for (i in 1:length(colors)) {
    
    varname <- paste("module", colors[i], sep = "_")
    ListModules <- c(ListModules, varname)
    x <- colnames(DatExpr)[ModuleColors==colors[i]]
    assign(varname, x)
}

# new column in GOI
GOI["ModuleMembership"] <- NA

# find module membership of GOI

# row-wise
for (i in 1:nrow(GOI)) {
    
    # go through length of vector
    for (j in 1:length(ListModules)) {
        
        module <- eval(parse(text=ListModules[j]))
        
        for (x in 1:length(module)) {
            
            if(as.character(GOI$EntrezID[i])==module[x]){ 
                GOI$ModuleMembership[i] <- str_replace(ListModules[j], "module_", "")
                #print(GOI$EntrezID[i])
                #print(ListModules[j])
            }
        }
    }
}

table(GOI$ModuleMembership)

write.table(GOI,
            file.path(resultsdir,
                      paste(scriptName,
                            "_GOI_ModuleMembership",
                            ".txt",
                            sep = "")),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)




#----------------------------------------------
# Get Modules and write gml files for cytoscape
#----------------------------------------------
vertexAttributes <- fData(EsetAnalysis)
# default bei adjacency corfnc=pearson 
adjacency <- adjacency(datExpr = DatExpr,
                       type = "signed",
                       power = beta)

write.table(adjacency, file = file.path(resultsdir,
                                        paste(scriptName,
                                              "_AdjacencyMatrix.txt",
                                              sep = "")),
            quote = F,
            sep = "\t")

graph <- CreateNetwork(adjacency,
                       colors = labels2colors(Net$colors),
                       weighted = TRUE,
                       vertex.attributes = vertexAttributes)



for (color in labels2colors(unique(Net$colors))){
    subgraph <- induced.subgraph(graph, V(graph)[V(graph)$colors == color])
    write.graph(subgraph, 
                format = "gml",
                file = file.path(resultsdir, paste(scriptName, "_Networks_", color, ".gml", sep = "")))
    print(paste(color, "finished"))
}

#--------------------------------------------------------------------
# save entrezIDs of modules of interest
interestingmodules <- list(module_black, 
                           module_blue,
                           module_brown,
                           module_cyan,
                           module_darkgreen,
                           module_darkgrey,
                           module_darkorange,
                           module_darkred,
                           module_darkturquoise,
                           module_green,
                           module_greenyellow,
                           module_grey,
                           module_grey60,
                           module_lightcyan,
                           module_lightgreen,
                           module_lightyellow,
                           module_magenta,
                           module_midnightblue,
                           module_orange,
                           module_pink,
                           module_purple,
                           module_red,
                           module_royalblue,
                           module_salmon,
                           module_tan,
                           module_turquoise,
                           module_white,
                           module_yellow)

interestingmodulesnames <- list("module_black", 
                                "module_blue",
                                "module_brown",
                                "module_cyan",
                                "module_darkgreen",
                                "module_darkgrey",
                                "module_darkorange",
                                "module_darkred",
                                "module_darkturquoise",
                                "module_green",
                                "module_greenyellow",
                                "module_grey",
                                "module_grey60",
                                "module_lightcyan",
                                "module_lightgreen",
                                "module_lightyellow",
                                "module_magenta",
                                "module_midnightblue",
                                "module_orange",
                                "module_pink",
                                "module_purple",
                                "module_red",
                                "module_royalblue",
                                "module_salmon",
                                "module_tan",
                                "module_turquoise",
                                "module_white",
                                "module_yellow")

for (i in 1:length(interestingmodules)) {
    write.table(x = interestingmodules[i],
                file.path(resultsdir,
                          paste(scriptName,
                                "_",
                                interestingmodulesnames[i],
                                ".txt",
                                sep = "")),
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE
                )

}

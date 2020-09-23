## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 26.TCGA-MESO_singleWGCNA_Analysing_Modules.R
##
## Purpose of script: Analysis of interesting modules of single WGCNA, query of drugability
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 21. September 2020
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
scriptName <- "26.TCGA-MESO_singleWGCNA_Analysing_Modules"
library(ggplot2)
library(clusterProfiler)
library(DESeq2)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)
library(rDGIdb)
library(igraph)
library(plyr)
library(reshape2) 

Eset <- readRDS(file = file.path(resultsdir,
                                 "15.TCGA-MESO_Stratification.ExpressionSet.RDS"))
pData <- pData(Eset)
fData <- fData(Eset)

################################################################################################
# GOTEA of brown module
################################################################################################
brown <- read.delim(file.path(resultsdir,
                              "25.TCGA-MESO_WGCNA_module_brown.txt"))
colnames(brown) <- "EntrezID"

keep <- match(brown$EntrezID, fData$EntrezID)
brown <- cbind(brown, fData[keep,])
brown <- brown[c(-1)]

brown_ego <- enrichGO(gene          = brown$EntrezID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENTREZID',
                      ont           = "BP",
                      pAdjustMethod = "hochberg",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable = TRUE)
head(brown_ego)
head(brown_ego@result$Description, n=20)

saveRDS(brown_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_modulebrown",
                        ".RDS",
                        sep = "")))

# gene concept network
# all edges and nodes
p_genenetwork <- enrichplot::cnetplot(brown_ego, 
                                      foldChange=brown$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-TCGA",
                                      subtitle = "GOTEA of Module Brown")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_brown",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(brown_ego, 
                                      foldChange=brown$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-TCGA",
                                      subtitle = "GOTEA of Module Brown")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_brown",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

################################################################################################
# GOTEA of blue module
################################################################################################
blue <- read.delim(file.path(resultsdir,
                             "25.TCGA-MESO_WGCNA_module_blue.txt"))
colnames(blue) <- "EntrezID"

keep <- match(blue$EntrezID, fData$EntrezID)
blue <- cbind(blue, fData[keep,])
blue <- blue[c(-1)]

blue_ego <- enrichGO(gene          = blue$EntrezID,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = 'ENTREZID',
                     ont           = "BP",
                     pAdjustMethod = "hochberg",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable = TRUE)
head(blue_ego)
head(blue_ego@result$Description, n=20)

saveRDS(blue_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_moduleblue",
                        ".RDS",
                        sep = "")))

# gene concept network
# all edges and nodes
p_genenetwork_blue <- enrichplot::cnetplot(blue_ego, 
                                           foldChange = blue$Symbol,
                                           colorEdge = TRUE,
                                           showCategory = 10,
)

p_genenetwork_blue <- p_genenetwork_blue + labs(title = "Single WGCNA MPM-TCGA",
                                                subtitle = "GOTEA of Module blue")
p_genenetwork_blue

ggsave(plot = p_genenetwork_blue,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_blue",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(blue_ego, 
                                      foldChange=blue$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-TCGA",
                                      subtitle = "GOTEA of Module blue")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_blue",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

################################################################################################
# GOTEA of black module
################################################################################################
black <- read.delim(file.path(resultsdir,
                              "25.TCGA-MESO_WGCNA_module_black.txt"))
colnames(black) <- "EntrezID"

keep <- match(black$EntrezID, fData$EntrezID)
black <- cbind(black, fData[keep,])
black <- black[c(-1)]

black_ego <- enrichGO(gene          = black$EntrezID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = 'ENTREZID',
                      ont           = "BP",
                      pAdjustMethod = "hochberg",
                      pvalueCutoff  = 0.05,
                      qvalueCutoff  = 0.05,
                      readable = TRUE)
head(black_ego)
head(black_ego@result$Description, n=20)

saveRDS(black_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_moduleblack",
                        ".RDS",
                        sep = "")))

# gene concept network
# all edges and nodes
p_genenetwork_black <- enrichplot::cnetplot(black_ego, 
                                            foldChange = black$Symbol,
                                            colorEdge = TRUE,
                                            showCategory = 10,
)

p_genenetwork_black <- p_genenetwork_black + labs(title = "Single WGCNA MPM-TCGA",
                                                  subtitle = "GOTEA of Module black")
p_genenetwork_black

ggsave(plot = p_genenetwork_black,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_black",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(black_ego, 
                                      foldChange=black$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-TCGA",
                                      subtitle = "GOTEA of Module black")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_black",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

################################################################################################
# GOTEA of grey60 module
################################################################################################
grey60 <- read.delim(file.path(resultsdir,
                               "25.TCGA-MESO_WGCNA_module_grey60.txt"))
colnames(grey60) <- "EntrezID"

keep <- match(grey60$EntrezID, fData$EntrezID)
grey60 <- cbind(grey60, fData[keep,])
grey60 <- grey60[c(-1)]

grey60_ego <- enrichGO(gene          = grey60$EntrezID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP",
                       pAdjustMethod = "hochberg",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)
head(grey60_ego)
head(grey60_ego@result$Description, n=20)

saveRDS(grey60_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_modulegrey60",
                        ".RDS",
                        sep = "")))


# gene concept network
# all edges and nodes
p_genenetwork_grey60 <- enrichplot::cnetplot(grey60_ego, 
                                             foldChange = grey60$Symbol,
                                             colorEdge = TRUE,
                                             showCategory = 10,
)

p_genenetwork_grey60 <- p_genenetwork_grey60 + labs(title = "Single WGCNA MPM-TCGA",
                                                    subtitle = "GOTEA of Module grey60")
p_genenetwork_grey60

ggsave(plot = p_genenetwork_grey60,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_grey60",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(grey60_ego, 
                                      foldChange=grey60$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-TCGA",
                                      subtitle = "GOTEA of Module grey60")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_grey60",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

#----------------------------------------------------------------
# get drug target check with db of top 10 hub genes for each module
# brown
brown <- read.graph(file.path(resultsdir,
                              "25.TCGA-MESO_WGCNA_Networks_brown.gml"),
                    format = "gml")
brown_vertattr <- vertex_attr(brown)
brown_vertattr$kWithinnorm[brown_vertattr$Symbol=="TERT"]

# black
black <- read.graph(file.path(resultsdir,
                              "25.TCGA-MESO_WGCNA_Networks_black.gml"),
                    format = "gml")
black_vertattr <- vertex_attr(black)  

# blue
blue <- read.graph(file.path(resultsdir,
                             "25.TCGA-MESO_WGCNA_Networks_blue.gml"),
                   format = "gml")
blue_vertattr <- vertex_attr(blue)

# define connectivity threshold
threshold <- 0.9


length(brown_vertattr$EntrezID[brown_vertattr$kWithinnorm>=threshold])

length(black_vertattr$EntrezID[black_vertattr$kWithinnorm>=threshold])

length(blue_vertattr$EntrezID[blue_vertattr$kWithinnorm>=threshold])

# table
df <- as.data.frame(matrix(data <- 0,
                           nrow = 30,
                           ncol = 6))
colnames(df) <- c("ENTREZID",
                  "SYMBOL",
                  "GENENAME",
                  "CONNECTIVITYWITHINMODULE",
                  "MODULE",
                  "DRUG")

# 
entrezids <- c(brown_vertattr$EntrezID[brown_vertattr$kWithinnorm>=threshold],
               black_vertattr$EntrezID[black_vertattr$kWithinnorm>=threshold],
               blue_vertattr$EntrezID[blue_vertattr$kWithinnorm>=threshold])

symbols <- c(brown_vertattr$Symbol[brown_vertattr$kWithinnorm>=threshold],
             black_vertattr$Symbol[black_vertattr$kWithinnorm>=threshold],
             blue_vertattr$SymbolL[blue_vertattr$kWithinnorm>=threshold])

# genenames <- c(brown_vertattr$GENENAME[brown_vertattr$kWithinnorm>=threshold],
#                black_vertattr$GENENAME[black_vertattr$kWithinnorm>=threshold],
#                blue_vertattr$GENENAME[blue_vertattr$kWithinnorm>=threshold])

kwithinnorms <- c(brown_vertattr$kWithinnorm[brown_vertattr$kWithinnorm>=threshold],
                  black_vertattr$kWithinnorm[black_vertattr$kWithinnorm>=threshold],
                  blue_vertattr$kWithinnorm[blue_vertattr$kWithinnorm>=threshold])
modulecolors <- c(rep("brown", 20),
                  rep("black", 5),
                  rep("blue", 5))


# get drug queries
results <- queryDGIdb(genes = entrezids)

for (i in 1:nrow(df)) {
        df$ENTREZID[i] <- entrezids[i]
        #df$SYMBOL[i] <- symbols[i]
        #df$GENENAME[i] <- genenames[i]
        df$CONNECTIVITYWITHINMODULE[i] <- kwithinnorms[i]
        df$MODULE[i] <- modulecolors[i]
        
        
}

df$DRUG[df$ENTREZID %in% results@matchedTerms] <- 1

# get complete fData
# feature data information
fdata <- CreateAnnotation("org.Hs.eg", 
                          keytype = "ENTREZID")
fData <- fdata[c(3,7,8)]

keep <- match(df$ENTREZID, fData$ENTREZID)
fData <- fData[keep,]

df$GENENAME <- fData$GENENAME
df <- df[order(df$MODULE, -df$CONNECTIVITYWITHINMODULE),]
write.table(x = df,
            file.path(resultsdir,
                      paste(scriptName,
                            "_hubgenes.txt",
                            sep = "")),
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# get kwithinnorm of YBX1 greenyellow
greenyellow <- read.graph(file.path(resultsdir,
                              "25.TCGA-MESO_WGCNA_Networks_greenyellow.gml"),
                    format = "gml")
greenyellow_vertattr <- vertex_attr(greenyellow)
greenyellow_vertattr$kWithinnorm[greenyellow_vertattr$Symbol=="YBX1"]


# comparision module brown to turquiose of BUENO
turquoise <- read.delim(file="../BUENO/results/25.BUENO-MESO_WGCNA_module_turquoise.txt",
                        header = FALSE)

colnames(turquoise) <- "EntrezID"

table(brown$EntrezID %in% turquoise$EntrezID)

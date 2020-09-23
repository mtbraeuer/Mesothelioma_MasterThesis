## ------------------------------------------------------
## Project: BUENO-MESO
## Script name: 26.BUENO-MESO_singleWGCNA_Analysing_Modules.R
##
## Purpose of script: Analysis of interesting modules of single WGCNA
## query of hub genes for durgability
## 
##
## Author: Marie-Theres Bräuer
##
## Version: 20. September 202
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
##
## ------------------------------------------------------
source("Initialize_BUENO-MESO.R")
#--------------------------------------------------
scriptName <- "26.BUENO-MESO_singleWGCNA_Analysing_Modules"
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
                                 "15.BUENO-MESO_stratification.ExpressionSet.RDS"))
pData <- pData(Eset)
fData <- fData(Eset)

################################################################################################
# GOTEA of turquoise module
################################################################################################
turquoise <- read.delim(file.path(resultsdir,
                                  "25.BUENO-MESO_WGCNA_module_turquoise.txt"),
                        header = FALSE)

colnames(turquoise) <- "EntrezID"

keep <- match(turquoise$EntrezID, fData$ENTREZID)
turquoise <- cbind(turquoise, fData[keep,])
turquoise <- turquoise[c(-1)]

turquoise_ego <- enrichGO(gene          = turquoise$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          keyType       = 'ENTREZID',
                          ont           = "BP",
                          pAdjustMethod = "hochberg",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05,
                          readable = TRUE)
head(turquoise_ego)
head(turquoise_ego@result$Description, n=20)

saveRDS(turquoise_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_moduleturquoise",
                        ".RDS",
                        sep = "")))


# gene concept network
# all edges and nodes
p_genenetwork <- enrichplot::cnetplot(turquoise_ego, 
                                      foldChange=turquoise$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module turquoise")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_turquoise",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(turquoise_ego, 
                                      foldChange=turquoise$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module turquoise")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_turquoise",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)


################################################################################################
# GOTEA of darkgrey module
################################################################################################
darkgrey <- read.delim(file.path(resultsdir,
                                 "25.BUENO-MESO_WGCNA_module_darkgrey.txt"),
                       header = FALSE)
colnames(darkgrey) <- "EntrezID"

keep <- match(darkgrey$EntrezID, fData$ENTREZID)
darkgrey <- cbind(darkgrey, fData[keep,])
darkgrey <- darkgrey[c(-1)]

darkgrey_ego <- enrichGO(gene          = darkgrey$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         keyType       = 'ENTREZID',
                         ont           = "BP",
                         pAdjustMethod = "hochberg",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.05,
                         readable = TRUE)
head(darkgrey_ego)
head(darkgrey_ego@result$Description, n=20)

saveRDS(darkgrey_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_moduledarkgrey",
                        ".RDS",
                        sep = "")))


# gene concept network
# all edges and nodes
p_genenetwork <- enrichplot::cnetplot(darkgrey_ego, 
                                      foldChange=darkgrey$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module darkgrey")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_darkgrey",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(darkgrey_ego, 
                                      foldChange=darkgrey$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module darkgrey")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_darkgrey",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

################################################################################################
# GOTEA of yellow module
################################################################################################
yellow <- read.delim(file.path(resultsdir,
                               "25.BUENO-MESO_WGCNA_module_yellow.txt"),
                     header = FALSE)
colnames(yellow) <- "EntrezID"

keep <- match(yellow$EntrezID, fData$ENTREZID)
yellow <- cbind(yellow, fData[keep,])
yellow <- yellow[c(-1)]

yellow_ego <- enrichGO(gene          = yellow$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       keyType       = 'ENTREZID',
                       ont           = "BP",
                       pAdjustMethod = "hochberg",
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.05,
                       readable = TRUE)
head(yellow_ego)
head(yellow_ego@result$Description, n=20)

saveRDS(yellow_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_moduleyellow",
                        ".RDS",
                        sep = "")))


# gene concept network
# all edges and nodes
p_genenetwork <- enrichplot::cnetplot(yellow_ego, 
                                      foldChange=yellow$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module yellow")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_yellow",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(yellow_ego, 
                                      foldChange=yellow$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module yellow")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_yellow",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

################################################################################################
# GOTEA of tan module
################################################################################################
tan <- read.delim(file.path(resultsdir,
                            "25.BUENO-MESO_WGCNA_module_tan.txt"),
                  header = FALSE)
colnames(tan) <- "EntrezID"

keep <- match(tan$EntrezID, fData$ENTREZID)
tan <- cbind(tan, fData[keep,])
tan <- tan[c(-1)]

tan_ego <- enrichGO(gene          = tan$ENTREZID,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENTREZID',
                    ont           = "BP",
                    pAdjustMethod = "hochberg",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.05,
                    readable = TRUE)
head(tan_ego)
head(tan_ego@result$Description, n=20)

saveRDS(tan_ego,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_results_moduletan",
                        ".RDS",
                        sep = "")))


# gene concept network
# all edges and nodes
p_genenetwork <- enrichplot::cnetplot(tan_ego, 
                                      foldChange=tan$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module tan")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot",
                       "_tan",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)
# only category in networkplot
p_genenetwork <- enrichplot::cnetplot(tan_ego, 
                                      foldChange=tan$Symbol,
                                      colorEdge = TRUE,
                                      showCategory = 10,
                                      node_label="category"
)

p_genenetwork <- p_genenetwork + labs(title = "Single WGCNA MPM-BUENO",
                                      subtitle = "GOTEA of Module tan")
p_genenetwork

ggsave(plot = p_genenetwork,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA",
                       "_networkplot_category_only",
                       "_tan",
                       ".pdf",
                       sep = "")),
       width = 15,
       height = 15)

#----------------------------------------------------------------
# get drug target check with db of top 10 hub genes for each module
# darkgrey
darkgrey <- read.graph(file.path(resultsdir,
                                 "25.BUENO-MESO_WGCNA_Networks_darkgrey.gml"),
                       format = "gml")
darkgrey_vertattr <- vertex_attr(darkgrey)

# turquoise
turquoise <- read.graph(file.path(resultsdir,
                                  "25.BUENO-MESO_WGCNA_Networks_turquoise.gml"),
                        format = "gml")
turquoise_vertattr <- vertex_attr(turquoise)  

#tan
tan <- read.graph(file.path(resultsdir,
                            "25.BUENO-MESO_WGCNA_Networks_tan.gml"),
                  format = "gml")
tan_vertattr <- vertex_attr(tan)

# define connectivity threshold
threshold <- 0.9


length(darkgrey_vertattr$ENTREZID[darkgrey_vertattr$kWithinnorm>=threshold])

length(turquoise_vertattr$ENTREZID[turquoise_vertattr$kWithinnorm>=threshold])

length(tan_vertattr$ENTREZID[tan_vertattr$kWithinnorm>=threshold])

# table
df <- as.data.frame(matrix(data <- 0,
                           nrow = 29,
                           ncol = 6))
colnames(df) <- c("ENTREZID",
                  "SYMBOL",
                  "GENENAME",
                  "CONNECTIVITYWITHINMODULE",
                  "MODULE",
                  "DRUG")

# 
entrezids <- c(darkgrey_vertattr$ENTREZID[darkgrey_vertattr$kWithinnorm>=threshold],
               turquoise_vertattr$ENTREZID[turquoise_vertattr$kWithinnorm>=threshold],
               tan_vertattr$ENTREZID[tan_vertattr$kWithinnorm>=threshold])

symbols <- c(darkgrey_vertattr$SYMBOL[darkgrey_vertattr$kWithinnorm>=threshold],
             turquoise_vertattr$SYMBOL[turquoise_vertattr$kWithinnorm>=threshold],
             tan_vertattr$SYMBOL[tan_vertattr$kWithinnorm>=threshold])

genenames <- c(darkgrey_vertattr$GENENAME[darkgrey_vertattr$kWithinnorm>=threshold],
               turquoise_vertattr$GENENAME[turquoise_vertattr$kWithinnorm>=threshold],
               tan_vertattr$GENENAME[tan_vertattr$kWithinnorm>=threshold])

kwithinnorms <- c(darkgrey_vertattr$kWithinnorm[darkgrey_vertattr$kWithinnorm>=threshold],
                  turquoise_vertattr$kWithinnorm[turquoise_vertattr$kWithinnorm>=threshold],
                  tan_vertattr$kWithinnorm[tan_vertattr$kWithinnorm>=threshold])
modulecolors <- c(rep("darkgrey", 3),
                  rep("turquoise", 20),
                  rep("tan", 6))


# get drug queries
results <- queryDGIdb(genes = entrezids)

for (i in 1:nrow(df)) {
        df$ENTREZID[i] <- entrezids[i]
        df$SYMBOL[i] <- symbols[i]
        df$GENENAME[i] <- genenames[i]
        df$CONNECTIVITYWITHINMODULE[i] <- kwithinnorms[i]
        df$MODULE[i] <- modulecolors[i]
        
        
}

df$DRUG[df$ENTREZID %in% results@matchedTerms] <- 1
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

# get kwithinnorm of GOI
turquoise_vertattr$kWithinnorm[turquoise_vertattr$SYMBOL=="TERT"]

#yellow
yellow <- read.graph(file.path(resultsdir,
                            "25.BUENO-MESO_WGCNA_Networks_yellow.gml"),
                  format = "gml")
yellow_vertattr <- vertex_attr(yellow)
yellow_vertattr$kWithinnorm[yellow_vertattr$SYMBOL=="TPH1"]

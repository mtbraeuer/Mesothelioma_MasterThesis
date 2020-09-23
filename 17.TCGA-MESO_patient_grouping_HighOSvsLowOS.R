## ------------------------------------------------------
## Project: TCGA-MESO
## Script name: 17.TCGA-MESO_patient_grouping_HighOSvsLowOS.R
##
## Purpose of script: divide patients in 2 groups high and low survivors, DEG, GOTEA
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
## good for patients
## CDKN2A high
## IL6R high
## NF2 high
## TERT low
## TPH1 high
## YBX1 low
##
## ------------------------------------------------------
source("Initialize_TCGA-MESO.R")
scriptName <- "17.TCGA-MESO_patient_grouping_HighOSvsLowOS"
library(VennDiagram)
library(EnhancedVolcano)
library(ggplot2)
library(clusterProfiler)
library(DESeq2)
library(enrichplot)
library(org.Hs.eg.db)


# filtered and mapped raw count countmatrix (not normalized)
cnts <- read.delim(file.path(resultsdir,
                             "05.TCGA-MESO_normalization_countmatrix_only_filteredmapped.txt"),
                   header = TRUE,
                   sep = "\t")

# convert numeric values in matrix to integer values
cnt <- mapply(cnts, FUN=as.integer)
str(cnt)
rownames(cnt) <- rownames(cnts)

# pData
Eset <- readRDS(file = file.path(resultsdir,
                                 "15.TCGA-MESO_Stratification.ExpressionSet.RDS"))
pData <- pData(Eset)
fData <- fData(Eset)

GOI <- read.delim(file.path("../GOI.txt"),
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

all(colnames(cnt)==pData$Barcode)
colnames(cnt) <- pData$PatientID
# patient grouping
colnames(pData)
sub <- pData[c(35,42,48,50,53,54)]
colnames(sub)

# high --> good
sub$strat_CDKN2A[sub$strat_CDKN2A=="high"] <- 1
sub$strat_CDKN2A[sub$strat_CDKN2A=="low"] <- 0

sub$strat_IL6R[sub$strat_IL6R=="high"] <- 1
sub$strat_IL6R[sub$strat_IL6R=="low"] <- 0

sub$strat_TPH1[sub$strat_TPH1=="high"] <- 1
sub$strat_TPH1[sub$strat_TPH1=="low"] <- 0

sub$strat_NF2[sub$strat_NF2=="high"] <- 1
sub$strat_NF2[sub$strat_NF2=="low"] <- 0

sub$strat_TERT[sub$strat_TERT=="high"] <- 0
sub$strat_TERT[sub$strat_TERT=="low"] <- 1

sub$strat_YBX1[sub$strat_YBX1=="high"] <- 0
sub$strat_YBX1[sub$strat_YBX1=="low"] <- 1

str(sub)
sub[,1:ncol(sub)] <- lapply(sub[,1:ncol(sub)], as.numeric)
str(sub)

sub[,"RowSum"] <-  apply(sub, 1, FUN=function(x) sum(x))
table(sub$RowSum)

sub["HighOSvsLowOS"] <- "low"
sub$HighOSvsLowOS[sub$RowSum>=4] <- "high"
table(sub$HighOSvsLowOS)

sub <- cbind(sub, pData$OS_days)

mean(sub$`pData$OS_days`[sub$HighOSvsLowOS=="high"], na.rm=TRUE)
mean(sub$`pData$OS_days`[sub$HighOSvsLowOS=="low"], na.rm=TRUE)


for (i in 1:5) {
    print(i)
    print(mean(sub$`pData$OS_days`[sub$RowSum==i], na.rm=TRUE))
}

pData<- cbind(pData, sub$HighOSvsLowOS)
colnames(pData)[ncol(pData)] <- "HighOSvsLowOS"

write.table(pData,
            file.path(resultsdir,
                      paste(scriptName,
                            "_pData_HighOSvsLowOS.txt",
                            sep = "")),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)


pData$HighOSvsLowOS <- as.factor(pData$HighOSvsLowOS)
gene <- "HighOSvsLowOS"

############################################################
# Survival curve of high vs low survivors
############################################################
colnames(pData)
sub <- pData[,c(1,3,4,55)]
colnames(sub)
sub$HighOSvsLowOS <- factor(sub$HighOSvsLowOS)

#surv_object <- Surv(time = sub$OS_days, event = sub$VitalStatus)
#surv_object

fit <-survminer::surv_fit( formula = survival::Surv(time = sub$OS_days) ~ sub$HighOSvsLowOS, 
                           data = sub)

names(fit$strata) <- c("high survivors", "low survivors")

pdf(file=file.path(plotsdir, 
                   paste(scriptName,
                         "SurvivalCurve_highvslowsurvivors",
                         "pdf",
                         sep=".")),
    #paper="a4",
    onefile = FALSE)

p <- ggsurvplot(fit,
                data = sub,
                risk.table = TRUE,
                pval = TRUE,
                conf.int = TRUE,
                xlim = c(0,2800),
                xlab = "Time in days",
                ggtheme = theme_bw() +
                    theme(panel.grid.major = element_blank(),
                          panel.grid.minor = element_blank(),
                          panel.background = element_blank()),
                break.time.by = 250,
                risk.table.y.text.col = T,
                risk.table.y.text = FALSE,
                title = paste("Survival Curve of high vs low survivors", "\n",
                              "data set: TCGA; patients: ", nrow(sub),
                              sep = ""))
print(p)
dev.off()

############################################################
# Differential expression analysis high vs low survivors
############################################################
dds <- DESeqDataSetFromMatrix(countData = cnt,
                              colData = pData,
                              design = ~HighOSvsLowOS)

# prefiltering necessary?
table(rowSums(counts(dds)) >= 10)

# dds data set creation
dds$HighOSvsLowOS <- factor(dds$HighOSvsLowOS, levels = c("high","low"))

#dds$HighOSvsLowOS <- relevel(dds$HighOSvsLowOS, ref = "high")

dds <- DESeq(dds)
res <- results(dds)
#res <- results(dds, contrast=c("HighOSvsLowOS","high", "low"))
head(res)
# 
saveRDS(res,
        file.path(resultsdir,
                  paste(scriptName,
                        "_DEG_results",
                        ".RDS",
                        sep = "")))


# write to data frame for further processing
res <- as.data.frame(res)

# keep only entries with entrez id
keep <- match(rownames(res), fData$EntrezID)
res$GeneSymbol <- fData$Symbol[keep]
res$EntrezID <- fData$EntrezID[keep]

dim(res)
res <- res[!is.na(res$padj),]
res <- res[!is.na(res$GeneSymbol),]
dim(res)

#save results of DEG
write.table(res,
            file.path(resultsdir,
                      paste(scriptName,
                            "_DEG_results",
                            ".txt",
                            sep = "")),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

#How many adjusted p-values were less than 0.5?
# get results low
g1low <- res$GeneSymbol[res$log2FoldChange<=-1.5 & res$padj <=0.05]
length(g1low)
keep <- match(g1low, res$GeneSymbol)
low <- res[keep,]

# save low
write.table(low,
            file.path(resultsdir,
                      paste(scriptName,
                            "_DEG_low_significant_results",
                            ".txt",
                            sep = "")),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# get results high
g1high <- res$GeneSymbol[res$log2FoldChange >= 1.5 & res$padj <=0.05]
length(g1high)
keep <- match(g1high, res$GeneSymbol)
high <- res[keep,]

# save high
write.table(high,
            file.path(resultsdir,
                      paste(scriptName,
                            "_DEG_high_significant_results",
                            ".txt",
                            sep = "")),
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = TRUE)

# get custom color labels for volcano plot
# up-regulated: red
# down-regulated: blue

GOI <- GOI[c(4, 11, 17, 19, 22,23),]
for (i in 1:nrow(GOI)) {
    print(GOI$GeneSymbol[i])
    print(res$padj[res$EntrezID==GOI$EntrezID[i]])
    print(res$log2FoldChange[res$EntrezID==GOI$EntrezID[i]])
    
}

keyvals <- ifelse(res$padj>0.05, 'grey',
                  ifelse(res$log2FoldChange < -1.5, 'royalblue',
                         ifelse(res$log2FoldChange > 1.5, 'red',
                                'grey40')))

table(keyvals)

keyvals[is.na(keyvals)] <- 'grey'

names(keyvals)[keyvals == 'red'] <- 'up-regulated'
names(keyvals)[keyvals == 'royalblue'] <- 'down-regulated'
names(keyvals)[keyvals == 'grey'] <- 'NS'
names(keyvals)[keyvals == 'grey40'] <- 'padj. > 0.05'

mylabs <- c("CDKN2A",
            "IL6R",
            "NF2",
            "TERT",
            "TPH1",
            "YBX1")

mylabs <- c(mylabs,
            res$GeneSymbol[res$log2FoldChange>=4 & res$padj<=0.01],
            res$GeneSymbol[res$log2FoldChange<=-4 & res$padj<=0.01],
            res$GeneSymbol[res$padj<=1e-9 & res$log2FoldChange<=-2],
            res$GeneSymbol[res$padj<=1e-9 & res$log2FoldChange>=1.5])


# plotting **** FEHLER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
p_volc <- EnhancedVolcano(res,
                          lab = res$GeneSymbol,
                          x = 'log2FoldChange',
                          y = 'padj',
                          xlim = c(-7.5,7.5),
                          ylim = c(0,25),
                          selectLab = mylabs, #res$GeneSymbol[which(names(keyvals) %in% c('up-regulated', 'down-regulated'))],
                          legendLabSize = 8,
                          title = "DEG MPM-TCGA",
                          subtitle = "High vs. Low Survivors",
                          pCutoff = 0.05,
                          FCcutoff = 1.5,
                          xlab = bquote(~Log[2]~ 'fold change'),
                          ylab = bquote(~-Log[10]~adjusted~italic(P)),
                          colCustom = keyvals,
                          drawConnectors = TRUE,
                          widthConnectors = 0.3,
                          colConnectors = 'black',
                          gridlines.major = FALSE,
                          gridlines.minor = FALSE,
                          border = 'partial',
                          borderWidth = 1.0,
                          borderColour = 'black',
                          labSize = 2.5,
)
p_volc
ggsave(file.path(plotsdir,
                 paste(scriptName,
                       "_DEG_volcanoplot",
                       ".pdf",
                       sep = "")),
       plot = p_volc)


################################################################################################
# GOTEA of DEG low
################################################################################################
egol <- enrichGO(gene          = low$EntrezID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "hochberg",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)
head(egol)
head(egol@result$Description, n=10)

saveRDS(egol,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_low_results",
                        ".RDS",
                        sep = "")))

# dot plot for displaying most significant enriched terms
p_dot_low <- clusterProfiler::dotplot(egol, showCategory=10)
p_dot_low <- p_dot_low + 
    labs(title = "GOTEA MPM-TCGA",
         subtitle = "down-regulated") +
    xlim(c(0, 0.3))

p_dot_low
ggsave(plot = p_dot_low,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA_low",
                       "_dotplot",
                       ".pdf",
                       sep = "")))

# gene concept network
# convert gene ID to Symbol
genelist_low <- low$log2FoldChange
names(genelist_low) <- low$EntrezID
genelist_low <- sort(genelist_low, decreasing = TRUE)

p_genenetwork_low <- enrichplot::cnetplot(egol, 
                                          foldChange=genelist_low,
                                          categorySize="pvalue",
                                          circular = TRUE,
                                          colorEdge = TRUE,
                                          showCategory = 10)

p_genenetwork_low <- p_genenetwork_low + labs(title = "GOTEA MPM-TCGA",
                                              subtitle = "down-regulated")
p_genenetwork_low

ggsave(plot = p_genenetwork_low,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA_low",
                       "_networkplot",
                       ".pdf",
                       sep = "")),
       width = 12,
       height = 9)

################################################################################################
# GOTEA of DEG high
################################################################################################
egoh <- enrichGO(gene         = high$EntrezID,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable = TRUE)
head(egoh)
head(egoh@result$Description)

saveRDS(egoh,
        file.path(resultsdir,
                  paste(scriptName,
                        "_GOETA_high_results",
                        ".RDS",
                        sep = "")))

# dot plot for displaying most significant enriched terms
p_dot_high <- clusterProfiler::dotplot(egoh, showCategory=10)
p_dot_high <- p_dot_high + labs(title = "GOTEA MPM-TCGA",
                                subtitle = "up-regulated")

p_dot_high
ggsave(plot = p_dot_high,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA_high",
                       "_dotplot",
                       ".pdf",
                       sep = "")))

# gene concept network
# convert gene ID to Symbol
genelist_high <- high$log2FoldChange
names(genelist_high) <- high$EntrezID
genelist_high <- sort(genelist_high, decreasing = TRUE)

p_genenetwork_high <- enrichplot::cnetplot(egoh, 
                                           foldChange=genelist_high,
                                           categorySize="pvalue",
                                           circular = TRUE,
                                           colorEdge = TRUE,
                                           showCategory = 10)

p_genenetwork_high <- p_genenetwork_high + labs(title = "GOTEA MPM-TCGA ",
                                                subtitle = "up-regulated")
p_genenetwork_high

ggsave(plot = p_genenetwork_high,
       file.path(plotsdir,
                 paste(scriptName,
                       "_GOTEA_high",
                       "_networkplot",
                       ".pdf",
                       sep = "")),
       width = 12,
       height = 9)


## ------------------------------------------------------
## Project: Mesothelioma Masterthesis
## Script name: function-StratificationLowHighExpressors_SurvivialCurves.R
##
## Purpose of script: StratificationLowHighExpressors_SurvivialCurves
## 
##
## Author: Marie-Theres Bräuer
##
## Version: October 2019
##
## Copyright (c) Marie-Theres Bräuer, 2020
## Email: mt_braeuer@hotmail.com
##
## ------------------------------------------------------
## Notes: 
## input
## expressionmatrix
## stratificationDF
## goiEntrezID
## goiName
## output
## returnDF stratificated DF
## pval
##
## http://r-addict.com/2016/11/21/Optimal-Cutpoint-maxstat.html
## https://rpkgs.datanovia.com/survminer/index.html
##
## ------------------------------------------------------
library(survminer)

StratificationLowHighExpressors_SurvivialCurves <- function(expressionMatrix,
                                                            stratificationDF,
                                                            goiEntrezID,
                                                            goiName){
    
    # add expression values to df
    stratificationDF <- cbind(stratificationDF, as.data.frame(expressionMatrix[goiEntrezID,]))
    colnames(stratificationDF)[ncol(stratificationDF)] <- goiName
    
    # Determining the optimal cutpoint for GeneX gene's expression
    gene_cut <- survminer::surv_cutpoint(data = stratificationDF,
                                         time = "OS_days",
                                         event = "VitalStatus",
                                         variables = c(goiName))
    print(summary(gene_cut))
    
    #Plot the cutpoint
    pdf(file=file.path(plotsdir,
                       paste(scriptName,
                             goiName,
                             "Distribution_MaxRankStat",
                             "pdf",
                             sep=".")),
        #paper="a4",
        onefile = FALSE)
    
    print(plot(gene_cut, goiName,
               palette = "npg"))
    dev.off()
    
    # Categorize GeneX variable into high and low expressors
    gene_cat <- survminer::surv_categorize(x=gene_cut)
    print(table(gene_cat[goiName]))
    
    # Fit and visualize Kaplan-Meier estimates of survival curves
    fit <- survminer::surv_fit(formula = survival::Surv(gene_cat$OS_days) ~ gene_cat[,goiName], 
                               data = gene_cat)
    # change name in fit$strata for plotting
    names(fit$strata) <- c(paste(goiName,
                                 " high", 
                                 sep = ""), 
                           paste(goiName,
                                 " low", 
                                 sep = ""))
    #print(summary(fit))
    #class(fit)
    pval <- survminer::surv_pvalue(fit)
    print(pval[[4]])
    
    pdf(file=file.path(plotsdir, 
                       paste(scriptName,
                             goiName,
                             "SurvivalCurve",
                             "pdf",
                             sep=".")),
        #paper="a4",
        onefile = FALSE)
    
    p <- ggsurvplot(fit,
                    data = stratificationDF,
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
                    title = paste("Survival curve of ", goiName, " stratification\n", 
                                  "data set: BUENO; patients: ", nrow(stratificationDF),
                                  sep = ""))
    print(p)
    dev.off()
    
    # median overall survival 
    low <- median(gene_cat$OS_days[gene_cat[,goiName]=="low"], na.rm = T)
    high <- median(gene_cat$OS_days[gene_cat[,goiName]=="high"], na.rm = T)
    
    print(paste("median overall survival high:", high))
    print(paste("median overall survival low:", low))
    
    # add stratification col to return data frame
    stratificationDF <- cbind(stratificationDF, gene_cat[,goiName])
    colnames(stratificationDF)[ncol(stratificationDF)] <- paste("stratification_",
                                                                goiName,
                                                                sep = "")
    returnDF <- subset(stratificationDF,
                       select = paste("stratification_",
                                      goiName,
                                      sep = ""))
    
    # write log file 
    printer = file(file.path(resultsdir,
                             paste(scriptName,
                                   goiName,
                                   "log",
                                   "txt",
                                   sep = ".")),
                   "w")
    writeLines(goiName, 
               con = printer)
    writeLines(paste("Number_High",table(gene_cat[goiName])[[1]], sep = "\t"), 
               con = printer)
    writeLines(paste("Number_Low",table(gene_cat[goiName])[[2]], sep = "\t"), 
               con = printer)
    writeLines(paste("Median_overall_survival_high", high, sep = "\t"),
               con = printer)
    writeLines(paste("Median_overall_survival_low", low, sep = "\t"),
               con = printer)
    close(printer)
    
    
    
    # return
    return(c(returnDF, pval$pval.txt))
} 
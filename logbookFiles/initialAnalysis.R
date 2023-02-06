# loading libraries for the workflow
library(devtools)
library(remotes)
remotes::install_github("b-klaus/maEndToEnd", ref="master")
suppressPackageStartupMessages({library("maEndToEnd")})


#General Bioconductor packages
library(Biobase)
library(oligoClasses)


#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hgu133plus2.db)


#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)


#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)


#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)


#Formatting/documentation packages
library(dplyr)
library(tidyr)


#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)


# downloading the raw data from ArrayExpress
raw_data_dir <- "E-GEOD-68848"
anno_AE <- getAE("E-GEOD-68848", path = raw_data_dir, type = "raw") # took > 1 hour to load


# importing the data into R
sdrf_location <- file.path(raw_data_dir, "E-GEOD-68848.sdrf.txt")
SDRF <- read.delim(sdrf_location)


rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)


raw_data <- oligo::read.celfiles(filesnames = file.path(raw_data_dir,
                                                       SDRF$Array.Data.File),
                                verbose = FALSE, phenoData = SDRF)


## looking at our columns of the whole table
head(Biobase::pData(raw_data)) # Table 1: Raw Data with important columns


## we subselect the columns important to us
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Source.Name", "FactorValue..tumor.grading.",
                                                        "FactorValue..disease.state.")]


# Quality control of the data


## checking the first few data values of the raw data
Biobase::exprs(raw_data)[1:5, 1:5] # Table 2: expression values


## performing PCA and plotting it. every point in the plot represents one sample
expr_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(expr_raw), scale. = FALSE) # 4:41PM to 5:45PM: ~1 hr to run


percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)


dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                    Disease = pData(raw_data)$FactorValue..disease.state.,
                    Grade = pData(raw_data)$FactorValue..tumor.grading.,
                    Individual = pData(raw_data)$Source.Name)


ggplot(dataGG, aes(PC1, PC2)) +
 geom_point(aes(shape = Disease, colour = Grade)) +
 ggtitle("PCA plot of the log-transformed raw expression data") +
 xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
 ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
 theme(plot.title = element_text(hjust = 0.5))+
 scale_shape_manual(values = c(4,5,6,7,8,10,13,15)) +
 scale_color_manual(values = c("antiquewhite4", "dodgerblue4", "blueviolet", "chartreuse", "darkgoldenrod1")) # Figure 1 in the graphs


oligo::boxplot(raw_data, target = "core",
              main = "Boxplot of log2-intensitites for the raw data") # Figure 2
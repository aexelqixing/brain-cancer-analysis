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
head(Biobase::pData(raw_data))

## we subselect the columns important to us
Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Source.Name", "FactorValue..tumor.grading.",
                                                         "FactorValue..disease.state.")]

# Quality control of the data

## checking the first few data values of the raw data
Biobase::exprs(raw_data)[1:5, 1:5]

## performing PCA and plotting it. every point in the plot represents one sample
expr_raw <- log2(Biobase::exprs(raw_data)) 
PCA_raw <- prcomp(t(expr_raw), scale. = FALSE) # 4:41PM to 5:45PM 

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
  scale_color_manual(values = c("antiquewhite4", "dodgerblue4", "blueviolet", "chartreuse", "darkgoldenrod1"))

oligo::boxplot(raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

## Filtering astrocytoma and oligodendroglioma 
filter <- colnames(raw_data)[raw_data@phenoData@data$`FactorValue..disease.state.`=="astrocytoma" | raw_data@phenoData@data$`FactorValue..disease.state.`=="oligodendroglioma"]
raw_data_ao <- raw_data[,filter]

filter <- colnames(raw_data_ao)[raw_data_ao@phenoData@data$`FactorValue..tumor.grading.`=="Grade 2" | raw_data_ao@phenoData@data$`FactorValue..tumor.grading.`=="Grade 3"]
raw_data_ao <- raw_data_ao[,filter]

### making the PCA plot for the astrocytoma and oligodendroglioma plot
expr_raw_ao <- log2(Biobase::exprs(raw_data_ao)) 
PCA_raw_ao <- prcomp(t(expr_raw_ao), scale. = FALSE) #10:40PM to 10:50PM

percentVar_ao <- round(100*PCA_raw_ao$sdev^2/sum(PCA_raw_ao$sdev^2),1)

dataGG_ao <- data.frame(PC1 = PCA_raw_ao$x[,1], PC2 = PCA_raw_ao$x[,2],
                     Disease = pData(raw_data_ao)$FactorValue..disease.state.,
                     Grade = pData(raw_data_ao)$FactorValue..tumor.grading.,
                     Individual = pData(raw_data_ao)$Source.Name)

ggplot(dataGG_ao, aes(PC1, PC2)) +
  geom_point(aes(shape = Grade, colour = Disease)) +
  ggtitle("PCA plot of the log-transformed raw expression data (astrocytoma and oligodendroglioma)") +
  xlab(paste0("PC1, VarExp: ", percentVar_ao[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar_ao[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

oligo::boxplot(raw_data_ao, target = "core",
               main = "Boxplot of log2-intensitites for the raw data (astrocytoma and oligodendroglioma) ")

## Filtering astrocytoma, oligodendroglioma and normal
filter <- colnames(raw_data)[raw_data@phenoData@data$`FactorValue..disease.state.`=="astrocytoma" | raw_data@phenoData@data$`FactorValue..disease.state.`=="oligodendroglioma" | raw_data@phenoData@data$`FactorValue..disease.state.`=="NON_TUMOR"]
raw_data_aon <- raw_data[,filter]

filter <- colnames(raw_data_aon)[raw_data_aon@phenoData@data$`FactorValue..tumor.grading.`=="Grade 3" | raw_data_aon@phenoData@data$`FactorValue..tumor.grading.`=="--"]
raw_data_aon <- raw_data_aon[,filter]

### making PCA plot for this 
expr_raw_aon <- log2(Biobase::exprs(raw_data_aon)) 
PCA_raw_aon <- prcomp(t(expr_raw_aon), scale. = FALSE) #from 11:02PM to 11:08PM

percentVar_aon <- round(100*PCA_raw_aon$sdev^2/sum(PCA_raw_aon$sdev^2),1)

dataGG_aon <- data.frame(PC1 = PCA_raw_aon$x[,1], PC2 = PCA_raw_aon$x[,2],
                        Disease = pData(raw_data_aon)$FactorValue..disease.state.,
                        Grade = pData(raw_data_aon)$FactorValue..tumor.grading.,
                        Individual = pData(raw_data_aon)$Source.Name)

ggplot(dataGG_aon, aes(PC1, PC2)) +
  geom_point(aes(shape = Grade, colour = Disease)) +
  ggtitle("PCA plot of the log-transformed raw expression data (astrocytoma, oligodendroglioma, and normal)") +
  xlab(paste0("PC1, VarExp: ", percentVar_aon[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar_aon[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("antiquewhite4", "blueviolet", "chartreuse"))

oligo::boxplot(raw_data_aon, target = "core",
               main = "Boxplot of log2-intensitites for the raw data (astrocytoma, oligodendroglioma, and normal) ")

## Filtering astrocytoma and normal cells
filter <- colnames(raw_data)[raw_data@phenoData@data$`FactorValue..disease.state.`=="astrocytoma" | raw_data@phenoData@data$`FactorValue..disease.state.`=="NON_TUMOR"]
raw_data_an <- raw_data[,filter]

filter <- colnames(raw_data_an)[raw_data_an@phenoData@data$`FactorValue..tumor.grading.`=="Grade 3" | raw_data_an@phenoData@data$`FactorValue..tumor.grading.`=="--"]
raw_data_an <- raw_data_an[,filter]

### making PCA plot for this 
expr_raw_an <- log2(Biobase::exprs(raw_data_an)) 
PCA_raw_an <- prcomp(t(expr_raw_an), scale. = FALSE) #from 9:11AM to 9:15AM

percentVar_an <- round(100*PCA_raw_an$sdev^2/sum(PCA_raw_an$sdev^2),1)

dataGG_an <- data.frame(PC1 = PCA_raw_an$x[,1], PC2 = PCA_raw_an$x[,2],
                         Disease = pData(raw_data_an)$FactorValue..disease.state.,
                         Grade = pData(raw_data_an)$FactorValue..tumor.grading.,
                         Individual = pData(raw_data_an)$Source.Name)

ggplot(dataGG_an, aes(PC1, PC2)) +
  geom_point(aes(shape = Grade, colour = Disease)) +
  ggtitle("PCA plot of the log-transformed raw expression data (astrocytoma and normal)") +
  xlab(paste0("PC1, VarExp: ", percentVar_an[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar_an[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("antiquewhite4", "blueviolet"))

oligo::boxplot(raw_data_aon, target = "core",
               main = "Boxplot of log2-intensitites for the raw data (astrocytoma and normal) ")

# RLE data quality analysis
palmieri_eset <- oligo::rma(raw_data_an, normalize=FALSE)

row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(palmieri_eset)))

RLE_data <- sweep(Biobase::exprs(palmieri_eset), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)

RLE_data_gathered <-
  tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

ggplot2::ggplot(RLE_data_gathered, aes(patient_array,
                                       log2_expression_deviation)) + 
  geom_boxplot(outlier.shape = NA) + 
  ylim(c(-2, 2)) + 
  theme(axis.text.x = element_text(colour = "aquamarine4",
                                   angle = 60, size = 6.5, hjust = 1 ,
                                   face = "bold"))

# RMA calibration of the data
palmieri_eset_norm <- oligo::rma(raw_data_an)

# Quality assessment of calibrated data
exp_palmieri <- Biobase::exprs(palmieri_eset_norm)
PCA <- prcomp(t(exp_palmieri), scale = FALSE) # 2 seconds to run

percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     Disease =
                       Biobase::pData(palmieri_eset_norm)$FactorValue..disease.state.,
                     Grade =
                       Biobase::pData(palmieri_eset_norm)$FactorValue..tumor.grading.)

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Grade, colour = Disease)) +
  ggtitle("PCA plot of the calibrated, summarized data (astrocytoma and normal)") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("antiquewhite4", "blueviolet"))

# plotting before quantile normalization
hist(raw_data,lwd=2,xlab='log intensity', which='pm',
     main="CEL file densities before quantile normalisation")

# plotting after quantile normalization
oligo_normalised <- normalize(raw_data,method='quantile',which='pm')

hist(oligo_normalised,lwd=2,xlab='log intensity', which='pm',
     main="CEL file densities after quantile normalisation")
library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]

pData <- pData(gse108474)

filter <- colnames(gse108474)[gse108474@phenoData@data$`disease:ch1`=="astrocytoma" | gse108474@phenoData@data$`disease:ch1`=="normal" | gse108474@phenoData@data$`disease:ch1`=="glioblastoma multiforme" | gse108474@phenoData@data$`disease:ch1`=="oligodendroglioma"]
gse108474 <- gse108474[,filter]

library(oligo)
library(oligoClasses)
library(tidyverse)
pd <- pData(gse108474)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)

# Loading the raw data
gse108474_celdata <- read.celfiles(pd$cel_file,phenoData=phenoData(gse108474))

# Showing the columns that matter
pData(gse108474_celdata)[,c("geo_accession", "disease:ch1")]

# Displaying the first few expression data
exprs(gse108474_celdata)[1:4, 1:4]

gse108474_eset <- rma(gse108474_celdata)

library(hgu133plus2.db)
library(AnnotationDbi)

AnnotationDbi::select(hgu133plus2.db,c("IDH1"),c("PROBEID","ENTREZID","GENENAME"),keytype="SYMBOL") 
expression_df <- exprs(gse108474_eset)
metadata <- pData(gse108474_eset)

expression_df <- as.data.frame(expression_df) # turn into dataframe instead of matrix
idh1_probes <- t(expression_df[c("1555037_a_at", "201193_at", "242956_at"),])

metadata[c("disease:ch1")]

idh1_probes <- cbind(idh1_probes, metadata[c("disease:ch1")])

astrocytoma <- idh1_probes[which(idh1_probes$`disease:ch1`=="astrocytoma"),]
astrocytoma[1:3] <- sapply(astrocytoma[1:3], as.numeric)

apply(astrocytoma[1:3], 2, mean)
apply(astrocytoma[1:3], 2, sd)

oligodendroglioma <- idh1_probes[which(idh1_probes$`disease:ch1`=="oligodendroglioma"),]

apply(oligodendroglioma[1:3], 2, mean)
apply(oligodendroglioma[1:3], 2, sd)

gbm <- idh1_probes[which(idh1_probes$`disease:ch1`=="glioblastoma multiforme"),]

apply(gbm[1:3], 2, mean)
apply(gbm[1:3], 2, sd)

normal <- idh1_probes[which(idh1_probes$`disease:ch1`=="normal"),]

apply(normal[1:3], 2, mean)
apply(normal[1:3], 2, sd)


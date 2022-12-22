library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]

# Filtering for only astrocytoma and normal samples
filter <- colnames(gse108474)[gse108474@phenoData@data$`disease:ch1`=="astrocytoma" | gse108474@phenoData@data$`disease:ch1`=="normal"]
gse108474 <- gse108474[,filter]

library(oligo)
library(oligoClasses)

library(tidyverse)
pd <- pData(gse108474)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)

gse108474_celdata <- read.celfiles(pd$cel_file,phenoData=phenoData(gse108474))
gse108474_eset <- rma(gse108474_celdata)
expression_df <- as.data.frame(exprs(gse108474_eset))
metadata <- pData(gse108474_eset)

library(hgu133plus2.db)
ls('package:hgu133plus2.db')

ref_genes = c("PCLO", "BSN", "ERC2", "NEFL", "NEFH", "SHANK2")
library(AnnotationDbi)
df <- AnnotationDbi::select(hgu133plus2.db,ref_genes,c("PROBEID","ENTREZID","GENENAME"),keytype="SYMBOL")

# PCLO's gene expression comparison
expression_df <- exprs(gse108474_eset)
metadata <- pData(gse108474_eset)

expression_df <- as.data.frame(expression_df) # turn into dataframe instead of matrix
top_gene_df <- t(expression_df[c("210650_s_at"),])
top_gene_df <- cbind(rownames(top_gene_df), top_gene_df)
colnames(top_gene_df)[1] = "geo_accession"

top_gene_df <- merge(top_gene_df, metadata, by.x = "geo_accession")

rownames(top_gene_df) <- top_gene_df$geo_accession

library(ggplot2)
ggplot(top_gene_df, aes(x = `disease:ch1`, y = `210650_s_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()









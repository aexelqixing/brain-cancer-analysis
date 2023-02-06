# Loading annotated data
library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]


library(oligo)
library(oligoClasses)


library(tidyverse)
varLabels(gse108474) # What different categories are there?


# Filtering for only astrocytoma and normal samples
filter <- colnames(gse108474)[gse108474@phenoData@data$`disease:ch1`=="astrocytoma" | gse108474@phenoData@data$`disease:ch1`=="normal"]


gse108474 <- gse108474[,filter]


pd <- pData(gse108474)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)


# Loading the raw data
gse108474_celdata <- read.celfiles(pd$cel_file,phenoData=phenoData(gse108474))


# Showing the columns that matter
pData(gse108474_celdata)[,c("geo_accession", "disease:ch1")]


# Displaying the first few expression data
exprs(gse108474_celdata)[1:4, 1:4]


hist(gse108474_celdata,lwd=2,xlab='log intensity', which='pm',
    main="CEL file densities before quantile normalization")                               


oligo_normalized <- normalize(gse108474_celdata, method='quantile', which='pm')
hist(oligo_normalized, lwd=2, xlab='log intensity', which='pm',
    main="CEL file densities after quantile normalization")


oligo_summarized <- rma(oligo_normalized, background=FALSE, normalize=FALSE)
hist(oligo_summarized)


gse108474_eset <- rma(gse108474_celdata)


# Doing differential analysis
library(limma)
design <- model.matrix(~ `disease:ch1` + 0, data = pd) # turning strings to categories
colnames(design) <- levels(as.factor(gse108474_eset[['disease:ch1']])) # rename columns


fit <- lmFit(gse108474_eset, design) # fit a linear model
fit <- eBayes(fit) # get rid of any errors


## Using contrasts to make linear model
contrast_matrix <- makeContrasts(astrocytoma - normal, levels=design)


contrasts_fit <- contrasts.fit(fit, contrast_matrix)
contrasts_fit <- eBayes(contrasts_fit)


stats_df <- topTable(contrasts_fit) # top 10 probes that are differentially expressed


summary(decideTests(contrasts_fit,lfc=1)) # downregulated and upregulated genes


## Finding gene names for top 10 probes
library(hgu133plus2.db)
library(AnnotationDbi)
ls('package:hgu133plus2.db')


### extra info about the packages
columns(hgu133plus2.db)
keytypes(hgu133plus2.db)
head(keys(hgu133plus2.db,keyType="PROBEID"))


ps <- rownames(topTable(contrasts_fit))
AnnotationDbi::select(hgu133plus2.db,rownames(stats_df),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")


# SERTM1's gene expression comparison
expression_df <- exprs(gse108474_eset)
metadata <- pData(gse108474_eset)


expression_df <- as.data.frame(expression_df) # turn into dataframe instead of matrix
top_gene_df <- t(expression_df[c("241672_at"),])
top_gene_df <- cbind(rownames(top_gene_df), top_gene_df)
colnames(top_gene_df)[1] = "geo_accession"


top_gene_df <- merge(top_gene_df, metadata, by.x = "geo_accession")


rownames(top_gene_df) <- top_gene_df$geo_accession


library(ggplot2)
ggplot(top_gene_df, aes(x = `disease:ch1`, y = `241672_at`, color = `disease:ch1`)) +
 geom_jitter(width = 0.2, height = 0) +
 theme_classic()


# Making volcano plot
volcanoplot(contrasts_fit)


interesting_genes <- topTable(contrasts_fit,number=Inf,p.value=0.05,lfc=2)


volcanoplot(contrasts_fit,main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')


# Making heatmap
library(RColorBrewer)


eset_of_interest <- gse108474_eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest),
       labCol=gse108474_eset[['disease:ch1']],labRow=NA,
       col=rev(brewer.pal(10,"RdBu")),
       distfun = function(x) as.dist(1-cor(t(x))))
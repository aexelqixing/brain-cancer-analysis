# This file is where I analyzes all of the genes that were enriched between astrocytoma and normal samples. It allowed me to see which genes were upregulated and downregulated between these samples and also do some enrichment analysis.

# Loading annotated data
library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]

library(oligo)
library(oligoClasses)

library(tidyverse)
varLabels(gse108474) # what different categories are there?

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
#--------------
stats_df <- topTable(fit)

expression_df <- exprs(gse108474_eset)
metadata <- pData(gse108474_eset)

expression_df <- as.data.frame(expression_df)
t <- t(expression_df)
top_gene_df <- expression_df[c("241672_at"),]
top_gene_df <- t(top_gene_df)
top_gene_df <- cbind(rownames(top_gene_df), top_gene_df)
colnames(top_gene_df)[1] = "geo_accession"

top_gene_df <- merge(top_gene_df, metadata, by.x = "geo_accession")

rownames(top_gene_df) <- top_gene_df$geo_accession

library(ggplot2)
ggplot(top_gene_df, aes(x = `disease:ch1`, y = `241672_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

fit2 <- contrasts_fit
summary(decideTests(fit2,lfc=1))

library(hgu133plus2.db)
ls('package:hgu133plus2.db')

fitted.ebayes <- fit2
topTable(fitted.ebayes)

ps <- rownames(topTable(fitted.ebayes))
ps

unlist(mget(ps,hgu133plus2SYMBOL))

columns(hgu133plus2.db)
keytypes(hgu133plus2.db)
head(keys(hgu133plus2.db,keyType="PROBEID"))

library(AnnotationDbi)
AnnotationDbi::select(hgu133plus2.db,ps,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

ps2 <- topTable(fitted.ebayes,number=Inf,p.value=0.05,lfc=1)
ps2_up <- rownames(ps2[ps2$logFC>0,])
df <- AnnotationDbi::select(hgu133plus2.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

dplyr::mutate(df,GENENAME=stringr::str_trunc(GENENAME,30))

volcanoplot(fitted.ebayes)

interesting_genes <- topTable(fitted.ebayes,number=Inf,p.value=0.05,lfc=2)

volcanoplot(fitted.ebayes,main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

eset_of_interest <- gse108474_eset[rownames(interesting_genes),]
heatmap(exprs(eset_of_interest))

write.csv(expression_df, "expressions_for_astrocytoma_vs_normal.csv")
write.csv(exprs(eset_of_interest), "interesting_genes_expressions_for_astrocytoma_vs_normal.csv")

write.csv(metadata, "metadata_for_all_astrocytoma_vs_normal.csv")
write.csv(pData(eset_of_interest), "interesting_metadata_for_astrocytoma_vs_normal.csv")

library(RColorBrewer)
heatmap(exprs(eset_of_interest),
        labCol=gse108474_eset[['disease:ch1']],labRow=NA,
        col=rev(brewer.pal(10,"RdBu")),
        distfun = function(x) as.dist(1-cor(t(x))))

pdf("heatmap of interesting genes")
dev.off()

# Create the data folder if it doesn't exist
if (!dir.exists("data")) {
  dir.create("data")
}

# Define the file path to the plots directory
plots_dir <- "plots"

# Create the plots folder if it doesn't exist
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# Define the file path to the results directory
results_dir <- "results"

# Create the results folder if it doesn't exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("msigdbr" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("msigdbr", update = FALSE)
}

if (!("org.Dr.eg.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.Dr.eg.db", update = FALSE)
}

# Attach the library
library(clusterProfiler)

# Package that contains MSigDB gene sets in tidy format
library(msigdbr)

# Zebrafish annotation package we'll use for gene identifier conversion
library(org.Dr.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

dge_df <- ps2
dge_df

dr_hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)
head(dr_hallmark_df)

dge_df = cbind(row.names(dge_df), dge_df)
colnames(dge_df)[1] <- "Gene"

dge_mapped_df <- data.frame(
  entrez_id = mapIds(
    hgu133plus2.db,
    keys = dge_df$Gene,
    keytype = "PROBEID",
    column = "ENTREZID",
    multiVals = "first"
  )
) %>%
  dplyr::filter(!is.na(entrez_id)) %>%
  tibble::rownames_to_column("ProbeID") %>%
  dplyr::inner_join(dge_df, by=c("ProbeID" = "Gene"))
head(dge_mapped_df)

any(duplicated(dge_mapped_df$entrez_id))

dup_entrez_ids <- dge_mapped_df %>%
  dplyr::filter(duplicated(entrez_id)) %>%
  dplyr::pull(entrez_id)

dup_entrez_ids

dge_mapped_df %>%
  dplyr::filter(entrez_id %in% dup_entrez_ids)

filtered_dge_mapped_df <- dge_mapped_df %>%
  dplyr::arrange(dplyr::desc(abs(t))) %>%
  dplyr::distinct(entrez_id, .keep_all = TRUE)

filtered_dge_mapped_df %>%
  dplyr::filter(entrez_id %in% dup_entrez_ids)

t_vector <- filtered_dge_mapped_df$t
names(t_vector) <- filtered_dge_mapped_df$entrez_id

t_vector <- sort(t_vector, decreasing = TRUE)
head(t_vector)

set.seed(2022)

gsea_results <- GSEA(
  geneList = t_vector,
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    dr_hallmark_df,
    gs_name,
    entrez_gene
  )
)

head(gsea_results@result)

gsea_result_df <- data.frame(gsea_results@result)

gsea_result_df %>%
  dplyr::slice_max(n=3,order_by=NES)

most_positive_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_HYPOXIA",
  title = "HALLMARK_HYPOXIA",
  color.line = "#0d76ff"
)

most_positive_nes_plot

gsea_result_df %>%
  dplyr::slice_min(n=3, order_by = NES)

most_negative_nes_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_E2F_TARGETS",
  title = "HALLMARK_E2F_TARGETS",
  color.line = "#0d76ff"
)
most_negative_nes_plot

readr::write_tsv(
  gsea_result_df,
  file.path(
    results_dir,
    "GSE108474_gsea_results.tsv"
  )
)

sessioninfo::session_info()
library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]

library(oligo)
library(oligoClasses)

library(tidyverse)

filter <- colnames(gse108474)[gse108474@phenoData@data$`disease:ch1`=="glioblastoma multiforme" | gse108474@phenoData@data$`disease:ch1`=="normal"]
gse108474 <- gse108474[,filter]

pd <- pData(gse108474)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)

gse108474_celdata <- read.celfiles(pd$cel_file,phenoData=phenoData(gse108474))

pData(gse108474_celdata)[,c("geo_accession", "disease:ch1")]

exprs(gse108474_celdata)[1:4, 1:4]
exprs <- exprs(gse108474_celdata)
rownames(exprs(gse108474_celdata))

hist(gse108474_celdata,lwd=2,xlab='log intensity', which='pm', 
     main="CEL file densities before quantile normalization") 

oligo_normalized <- normalize(gse108474_celdata, method='quantile', which='pm')
hist(oligo_normalized, lwd=2, xlab='log intensity', which='pm', 
     main="CEL file densities after quantile normalization")

gse108474_eset <- rma(gse108474_celdata)

library(limma)
design <- model.matrix(~ `disease:ch1` + 0, data = pd)
colnames(design) <- levels(as.factor(gse108474_eset[['disease:ch1']]))

colnames(design)[1] <- "GBM"

fit <- lmFit(gse108474_eset, design)
fit <- eBayes(fit)

contrast_matrix <- makeContrasts(GBM - normal, levels=design)
contrast_matrix

contrasts_fit <- contrasts.fit(fit, contrast_matrix)

contrasts_fit <- eBayes(contrasts_fit)

stats_df <- topTable(contrasts_fit)

expression_df <- as.data.frame(exprs(gse108474_eset))
metadata <- pData(gse108474_eset)

top_gene_df <- expression_df[c("230765_at"),]
top_gene_df <- t(top_gene_df)
top_gene_df <- cbind(rownames(top_gene_df), top_gene_df)
colnames(top_gene_df)[1] = "geo_accession"
top_gene_df <- as.data.frame(top_gene_df)

top_gene_df <- merge(top_gene_df, metadata, by.x = "geo_accession")
rownames(top_gene_df) <- top_gene_df$geo_accession

library(hgu133plus2.db)
library(AnnotationDbi)
AnnotationDbi::select(hgu133plus2.db,rownames(stats_df),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID") 

## NWD2
ggplot(top_gene_df, aes(x = `disease:ch1`, y = `230765_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

summary(decideTests(contrasts_fit,lfc=1))

tops <- topTable(contrasts_fit, number=Inf)
upregulated <- tops[which(tops$logFC > 0), ][1:1814,]

write.csv(upregulated, "upregulated_genes_GBM.csv")

ps2 <- topTable(contrasts_fit,number=Inf,p.value=0.05,lfc=1)
ps2_up <- rownames(ps2)
df <- AnnotationDbi::select(hgu133plus2.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

df[df$PROBEID == "219752_at", ]

volcanoplot(contrasts_fit)
interesting_genes <- topTable(contrasts_fit,number=Inf,p.value=0.05,lfc=2)

upregulated <- interesting_genes[which(interesting_genes$logFC > 0), ]

write.csv(interesting_genes, "GBM_limma_results.csv")

volcanoplot(contrasts_fit,main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

eset_of_interest <- gse108474_eset[rownames(interesting_genes),]

expression_df <- as.data.frame(exprs(eset_of_interest))
metadata <- pData(eset_of_interest)

interesting_probes <- AnnotationDbi::select(hgu133plus2.db,rownames(interesting_genes),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

interesting_probes$SYMBOL

aamp = c("AAMP")
aamp_probes <- AnnotationDbi::select(hgu133plus2.db, aamp, c("GENENAME", "PROBEID", "ENTREZID"), keytype="SYMBOL")

aamp_df <- expression_df[c("201511_at"),]
aamp_df <- t(aamp_df)
aamp_df <- cbind(rownames(aamp_df), aamp_df)
colnames(aamp_df)[1] = "geo_accession"
aamp_df <- as.data.frame(aamp_df)

aamp_df <- merge(aamp_df, metadata, by.x = "geo_accession")
rownames(aamp_df) <- aamp_df$geo_accession

ggplot(aamp_df, aes(x = `disease:ch1`, y = `2201511_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

srpx = c("SRPX")
srpx_probes <- AnnotationDbi::select(hgu133plus2.db, srpx, c("GENENAME", "PROBEID", "ENTREZID"), keytype="SYMBOL")

srpx_df <- expression_df[c("204955_at"),]
srpx_df <- t(srpx_df)
srpx_df <- cbind(rownames(srpx_df), srpx_df)
colnames(srpx_df)[1] = "geo_accession"
srpx_df <- as.data.frame(srpx_df)

srpx_df <- merge(srpx_df, metadata, by.x = "geo_accession")
rownames(srpx_df) <- srpx_df$geo_accession

ggplot(srpx_df, aes(x = `disease:ch1`, y = `204955_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

pak1 = c("PAK1")
pak1_probes <- AnnotationDbi::select(hgu133plus2.db, pak1, c("GENENAME", "PROBEID", "ENTREZID"), keytype="SYMBOL")

pak1_df <- expression_df[c("209615_s_at", "226507_at", "230100_x_at"),]
pak1_df <- t(pak1_df)
pak1_df <- cbind(rownames(pak1_df), pak1_df)
colnames(pak1_df)[1] = "geo_accession"
pak1_df <- as.data.frame(pak1_df)

pak1_df <- merge(pak1_df, metadata, by.x = "geo_accession")
rownames(pak1_df) <- pak1_df$geo_accession

ggplot(pak1_df, aes(x = `disease:ch1`, y = `230100_x_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

common_genes = c("KCTD4", "KCNV1", "SPARC", "ARPP21", "AKAP5", "KCTD16", "CACNG3", "SST", "SLITRK4", "PPP4R4", "KCNJ6", "DLGAP2", "PAK1", "SLC35F3", "GABRA1", "GRM1", "SLC32A1", "MFSD4A", "CDC42", "HCN1", "NES", "SMIM10L2A", "NEGR1", "LRRC8B", "PAK5", "CD44", "EPHX4", "RAB3C", "CUX2", "MYRIP", "RBFOX1", "SNCA", "FRRS1L", "GNG3", "KCNMA1", "HPCAL4", "NPY1R", "HPCA", "RALYL", "ATP2B1", "RAP1GAP2", "SYNGR3", "GABRB3", "KCNQ5", "GALNT17", "CDH18", "HSPA12A", "GABRB2", "LPL", "CELF5", "SULT4A1", "AMPH", "CCK", "DNM1", "PHYHIP", "TMEM130", "DCLK1", "RAPGEF4", "CPLX2", "KCNQ2", "SRPX", "ZNF365", "GRM3", "GFAP", "DIRAS2", "RIMS3", "SNAP91", "NEFM", "CPLX1", "CAMK2B", "FXYD7", "TSPYL5", "SOX4", "TIMP4", "SLIT2")
  
common_geneNames <- AnnotationDbi::select(hgu133plus2.db,common_genes,c("GENENAME"),keytype="SYMBOL")

write.csv(common_geneNames, "common_genes_of_astrocytoma_and_GBM.csv")
write.csv(expression_df, "interesting_genes_expressions_for_GBM_vs_normal.csv")
write.csv(metadata, "interesting_metadata_for_GBM_vs_normal.csv")



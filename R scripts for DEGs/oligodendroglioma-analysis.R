library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]

library(oligo)
library(oligoClasses)

library(tidyverse)

filter <- colnames(gse108474)[gse108474@phenoData@data$`disease:ch1`=="oligodendroglioma" | gse108474@phenoData@data$`disease:ch1`=="normal"]
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

fit <- lmFit(gse108474_eset, design)
fit <- eBayes(fit)

contrast_matrix <- makeContrasts(oligodendroglioma - normal, levels=design)
contrast_matrix

contrasts_fit <- contrasts.fit(fit, contrast_matrix)

contrasts_fit <- eBayes(contrasts_fit)

stats_df <- topTable(contrasts_fit)

expression_df <- as.data.frame(exprs(gse108474_eset))
metadata <- pData(gse108474_eset)

top_gene_df <- expression_df[c("219752_at"),]
top_gene_df <- t(top_gene_df)
top_gene_df <- cbind(rownames(top_gene_df), top_gene_df)
colnames(top_gene_df)[1] = "geo_accession"
top_gene_df <- as.data.frame(top_gene_df)

top_gene_df <- merge(top_gene_df, metadata, by.x = "geo_accession")
rownames(top_gene_df) <- top_gene_df$geo_accession

library(hgu133plus2.db)
library(AnnotationDbi)
AnnotationDbi::select(hgu133plus2.db,rownames(stats_df),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID") 

## RASAL1
ggplot(top_gene_df, aes(x = `disease:ch1`, y = `219752_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

summary(decideTests(contrasts_fit,lfc=1))

ps2 <- topTable(contrasts_fit,number=Inf,p.value=0.05,lfc=1)
ps2_up <- rownames(ps2)
df <- AnnotationDbi::select(hgu133plus2.db,ps2_up,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

df[df$PROBEID == "219752_at", ]

volcanoplot(contrasts_fit)
interesting_genes <- topTable(contrasts_fit,number=Inf,p.value=0.05,lfc=2)

volcanoplot(contrasts_fit,main=sprintf("%d features pass our cutoffs",nrow(interesting_genes)))
points(interesting_genes[['logFC']],-log10(interesting_genes[['P.Value']]),col='red')

eset_of_interest <- gse108474_eset[rownames(interesting_genes),]

expression_df <- as.data.frame(exprs(eset_of_interest))
metadata <- pData(eset_of_interest)

interesting_probes <- AnnotationDbi::select(hgu133plus2.db,rownames(interesting_genes),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

interesting_probes$SYMBOL

common_genes = c("STYK1", "LY86-AS1", "KCTD4", "SPARC", "SEPTIN11", "HTR2C", "CACNG3", "LRRC7", "SLC32A1", "KCTD16", "CDC42", "STX1A", "PWWP3B", "KCNJ6", "MFSD4A", "RBP4", "CDKL5", "DLGAP2", "SYN2", "VANGL2", "ARPP21", "ERC2", "EPHX4", "MICAL2", "HCN1", "RBFOX1", "CRYM", "BSN", "GABRA1", "ARHGEF7", "HPCA", "FAM153B", "FAM153CP", "MEF2C", "LRFN5", "RAB3C", "PDYN", "GALNT17", "MCTP1", "NEGR1", "CPEB3", "LPL", "SULT4A1", "CREG2", "SLC8A1", "NEFL", "DCLK1", "SYNGR3", "SMIM10L2A", "SOX4", "MYRIP", "SLC12A5", "HPCAL4", "DNM1", "GRM3", "PGM2L1", "KCNJ3", "PRKCZ", "CAMK2B", "TSPYL5", "RIMS3", "VSNL1", "SSTR1", "ITPR1", "NEFM", "SERPINI1", "SRPX", "TAC1", "LINC01088", "ETV1", "NRGN", "DIRAS2", "SOX11", "RAPGEF4", "PEX5L", "HHIP", "SPOCK3")

common_geneNames <- AnnotationDbi::select(hgu133plus2.db,common_genes,c("GENENAME"),keytype="SYMBOL")

write.csv(common_geneNames, "common_genes_of_astrocytoma_and_oligodendroglioma.csv")
write.csv(expression_df, "interesting_genes_expressions_for_oligodendroglioma_vs_normal.csv")
write.csv(metadata, "interesting_metadata_for_oligodendroglioma_vs_normal.csv")



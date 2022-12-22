library(hgu133plus2.db)
ls('package:hgu133plus2.db')

# The reference genes we want to look at
ref_genes = c("HPRT1", "GAPDH", "TBP", "B2M", "RPL13A", "RN18S1")

# Grab the probes that map to the ref genes
library(AnnotationDbi)
df <- AnnotationDbi::select(hgu133plus2.db,ref_genes,c("PROBEID","ENTREZID","GENENAME"),keytype="SYMBOL")

# The actual probes
ref_gene_probes <- df['PROBEID']

# Load the annotation data from REMBRANDT
library(GEOquery)
gse108474 <- getGEO('GSE108474')
gse108474 <- gse108474[[1]]

# Filter for only astrocytoma and normal
library(tidyverse)
filter <- colnames(gse108474)[gse108474@phenoData@data$`disease:ch1`=="astrocytoma" | gse108474@phenoData@data$`disease:ch1`=="normal"]
gse108474 <- gse108474[,filter]

# Read in the raw data
library(oligo)
library(oligoClasses)
pd <- pData(gse108474)
pd['cel_file'] <- str_split(pd$supplementary_file,"/") %>% map_chr(tail,1)
gse108474_celdata <- read.celfiles(pd$cel_file,phenoData=phenoData(gse108474))

# Normalize using RMA
gse108474_eset <- rma(gse108474_celdata)

# Do differential gene expression between normal and astrocytoma
library(limma)
design <- model.matrix(~ `disease:ch1` + 0, data = pd)
colnames(design) <- levels(as.factor(gse108474_eset[['disease:ch1']]))

# Get only the expression for the reference genes
expression_df <- exprs(gse108474_eset)
idx <- rownames(expression_df) %in% ref_gene_probes$PROBEID
expression_df <- expression_df[idx,]

metadata <- pData(gse108474_eset)
expression_df <- as.data.frame(expression_df)

expression_df <- expression_df %>%
  dplyr::select(metadata$geo_accession)

fit <- lmFit(expression_df, design)
fit <- eBayes(fit)

contrast_matrix <- makeContrasts(astrocytoma - normal, levels=design)

contrasts_fit <- contrasts.fit(fit, contrast_matrix)
contrasts_fit <- eBayes(contrasts_fit)

stats_df <- topTable(contrasts_fit, number = nrow(expression_df)) # order the probes by most to least differentially expressed

# Look at the top probe expression and plot it
top_gene_df <- expression_df %>%
  dplyr::filter(rownames(.) == "202854_at") %>%
  t() %>%
  data.frame %>%
  tibble::rownames_to_column("geo_accession") %>%
  dplyr::inner_join(dplyr::select(
    metadata,
    geo_accession,
    `disease:ch1`
  ))

top_gene <- c("202854_at")
AnnotationDbi::select(hgu133plus2.db,top_gene,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID") # Match it to a gene

ggplot(top_gene_df, aes(x = `disease:ch1`, y = `X202854_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

# Look at the lowest probe expression and plot it
low_gene_df <- expression_df %>%
  dplyr::filter(rownames(.) == "1565446_at") %>%
  t() %>%
  data.frame %>%
  tibble::rownames_to_column("geo_accession") %>%
  dplyr::inner_join(dplyr::select(
    metadata,
    geo_accession,
    `disease:ch1`
  ))

low_gene <- c("1565446_at")
AnnotationDbi::select(hgu133plus2.db,low_gene,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID") # Match it to a gene

ggplot(low_gene_df, aes(x = `disease:ch1`, y = `X1565446_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()

# Put all the information in a single table
probe_ids <- rownames(stats_df)
probe_data <- AnnotationDbi::select(hgu133plus2.db,probe_ids,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

stats_df <- cbind(rownames(stats_df), stats_df)
colnames(stats_df)[1] = "PROBEID"

full_housekeeping <- merge(stats_df, probe_data, by.x="PROBEID") # full table (modified)



write.csv(full_housekeeping, "full_housekeeping.csv")

# ref_genes = c("A4GALT", "ACTB", "B2M", "CCK", "CSNK1G2", "DECR1", "DIMT1", 
#   "EEF1A1", "FARP1", "FPGS", "GAPDH", "GINS2", "GUSB", "HMBS", "HPRT1", 
#   "HSP90AB1", "MAPRE2", "PEX16", "PGK1", "POLR2A", "PPIA", "PPIB", 
#   "PUM1", "RPL4", "RPLP2", "SDHA", "SRSF4", "TBP", "TFRC", "TRAP1", 
#   "UBC", "YWHAG");

write.csv(df, file="genes_to_probes")

featureData(gse108474)[["ID"]]
idx <- featureData(gse108474)[["ID"]] %in% ref_gene_probes$PROBEID
gse108474_ref_genes <- gse108474[idx,]

write.csv(stats_df, "housekeeping_genes.csv")

AnnotationDbi::select(hgu133plus2.db,c("GAPDH"),c("PROBEID","ENTREZID","GENENAME"),keytype="SYMBOL")
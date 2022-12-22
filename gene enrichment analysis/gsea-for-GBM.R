library(clusterProfiler)
library(msigdbr)
library(magrittr)

dge_df <- as.data.frame(readr::read_csv("GBM_limma_results.csv"))
colnames(dge_df)[1] <- "Gene"

dr_hallmark_df <- msigdbr(
  species = "Homo sapiens",
)

library(hgu133plus2.db)
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
  dplyr::inner_join(dge_df, by = c("ProbeID" = "Gene"))

any(duplicated(dge_mapped_df$entrez_id))                            

dup_entrez_ids <- dge_mapped_df %>%
  dplyr::filter(duplicated(entrez_id)) %>%
  dplyr::pull(entrez_id)

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

gsea_result_df <- data.frame(gsea_results@result)

EMT <- gsea_result_df[which(gsea_result_df$ID=="HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),]
  
EMT$core_enrichment

df <- AnnotationDbi::select(hgu133plus2.db,c("6678", "3371", "2316", "960", "5396", "3485", "2335", "1284", "1282", "1278", "800", "3486", "2014", "1277", "1281", "4017", "7076", "7045", "7422", "5054", "7431", "5806", "4837", "1292", "4015", "10631", "4060", "10085"),c("PROBEID","SYMBOL","GENENAME"),keytype="ENTREZID") 

genenames <- AnnotationDbi::select(hgu133plus2.db, c("6678", "3371", "2316", "960", "5396", "3485", "2335", "1284", "1282", "1278", "800", "3486", "2014", "1277", "1281", "4017", "7076", "7045", "7422", "5054", "7431", "5806", "4837", "1292", "4015", "10631", "4060", "10085"), c("SYMBOL", "GENENAME"), keytype="ENTREZID")

genes <- gsea_results@result$core_enrichment

enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  title = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  color.lin = "#0d76fe"
)

data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]

entrez_ids <- filtered_dge_mapped_df$entrez_id
ggo <- groupGO(
  gene = entrez_ids,
  OrgDb = hgu133plus2.db,
  ont = "CC",
  level = 3,
  readable = TRUE
)                            
head(ggo)

ego <- enrichGO(
  gene = entrez_ids,
  OrgDb = hgu133plus2.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)
head(ego)  

gene.df <- bitr(entrez_ids, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = hgu133plus2.db)  

goplot(ego)

AnnotationDbi::select(hgu133plus2.db,c("200665_s_at"),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
up_rows <- rownames(upregulated)

up_rows <- AnnotationDbi::select(hgu133plus2.db,up_rows,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")

up_rows <- cbind(AnnotationDbi::select(hgu133plus2.db,up_rows,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID"), up_rows)

up_rows$SYMBOL

sparc_df <- expression_df[c("200665_s_at"),]
sparc_df <- t(sparc_df)
sparc_df <- cbind(rownames(sparc_df), sparc_df)
colnames(sparc_df)[1] = "geo_accession"
sparc_df <- as.data.frame(sparc_df)

sparc_df <- merge(sparc_df, metadata, by.x = "geo_accession")
rownames(sparc_df) <- sparc_df$geo_accession

ggplot(sparc_df, aes(x = `disease:ch1`, y = `200665_s_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()
  
vegfa = c("VEGFA")
vegfa_probes <- AnnotationDbi::select(hgu133plus2.db, vegfa, c("GENENAME", "PROBEID", "ENTREZID"), keytype="SYMBOL")

vegfa_df <- expression_df[c("210512_s_at", "210513_s_at", "211527_x_at", "212171_x_at"),]
vegfa_df <- t(vegfa_df)
vegfa_df <- cbind(rownames(vegfa_df), vegfa_df)
colnames(vegfa_df)[1] = "geo_accession"
vegfa_df <- as.data.frame(vegfa_df)

vegfa_df <- merge(vegfa_df, metadata, by.x = "geo_accession")
rownames(vegfa_df) <- vegfa_df$geo_accession

ggplot(vegfa_df, aes(x = `disease:ch1`, y = `212171_x_at`, color = `disease:ch1`)) +
  geom_jitter(width = 0.2, height = 0) +
  theme_classic()
  
  
  
  
  
  
library(clusterProfiler)
library(msigdbr)
library(magrittr)


dge_df <- as.data.frame(readr::read_csv("GBM_limma_results.csv"))
colnames(dge_df)[1] <- "Gene"


dr_hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")


library(hgu133plus2.db)
dge_mapped_df <- data.frame(entrez_id = mapIds(hgu133plus2.db, keys = dge_df$Gene, keytype = "PROBEID", column = "ENTREZID", multiVals = "first")) %>%
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


set.seed(2022)


gsea_results <- GSEA(geneList = t_vector, minGSSize = 25, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, seed = TRUE, pAdjustMethod = "BH", TERM2GENE = dplyr::select(dr_hallmark_df, gs_name, entrez_gene))


gsea_result_df <- data.frame(gsea_results@result)
enrichplot::gseaplot(gsea_results, geneSetID = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", title = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", color.lin = "#0d76fe")


data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]


entrez_ids <- filtered_dge_mapped_df$entrez_id
ggo <- groupGO(gene = entrez_ids, OrgDb = hgu133plus2.db, ont = "CC", level = 3, readable = TRUE)                           


ego <- enrichGO(gene = entrez_ids, OrgDb = hgu133plus2.db, ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05, readable = TRUE)


gene.df <- bitr(entrez_ids, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"),
               OrgDb = hgu133plus2.db) 


goplot(ego)


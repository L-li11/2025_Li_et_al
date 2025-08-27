setwd("~/Desktop/transdiff_GRN")

library(dplyr)
library(org.Mm.eg.db)
library(clusterProfiler)

gene_ids <- read.csv("data/GSE208199/Mm_ID.csv", sep = ",")
gene_ids <- data.frame(gene_ids)
gene_ids <- gene_ids[, c(2, 3)]

mdl <- read.csv("results/GSE208199/rank_community_Ngn2n_t5.txt", sep = "\t")
mdl <- unique(data.frame(mdl))
entrezmdl <- left_join(mdl, gene_ids, by = c("geneid" = "Gene.name"))
entrezmdl <- na.omit(entrezmdl)
entrezmdl$thmdls <- paste(entrezmdl$run, entrezmdl$nodes,
                          entrezmdl$community, sep = "_")

mdls <- unique(entrezmdl$thmdls)

for (m in mdls){
  mdlid <- entrezmdl[entrezmdl$thmdls == m, ]$NCBI.gene..formerly.Entrezgene..ID
  if (length(mdlid) >= 20) {
    result <- enrichGO(mdlid, OrgDb = "org.Mm.eg.db",
                       keyType = "ENTREZID",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       qvalueCutoff = 0.05,
                       minGSSize = 10,
                       maxGSSize = 2000,
                       readable = FALSE)
    resultframe <- data.frame(result)
    if (nrow(resultframe) >= 1) {
        rownames(resultframe) <- NULL
        write.csv(resultframe, row.names = FALSE,
                  paste("results/GSE208199/GO/Ngn2n_", as.character(m), "_t5.csv", sep = "")
                  )
    }
  }
}

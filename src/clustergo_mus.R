setwd("~/Desktop/transdiff_GRN")

library(dplyr)
library(org.Mm.eg.db)
library(clusterProfiler)

gene_ids <- read.csv("data/GSE208199/Mm_ID.csv", sep = ",")
gene_ids <- data.frame(gene_ids)
gene_ids <- gene_ids[, c(2, 3)]

mdl <- read.csv("results/temp/community.csv")
entrezmdl <- left_join(mdl, gene_ids, by = c("X" = "Gene.name"))
entrezmdl <- na.omit(entrezmdl)
mdls <- unique(entrezmdl$community)

for (m in mdls){
  mdlid <- entrezmdl[entrezmdl$community == m, ]$NCBI.gene..formerly.Entrezgene..ID
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
    } else {
      resultframe <- data.frame(0)
    }
  } else {
    resultframe <- data.frame(0)
  }
  write.csv(resultframe, row.names = FALSE,
            paste("results/temp/", as.character(m), ".csv", sep = ""))
}
write.csv(entrezmdl, row.names = FALSE, paste("results/temp/community.csv"))
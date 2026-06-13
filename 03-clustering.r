#“Resolution was selected using a composite metric integrating cluster stability (ARI), 
# sample mixing (LISI), entropy, and cluster size constraints.”
# LISI : = 1/ sum_i(1/p_i^2) p_i : probability of neighbors from batch b


input_info="data/sample_input.csv"
input="bigdata/seurat5/s13_integrated.rds"
output="bigdata/seurat5/s13_clustered.rds"

library(Seurat)
d=readRDS(input)


d <- FindNeighbors(d, reduction = "integrated.cca", dims = 1:30)

results <- list()

for (res in seq(0.2, 1.2, by = 0.2)) {
  d <- FindClusters(d, resolution = res, verbose = FALSE)
  cluster_col <- paste0("RNA_snn_res.", res)
  clusters <- d@meta.data[[cluster_col]]
  emb <- Embeddings(d, "umap")
  samples <- d$orig.ident

  library(FNN)
  k <- 30
  nn <- get.knn(emb, k = k)

  lisi_scores <- sapply(1:nrow(emb), function(i) {
    neigh_samples <- samples[nn$nn.index[i, ]]
    p <- table(neigh_samples) / k
    1 / sum(p^2)
  })

  d$lisi_tmp <- lisi_scores
  meta <- d@meta.data
  meta$cluster_tmp <- clusters

  library(dplyr)
  summary_stats <- meta %>%
    group_by(cluster_tmp) %>%
    summarise( mean_lisi = mean(lisi_tmp), entropy = {
        p <- table(orig.ident) / n(); -sum(p * log(p + 1e-10))
      }, n_cells = n(), .groups = "drop")

  results[[as.character(res)]] <- list( n_clusters = length(unique(clusters)), stats = summary_stats)
}


resolution_summary <- data.frame()

for (res in names(results)) {
  stats <- results[[res]]$stats
  resolution_summary <- rbind( resolution_summary,
    data.frame( resolution = as.numeric(res), n_clusters = results[[res]]$n_clusters, mean_lisi = mean(stats$mean_lisi),
      mean_entropy = mean(stats$entropy), min_cluster_size = min(stats$n_cells)))
}
resolution_summary
#resolution_summary
#  resolution n_clusters mean_lisi mean_entropy min_cluster_size
#1        0.2         12  4.393888     1.915863              136
#2        0.4         19  4.133513     1.804466               76
#3        0.6         22  4.407481     1.893463              104
#4        0.8         26  4.411981     1.866129              104
#5        1.0         32  4.296710     1.803040               76
#6        1.2         35  4.298096     1.805409               77
# Choose resolution where:
#- LISI is near maximum
#- Entropy not decreasing
#- Minimum cluster size not collapsing
#- Before cluster number accelerates

library(clustree)
p=clustree(d@meta.data, prefix = "RNA_snn_res.")
ggsave(p,file="figures/clustree.svg");

best_res=0.6
d <- FindClusters(d, resolution = best_res, verbose = FALSE)
saveRDS(d,file=output)

p=DimPlot(d, reduction = "umap", group.by = c("seurat_clusters"),label=T,repel=T) + theme_classic()
ggsave(p,file="figures/umap_cluster.svg",width=7,height=7)

## redo after celltype identification
#library(mclust)
#res_cols <- grep("RNA_snn_res", colnames(d@meta.data), value = TRUE)
#ari_matrix <- matrix(NA, nrow = length(res_cols), ncol = length(res_cols))
#rownames(ari_matrix) <- res_cols
#colnames(ari_matrix) <- res_cols
#
#for (i in 1:length(res_cols)) {
#  for (j in 1:length(res_cols)) {
#    ari_matrix[i, j] <- adjustedRandIndex(
#      d@meta.data[[res_cols[i]]],
#      d@meta.data[[res_cols[j]]]
#    )
#  }
#}
#
#ari_matrix
#

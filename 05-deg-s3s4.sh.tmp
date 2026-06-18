library(data.table)
library(Seurat)
library(edgeR)

input="bigdata/s13_celltyped.rds")
d=readRDS(input)
d1=subset(d,seurat_clusters %in% c("s3","s4"))

results <- list()
Idents(d1) <- "seurat_clusters"
for (cl in levels(d1)) {
  obj <- subset(d1, idents = cl)
  tab <- table(obj$selfrenewal)
  if (all(c("high_selfrenewal", "low_selfrenewal") %in% names(tab)) && all(tab[c("high_selfrenewal", "low_selfrenewal")] >= 3)) {
    Idents(obj) <- obj$selfrenewal
    results[[cl]] <- FindMarkers( obj, ident.1 = "high_selfrenewal", ident.2 = "low_selfrenewal", test.use = "wilcox")
  }
}

all_res <- rbindlist(lapply(names(results), function(cl) {
  res <- as.data.table(results[[cl]], keep.rownames = "gene")
  res[, cluster := cl]
  setnames(res, old = c("avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2"),
    new = c("logFC", "pval", "padj", "pct_good", "pct_bad"))
  res
}), fill = TRUE)
setorder(all_res, padj)
fwrite(all_res, "data/s3s4_selfrenewal_DE.csv")
all_res=fread("data/s3s4_selfrenewal_DE.csv")

j=c(0, 1, 2, 14, 15, 11, 17)
res_filt <- all_res[ cluster %in% j & padj < 0.01 & abs(logFC) > 1.5 & pmax(pct_good,pct_bad) > 0.5 ]
gene_score <- res_filt[, .(score = mean(abs(logFC),na.rm = TRUE)), by = gene]

core_genes <- gene_score[order(-score)][1:200, gene]
targ_genes = c("IER2","ZFP36","FOSB","BTG2","CEBPD","CD74")
intersect(core_genes,targ_genes)

## include non-sig but core expression per cluster
mat_dt <- all_res[cluster %in% j & gene %in% core_genes, .(logFC = mean(logFC, na.rm = TRUE)), by = .(gene, cluster)]
mat <- dcast(mat_dt, cluster ~ gene, value.var = "logFC", fill = 0)

n <- mat$cluster; mat$cluster <- NULL
mat <- as.matrix(mat)
rownames(mat) <- n
colors <- ifelse(colnames(mat) %in% targ_genes, "red", "black")
library(ComplexHeatmap)
svg("figures/heatmap_core_genes_epionly.svg", width = 30, height = 5)
Heatmap(mat,column_names_gp = gpar(col=colors))
dev.off()

pdf("figures/heatmap_core_genes_epionly.pdf", width = 30, height = 5)
Heatmap(mat,column_names_gp = gpar(col=colors))
dev.off()

fwrite(all_res[gene%in% core_genes & cluster %in% j], "data/s3s4_selfrenewal_DE_core_epi.csv")



dt <- copy(all_res[ cluster %in% j] )
mat_dt <- dt[gene %in% core_genes, .(logFC = mean(logFC, na.rm = TRUE)), by = .(gene, cluster)]
mat <- dcast(mat_dt, gene ~ cluster, value.var = "logFC", fill = 0)

gene_names <- mat$gene
mat$gene <- NULL
mat <- as.matrix(mat)
rownames(mat) <- gene_names
library(ComplexHeatmap)
svg("figures/heatmap_exclusive_genes_epionly.svg", width = 30, height = 5)
Heatmap(t(mat))
dev.off()

pdf("figures/heatmap_exclusive_genes_epionly.pdf", width = 30, height = 5)
Heatmap(t(mat))
dev.off()






## gprofiler
library(gprofiler2)
res <- gost( query = core_genes, organism = "hsapiens", correction_method = "fdr")

sig <- df[df$p_value < 0.05 & df$intersection_size >= 3, ]

res_filt[, direction := ifelse(logFC > 0, "good", "bad")]
core_dt <- res_filt[gene %in% core_genes]

library(gprofiler2)
gp_core <- lapply(split(core_dt$gene, core_dt$direction), function(g) {
  gost( query = unique(g), organism = "hsapiens", correction_method = "fdr", significant = TRUE)
})

gp_dt <- rbindlist(lapply(names(gp_core), function(nm) {
  gp <- gp_core[[nm]]
  if (is.null(gp$result)) return(NULL)
  data.table( direction = nm, pathway = gp$result$term_name, source = gp$result$source, pval = gp$result$p_value)
}), fill = TRUE)
gp_top <- gp_dt[ order(pval), head(.SD, 20), by = direction ]
library(ggplot2)

ggplot(gp_top, aes(x = direction, y = pathway)) +
  geom_point(aes( size = -log10(pval), color = -log10(pval))) +
  scale_color_viridis_c() +
  theme_classic() +
  labs( x = "Repair direction", y = "Core pathways")


top_genes <- all_res[gene %in% core_genes][order(padj), head(.SD, 20), by = cluster]

top_markers <- top_genes[ order(-abs(logFC)), head(.SD, 5), by = cluster ]
library(ggplot2)
ggplot(top_genes, aes(x = cluster, y = gene, color = logFC)) +
  geom_point(aes(size = pct_good)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red") +
  theme_classic()



### gprofiler2
library(gprofiler2)
library(data.table)

genes <- unique(top_genes[sig == TRUE]$gene)
gp <- gost( query = genes, organism = "hsapiens", correction_method = "fdr", significant = TRUE)

res <- gp$result

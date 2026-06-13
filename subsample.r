library(Seurat)
library(data.table)
input="bigdata/s13_celltyped.rds"
#tt=readRDS("bigdata/s13_clustered.rds")
#d$orig.ident <- tt$orig.ident[names(d$orig.ident)]
#saveRDS(d,file="bigdata/s13_celltyped.rds")
d=readRDS(input)
d1 <- subset(d, subset = orig.ident %in% c("s3", "s4"))

library(ggplot2)
p <- DimPlot(d1, reduction="umap",group.by = "seurat_clusters", split.by = "repairment",label=T)
ggsave("figures/umap_s3s4.svg", plot = p , width = 12 ,height = 6) 
ggsave("figures/umap_s3s4.pdf", plot = p , width = 12, height = 6)

p <- DimPlot(d1, reduction="umap", group.by = "celltype", split.by = "repairment",label=T)
ggsave("figures/umap_s3s4_celltype.svg", plot = p , width = 12 ,height = 6) 
ggsave("figures/umap_s3s4_celltype.pdf", plot = p , width = 12, height = 6)


library(Seurat)
library(ggplot2)
library(cowplot)

features= unlist(strsplit(m$gene,","))
a=VlnPlot(d[d$repairement=="good",], features,group.by="seurat_clusters", stack = TRUE, sort = TRUE, flip = F) +
     theme(legend.position = "none") + ggtitle("Target Genes by Cluster")

b=VlnPlot(d[d$repairment=="bad", features,group.by="seurat_clusters", stack = TRUE, sort = TRUE, flip = F) +
     theme(legend.position = "none") + ggtitle("Target Genes by Cluster")


library(Seurat)
library(ggplot2)

# Subset good repairment
a <- VlnPlot( subset(d, subset = repairment == "good_repair"),
  features = features, group.by = "seurat_clusters", stack = TRUE, flip = FALSE
) + theme(legend.position = "none") + ggtitle("Good Repairment")

# Subset bad repairment
b <- VlnPlot( subset(d, subset = repairment == "bad_repair"),
  features = features, group.by = "seurat_clusters", stack = TRUE, flip = FALSE
) + theme(legend.position = "none") + ggtitle("Bad Repairment")

ggsave(a+b,file="figures/vln_target.svg")

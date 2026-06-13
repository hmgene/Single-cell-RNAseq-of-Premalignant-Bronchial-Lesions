library(Seurat)
library(data.table)
input="bigdata/s13_celltyped.rds"

d=readRDS(input)
m=fread("data/target.tsv")

features= unlist(strsplit(m$gene,","))
p <- VlnPlot( d, features = features, group.by = "seurat_clusters", split.by = "repairment",  stack = TRUE, sort = TRUE, flip = FALSE) +
  ggtitle("Target Genes by Cluster (Good vs Bad Repair)")

ggsave(p,file="figures/vln_target.svg")
ggsave(p,file="figures/vln_target.pdf")

d1=subset(d, orig.ident %in% c("s3","s4"))
p <- VlnPlot( d1, features = features, group.by = "seurat_clusters", split.by = "repairment",  stack = TRUE, sort = TRUE, flip = FALSE) +
  ggtitle("Target Genes by Cluster (Good vs Bad Repair)")
ggsave(p,file="figures/vln_target_s3s4.svg")
ggsave(p,file="figures/vln_target_s3s4.pdf")


j=c(0, 1, 2, 14, 15, 11, 17)
d2=subset(d1, orig.ident %in% c("s3","s4") & seurat_clusters %in% j)

p <- VlnPlot( d2, features = features, group.by = "seurat_clusters", split.by = "repairment",  stack = TRUE, sort = TRUE, flip = FALSE) +
  ggtitle(paste0("Target Genes by Cluster,",paste(j,collapse=" ")," (Good vs Bad Repair)"))

ggsave(p,file="figures/vln_target_s3s4_epi.svg")
ggsave(p,file="figures/vln_target_s3s4_epi.pdf")


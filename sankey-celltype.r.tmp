
library(data.table)
library(Seurat)
m=fread("data/s13-meta.data.csv",); rownames(m)=m[,V1]; m[,V1:=NULL];
tt=readRDS("bigdata/s13-harmony.rds")

m=m[,cluster0 := tt@meta.data[row.names(m),]$seurat_clusters]

links1 <- m[, .N, by = .(source = orig.ident, target = seurat_clusters)]
links2 <- m[, .N, by = .(source = seurat_clusters, target = celltype)]
## sankey write
links2[order(target),
       cat(paste0( source, " [", N, "] ", target, "\n"))]
links1[order(source),
       cat(paste0( source, " [", N, "] ", target, "\n"))]




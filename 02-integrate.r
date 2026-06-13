input="bigdata/seurat5/s13_merged.rds"
input_info="data/sample_input.csv"
output="bigdata/seurat5/s13_integrated.rds"


library(data.table)
d=readRDS(input);
m=fread(input_info)

d = subset(d, subset = nFeature_RNA > 200 & percent.mt < 25 )
#dim(d)
#[1] 45068 48347
# table(d$orig.ident)
#  s1  s10  s11  s12  s13   s2   s3   s4   s5   s6   s7   s8   s9 
# 495 3877 9024 1600 2202 1063 1431 2740 1707 1192 8167 6673 8176 


#d[["RNA"]] <- split(d[["RNA"]], f = d$orig.ident)
d <- NormalizeData(d)
d <- FindVariableFeatures(d)
d <- ScaleData(d)
d <- RunPCA(d)
#d <- FindNeighbors(d, dims = 1:30, reduction = "pca")
#d <- FindClusters(d, resolution = 0.4, cluster.name = "unintegrated_clusters")
d = RunUMAP(d,dims=1:30,reduction = "pca", reduction.name = "umap.unintegrated")
library(ggplot2)
p=DimPlot(d, reduction = "umap.unintegrated", group.by = c("orig.ident", "repairment"))
ggsave(p,file="figures/umap_before_integration.svg",width=14,height=7)


#
# integration
# 
d <- IntegrateLayers(object = d, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = FALSE)
d[["RNA"]] <- JoinLayers(d[["RNA"]])
#d <- FindNeighbors(d, reduction = "integrated.cca", dims = 1:30)
#d <- FindClusters(d, resolution = 0.8)
d=RunUMAP(d, dims = 1:30, reduction = "integrated.cca")
p=DimPlot(d, reduction = "umap", group.by = c("orig.ident", "repairment"))
ggsave(p,file="figures/umap_after_integration.svg",width=14,height=7)

saveRDS(d,file=output)


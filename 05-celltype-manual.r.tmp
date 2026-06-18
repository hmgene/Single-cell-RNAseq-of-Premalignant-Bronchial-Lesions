input="bigdata/s13_celltyped.rds"
input_info="data/sample_input.csv"
input_marker="data/Marker Genes for Annotation.xlsx";
input_pathology="data/Pathology and sample ID.xlsx";
output="bigdata/s13_celltyped.rds"

library(Seurat)
library(Matrix)
library(data.table)
library(ggplot2)
library(readxl)

d=readRDS(input);


d$celltype <- "Unknown"

# Epithelial
d$celltype[d$seurat_clusters %in% c(0,1,2,14)] <- "Basal"
d$celltype[d$seurat_clusters == 15] <- "Proliferating"
d$celltype[d$seurat_clusters == 11] <- "Secretory"
d$celltype[d$seurat_clusters == 7] <- "Ciliated"
d$celltype[d$seurat_clusters == 17] <- "Deuterosomal"

# Immune
d$celltype[d$seurat_clusters == 3]  <- "Plasma"
d$celltype[d$seurat_clusters == 18] <- "B cells"
d$celltype[d$seurat_clusters == 8]  <- "CD4 T"
d$celltype[d$seurat_clusters == 9]  <- "CD8 T"
d$celltype[d$seurat_clusters == 16] <- "Macrophage"
d$celltype[d$seurat_clusters == 20] <- "cDC1"
d$celltype[d$seurat_clusters == 21] <- "cDC2"

# Endothelial
d$celltype[d$seurat_clusters %in% c(19,4,6)] <- "Endothelial"

# Fibroblasts
d$celltype[d$seurat_clusters %in% c(5,10)] <- "Fibroblast"

# VSMC / Pericyte
d$celltype[d$seurat_clusters %in% c(12,13)] <- "VSMC/Pericyte"
celltype_colors <- c(
  # Epithelial (orange tones)
  "Basal" = "#f4a261",
  "Proliferating" = "#f4a261",
  "Secretory" = "#f4a261",
  "Ciliated" = "#f4a261",
  "Deuterosomal" = "#f4a261",

  # Immune (green tones)
  "Plasma" = "#2ecc71",
  "B cells" = "#2ecc71",
  "CD4 T" = "#2ecc71",
  "CD8 T" = "#2ecc71",
  "Macrophage" = "#2ecc71",
  "cDC1" = "#2ecc71",
  "cDC2" = "#2ecc71",

  # Endothelial (purple/pink)
  "Endothelial" = "#e377c2",

  # Fibroblast (yellow)
  "Fibroblast" = "#f1c40f",

  # VSMC / Pericyte (red)
  "VSMC/Pericyte" = "#e74c3c",

  "Unknown" = "grey80"
)

d$celltype <- factor(d$celltype, levels = names(celltype_colors))


p=DimPlot(d,group.by="celltype",reduction="umap",label = TRUE) + ggtitle("UMAP by Cell Type")
ggsave(p,file="figures/umap_celltype.svg")
ggsave(p,file="figures/umap_celltype.pdf")

### convex hall

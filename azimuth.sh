library(Seurat)
library(Azimuth)
setwd("/w")
d=readRDS("bigdata/s13_clustered.rds")
r=LoadReference("./") #ref.Rds");

d <- SCTransform(
  object = d, assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = r$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi', do.correct.umi = FALSE, do.scale = FALSE, do.center = TRUE)


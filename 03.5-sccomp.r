library(Seurat)
library(tidyr)
library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)
library(tibble)
library(data.table)
library(cmdstanr)
cmdstanr::set_cmdstan_path("~/.cmdstan/cmdstan-2.38.0/")


data("seurat_obj")
data("sce_obj")
data("counts_obj")

input="bigdata/s13_clustered.rds"
tt=readRDS(input)

p=fread("data/pathology.csv")
library(data.table)
tt@meta.data$dysplasia <- p$description[match(tt@meta.data$orig.ident, p$id)]
tt@meta.data$patient <- p$patient[match(tt@meta.data$orig.ident, p$id)]

x <- as.data.table(tt@meta.data)
x <- x[, .N, by = .(seurat_clusters,patient,selfrenewal,dysplasia,orig.ident)]
setnames(x, "N", "count")
setnames(x, "orig.ident", "sample")
x[, selfrenewal := factor(selfrenewal, levels = c("low_selfrenewal","high_selfrenewal"))]
x[, dysplasia := factor(dysplasia, levels = c("Non-dysplastic", "Dysplastic"))]

fit <- x |> sccomp_estimate( formula_composition = ~ selfrenewal * dysplasia + ( 1 | patient),
  .sample = "sample", .cell_group = "seurat_clusters", abundance = "count", variational_inference = F,core=8) 
res= fit |> sccomp_test()
plots =  res |> plot()
ggsave("figures/credible_interval_all.svg", plots$credible_intervals_1D)

ggsave("figures/credible_interval_selfrenewal.svg",
plots$credible_intervals_1D[[1]] + labs( title = "Repair State Credible Intervals", y = "Cluster ID") + theme_minimal()  )

ggsave("figures/credible_interval_dysplasia.svg",
plots$credible_intervals_1D[[2]] + labs( title = "Dysplasia State Credible Intervals", y = "Cluster ID") + theme_minimal()  )

# sanity check
library(ggplot2)
p=ggplot(x, aes(x = selfrenewal, y = count, fill = dysplasia)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  facet_wrap(~ seurat_clusters)
ggsave("figures/barplot_selfrenewal_dysplasia.svg",p)

y=subset(tt,subset= tt$seurat_clusters==4 & tt$dysplasia=="Dysplastic" & tt$selfrenewal=="high_selfrenewal")
table(y$orig.ident)
#s4 
#350 

library(data.table)
library(pheatmap)

x[selfrenewal == "high_selfrenewal", selfrenewal := "high_selfrenewal"]
x[selfrenewal == "low_selfrenewal",  selfrenewal := "low_selfrenewal"]
x[, selfrenewal := factor(selfrenewal)]
dt <- as.data.table(x)
x[, state := paste(selfrenewal, dysplasia, sep = " / ")]

mat_dt <- dt[ , .(count = sum(count)), by = .(patient, state) ]
mat <- dcast( mat_dt, state ~ patient, value.var = "count", fill = 0)

mat <- as.data.frame(mat)
rownames(mat) <- mat$state
mat$state <- NULL
mat <- as.matrix(mat)


pheatmap( mat, color = colorRampPalette(c("white", "red"))(100),
  border_color = NA, cluster_rows = TRUE, cluster_cols = FALSE,
  main = "Patient-wise counts of 4 classes", filename = "figures/heatmap_state_patient.pdf",
  width = 6, height = 3)


library(data.table)
library(pheatmap)

dt_agg <- dt[ , .(count = sum(count)), by = .(seurat_clusters, state) ]
mat <- dcast( dt_agg, state ~ seurat_clusters, value.var = "count", fill = 0)

mat <- as.data.frame(mat)
rownames(mat) <- mat$state
mat$state <- NULL
mat <- as.matrix(mat)
cluster_order <- sort(as.numeric(colnames(mat)))
mat <- mat[, as.character(cluster_order), drop = FALSE]

pheatmap( mat, color = colorRampPalette(c("white", "red"))(100), filename="figures/heatmap_state_cluster.pdf",
  border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE, main = "States × Clusters (counts)",
  width = 8, height = 3)

pheatmap( mat, color = colorRampPalette(c("white", "red"))(100), #filename="figures/heatmap_state_cluster.pdf",
  border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE, main = "States × Clusters (counts)",
  width = 8, height = 3)


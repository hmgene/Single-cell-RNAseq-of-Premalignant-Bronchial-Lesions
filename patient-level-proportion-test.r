library(Seurat)
library(dplyr)
library(sccomp)
library(ggplot2)
library(forcats)
library(tidyr)
library(data.table)
data("seurat_obj")
data("sce_obj")
data("counts_obj")

tt=readRDS("bigdata/s13-integrate-by-selfrenewal.rds")

#Perform differential composition analysis (e.g., sccomp), adjusting for histology and smoking.
p=fread("data/pathology.csv")

library(data.table)
md <- as.data.table(tt@meta.data)
counts <- md[, .N, by = .(seurat_clusters, selfrenewal,orig.ident)]
setnames(counts, "N", "count")

sccomp_result <- sccomp_glm( counts[seurat_clusters==0],
  formula_composition = ~ selfrenewal,
  .sample = sample,
  .cell_group = seurat_clusters,
  .count = count
) |> test_contrasts("low_selfrenewal - high_selfrenewal ") 


perm_test_cluster_patient <- function(dt_cluster, n_perm = 10000) {
  # patient-level totals
  pt <- dt_cluster[
    , .(count = sum(count)),
    by = .(orig.ident, selfrenewal)
  ]

  obs <- with(pt,
    sum(count[selfrenewal == "low_selfrenewal"]) -
    sum(count[selfrenewal == "high_selfrenewal"])
  )

  perm_stats <- replicate(n_perm, {
    perm_rep <- sample(pt$selfrenewal)
    sum(pt$count[perm_rep == "low_selfrenewal"]) -
      sum(pt$count[perm_rep == "high_selfrenewal"])
  })

  mean(abs(perm_stats) >= abs(obs))
}

prop_results <- dt_prop[
  , {
      x <- prop[selfrenewal == "low_selfrenewal"]
      y <- prop[selfrenewal == "high_selfrenewal"]

      if (length(x) >= 2 && length(y) >= 2) {
        list(p_value = wilcox.test(x, y)$p.value)
      } else {
        list(p_value = NA_real_)
      }
    },
  by = seurat_clusters
]

perm_test_cluster_patient <- function(dt_cluster, n_perm = 10000) {

  # patient-level totals
  pt <- dt_cluster[
    , .(count = sum(count)),
    by = .(orig.ident, selfrenewal)
  ]

  obs <- with(pt,
    sum(count[selfrenewal == "low_selfrenewal"]) -
    sum(count[selfrenewal == "high_selfrenewal"])
  )

  perm_stats <- replicate(n_perm, {
    perm_rep <- sample(pt$selfrenewal)
    sum(pt$count[perm_rep == "low_selfrenewal"]) -
      sum(pt$count[perm_rep == "high_selfrenewal"])
  })

  mean(abs(perm_stats) >= abs(obs))
}

results <- counts[
  ,
  {
    pval <- perm_test_cluster_patient(.SD, n_perm = 10000)

    list(
      p_value = pval,
      bad_total  = sum(count[selfrenewal == "low_selfrenewal"]),
      good_total = sum(count[selfrenewal == "high_selfrenewal"])
    )
  },
  by = seurat_clusters
]

fwrite(results,"data/cluster_composition_test_result.csv")


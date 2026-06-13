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

tt=readRDS("bigdata/s13-integrate-by-repairment.rds")

#Perform differential composition analysis (e.g., sccomp), adjusting for histology and smoking.
p=fread("data/pathology.csv")

library(data.table)
md <- as.data.table(tt@meta.data)
counts <- md[, .N, by = .(seurat_clusters, repairment,orig.ident)]
setnames(counts, "N", "count")

sccomp_result <- sccomp_glm( counts[seurat_clusters==0],
  formula_composition = ~ repairment,
  .sample = sample,
  .cell_group = seurat_clusters,
  .count = count
) |> test_contrasts("bad_repair - good_repair ") 


perm_test_cluster_patient <- function(dt_cluster, n_perm = 10000) {
  # patient-level totals
  pt <- dt_cluster[
    , .(count = sum(count)),
    by = .(orig.ident, repairment)
  ]

  obs <- with(pt,
    sum(count[repairment == "bad_repair"]) -
    sum(count[repairment == "good_repair"])
  )

  perm_stats <- replicate(n_perm, {
    perm_rep <- sample(pt$repairment)
    sum(pt$count[perm_rep == "bad_repair"]) -
      sum(pt$count[perm_rep == "good_repair"])
  })

  mean(abs(perm_stats) >= abs(obs))
}

prop_results <- dt_prop[
  , {
      x <- prop[repairment == "bad_repair"]
      y <- prop[repairment == "good_repair"]

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
    by = .(orig.ident, repairment)
  ]

  obs <- with(pt,
    sum(count[repairment == "bad_repair"]) -
    sum(count[repairment == "good_repair"])
  )

  perm_stats <- replicate(n_perm, {
    perm_rep <- sample(pt$repairment)
    sum(pt$count[perm_rep == "bad_repair"]) -
      sum(pt$count[perm_rep == "good_repair"])
  })

  mean(abs(perm_stats) >= abs(obs))
}

results <- counts[
  ,
  {
    pval <- perm_test_cluster_patient(.SD, n_perm = 10000)

    list(
      p_value = pval,
      bad_total  = sum(count[repairment == "bad_repair"]),
      good_total = sum(count[repairment == "good_repair"])
    )
  },
  by = seurat_clusters
]

fwrite(results,"data/cluster_composition_test_result.csv")


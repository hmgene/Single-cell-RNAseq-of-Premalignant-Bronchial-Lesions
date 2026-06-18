library(Seurat)
library(data.table)

add_modules_named <- function(seurat_obj, mk_table, prefix) {
  
  # Split features by type
  modules <- split(mk_table$feature, mk_table$type)
  
  # Keep only genes that exist in Seurat object and are unique
  modules <- lapply(modules, function(g) {
    g <- unique(g)
    g[g %in% rownames(seurat_obj)]
  })
  
  modules <- modules[sapply(modules, length) > 0]
  if (length(modules) == 0) {
    warning("No modules with 0 genes found in Seurat object.")
    return(seurat_obj)
  }
  
  type_names <- names(modules)
  
  # Add module scores (temporary numeric suffix)
  seurat_obj <- AddModuleScore(
    seurat_obj,
    features = modules,
    name = paste0(prefix, "_tmp")
  )
  
  # Find newly created columns
  new_cols <- grep(paste0("^", prefix, "_tmp"), 
                   colnames(seurat_obj@meta.data),
                   value = TRUE)
  
  # Rename columns to prefix_type
  for (i in seq_along(type_names)) {
    clean_name <- gsub("[^A-Za-z0-9]+", "_", type_names[i])
    new_name <- paste0(prefix, "_", clean_name)
    
    colnames(seurat_obj@meta.data)[
      colnames(seurat_obj@meta.data) == new_cols[i]
    ] <- new_name
  }
  
  # Return the updated Seurat object
  return(seurat_obj)
}


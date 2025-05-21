run_ProjecTILs <- function(
    ref_path,
    query_path,
    save.loc,
    query_assay = "RNA",
    split_by = "sampleSource",
    reduction = "umap",
    k = 20,
    filter_cell = TRUE,
    skip_normalize = TRUE,
    metadata_column = "predictedCluster",
    return.seu.obj = TRUE
) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(ProjecTILs)
  })
  
  message("[MSG] Loading reference and query...")
  ref <- readRDS(ref_path)
  query <- readRDS(query_path)
  
  # Convert `query` assay to Seurat v4 assay as v5 is incompatible for now
  DefaultAssay(query) <- query_assay
  query[[query_assay]] <- as(query[[query_assay]], Class = "Assay")
  
  message("[MSG] Running ProjecTILs projection...")
  projected <- Run.ProjecTILs(
    query,
    ref = ref,
    split.by = split_by,
    reduction = reduction,
    filter.cell = filter_cell,
    skip.normalize = skip_normalize,
    k = k
  )
  
  message("[MSG] Saving predicted cluster metadata...")
  cluster_annotations <- setNames(projected$functional.cluster, colnames(projected))
  dir.create(save.loc, showWarnings = FALSE, recursive = TRUE)
  saveRDS(cluster_annotations, file.path(save.loc, paste0(query@project.name, ".ProjecTILs.cell_annotations.rds")))
  
  if (isTRUE(return.seu.obj)) {
    message("[MSG] Adding predicted cluster annotations to query Seurat object...")
    query <- AddMetaData(query, cluster_annotations, col.name = metadata_column)
    return(query)
  }
}

# # save legend only
# library(patchwork)
# library(cowplot)
# p <- DimPlot(query, label = T, repel = T, label.size = 4, group.by = metadata_column, reduction = "umap", cols = ref_col) +
#   DimPlot(query, label = T, repel = T, label.size = 4, group.by = "seuratClusters", reduction = "umap", cols = cluster_col) & theme(aspect.ratio = 1)
# 
# ggsave(file.path(save.loc,"predictedCluster.legend_only.pdf"),
#        as_ggplot(get_legend(p)), device = "pdf",
#        width = 4, height = 5)
# 

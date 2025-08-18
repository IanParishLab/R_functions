run_ProjecTILs <- function(
  ref_path,
  query_path,
  classifier_mode = TRUE,
  save.loc,
  query_assay = "RNA",
  split_by = "sampleSource",
  reduction = "pca",
  k = 20,
  filter_cell = TRUE,
  skip_normalize = TRUE,
  metadata_column = "predicted_cluster",
  return.seu.obj = TRUE
) {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(ProjecTILs)
  })
  
  # Create output directory
  dir.create(file.path(save.loc), showWarnings = FALSE, recursive = TRUE)

  
  message("[MSG] Loading reference and query objects...")
  ref <- readRDS(ref_path)
  query <- readRDS(query_path)
  
  # Ensure assay compatibility by converting to Seurat v4 Assay class
  DefaultAssay(query) <- query_assay
  query[[query_assay]] <- as(query[[query_assay]], Class = "Assay")
  
  if (classifier_mode) {
    message("[MSG] Using ProjecTILs.classifier...")
    projected <- ProjecTILs.classifier(
      query,
      ref = ref,
      split.by = split_by,
      filter.cells = filter_cell,
      reduction = reduction,
      ndim = NULL,
      k = k,
      labels.col = metadata_column
    )
  } else {
    message("[MSG] Using Run.ProjecTILs...")
    projected <- Run.ProjecTILs(
      query,
      ref = ref,
      split.by = split_by,
      reduction = reduction,
      filter.cells = filter_cell,
      skip.normalize = skip_normalize,
      k = k
    )
  }
  
  message("[MSG] Extracting and saving predicted cluster annotations...")
  cluster_annotations <- setNames(projected$functional.cluster, colnames(projected))
  saveRDS(cluster_annotations, file.path(save.loc, paste0(query@project.name, ".ProjecTILs.cell_annotations.rds")))
  
  if (return.seu.obj) {
    message("[MSG] Adding predicted cluster annotations to query Seurat object...")
    query <- AddMetaData(query, cluster_annotations, col.name = metadata_column)
    return(query)
  }
}

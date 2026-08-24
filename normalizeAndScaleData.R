normalizeAndScaleData <- function(seu.obj,
                                  assay = "RNA",
                                  sample_name,
                                  min.pc,
                                  seed.use,
                                  save.loc,
                                  vars.to.regress = "percent.mito",
                                  dry_run = TRUE,
                                  bk.list = NULL,
                                  verbose = FALSE,
                                  ...) {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(scGate)
    library(tidyverse)
  })

  # Create output directory for plots
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  message("[MSG] Setting default assay to RNA and starting normalization...")
  DefaultAssay(seu.obj) <- assay
  
  seu.obj <- NormalizeData(object = seu.obj, verbose = verbose)
  seu.obj <- FindVariableFeatures(object = seu.obj, verbose = verbose)
  
  if (!is.null(bk.list)) {
    message("[MSG] Removing blacklist features from variable features...")
    VariableFeatures(seu.obj) <- unique(setdiff(VariableFeatures(seu.obj), bk.list))
  }

  message("[MSG] Scaling data and running PCA...")
  seu.obj <- ScaleData(seu.obj, vars.to.regress = vars.to.regress, verbose = verbose, ...)
  seu.obj <- RunPCA(seu.obj, verbose = verbose)

  # Save ElbowPlot
  elbow_path <- file.path(save.loc, "plots", paste0(sample_name, ".ElbowPlot.pdf", collapse = "_"))
  message("[MSG] Saving ElbowPlot to: ", elbow_path)
  ggsave(Seurat::ElbowPlot(seu.obj), filename = elbow_path, device = "pdf", width = 5, height = 5)

  if (isFALSE(dry_run)) {
    message("[MSG] Running UMAP, neighbor graph, and clustering...")
    seu.obj <- RunUMAP(seu.obj, dims = 1:min.pc, return.model = TRUE, seed.use = seed.use, verbose = verbose)
    seu.obj <- FindNeighbors(seu.obj, dims = 1:min.pc, assay = "RNA", reduction = "pca", verbose = verbose)
    seu.obj <- FindClusters(seu.obj, resolution = 1, verbose = verbose)
  } else {
    message("[MSG] Dry run enabled — skipping UMAP, neighbors, and clustering.")
  }

  return(seu.obj)
}

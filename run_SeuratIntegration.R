run_SeuratIntegration <- function(seu.obj.list, merged.seu.obj.name, min.pc, seed.use, 
                                  save.loc = "integration/", integration_method, integration_method_name,
                                  vars.to.regress = "percent.mito",
                                  conda_env_path = "/home/nsaw/miniconda3/envs/scvi-env/", 
                                  DimPlot.ident = "HTO_maxID", bk.list = NULL, add.cell.ids = NULL,
                                  cols = list(HTO_col, capture_col),
                                  verbose = FALSE) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratWrappers)
    library(reticulate)
    library(tidyverse)
  })
  
  # Create output directory
  dir.create(file.path(save.loc), showWarnings = FALSE, recursive = TRUE)

  # get color palette if not provided
  if (is.null(capture_col) && exists("capture_col", envir = .GlobalEnv)) {
    capture_col <- get("capture_col", envir = .GlobalEnv)
  }
  if (is.null(HTO_col) && exists("HTO_col", envir = .GlobalEnv)) {
    HTO_col <- get("HTO_col", envir = .GlobalEnv)
  }
  
  reduction <- paste0("integrated.", integration_method_name)
  reduction.name <- paste0("umap.", integration_method_name)
  
  # Validate seu.obj.list
  if (length(seu.obj.list) < 2) {
    stop("[ERROR] seu.obj.list must contain at least two Seurat objects.")
  }
  
  # Merge scRNAseq objects
  message("[MSG] Merging Seurat objects ...")
  if (is.null(add.cell.ids)) {
    add.cell.ids <- sapply(seu.obj.list, \(seu) seu@project.name)
  }
  merged.seu.obj <- merge(
    seu.obj.list[[1]],
    y = seu.obj.list[2:length(seu.obj.list)],
    add.cell.ids = if (!identical(add.cell.ids, "none")) add.cell.ids else NULL
  )
  
  message("[MSG] Normalizing and scaling data ...")
  merged.seu.obj <- normalizeAndScaleData(
    merged.seu.obj, 
    dry_run = FALSE,
    sample_name = merged.seu.obj.name, 
    min.pc = min.pc, 
    seed.use = seed.use, 
    save.loc = save.loc,
    vars.to.regress = vars.to.regress,
    bk.list = bk.list, 
    verbose = verbose
  )
  
  # Run UMAP before integration
  merged.seu.obj <- RunUMAP(
    merged.seu.obj, dims = 1:min.pc, seed.use = seed.use, 
    reduction = "pca", reduction.name = "umap.unintegrated", verbose = verbose
  )
  
  theme <- theme(aspect.ratio = 1) & NoAxes()
  grouped_by <- if (!is.null(DimPlot.ident)) c(DimPlot.ident, "orig.ident") else "orig.ident"
  
  before_int_umap <- DimPlot(
    merged.seu.obj, reduction = "umap.unintegrated", 
    group.by = grouped_by, cols = cols
  ) & theme
  
  # Integrate layers
  message("[MSG] Running IntegrateLayers ...")
  combined <- IntegrateLayers(
    object = merged.seu.obj,
    method = integration_method,
    orig.reduction = "pca",
    new.reduction = reduction,
    conda_env = conda_env_path,
    verbose = TRUE
  )
  
  message("[MSG] Running JoinLayers ...")
  combined <- JoinLayers(combined)
  
  combined <- RunUMAP(
    combined, dims = 1:min.pc, seed.use = seed.use,
    reduction = reduction, reduction.name = reduction.name, 
    verbose = FALSE
  )
  
  after_int_umap <- DimPlot(
    combined, reduction = reduction.name, 
    group.by = grouped_by, cols = cols
  ) & theme
  
  message("[MSG] Saving integrated object ...")
  saveRDS(combined, file.path(save.loc, paste0(merged.seu.obj.name, ".rds")))
  
  ggsave(
    before_int_umap / after_int_umap,
    filename = file.path(save.loc, paste0("before_after_integration_umap.", merged.seu.obj.name, ".pdf")),
    width = 10, height = 10, device = "pdf"
  )
}

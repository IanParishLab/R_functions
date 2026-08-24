run_SeuratIntegration <- function(seu.obj.list, sample_name, min.pc, seed.use, 
                                  save.loc = "integration/", integration_method, integration_method_name,
                                  vars.to.regress = "percent.mito",
                                  conda_env_path = "/home/nsaw/.conda/envs/scvi-tools/", 
                                  bk.list = NULL, add.cell.ids = NULL,
                                  DimPlot.ident = "HTO_maxID", cols = NULL,
                                  verbose = FALSE) {
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratWrappers)
    library(reticulate)
    library(tidyverse)
    library(qs2)
  })
  
  # Create output directory
  dir.create(file.path(save.loc), showWarnings = FALSE, recursive = TRUE)

  reduction <- paste0("integrated.", integration_method_name)
  reduction.name <- paste0("umap.", integration_method_name)
  
  # Validate seu.obj.list
  if (length(seu.obj.list) < 2) {
    stop("[ERROR] seu.obj.list must contain at least two Seurat objects.")
  }

  # Merge scRNAseq objects
  message("[MSG] Merging Seurat objects ...")
  merged.seu.obj <- merge(
    seu.obj.list[[1]],
    y = seu.obj.list[2:length(seu.obj.list)],
    add.cell.ids = add.cell.ids
  )

  # # add chunk for covariates
  # merged.seu.obj[["RNA"]] <- split(merged.seu.obj[["RNA"]], f = merged.seu.obj$covariates)
  
  message("[MSG] Normalizing and scaling data ...")
  merged.seu.obj <- normalizeAndScaleData(
    merged.seu.obj, 
    dry_run = FALSE,
    sample_name = sample_name, 
    min.pc = min.pc, 
    seed.use = seed.use, 
    save.loc = save.loc,
    vars.to.regress = vars.to.regress,
    bk.list = bk.list, 
    verbose = verbose
  )

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
  
  message("[MSG] Saving integrated object in", save.loc, "...")
  dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)
  qs_save(combined, file.path(save.loc, "int_obj", paste0(sample_name,".",integration_method_name,".integrated.qs")))
  
  # # plot before after umap
  # theme <- theme(aspect.ratio = 1) & NoAxes()
  # grouped_by <- if (!is.null(DimPlot.ident)) c(DimPlot.ident, "orig.ident") else "orig.ident"
  # 
  # before_int_umap <- DimPlot(
  #   merged.seu.obj, reduction = "umap",
  #   group.by = grouped_by, cols = cols
  # ) & theme
  # 
  # after_int_umap <- DimPlot(
  #   combined, reduction = reduction.name,
  #   group.by = grouped_by, cols = cols
  # ) & theme
  # 
  # ggsave(
  #   before_int_umap / after_int_umap,
  #   filename = file.path(save.loc, paste0("before_after_integration_umap.", sample_name, ".pdf")),
  #   width = 10, height = 10, device = "pdf"
  # )
}

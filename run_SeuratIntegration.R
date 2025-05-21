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
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc, " directory exists"))
  
  reduction <- paste0("integrated.", integration_method_name)
  reduction.name <- paste0("umap.", integration_method_name)
  
  # Merge scRNAseq into one dataset
  if(is.null(add.cell.ids)){
  print("[MSG] Merging data & add.cell.ids with `seu.obj@project.name` ...")
    add.cell.ids <- sapply(seu.obj.list, function(seu) { seu@project.name }) %>% unlist %>% unname
    merged.seu.obj <- merge(seu.obj.list[[1]], y = seu.obj.list[2:length(seu.obj.list)], add.cell.ids = add.cell.ids)
  } else if(!is.null(add.cell.ids)) {
    print("[MSG] Merging data & add provided add.cell.ids ...")
    merged.seu.obj <- merge(seu.obj.list[[1]], y = seu.obj.list[2:length(seu.obj.list)], add.cell.ids = add.cell.ids)
  }else if(add.cell.ids == "none") {
    print("[MSG] Merging data ...")
    merged.seu.obj <- merge(seu.obj.list[[1]], y = seu.obj.list[2:length(seu.obj.list)])
  }
  
  print("[MSG] normalizeAndScaleData ...")
  merged.seu.obj <- normalizeAndScaleData(merged.seu.obj, 
                                          dry_run = FALSE,
                                          sampleName = merged.seu.obj.name, 
                                          min.pc = min.pc, 
                                          seed.use = seed.use, 
                                          save.loc = save.loc,
                                          vars.to.regress = vars.to.regress,
                                          bk.list = bk.list, 
                                          verbose = verbose)
  
  merged.seu.obj <- RunUMAP(merged.seu.obj, dims = 1:min.pc, seed.use = seed.use, 
                            reduction = "pca", reduction.name = "umap.unintegrated", verbose = verbose)
  # set theme
  theme = theme(aspect.ratio = 1) & 
    NoAxes()
  
  # plot
  grouped_by <- if(!is.null(DimPlot.ident)){
    c(DimPlot.ident, "orig.ident")
  } else {
    "orig.ident"
  }
  before_int_umap <- DimPlot(merged.seu.obj, reduction = "umap.unintegrated", 
                             group.by = grouped_by, 
                             cols = cols) & theme
  
  print("[MSG] IntegrateLayers ...")
  combined <- IntegrateLayers(object = merged.seu.obj,
                              method = integration_method,
                              orig.reduction = "pca",
                              new.reduction = reduction,
                              conda_env = conda_env_path,
                              verbose = TRUE)
  
  print("[MSG] JoinLayers ...")
  combined <- JoinLayers(combined)
  
  combined <- RunUMAP(combined, dims = 1:min.pc, seed.use = seed.use,
                      reduction = reduction, reduction.name = reduction.name, 
                      verbose = FALSE)
  
  # plot
  after_int_umap <- DimPlot(combined, 
                            reduction = reduction.name, 
                            group.by = grouped_by, 
                            cols = cols) & theme
  
  print("[MSG] Saving integrated object ...")
  saveRDS(combined, file.path(save.loc, paste0(merged.seu.obj.name, ".rds")))
  
  # Save plot
  ggsave(before_int_umap / after_int_umap,
         filename = file.path(save.loc, paste0("before_after_integratelayers_umap.", merged.seu.obj.name, ".pdf")),
         width = 10, height = 10, device = 'pdf')
}
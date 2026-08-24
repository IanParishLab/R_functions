run_scMerge2Integration <- function(seu_list, group.by = "PatientID", ctl,
                                    sample_name, min.pc, seed.use, vars.to.regress = "percent.mito",
                                    dry_run = TRUE, verbose = FALSE, save.loc = "Integration") {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(qs2)
    library(scater)
    library(SingleCellExperiment)
    library(scMerge)
    library(tidyverse)
    library(patchwork)
  })

    if (!dry_run) {
    dir.create(file.path(save.loc, "int_obj"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(save.loc, "plots"), recursive = TRUE, showWarnings = FALSE)
  }
  
  ## Rename cells / harmonise metadata
  add.cell.ids <- sapply(seu_list, function(seu) seu@project.name) %>% unlist()
  
  for (i in seq_along(seu_list)) {
    object <- seu_list[[i]]
    
    if (verbose) message("[MSG] Processing: ", object@project.name)
    
    if (is.null(object@meta.data$sample)) {
      object$sample <- object@project.name
    }
    if (is.null(object@meta.data$cellBarcode)) {
      object$cellBarcode <- colnames(object)
    }
    
    object <- RenameCells(object, new.names = paste0(object@meta.data$sample, "_", object@meta.data$cellBarcode))
    
    DefaultAssay(object) <- "RNA"
    
    seu_list[[i]] <- object
  }
  
  ## Metadata ----
  metadata <- lapply(seu_list, function(seu){ seu@meta.data }) %>%
    plyr::rbind.fill() %>%
    unite("cell_bc", sample, cellBarcode, sep = "_", remove = FALSE)
  
  ## Convert to SCE and merge ----
  merged_sce <- lapply(seu_list, function(object) {
    object[["RNA"]] <- as(object[["RNA"]], Class = "Assay")
    return(as.SingleCellExperiment(object))
  }) %>%
    sce_cbind(batch_names = add.cell.ids)
  
  # ensure cell names are consistent
  stopifnot(ncol(merged_sce) == nrow(metadata))
  colnames(merged_sce) <- metadata$cell_bc
  
  merged_sce <- logNormCounts(merged_sce)
  merged_sce <- runPCA(merged_sce, exprs_values = "logcounts")
  
  # sanity check ----
  stopifnot(
    identical(colnames(merged_sce), metadata$cell_bc),
    identical(colnames(merged_sce), colnames(logcounts(merged_sce)))
  )
  
  colData(merged_sce)[[group.by]] <- metadata[[group.by]]
  
  # plot before scMerge integration
  p1 <- plotPCA(merged_sce, colour_by = group.by)
  
  ## ---------------------------
  ## scMerge2
  ## ---------------------------
  message("scMerge2...")
  
  scMerge2_res <- scMerge2(
    exprsMat = logcounts(merged_sce),
    batch = merged_sce$batch,
    ctl = ctl,
    verbose = verbose
  )
  
  assay(merged_sce, "scMerge2") <- scMerge2_res$newY
  merged_sce <- runPCA(merged_sce, exprs_values = "scMerge2")
  
  # plot after scMerge integration
  p2 <- plotPCA(merged_sce, colour_by = group.by)
  
  if (!dry_run) {
    message("Saving `merged_sce` & its `metadata`...")
    qs_save(merged_sce, file.path(save.loc, "int_obj", paste0(sample_name,".scMerge.sce.qs")))
    qs_save(metadata, file.path(save.loc, "int_obj", paste0(sample_name,".scMerge.sce.metadata.qs")))
  }
  
  ## ---------------------------
  ## Convert back to Seurat
  ## ---------------------------
  message("CreateSeuratObject...")
  seu <- CreateSeuratObject(
    counts = counts(merged_sce),
    meta.data = metadata,
    min.cells = 1,
    project = sample_name
  )
  
seu <- SetAssayData(
    seu,
    assay = "RNA",
    slot = "data",
    new.data = assay(merged_sce, "scMerge2")
  )

  seu <- FindVariableFeatures(seu, verbose = FALSE)
  
  bk.list <- get_blacklist(seu)
  if (!is.null(bk.list)) {
    message("[MSG] Removing blacklist features from variable features...")
    VariableFeatures(seu) <- unique(setdiff(VariableFeatures(seu), bk.list))
  }
  
  seu <- ScaleData(seu, vars.to.regress = vars.to.regress, verbose = FALSE)
  seu <- RunPCA(seu, npcs = 50, verbose = FALSE)
  
  ## ---------------------------
  ## Diagnostics & clustering
  ## ---------------------------
  elbow_path <- file.path(save.loc,"plots", paste0(sample_name, ".ElbowPlot.pdf"))
  
  if (!dry_run) {
    ggsave(ElbowPlot(seu),filename = elbow_path, device = "pdf", width = 5, height = 5)
  } else {
    print(ElbowPlot(seu))
  }
  
  message("RunUMAP...")
  seu <- RunUMAP(seu,dims = 1:min.pc, return.model = TRUE,seed.use = seed.use,verbose = FALSE)
  seu <- FindNeighbors(seu,dims = 1:min.pc,assay = "RNA", reduction = "pca",verbose = FALSE)
  seu <- FindClusters(seu,resolution = 1,verbose = FALSE)
  
  if (!dry_run) {
    message("Saving normalised scMerge-integrated Seurat object...")
    ggsave((p1|p2),filename = file.path(save.loc, "plots", paste0(sample_name,".before_after_scMerge.pca.pdf")), 
           device = "pdf", width = 12, height = 8)
    qs_save(seu, file.path(save.loc, "int_obj", paste0(sample_name,".scMerge.seu.normalised.qs")))
  }
  
  return(seu)
}

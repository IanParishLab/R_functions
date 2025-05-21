normalizeAndScaleData <- function(seu.obj, sampleName, min.pc, seed.use = 800, 
                                  save.loc, 
                                  vars.to.regress = "percent.mito",
                                  dry_run = TRUE, 
                                  bk.list = NULL,
                                  verbose = FALSE, 
                                  ...){
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(scGate)
    library(tidyverse)
  })
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc, "plots")),
         dir.create(file.path(save.loc, "plots"), recursive = TRUE), paste0(save.loc," directory exists"))
  
  # normalize
  DefaultAssay(seu.obj) <- "RNA"
  # print("[MSG] NormalizeData...")
  seu.obj <- NormalizeData(object = seu.obj, verbose = verbose)
  # print("[MSG] FindVariableFeatures...")
  seu.obj <- FindVariableFeatures(object = seu.obj, verbose = verbose)
  if(!is.null(bk.list)){
    VariableFeatures(seu.obj) <- unique(setdiff(VariableFeatures(seu.obj), bk.list))
  }
  # print("[MSG] ScaleData...")
  seu.obj <- ScaleData(seu.obj, vars.to.regress = vars.to.regress, verbose = verbose, ...)
  # print("[MSG] RunPCA...")
  seu.obj <- RunPCA(seu.obj, npcs = 50, verbose = verbose)
  
  # ElbowPlot
  ggsave(ElbowPlot(seu.obj, ndims = 50), 
         filename = file.path(save.loc, "plots", paste0(sampleName,".ElbowPlot.pdf")), 
         device = "pdf", width = 5, height = 5)
  
  if(isFALSE(dry_run)){
    print("[MSG] RunUMAP...")
    seu.obj <- RunUMAP(seu.obj, dims = 1:min.pc, return.model = TRUE, seed.use = seed.use, verbose = verbose)
    print("[MSG] FindNeighbors...")
    seu.obj <- FindNeighbors(seu.obj, dims = 1: min.pc, assay = "RNA", reduction = "pca", verbose = verbose)
    print("[MSG] FindClusters...")
    seu.obj <- FindClusters(seu.obj, resolution = 1, verbose = verbose)
  } else if(isTRUE(dry_run)){
    seu.obj <- seu.obj
  }
  
  return(seu.obj)
}


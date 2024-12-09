normalizeAndScaleData <- function(seu.obj, sampleName, min.pc, seed.use = 800, 
                                  save.loc = "initial_filtering/", 
                                  vars.to.regress = c("percent.mito","CC.Difference","RP.Score1"),
                                  dry_run = TRUE, 
                                  bk.list = scGate::genes.blacklist.default$Mm$Ribo){
  
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
  seu.obj <- NormalizeData(object = seu.obj)
  seu.obj <- FindVariableFeatures(object = seu.obj)
  VariableFeatures(seu.obj) <- unique(setdiff(VariableFeatures(seu.obj), bk.list))
  seu.obj <- ScaleData(seu.obj, verbose = FALSE, vars.to.regress = vars.to.regress)
  seu.obj <- RunPCA(seu.obj, npcs = 50, verbose = FALSE)
  
  # ElbowPlot
  ggsave(ElbowPlot(seu.obj, ndims = 50), 
         filename = file.path(save.loc, "plots", paste0(sampleName,".ElbowPlot.pdf")), 
         device = "pdf", width = 5, height = 5)
  
  if(isFALSE(dry_run)){
    seu.obj <- RunUMAP(seu.obj, dims = 1:min.pc, return.model = TRUE, seed.use = seed.use)
    seu.obj <- FindNeighbors(seu.obj, dims = 1: min.pc, assay = "RNA", reduction = "pca")
    seu.obj <- FindClusters(seu.obj, resolution = 1)
  } else if(isTRUE(dry_run)){
    seu.obj <- seu.obj
  }
  
  return(seu.obj)
}


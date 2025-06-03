plot_CellCycleRegression <- function(seu.obj, sample_name, ccgenes, save.loc , npcs = 30){
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
  })
  
  # make output directory
dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)

  # see if cell cycle needs to be regressed
  seu.obj <- NormalizeData(object = seu.obj)
  seu.obj <- CellCycleScoring(seu.obj, 
                              g2m.features = ccgenes$g2m.genes, 
                              s.features = ccgenes$s.genes)
  seu.obj$CC.Difference <- seu.obj$S.Score - seu.obj$G2M.Score # calculate cell cycle difference 
  seu.obj <- ScaleData(seu.obj, verbose = FALSE)
  seu.obj <- RunPCA(seu.obj, npcs = npcs, features = unlist(ccgenes), verbose = FALSE)
  
  # save plot
  ggsave(DimPlot(seu.obj, reduction = "pca", group.by = "Phase"), 
         filename = file.path(save.loc,"plots", paste0(sample_name, ".CellCyclePhase.pdf")), 
         device = "pdf", width = 5, height = 5)
  
  return(seu.obj)
}

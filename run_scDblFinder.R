run_scDblFinder <- function(seu.obj, sampleName, plot = FALSE, save.loc, samples = NULL, clusters = NULL){
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(SingleCellExperiment)
    require(scDblFinder)
    require(patchwork)
    require(tidyverse)
  })
  
  sce <- as.SingleCellExperiment(seu.obj, assay = "RNA")
  sce.dbl <- scDblFinder(sce, clusters = clusters, samples = samples)
  seu.obj[["scDblFinder.class"]] <- sce.dbl[["scDblFinder.class"]]
  seu.obj[["scDblFinder.score"]] <- sce.dbl[["scDblFinder.score"]]
  seu.obj@meta.data[["scDblFinder.score"]] <- as.numeric(seu.obj@meta.data[["scDblFinder.score"]])
  
  if(isTRUE(clusters)){
    seu.obj[["scDblFinder.cluster"]] <- sce.dbl[["scDblFinder.cluster"]]
  }
  
  if(isTRUE(plot)){
    if(isTRUE(clusters)){
      group_by.clusters <- "scDblFinder.cluster"
    } else {
      group_by.clusters <- NULL
    }
    
    # make output directory
    ifelse(!dir.exists(file.path(save.loc, "plots")),
           dir.create(file.path(save.loc, "plots"), recursive = TRUE), paste0(save.loc," directory exists"))
    
    p1 <- DimPlot(seu.obj, reduction = "umap", group.by = "scDblFinder.class", order = TRUE) + theme(aspect.ratio = 1) 
    p2 <- DimPlot(seu.obj, reduction = "umap", group.by = group_by.clusters, order = TRUE) + theme(aspect.ratio = 1) 
    p3 <- FeaturePlot(seu.obj, reduction = "umap", "scDblFinder.score", order = TRUE, repel = TRUE) +
      theme(aspect.ratio = 1) 
    p4 <- VlnPlot(seu.obj, group.by = "orig.ident", 
                  features = c("nCount_RNA","nFeature_RNA","percent.mito"), 
                  split.by = "scDblFinder.class", ncol = 3) &
      theme(text = element_text(size = 7), axis.text = element_text(size = 7))
    
    # scDblFinder plot
    pdf(file.path(save.loc, "plots", paste0(sampleName,".scDblFinder.pdf")), width = 12, height = 5)
    print(p1+p2+p3)
    print(p4)
    dev.off()
  }
  
  # return obj
  return(seu.obj)
  
}

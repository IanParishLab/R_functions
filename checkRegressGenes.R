checkRegressGenes <- function(seu.obj, sample_name, pca_npcs = 50, ccgenes, save.loc, regress_genes, label_genes = NULL){        
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(UCell)
    require(patchwork)
    require(tidyverse)
    require(RColorBrewer)
  })
  
  # make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)

  seu.obj <- NormalizeData(seu.obj, verbose = FALSE)
  seu.obj <- FindVariableFeatures(seu.obj,selection.method = "vst", nfeatures = 2000)
  seu.obj <- ScaleData(seu.obj, verbose = FALSE)
  seu.obj <- RunPCA(seu.obj, npcs = pca_npcs, verbose = FALSE)
  
  # get gene lists
  if(!is.null(regress_genes)){
    g <- setNames(c(ccgenes, regress_genes), c("S","G2M", names(regress_genes)))
    g <- lapply(g, function(g){
      g[g %in% rownames(seu.obj)]
    })
  } else {
    g <- setNames(ccgenes, c("S","G2M"))
    g <- lapply(g, function(g){
      g[g %in% rownames(seu.obj)]
    })
  }
  
  # add module scores
  seu.obj <- AddModuleScore_UCell(seu.obj, features = g)
  signature.names <- paste0(names(g), "_UCell")
  seu.obj <- SmoothKNN(seu.obj,
                       signature.names = signature.names,
                       reduction = "pca")
  # plot module scores
  moduleScorePlots <- FeaturePlot(seu.obj, reduction = "pca", features = c(paste0(signature.names,"_kNN")), ncol = 2) &
    theme(aspect.ratio = 1, axis.text = element_text(size = 7), axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  # https://github.com/satijalab/seurat/issues/147#issuecomment-360896778
  # VariableFeaturePlot here shows that some Rp genes are highly variable amongst cells, 
  # so ideally we wouldn't filter these genes out in future analyses and only consider regressing them out instead.
  
  if (!is.null(regress_genes)) {
    if (is.list(regress_genes) && !is.null(names(regress_genes))) {
      # Handle named lists
      for (feature_name in names(regress_genes)) {
        seu.obj <- AddModuleScore(
          object = seu.obj, 
          features = list(regress_genes[[feature_name]]), 
          name = paste0(feature_name,".Score"), 
          assay = "RNA"
        )
      }
    } else {
      # Handle unnamed lists
      seu.obj <- AddModuleScore(
        object = seu.obj, 
        features = regress_genes, 
        name = names(regress_genes), 
        assay = "RNA"
      )
    }
    VariableFeatures(seu.obj) <- setdiff(VariableFeatures(seu.obj), unlist(regress_genes))
    label_genes <- if (is.null(label_genes)) {
      unlist(regress_genes)[unlist(regress_genes) %in% rownames(seu.obj)]
    } else {
      label_genes
    }
  } else {
    label_genes <- unlist(ccgenes)[unlist(ccgenes) %in% rownames(seu.obj)]
  }

    # plot HVGs
  HVGLabelledPlot <- LabelPoints(plot = VariableFeaturePlot(seu.obj), repel = TRUE, points = label_genes) +
    theme(aspect.ratio = 1)
  
  # save plots
  pdf(file.path(save.loc, "plots", paste0(sample_name, ".checkRegressGenes_plots.pdf")), width = 6, height = 6)
  print(moduleScorePlots)
  print(HVGLabelledPlot)
  dev.off()
  
  
  return(seu.obj)
}

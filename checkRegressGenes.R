checkRegressGenes <- function(seu.obj, sampleName, pca_npcs = 50, ccgenes, save.loc, regress_genes, ribo_name = "Ribo", label_genes = NULL){        
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(UCell)
    require(patchwork)
    require(tidyverse)
    require(RColorBrewer)
  })
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc, "plots")),
         dir.create(file.path(save.loc, "plots"), recursive = TRUE), paste0(save.loc," directory exists"))
  
  seu.obj <- NormalizeData(seu.obj)
  seu.obj <- FindVariableFeatures(seu.obj,selection.method = "vst", nfeatures = 2000)
  seu.obj <- ScaleData(seu.obj, verbose = FALSE)
  seu.obj <- RunPCA(seu.obj, npcs = pca_npcs, verbose = FALSE)
  
  # get gene lists
  g <- setNames(c(ccgenes, regress_genes), c("S","G2M", names(regress_genes)))
  g <- lapply(g, function(g){
    g[g %in% rownames(seu.obj)]
  })
  
  # add module scores
  seu.obj <- AddModuleScore_UCell(seu.obj, features = g)
  signature.names <- paste0(names(g), "_UCell")
  seu.obj <- SmoothKNN(seu.obj,
                       signature.names = signature.names,
                       reduction = "pca")
  # plot module scores
  moduleScorePlots <- FeaturePlot(seu.obj, reduction = "pca", features = c(paste0(signature.names,"_kNN")), ncol = 2) &
    theme(aspect.ratio = 1, axis.title.x = element_text(size = 7), axis.title.y = element_text(size = 7)) &
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  
  # https://github.com/satijalab/seurat/issues/147#issuecomment-360896778
  # VariableFeaturePlot here shows that some Rp genes are highly variable amongst cells, so ideally we wouldn't filter these genes out in future analyses and only consider regressing them out instead.
  
  if(is.null(label_genes)) { 
    label_genes = unlist(regress_genes)[unlist(regress_genes) %in% rownames(seu.obj)]
  }
  
  HVGLabelledPlot <- LabelPoints(plot = VariableFeaturePlot(seu.obj), repel = TRUE, points = label_genes) +
    theme(aspect.ratio = 1)
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  # save plots
  pdf(file.path(save.loc, "plots", paste0(sampleName, ".checkRegressGenes_plots.pdf")), width = 6, height = 6)
  print(moduleScorePlots)
  print(HVGLabelledPlot)
  dev.off()
  
  # aim to regress out RP genes, remove RP genes from VariableFeatures(seu.obj)
  seu.obj <- AddModuleScore(object = seu.obj,
                            features = g[ribo_name],
                            name = 'RP.Score',
                            assay = "RNA")
  # VariableFeatures(seu.obj) <- setdiff(VariableFeatures(seu.obj), grep(rp_gene_pattern, VariableFeatures(seu.obj), value = TRUE))
  VariableFeatures(seu.obj) <- setdiff(VariableFeatures(seu.obj), g[ribo_name])
  
  return(seu.obj)
}

find_clusters_to_remove_multiome <- function(so, 
                                             capture, 
                                             assays, 
                                             resolutions = seq(0.5, 2, by = 0.1), 
                                             remove_clusters = NULL, 
                                             remove_clusters_res = NULL, 
                                             npcs,
                                             bk.list, 
                                             vars.to.regress,
                                             prefix,
                                             raster = TRUE,
                                             verbose = FALSE,
                                             process_wsnn = TRUE, 
                                             save.loc,
                                             seed.use,
                                             normalize = TRUE,
                                             scale = TRUE) {
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
    require(clustree)
    require(cowplot)
    require(Signac)
  })
  
  vln_features <- if ("RNA" %in% assays) {
    c("nCount_RNA", "nFeature_RNA", "percent.mito")
  } else if ("ATAC" %in% assays) {
    c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "percent.mito")
  }
  
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Subset out clusters if specified
  if(!is.null(remove_clusters)){
    remove.cells <- colnames(so)[so@meta.data[[remove_clusters_res]] %in% remove_clusters]
    so <- subset(so, cells = remove.cells, invert = TRUE)
  }
  
  if ("RNA" %in% assays) {
    assay <- "RNA"
    message("processing RNA assay...")
    DefaultAssay(so) <- assay
    
    # Normalization and scaling
    if (isTRUE(normalize)) {
      so <- NormalizeData(so)
    }
    
    so <- FindVariableFeatures(so)
    VariableFeatures(so) <- setdiff(VariableFeatures(so), bk.list)
    if (!is.null(features)) {
      VariableFeatures(so) <- unique(c(VariableFeatures(so), features))
    }
    
    if (isTRUE(scale)) {
      so <- ScaleData(so, vars.to.regress = vars.to.regress)
    } 
    
    so <- RunPCA(so, npcs = npcs, verbose = verbose)
    so <- RunUMAP(so, dims = 1:npcs, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = verbose, seed.use = seed.use)
    message("seed.use = ", seed.use)
    
    graph.name <- paste0(assay, "_snn")
    so <- FindNeighbors(so, graph.name = graph.name, dims = 1:npcs, verbose = verbose)
    for (res in resolutions) {
      so <- FindClusters(so, graph.name = graph.name, resolution = res, verbose = verbose)
    }
    
    pdf(file.path(save.loc, "plots", paste0(c(prefix, capture, assay, "snn_res_VlnPlot.pdf"), collapse = "_")), height = 9, width = 12)
    plots <- lapply(resolutions, function(res) {
      so <- SetIdent(so, value = paste0(assay, "_snn_res.", res))
      vln <- VlnPlot(so, features = vln_features, pt.size = 0, ncol = 1, raster = raster)
      title <- ggdraw() + draw_label(paste0("resolution = ", res), fontface = 'bold')
      cowplot::plot_grid(title, vln, ncol = 1, rel_heights = c(0.1, 1))
    })
    invisible(lapply(plots, print))
    dev.off()
    
    pdf(file.path(save.loc, "plots", paste0(c(prefix, capture, "clustree.pdf"), collapse = "_")), height = 12, width = 10)
    print(clustree(so, prefix = paste0(assay, "_snn_res.")))
    dev.off()
  }
  
  if ("ATAC" %in% assays) {
    assay <- "ATAC"
    message("processing ATAC assay...")
    DefaultAssay(so) <- assay
    
    so <- RunTFIDF(so)
    so <- FindTopFeatures(so, min.cutoff = 'q0')
    so <- RunSVD(so)
    so <- RunUMAP(so, reduction = 'lsi', dims = 2:npcs, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = verbose, seed.use = seed.use)
    message("seed.use = ", seed.use)
    
    graph.name <- paste0(assay, "_snn")
    so <- FindNeighbors(so, dims = 2:npcs, assay = assay, reduction = "lsi", graph.name = graph.name, verbose = verbose)
    for (res in resolutions) {
      so <- FindClusters(so, graph.name = graph.name, resolution = res, algorithm = 3, verbose = verbose)
    }
    
    pdf(file.path(save.loc, "plots", paste0(c(prefix, capture, assay, "snn_res_VlnPlot.pdf"), collapse = "_")), height = 9, width = 12)
    plots <- lapply(resolutions, function(res) {
      so <- SetIdent(so, value = paste0(assay, "_snn_res.", res))
      vln <- VlnPlot(so, features = vln_features, pt.size = 0, ncol = 1)
      title <- ggdraw() + draw_label(paste0("resolution = ", res), fontface = 'bold')
      cowplot::plot_grid(title, vln, ncol = 1, rel_heights = c(0.1, 1))
    })
    invisible(lapply(plots, print))
    dev.off()
    
    pdf(file.path(save.loc, "plots", paste0(c(prefix, capture, assay, "clustree.pdf"), collapse = "_")), height = 12, width = 10)
    print(clustree(so, prefix = paste0(assay, "_snn_res.")))
    dev.off()
  }
  
  if (process_wsnn) {
    print("processing MultiModal information and generating weighted nearest neighbors (WNN) clustree...")
    so <- FindMultiModalNeighbors(so, reduction.list = list("pca", "lsi"), dims.list = list(1:npcs, 2:npcs))
    so <- RunUMAP(so, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
    for (res in resolutions) {
      so <- FindClusters(so, graph.name = "wsnn", algorithm = 3, resolution = res, verbose = verbose)
    }
    
    pdf(file.path(save.loc, "plots", paste0(c(prefix, capture, "wsnn_clustree.pdf"), collapse = "_")), height = 12, width = 10)
    print(clustree(so, prefix = "wsnn_res."))
    dev.off()
  }
  
  return(so)
}

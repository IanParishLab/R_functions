find_clusters_to_remove <- function(so,
                                    capture,
                                    assay,
                                    reduction,
                                    reduction.name,
                                    resolutions = seq(0.5, 2, by = 0.1),
                                    remove_clusters = NULL,
                                    remove_clusters_res = NULL,
                                    npcs,
                                    bk.list,
                                    vars.to.regress,
                                    prefix,
                                    features,
                                    save.loc,
                                    verbose = FALSE,
                                    raster = TRUE,
                                    seed.use) {
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
    require(clustree)
    require(cowplot)
    require(scGate)
  })
  
  # Create output directories if they don't exist
  dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  DefaultAssay(so) <- assay
  message("[MSG] assay = ", assay)
  
  # Subset out clusters if specified
  if(!is.null(remove_clusters)){
    remove.cells <- colnames(so)[so@meta.data[[remove_clusters_res]] %in% remove_clusters]
    so <- subset(so, cells = remove.cells, invert = TRUE)
  }
  
  # Normalization and scaling
  if (assay == "RNA") {
    so <- NormalizeData(so)
    so <- FindVariableFeatures(so)
    VariableFeatures(so) <- setdiff(VariableFeatures(so), bk.list)
    if (!is.null(features)) {
      VariableFeatures(so) <- unique(c(VariableFeatures(so), features))
    }
    so <- ScaleData(so, vars.to.regress = vars.to.regress)
  } else if (assay == "integrated") {
    so <- ScaleData(so, vars.to.regress = vars.to.regress)
  } else if (assay == "SCT") {
    so <- SCTransform(so, vst.flavor = "v2", vars.to.regress = vars.to.regress, verbose = verbose, seed.use = seed.use)
  }
  
  # Dimensionality reduction
  so <- RunPCA(so, npcs = 50, verbose = verbose)
  so <- RunUMAP(so, dims = 1:npcs, reduction = reduction, reduction.name = reduction.name, verbose = verbose, seed.use = seed.use)
  message("[MSG] reduction = ", reduction, "; reduction.name = ", reduction.name, "; seed.use = ", seed.use)
  
  # Clustering
  graph.name <- paste0(assay, "_snn")
  so <- FindNeighbors(so, graph.name = graph.name, dims = 1:npcs, verbose = verbose)
  for (res in resolutions) {
    so <- FindClusters(so, graph.name = graph.name, resolution = res, verbose = verbose)
  }
  
  # Save violin plots for each resolution
  pdf(file.path(save.loc, "plots", paste0(prefix, "_", capture, "_", assay, "_snn_res_VlnPlot.pdf")), height = 9, width = 12)
  plots <-  lapply(resolutions, function(res) {
    so <- SetIdent(so, value = paste0(assay, "_snn_res.", res))
    vln <- VlnPlot(so, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0, ncol = 1, raster = raster)
    title <- ggdraw() + draw_label(paste0("resolution = ", res), fontface = 'bold')
    cowplot::plot_grid(title, vln, ncol = 1, rel_heights = c(0.1, 1))
  })
  invisible(lapply(plots, print))
  dev.off()
  
  # Save clustree plot
  pdf(file.path(save.loc, "plots", paste0(prefix, "_", capture, "_clustree.pdf")), height = 12, width = 10)
  print(clustree(so, prefix = paste0(assay, "_snn_res.")))
  dev.off()
  
  return(so)
}

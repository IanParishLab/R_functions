run_scDblFinder <- function(seu.obj, sample_name, plot = FALSE, save.loc, samples = NULL, clusters = NULL) {
  
  suppressPackageStartupMessages({
    require(Seurat)
    require(SingleCellExperiment)
    require(scDblFinder)
    require(patchwork)
    require(tidyverse)
  })

  # make output directory
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)

  # Convert and run scDblFinder
  message("[MSG] Running scDblFinder...")
  sce <- as.SingleCellExperiment(seu.obj, assay = "RNA")
  sce.dbl <- scDblFinder(sce, clusters = clusters, samples = samples)

  # Transfer results to Seurat object
  seu.obj[["scDblFinder.class"]] <- sce.dbl[["scDblFinder.class"]]
  seu.obj[["scDblFinder.score"]] <- as.numeric(sce.dbl[["scDblFinder.score"]])

  if (clusters) {
    seu.obj[["scDblFinder.cluster"]] <- sce.dbl[["scDblFinder.cluster"]]
  }

  # Optional plotting
  if (plot) {
    message("[MSG] Generating scDblFinder plots...")
    group.by.clusters <- if (clusters) "scDblFinder.cluster" else NULL

    p1 <- DimPlot(seu.obj, reduction = "umap", group.by = "scDblFinder.class", order = TRUE) + 
      theme(aspect.ratio = 1)
    p2 <- DimPlot(seu.obj, reduction = "umap", group.by = group.by.clusters, order = TRUE) + 
      theme(aspect.ratio = 1)
    p3 <- FeaturePlot(seu.obj, reduction = "umap", features = "scDblFinder.score", order = TRUE, repel = TRUE) +
      theme(aspect.ratio = 1)
    p4 <- VlnPlot(seu.obj, group.by = "orig.ident", 
                  features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
                  split.by = "scDblFinder.class", ncol = 3) &
      theme(text = element_text(size = 7), axis.text = element_text(size = 7))

    # Save plots
    pdf(file.path(file.path(save.loc, "plots"), paste0(sample_name, "_scDblFinder.pdf")), width = 12, height = 5)
    print(p1 + p2 + p3)
    print(p4)
    dev.off()
    message("[MSG] Plots saved to: ", file.path(file.path(save.loc, "plots"), paste0(sample_name, "_scDblFinder.pdf")))
  }

  return(seu.obj)
}

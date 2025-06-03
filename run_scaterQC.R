run_scaterQC <- function(
    sce, 
    mito_genes = scGate::genes.blacklist.default$Mm$Mito,
    sample_name,
    nmads = c(low = 2, high = 3),
    save.loc = "QC",
    plot.width = 6,
    plot.height = 4,
    min.cells = 1,
    clusters = NULL, 
    samples = NULL,
    dry_run = TRUE
) {
  suppressPackageStartupMessages({
    library(celda)
    library(scater)
    library(scuttle)
    library(cowplot)
    library(SingleCellExperiment)
    library(Seurat)
    library(tidyverse)
    library(scDblFinder)
  })
  
  # Make output directory
  dir.create(file.path(save.loc, "metrics"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # Backup raw counts
  assay(sce, "raw_counts") <- counts(sce)
  counts(sce) <- decontXcounts(sce)
  
  # Add per-cell QC metrics
  message("[MSG] Adding QC metrics...")
  is_mito <- rownames(sce) %in% mito_genes
  sce <- addPerCellQCMetrics(sce, subsets = list(mito = is_mito))
  
  # Flag outliers
  sce$low_lib_size <- isOutlier(sce$sum, nmads = nmads["low"], type = "lower")
  sce$high_lib_size <- isOutlier(sce$sum, nmads = nmads["high"], type = "higher")
  sce$lib_size <- sce$low_lib_size | sce$high_lib_size
  sce$low_n_features <- isOutlier(sce$detected, nmads = nmads["low"], type = "lower")
  sce$high_n_features <- isOutlier(sce$detected, nmads = nmads["high"], type = "higher")
  sce$n_features <- sce$low_n_features | sce$high_n_features
  sce$high_subsets_mito_percent <- isOutlier(sce$subsets_mito_percent, nmads = nmads["high"], type = "higher")
  sce$discard <- sce$lib_size | sce$n_features | sce$high_subsets_mito_percent
  sce$percent.mito <- sce$subsets_mito_percent
  sce$zero_pct <- colMeans(counts(sce) == 0) * 100
  sce$sample <- sample_name
  
  # Plot pre-filter QC
  message("[MSG] Generating QC plots...")
  qc_theme <- theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
  
  p1_prefilter <- plot_grid(
    plotColData(sce, y = "sum", x = "sample", colour_by = "lib_size") +
      scale_y_log10() + annotation_logticks(sides = "l") + qc_theme,
    plotColData(sce, y = "detected", x = "sample", colour_by = "n_features") + qc_theme,
    plotColData(sce, y = "subsets_mito_percent", x = "sample", colour_by = "high_subsets_mito_percent") + qc_theme,
    plotColData(sce, y = "zero_pct", x = "sample", colour_by = "sample") + qc_theme,
    ncol = 4
  ) + labs(title = "Pre-filter")
  
  p2_plotdiscard <- plot_grid(
    plotColData(sce, x = "sum", y = "subsets_mito_percent", colour_by = "discard"),
    plotColData(sce, x = "sum", y = "detected", colour_by = "discard"),
    ncol = 2
  )
  
  p3_plotHighestExprs <- plotHighestExprs(sce, exprs_values = "counts", colour_cells_by = "detected") + qc_theme
  
  # Backup before filtering
  sce_dry_run <- sce
  
  # Filter cells and genes
  message("[MSG] Filtering mito genes and discarded cells...")
  sce <- sce[!is_mito, !sce$discard]
  
  # Detect doublets
  message("[MSG] Running scDblFinder...")
  sce <- sce %>%
    logNormCounts() %>%
    runPCA() %>%
    runUMAP()
  sce <- scDblFinder(sce, clusters = clusters, samples = samples)
  
  p4_postfilter <- plot_grid(
    plotColData(sce, y = "sum", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") +
      scale_y_log10() + annotation_logticks(sides = "l") + qc_theme,
    plotColData(sce, y = "detected", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + qc_theme,
    plotColData(sce, y = "subsets_mito_percent", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + qc_theme,
    plotColData(sce, y = "zero_pct", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + qc_theme,
    ncol = 4
  ) + labs(title = "Post-filter")
  
  p5_scDblFinder <- plot_grid(
    plotColData(sce, y = "scDblFinder.score", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + qc_theme,
    plotUMAP(sce, colour_by = "scDblFinder.class"),
    ncol = 2
  ) + labs(title = "scDblFinder")
  
  # Save plots
  message("[MSG] Saving plots...")
  pdf(file.path(save.loc, "plots", paste0("addPerCellQCMetrics.scDblFinder.", sample_name, ".pdf")), 
      width = plot.width, height = plot.height)
  print(p1_prefilter)
  print(p2_plotdiscard)
  print(p3_plotHighestExprs)
  print(p4_postfilter)
  print(p5_scDblFinder)
  dev.off()
  
  # Convert to Seurat and save
  message("[MSG] Converting to Seurat object...")
  seu.obj <- CreateSeuratObject(counts(sce),
                                meta.data = as.data.frame(colData(sce)),
                                min.cells = min.cells,
                                project = sample_name)
  
  saveRDS(seu.obj, file.path(save.loc, "int_obj", paste0(seu.obj@project.name, ".preprocessed.seu.obj.rds")))
  
  # Return object
  return(if (dry_run) sce_dry_run else sce)
}

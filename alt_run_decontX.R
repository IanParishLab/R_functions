run_decontX <- function(
  raw = NULL, 
  filt,
  sampleName,
  Read10X = TRUE,
  min.cells = 1,
  rna_density_plot_markers = NULL,
  use_decontPro = TRUE,
  save.loc = "1_decontX",
  seed.use = 42
) {
  suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(decontX)
    library(Seurat)
    library(Matrix)
    library(tidyverse)
  })

  # ------------------------- #
  # Directory setup
  # ------------------------- #
  dir.create(save.loc, showWarnings = FALSE, recursive = TRUE)
  lapply(c("int_obj", "plots"), function(d) {
    dir.create(file.path(save.loc, d), showWarnings = FALSE, recursive = TRUE)
  })

  message(glue::glue("[INFO] Processing {sampleName}"))

  # ------------------------- #
  # Read Data
  # ------------------------- #
  if (Read10X) {
    message("[INFO] Reading matrices with Read10X...")
    if (!is.null(raw)) raw <- Read10X(raw)
    filt <- Read10X(filt)
  }

  if (use_decontPro) {
    raw.counts <- if (!is.null(raw)) raw$`Gene Expression` else NULL
    filt.counts <- filt$`Gene Expression`
  } else {
    raw.counts <- raw
    filt.counts <- filt
  }

  # ------------------------- #
  # Filter low-quality cells
  # ------------------------- #
  message("[INFO] Filtering out genes with low counts...")
  filt.counts <- filt.counts[rowSums(filt.counts > 0) >= min.cells, ]
  filt.counts <- filt.counts[, colSums(filt.counts) > 0]

  if (!is.null(raw.counts)) {
    raw.counts <- raw.counts[rownames(filt.counts), colnames(filt.counts)]
  }

  # ------------------------- #
  # Run decontX
  # ------------------------- #
  message("[INFO] Running decontX...")
  sce <- SingleCellExperiment(list(counts = filt.counts))
  if (!is.null(raw.counts)) {
    sce.raw <- SingleCellExperiment(list(counts = raw.counts))
    sce <- decontX::decontX(sce, background = sce.raw)
  } else {
    sce <- decontX::decontX(sce)
  }

  saveRDS(sce, file.path(save.loc, "int_obj", paste0(sampleName, ".RNA_decontX_sce.rds")))
  message("[INFO] Saved decontaminated SCE.")

  # ------------------------- #
  # Plotting
  # ------------------------- #
  if (is.null(rna_density_plot_markers)) {
    message("[INFO] Marker genes not provided; using first two genes.")
    rna_density_plot_markers <- head(rownames(sce), 2)
  }

  message("[INFO] Plotting contamination and gene density...")
  original.filt.counts <- if (use_decontPro) filt$`Gene Expression` else filt

  p1 <- plotDecontXContamination(sce) + ggtitle(paste(ncol(sce), "cells"))
  p2 <- plotDensity(original.filt.counts,
                    round(decontXcounts(sce)),
                    rna_density_plot_markers)

  pdf(file.path(save.loc, "plots", paste0(sampleName, ".rna.umap.pdf")), height = 6, width = 6)
  print(p1)
  print(p2)
  dev.off()
  message("[INFO] Plots saved.")

  invisible(sce)
}

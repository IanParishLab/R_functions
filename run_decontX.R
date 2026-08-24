# decontx 
run_decontX <- function(
    raw = NULL, 
    filt,
    sample_name,
    Read10X = TRUE,
    min.cells = 1,
    rna_density_plot_markers = NULL,
    save.loc = "1_decontX",
    seed.use = seed.use
){
  
  suppressPackageStartupMessages({
    # library(celda)
    library(decontX)
    library(SingleCellExperiment)
    library(Seurat)
    library(tidyverse)
  })
  
  # make output directory
  dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  ############################################################
  # i=length(capture)
  # raw = raw.list[[i]]
  # filt = filt.list[[i]]
  # sample_name = capture[i]
  # Read10X = TRUE
  # min.cells = 1
  # rna_density_plot_markers = c("CD8A","CD8B")
  # save.loc = here("CD8","1_decontX")
  # seed.use = seed.use
  ############################################################
  
  # read files
  if (Read10X) {
    message("[INFO] Reading matrices with Read10X...")
    if (!is.null(raw)) raw <- Read10X(raw)
    filt <- Read10X(filt)
  }
  
  raw.counts <- if (!is.null(raw) & "Gene Expression" %in% names(raw)) raw$`Gene Expression` else raw
  filt.counts <- if ("Gene Expression" %in% names(filt)) filt$`Gene Expression` else filt
  
  # filter GEX matrix
  message("[MSG] Filtering out genes with low counts...")
  filt.counts <- filt.counts[which(rowSums(filt.counts > 0) >= min.cells),]
  filt.counts <- filt.counts[, which(colSums(filt.counts) > 0) ]
  
  if (!is.null(raw.counts)) {
    raw.counts <- raw.counts[rownames(filt.counts), colnames(filt.counts)]
  }
  
  # decontX
  sce <- SingleCellExperiment(list(counts = filt.counts))
  
  if (!is.null(raw.counts)) {
    message("[INFO] Running decontX...")
    sce.raw <- SingleCellExperiment(list(counts = raw.counts))
    sce <- decontX::decontX(sce, background = sce.raw)
  } else {
    message("[INFO] Running decontX without raw counts as background..")
    sce <- decontX::decontX(sce)
  }
  
  # save decontaminated RNA counts as SCE
  message("[MSG] Save decontaminated SCE ...")
  qs_save(sce, file.path(save.loc, "int_obj", paste0(sample_name, ".RNA_decontX_sce.qs")))
  
  # plot decontX contamination
  message("[INFO] Plotting contamination and gene density...")
  original.filt.counts <- if ("Gene Expression" %in% names(filt)) filt$`Gene Expression` else filt
  
  if (is.null(rna_density_plot_markers)) {
    message("[INFO] Marker genes not provided; using first two genes.")
    rna_density_plot_markers <- head(rownames(sce), 2)
  }
  
  p1 <- plotDecontXContamination(sce) + ggtitle(paste(ncol(sce), "cells"))
  p2 <- plotDensity(original.filt.counts,
                    round(decontXcounts(sce)),
                    rna_density_plot_markers)
  
  pdf(file.path(save.loc, "plots", paste0(sample_name, ".rna.umap.pdf")), height = 6, width = 6)
  print(p1)
  print(p2)
  dev.off()
  
  invisible(sce)
  
  
  # decontPro
  # if (isTRUE(use_decontPro)){
  # print(paste0("Prep decontPro..."))
  # if(!is.null(raw)){
  #   raw.hto <- raw$`Antibody Capture`[,colnames(sce)]
  #   raw.hto <- raw.hto[,which(colSums(raw.hto) > 0)]
  # }
  # filt.hto <- filt$`Antibody Capture`[,colnames(sce)]
  # filt.hto <- filt.hto[,which(colSums(filt.hto) > 0)]
  #
  # hto.obj <- CreateSeuratObject(filt.hto, assay = "HTO")
  #
  # npc = nrow(hto.obj) - 1
  # hto.obj <- NormalizeData(hto.obj, normalization.method = "CLR", margin = 2) %>%
  #   ScaleData(assay = "HTO") %>%
  #   RunPCA(assay = "HTO", features = rownames(hto.obj), npcs = npc, reduction.name = "pca_hto") %>%
  #   FindNeighbors(dims = 1:npc, assay = "HTO", reduction = "pca_hto") %>%
  #   FindClusters(resolution = 0.5)
  #
  # hto.obj <- RunUMAP(hto.obj,
  #                    dims = 1:npc,
  #                    assay = "HTO",
  #                    reduction = "pca_hto",
  #                    reduction.name = "umap_hto",
  #                    seed.use = seed.use,
  #                    verbose = FALSE)
  #
  # # decontPro
  # clusters <- as.integer(Idents(hto.obj))
  # counts <- as.matrix(filt.hto)
  #
  # if(!is.null(raw)){
  #   print("[MSG] decontPro...")
  #   decont.HTO <- decontPro(counts, clusters, ambient_counts = raw.hto)
  # } else {
  #   print("[MSG] decontPro without raw HTO counts as background...")
  #   decont.HTO <- decontPro(counts, clusters, ambient_counts = NULL)
  # }
  # # save HTO count object
  # print("[MSG] Saving decontPro HTO sce ...")
  # saveRDS(decont.HTO, file.path(save.loc, "int_obj", paste0(sample_name,".HTO_decontPro_sce.rds")))
  #
  # # plots
  # p1 <- DimPlot(hto.obj, reduction = "umap_hto", label = TRUE)
  # # before decontaminating cell barcodes
  # p2 <- plotDensity(counts,
  #                   decont.HTO$decontaminated_counts,
  #                   rownames(counts))
  # pdf(file.path(save.loc, "plots", paste0(sample_name, ".hto.density.pdf")), height = 6, width = 6)
  # print(p1)
  # print(p2)
  # dev.off()
  #
  # # subset for `hto.obj` cells too if decontpro is used
  # sce <- sce[,colnames(hto.obj)]
  # }
  
}

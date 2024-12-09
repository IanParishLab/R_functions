run_scaterQC <- function(
    sce, 
    mito_genes = scGate::genes.blacklist.default$Mm$Mito,
    sampleName,
    nmads = c(low = 2,high = 3),
    save.loc = "QC",
    plot.width = 6,
    plot.height = 4,
    seqType = "10X 5'scRNAseq",
    dropletType = "NextGEM",
    min.cells = 1,
    clusters, samples = NULL,
    dry_run = TRUE
){
  
  suppressPackageStartupMessages({
    require(celda)
    require(scater)
    require(scuttle)
    require(cowplot)
    require(SingleCellExperiment)
    require(Seurat)
    require(tidyverse)
    require(scDblFinder)
  })
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  lapply(c("metrics","int_obj","plots"),function(d){
    ifelse(!dir.exists(file.path(save.loc,d)),
           dir.create(file.path(save.loc,d), recursive = TRUE), paste0(d," directory exists"))
  })
  # SingleCellExperiment processing
  print(paste0("[MSG] processing SingleCellExperiment object..."))
  # create assay
  assay(sce, "raw_counts") <- counts(sce)
  counts(sce) <- decontXcounts(sce)
  
  print(paste0("[MSG] addPerCellQCMetrics: `sum`, `detected`, `zero_pct`, `subsets_mito_percent`..."))
  # add is_mito
  is_mito <- rownames(sce) %in% mito_genes
  sce <- addPerCellQCMetrics(sce, subsets = list(mito = is_mito))
  
  sce$low_lib_size <- isOutlier(sce$sum, nmads = nmads['low'], type = "lower")
  sce$high_lib_size <- isOutlier(sce$sum, nmads = nmads['high'], type = "higher")
  sce$lib_size <- (sce$low_lib_size|sce$high_lib_size) == TRUE
  sce$low_n_features <- isOutlier(sce$detected, nmads = nmads['low'], type = "lower")
  sce$high_n_features <- isOutlier(sce$detected, nmads = nmads['high'], type = "higher")
  sce$n_features <- (sce$low_n_features|sce$high_n_features) == TRUE
  sce$high_subsets_mito_percent <- isOutlier(sce$subsets_mito_percent, nmads = nmads['high'], type = "higher")
  sce$discard <- (sce$lib_size|sce$n_features|sce$high_subsets_mito_percent) == TRUE
  
  # add percent.mito column for streamlined pipeline use
  sce$percent.mito <- sce$subsets_mito_percent
  
  # calculate %cells with zero counts for each gene.
  sce$zero_pct <- colMeans(counts(sce) == 0)*100
  
  # add sample name
  sce$sample <- sampleName
  
  # plotColData
  print(paste0("[MSG] plotColData..."))
  # set theme
  theme = theme(legend.position = "none",
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
  
  #make plots
  p1_prefilter <- plot_grid(
    plotColData(sce, y = "sum", x = "sample", colour_by = "lib_size") +
      scale_y_log10() +
      annotation_logticks(sides = "l", short = unit(0.03, "cm"), mid = unit(0.06, "cm"), long = unit(0.09, "cm")) + theme,
    plotColData(sce, y = "detected", x = "sample", colour_by = "n_features") + theme,
    plotColData(sce, y = "subsets_mito_percent", x = "sample", colour_by = "high_subsets_mito_percent") + theme,
    plotColData(sce, y = "zero_pct", x = "sample", colour_by = "sample") + theme,
    ncol = 4,
    align = "hv"
  ) + labs(title = "Pre-filter")
  
  # plot discard cells
  p2_plotdiscard <- plot_grid(
    plotColData(sce, x="sum", y="subsets_mito_percent", colour_by="discard"),
    plotColData(sce, x="sum", y="detected", colour_by="discard"),
    ncol=2
  )
  
  p3_plotHighestExprs <- plotHighestExprs(sce, 
                                          exprs_values = "counts", 
                                          feature_names_to_plot = NULL, # rownames(sce) is default
                                          colour_cells_by="detected") + theme
  
  # filter discard cells from sce
  print(paste0("[MSG] Filter `mito_genes` & `discard` cells..."))
  # sce <- sce[!rowData(sce)$discard,!colData(sce)$discard]
  
  # back up sce
  sce_dry_run <- sce
  
  # filter mito genes, discard low quality cells
  sce <- sce[!is_mito, !colData(sce)$discard]
  
  # # detect doublets
  print(paste0("[MSG] Detect scDblFinder doublets..."))
  sce <- logNormCounts(sce) %>%
    runPCA() %>%
    runUMAP()
  sce <- scDblFinder(sce, clusters = clusters, samples = samples)
  
  p4_postfilter <- plot_grid(
    plotColData(sce, y = "sum", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") +
      scale_y_log10() +
      annotation_logticks(sides = "l", short = unit(0.03, "cm"), mid = unit(0.06, "cm"), long = unit(0.09, "cm")) + theme,
    plotColData(sce, y = "detected", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    plotColData(sce, y = "subsets_mito_percent", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    plotColData(sce, y = "zero_pct", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    ncol = 4, 
    align = "hv"
  ) + labs(title = "Post-filter")
  
  p5_scDblFinder <- plot_grid(
    plotColData(sce, y = "scDblFinder.score", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    plotUMAP(sce, colour_by = "scDblFinder.class"),
    ncol = 2, 
    align = "hv"
  ) + labs(title = "scDblFinder")
  
  # save plots
  print(paste0("[MSG] Save plots..."))
  pdf(file.path(save.loc, "plots", paste0("addPerCellQCMetrics.scDblFinder.",sampleName,".pdf")), width = plot.width, height = plot.height)
  p1_prefilter %>% print
  p2_plotdiscard %>% print
  p3_plotHighestExprs %>% print
  p4_postfilter %>% print
  p5_scDblFinder %>% print
  dev.off()
  
  # # save preprocessed sce
  # print(paste0("[MSG] Save preprocessed sce..."))
  # saveRDS(sce,  file.path(save.loc, "int_obj", paste0("preprocessed.sce.",sampleName,".rds")))
  
  print("[MSG] Create Seurat Object... ")
  seu.obj <- CreateSeuratObject(counts(sce),
                                meta.data = as.data.frame(colData(sce)),
                                min.cells = min.cells,
                                project = sampleName)
  seu.obj[["seqType"]] <- seqType
  seu.obj[["dropletType"]] <- dropletType
  
  # Save preprocessed objects
  print("Saving preprocessed Seurat object...")
  saveRDS(seu.obj,  file.path(save.loc, "int_obj", paste0("preprocessed.seu.obj.",sampleName,".rds")))
  
  if(isTRUE(dry_run)){
    return(sce_dry_run)
  } else {
    return(sce)
  }
}

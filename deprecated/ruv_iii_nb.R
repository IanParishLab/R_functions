run_scaterQC <- function(
    sce, 
    mito_genes = scGate::genes.blacklist.default$Mm$Mito,
    sampleName,
    nmads = c(low = 2,high = 3),
    save.loc = "QC",
    plot.width = 6,
    plot.height = 4,
    min.cells = 1,
    clusters, samples = NULL,
    dry_run = TRUE,
    genome = "mouse",
    low_abundant_genes_limit = -5
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
    require(ruvIIInb)
    require(edgeR)
    require(scMerge)
    require(igraph)
    require(DelayedArray)
    require(dittoSeq)
  })
  
  # QC Directory structure
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  lapply(c("metrics","int_obj","plots"),function(d){
    ifelse(!dir.exists(file.path(save.loc,d)),
           dir.create(file.path(save.loc,d), recursive = TRUE), paste0(d," directory exists"))
  })
  
  print(paste0("[MSG] processing SingleCellExperiment object..."))
  assay(sce, "raw_counts") <- counts(sce)
  counts(sce) <- decontXcounts(sce)
  
  print(paste0("[MSG] addPerCellQCMetrics..."))
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
  
  sce$percent.mito <- sce$subsets_mito_percent
  sce$zero_pct <- colMeans(counts(sce) == 0)*100
  sce$sample <- sampleName
  
  print(paste0("[MSG] plotColData..."))
  theme = theme(legend.position = "none",
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))
  
  p1_prefilter <- plot_grid(
    plotColData(sce, y = "sum", x = "sample", colour_by = "lib_size") + scale_y_log10() + annotation_logticks(sides = "l") + theme,
    plotColData(sce, y = "detected", x = "sample", colour_by = "n_features") + theme,
    plotColData(sce, y = "subsets_mito_percent", x = "sample", colour_by = "high_subsets_mito_percent") + theme,
    plotColData(sce, y = "zero_pct", x = "sample", colour_by = "sample") + theme,
    ncol = 4
  ) + labs(title = "Pre-filter")
  
  p2_plotdiscard <- plot_grid(
    plotColData(sce, x="sum", y="subsets_mito_percent", colour_by="discard"),
    plotColData(sce, x="sum", y="detected", colour_by="discard"),
    ncol=2
  )
  
  p3_plotHighestExprs <- plotHighestExprs(sce, exprs_values = "counts", colour_cells_by="detected") + theme
  
  sce_dry_run <- sce
  sce <- sce[!is_mito, !colData(sce)$discard]
  
  print(paste0("[MSG] Detect scDblFinder doublets..."))
  sce <- logNormCounts(sce) %>% runPCA() %>% runUMAP()
  sce <- scDblFinder(sce, clusters = clusters, samples = samples)
  
  p4_postfilter <- plot_grid(
    plotColData(sce, y = "sum", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + scale_y_log10() + annotation_logticks(sides = "l") + theme,
    plotColData(sce, y = "detected", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    plotColData(sce, y = "subsets_mito_percent", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    plotColData(sce, y = "zero_pct", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    ncol = 4
  ) + labs(title = "Post-filter")
  
  p5_scDblFinder <- plot_grid(
    plotColData(sce, y = "scDblFinder.score", x = "sample", colour_by = "scDblFinder.class", shape_by = "scDblFinder.class") + theme,
    plotUMAP(sce, colour_by = "scDblFinder.class"),
    ncol = 2
  ) + labs(title = "scDblFinder")
  
  print(paste0("[MSG] Save plots..."))
  pdf(file.path(save.loc, "plots", paste0("addPerCellQCMetrics.scDblFinder.",sampleName,".pdf")), width = plot.width, height = plot.height)
  p1_prefilter %>% print
  p2_plotdiscard %>% print
  p3_plotHighestExprs %>% print
  p4_postfilter %>% print
  p5_scDblFinder %>% print
  dev.off()
  
  print("[MSG] Create Seurat Object...")
  seu.obj <- CreateSeuratObject(counts(sce), meta.data = as.data.frame(colData(sce)), min.cells = min.cells, project = sampleName)
  
  print("Saving preprocessed Seurat object...")
  saveRDS(seu.obj,  file.path(save.loc, "int_obj", paste0("preprocessed.seu.obj.",sampleName,".rds")))
  
  print("[MSG] Running RUV-III-NB normalization...")
  sce <- addPerFeatureQCMetrics(sce)
  sce <- subset(sce, rowData(sce)$mean > 0 )
  lowcount_drop <- log(rowData(sce)$mean) < low_abundant_genes_limit
  sce <- sce[!lowcount_drop, ]
  
  data(segList)
  ctl <- if (genome == "mouse") segList$mouse$mouse_scSEG else segList$human$human_scSEG
  rowData(sce)$ctlLogical <- rownames(assays(sce)$counts) %in% ctl
  
  sce <- computeSumFactors(sce,assay.type="counts")
  data_norm_pre <- sweep(assays(sce)$counts,2,sce$sizeFactor,'/')
  assays(sce, withDimnames=FALSE)$lognormcounts<- log(data_norm_pre+1)
  snn_gr_init <- buildSNNGraph(sce, assay.type = "lognormcounts")
  clusters_init <- igraph::cluster_louvain(snn_gr_init)
  sce$cluster_init <- factor(clusters_init$membership)
  
  sce <- runUMAP(sce,exprs_values = "lognormcounts")
  sce <- runPCA(sce,exprs_values = "lognormcounts")
  
  png(file.path(save.loc, "plots", paste0("dittoDimPlot_cluster_init.", sampleName, ".png")), width=10, height=10, units="in", res=300)
  print(dittoDimPlot(sce, "cluster_init"))
  dev.off()
  
  sce$logLS <- log(colSums(assays(sce)$counts))
  png(file.path(save.loc, "plots", paste0("plotUMAP_cluster_init_logLS.", sampleName, ".png")), width=10, height=10, units="in", res=300)
  plotUMAP(sce,colour_by='logLS')
  dev.off()
  
  M <- matrix(0,ncol(assays(sce)$counts),length(unique(sce$cluster_init)))
  cl <- sort(unique(as.numeric(unique(sce$cluster_init))))
  for(CL in cl){ M[which(as.numeric(sce$cluster_init)==CL),CL] <- 1 }
  
  ruv3nb_out <- fastruvIII.nb(Y=DelayedArray(assays(sce)$counts), M=M, ctl=rowData(sce)$ctlLogical, k=2, ncores = 6)
  sce_ruv3nb <- makeSCE(ruv3nb_out,cData=colData(sce))
  seurat_ruv3nb <- as.Seurat(sce_ruv3nb, counts = "counts", data ="logPAC")
  
  saveRDS(ruv3nb_out, file.path(save.loc, "int_obj", paste0("ruv3nb_out.", sampleName, ".rds")))
  saveRDS(sce_ruv3nb, file.path(save.loc, "int_obj", paste0("sce_ruv3nb.", sampleName, ".rds")))
  saveRDS(seurat_ruv3nb, file.path(save.loc, "int_obj", paste0("seurat_ruv3nb.", sampleName, ".rds")))
  
  if(isTRUE(dry_run)) return(sce_dry_run) else return(sce_ruv3nb)
}

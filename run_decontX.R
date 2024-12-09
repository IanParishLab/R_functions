# decontx 
run_decontX <- function(
    raw = NULL, 
    filt,
    sampleName,
    Read10X = TRUE,
    min.cells = 1,
    rna_density_plot_markers = NULL,
    use_decontPro = TRUE,
    save.loc = "decontX",
    seed.use = seed.use
){
  
  suppressPackageStartupMessages({
    require(celda)
    require(decontX)
    require(SingleCellExperiment)
    require(Seurat)
    require(tidyverse)
  })
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  lapply(c("int_obj","plots"),function(d){
    ifelse(!dir.exists(file.path(save.loc,d)),
           dir.create(file.path(save.loc,d), recursive = TRUE), paste0(d," directory exists"))
  })
  
  # read files
  if (Read10X == TRUE){
    print("[MSG] Using Read10X...")
    
    if(!is.null(raw)){
      raw <- Read10X(raw)
    }
    
    filt <- Read10X(filt)
  }
  
  if (isTRUE(use_decontPro)){
    if(!is.null(raw)){
      raw.counts <- raw$`Gene Expression`
    }
    filt.counts <- filt$`Gene Expression`
  } else if(isFALSE(use_decontPro)){
    if(!is.null(raw)){
      raw.counts <- raw
    }
    filt.counts <- filt
  }
  
  # filter GEX matrix
  print("[MSG] Filter cells with zero counts across all genes...")
  filt.counts <- filt.counts[which(rowSums(filt.counts > 0) >= min.cells), ]
  filt.counts <- filt.counts[,which(colSums(filt.counts) > 0)]
  if(!is.null(raw)){
    raw.counts <- raw.counts[rownames(filt.counts), colnames(filt.counts)]
  }
  
  # create sce object
  sce <- SingleCellExperiment(list(counts = filt.counts))
  
  # decontX
  if(!is.null(raw)){
    sce.raw <- SingleCellExperiment(list(counts = raw.counts))
    print("[MSG] decontX...")
    sce <- decontX(sce, background = sce.raw)
  } else {
    print("[MSG] decontX without raw counts as background...")
    sce <- decontX(sce)
  }
  
  # save decontaminated RNA counts as SCE
  saveRDS(sce, file.path(save.loc, "int_obj", paste0(sampleName,".RNA_decontX_sce.rds")))
  
  # plot decontX contamination
  if(is.null(rna_density_plot_markers)) {
    print(paste0("[MSG] Setting first 2 markers to plot: ", 
                 paste0(rownames(sce)[1:2],collapse = ",")
    ))
    
    rna_density_plot_markers <- rownames(sce)[1:2]
  } 
  
  print(paste0("[MSG] plot DecontX..."))
  if(isTRUE(use_decontPro)){
    original.filt.counts <- filt$`Gene Expression`
  } else if (isFALSE(use_decontPro)){
    original.filt.counts <- filt
  }
  
  p1 <- plotDecontXContamination(sce) +
    ggtitle(paste0(ncol(sce), " cells"))
  p2 <- plotDensity(original.filt.counts,
                    round(decontXcounts(sce)),
                    rna_density_plot_markers 
  )
  
  pdf(file.path(save.loc, "plots", paste0(sampleName, ".rna.umap.pdf")), height = 6, width = 6)
  print(p1)
  print(p2)
  dev.off()
  
  # decontPro
  if (isTRUE(use_decontPro)){
    print(paste0("Prep decontPro..."))
    if(!is.null(raw)){
      raw.hto <- raw$`Antibody Capture`[,colnames(sce)]
      raw.hto <- raw.hto[,which(colSums(raw.hto) > 0)]
    }
    filt.hto <- filt$`Antibody Capture`[,colnames(sce)]
    filt.hto <- filt.hto[,which(colSums(filt.hto) > 0)]
    
    hto.obj <- CreateSeuratObject(filt.hto, assay = "HTO")
    
    npc = nrow(hto.obj) - 1
    hto.obj <- NormalizeData(hto.obj, normalization.method = "CLR", margin = 2) %>%
      ScaleData(assay = "HTO") %>%
      RunPCA(assay = "HTO", features = rownames(hto.obj), npcs = npc, reduction.name = "pca_hto") %>%
      FindNeighbors(dims = 1:npc, assay = "HTO", reduction = "pca_hto") %>%
      FindClusters(resolution = 0.5)
    
    hto.obj <- RunUMAP(hto.obj,
                       dims = 1:npc,
                       assay = "HTO",
                       reduction = "pca_hto",
                       reduction.name = "umap_hto",
                       seed.use = seed.use,
                       verbose = FALSE)
    
    # decontPro
    clusters <- as.integer(Idents(hto.obj))
    counts <- as.matrix(filt.hto)
    
    if(!is.null(raw)){
      print("[MSG] decontPro...")
      decont.HTO <- decontPro(counts, clusters, ambient_counts = raw.hto)
    } else {
      print("[MSG] decontPro without raw HTO counts as background...")
      decont.HTO <- decontPro(counts, clusters, ambient_counts = NULL)
    }
    # save HTO count object
    saveRDS(decont.HTO, file.path(save.loc, "int_obj", paste0(sampleName,".HTO_decontPro_sce.rds")))
    
    # plots
    p1 <- DimPlot(hto.obj, reduction = "umap_hto", label = TRUE)
    # before decontaminating cell barcodes 
    p2 <- plotDensity(counts,
                      decont.HTO$decontaminated_counts,
                      rownames(counts))
    pdf(file.path(save.loc, "plots", paste0(sampleName, ".hto.density.pdf")), height = 6, width = 6)
    print(p1)
    print(p2)
    dev.off()
  }
  
  sce <- sce[,colnames(hto.obj)]
  
  # save object
  print("[MSG] Saving final object for run_decontX...")
  saveRDS(sce, file.path(save.loc, "int_obj", paste0(sampleName,".run_decontX.sce.rds")))
}

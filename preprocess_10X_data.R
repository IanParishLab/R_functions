# --- load_so_from_stage function  ------------------------------------------------------------- 
load_so_from_stage <- function(capture, stage, save.loc, pattern) {
  message("[INFO] Loading Seurat objects from: ", stage)
  
  lapply(capture, function(name) {
    
    # build base path (without extension)
    base_path <- here(save.loc, stage, "int_obj", paste(c(name, pattern), collapse = "."))
    
    qs_path  <- paste0(base_path, ".qs")
    rds_path <- paste0(base_path, ".rds")
    
    if (file.exists(qs_path)) {
      qs2::qs_read(qs_path)
    } else if (file.exists(rds_path)) {
      readRDS(rds_path)
    } else {
      stop("[ERROR] No `.qs` or `.rds` file found for: ", name)
    }
    
  }) %>% setNames(., capture)
}

# --- preprocess_10X_data function  ------------------------------------------------------------- 
preprocess_10X_data <- function(
    project.loc,
    species = "human",
    seed.use, min.pc, bk.list = NULL, vars.to.regress = "percent.mito",
    exclude_HTOs = NULL, save.loc = project.loc, dry_run = TRUE,
    overwrite = FALSE, parallel = FALSE, workers = 4, 
    steps = list(decontX = FALSE, normalise = FALSE, demux = FALSE,
                 filtering_singlets = FALSE, cluster_qc = FALSE)
) {
  
  # large function, basically includes run_decontX, run_scaterQC, prep_HTODemux, plot_HTODemux, find_clusters_to_remove
  
  # --- load packages --------------------------------------------------------------------- 
  setwd(project.loc)
  suppressPackageStartupMessages({
    library(BiocParallel)
    library(here)
    library(SeuratExtend)
    library(qs2)
  })
  source(here("0_pipeline_setup.R"))
  set.seed(seed.use)
  message("seed.use = ", seed.use)
  
  # --- read in Cell Ranger outputs --------------------------------------------------------------------- 
  message("[INFO] Obtaining capture/sample_name from named vector `folderName`")
  if (all(names(raw.list) == names(filt.list))) {
    capture <- names(raw.list)
    message("Processing captures: ", paste0(capture, collapse = ","))
  } else {
    stop("`raw.list` and `filt.list` must be named and must be in the same order")
  }
  
  # --- BiocParallel --------------------------------------------------------------------- 
  bpparam <- if (parallel) MulticoreParam(workers = workers) else SerialParam()
  
  # --- decontX --------------------------------------------------------------------------
  if (!steps$decontX) {
    message("[SKIP] Skipping decontX step...")
  } else {
    message("[INFO] Running decontX and recording pre-decontX cell counts...")
    bplapply(seq_along(capture), function(i) {
      filt <- Read10X(filt.list[[i]])
      filt.counts <- if ("Gene Expression" %in% names(filt)) filt$`Gene Expression` else filt
      
      get_cell_counts(data = filt.counts, sample_name = capture[i], stage = "Original filtered GEX matrix", save.loc = save.loc) # cell counts
      
      run_decontX(
        raw = raw.list[[i]],
        filt = filt.list[[i]],
        sample_name = capture[i],
        Read10X = TRUE,
        min.cells = 1,
        rna_density_plot_markers = NULL,
        save.loc = here(save.loc, "1_decontX"),
        seed.use = seed.use
      )
    }, BPPARAM = bpparam)
  }
  
  # --- scater & normalise --------------------------------------------------------------------------
  if (!steps$normalise) {
    message("[SKIP] Skipping normalise step...")
  } else {
    message("[INFO] Running scaterQC, cell cycle, normalization, and scaling...")
    so <- lapply(seq_along(capture), function(i) {
      
      base_path <- here(save.loc,"1_decontX","int_obj",paste0(capture[i], ".RNA_decontX_sce"))
      qs_path  <- paste0(base_path, ".qs")
      rds_path <- paste0(base_path, ".rds")
      
      if (file.exists(qs_path)) {
        sce <- qs2::qs_read(qs_path)
      } else if (file.exists(rds_path)) {
        sce <- readRDS(rds_path)
      } else {
        stop("[ERROR] No `.qs` or `.rds` file found for capture: ",capture[i])
      }
    
    get_cell_counts(data = sce, sample_name = capture[i], stage = "Filtered GEX matrix", save.loc = save.loc) # cell counts
    
    if(species == "human") {
      mito_genes <- scGate::genes.blacklist.default$Hs$Mito 
    } else if(species == "mouse") {
      mito_genes <- scGate::genes.blacklist.default$Mm$Mito
    } else {
      stop("Human or mouse `mito_genes` only.")
    }
    
    sce_filt <- run_scaterQC(
      sce = sce,
      nmads = c(low = 2, high = 3),
      mito_genes = mito_genes,
      sample_name = capture[i],
      save.loc = here(save.loc, "2_QC"),
      plot.width = 12,
      plot.height = 4,
      min.cells = 3,
      dry_run = dry_run
    )
    get_cell_counts(data = sce_filt, sample_name = capture[i], stage = "Remove library outliers", save.loc = save.loc) # cell counts
    
    # read in preprocessed object
    seu <- qs_read(here(save.loc, "2_QC", "int_obj", paste0(capture[i],".preprocessed.seu.obj.qs")))
    
    # get species specific cc genes
    if(species == "mouse"){
      ccgenes <- lapply(Seurat::cc.genes.updated.2019, function(genes) {
        HumanToMouseGenesymbol(genes, keep.seq = TRUE) %>% unname %>% na.omit %>% as.vector
      })
    } else if (species == "human"){
      ccgenes <- Seurat::cc.genes.updated.2019
    }

    message("[INFO] plot_CellCycleRegression...")
    seu <- plot_CellCycleRegression(seu, sample_name = capture[i], ccgenes = ccgenes, save.loc = here(save.loc,"2_QC"), npcs = 50)

    message("[INFO] warning: skipping checkRegressGenes...")
    # regress_genes <- get_blacklist(seu, species = species, unlist = FALSE)
    # regress_genes <- regress_genes[lengths(regress_genes) > 0]
    # seu <- checkRegressGenes(seu, sample_name = capture[i], pca_npcs = 50, ccgenes = ccgenes,
    #                          regress_genes = regress_genes, label_genes = NULL, save.loc = here(save.loc,"2_QC"))

    # get species specific `bk.list` genes
    bk.list <- if (is.null(bk.list)) get_blacklist(seu, species = species) %>% unlist() else NULL

    message("[INFO] normalizeAndScaleData...")
    seu <- normalizeAndScaleData(seu, dry_run = dry_run, sample_name = capture[i], min.pc = min.pc,
                                 seed.use = seed.use, save.loc = here(save.loc,"2_QC"),
                                 vars.to.regress = vars.to.regress, bk.list = bk.list, verbose = FALSE)
    
    return(seu)
  }) %>% setNames(.,capture)

# save objects
lapply(seq_along(capture), function(i) {
  message("[INFO] save normalised objects...")
  qs_save(so[[capture[i]]], here(save.loc, "2_QC","int_obj", paste0(so[[capture[i]]]@project.name, ".normalised.qs")))
})
}

# --- demux --------------------------------------------------------------------------
if (!steps$demux) {
  message("[SKIP] Skipping demux step...")
} else {
  
  if (!exists("so") || is.null(so) || length(so) == 0) {
    so <- load_so_from_stage(capture, "2_QC", save.loc, pattern = "preprocessed.seu.obj")
  }
  
  message("[INFO] Running HTO demultiplexing...")
  so <- lapply(seq_along(capture), function(i) {
    if (all(names(exclude_HTOs) == capture)) {
      message("Processing capture(s): ", paste0(capture[i], collapse = ","))
    } else {
      stop("`exclude_HTOs` must be named and must be in the same order as `capture`, or `names(raw.list)`")
    }
    
    seu <- prep_HTODemux(
      seu.obj = so[[capture[i]]], HTO.counts.path = filt.list[[i]], sample_name = capture[i], 
      HTO_AssayName = "HTO", RNA_AssayName = "RNA", 
      exclude_samples = exclude_HTOs[[i]], demux = TRUE
    )
    
    plot_HTODemux(
      seu, sample_name = capture[i], save.loc = here(save.loc, "3_HTODemux"),
      HTO_AssayName = "HTO", RNA_AssayName = "RNA", 
      hash_class_column = "HTO_classification.global", 
      hash.ident = "HTO_maxID", ncol = 4, nrow = 3,
      plot.width = 12, plot.height = 6
    )
    return(seu)
  }) %>% setNames(.,capture)
  
  # save objects
  lapply(seq_along(capture), function(i) {
    message("[INFO] save demux'd objects...")
    so[[capture[i]]]@project.name <- capture[i]
    qs_save(so[[capture[i]]], here(save.loc, "3_HTODemux","int_obj", paste0(so[[capture[i]]]@project.name, ".demux.qs")))
  })
}

# --- filtering_singlets --------------------------------------------------------------------------
if (!steps$filtering_singlets) {
  message("[SKIP] Skipping singlet filtering step...")
} else {
  message("[INFO] Filtering scDblFinder and HTO empty droplets & doublets...")
  
  if (!exists("so") || is.null(so) || length(so) == 0) {
    so <- load_so_from_stage(capture, "3_HTODemux", save.loc, pattern = "demux")
  }
  
  so <- lapply(seq_along(capture), function(i) {
    so_filt1 <- subset(so[[capture[i]]], scDblFinder.class == "singlet")
    so_filt2 <- subset(so_filt1, HTO_classification.global == "Singlet")
    get_cell_counts(data = so_filt1, sample_name = capture[i], stage = "Retain only GEX singlets", save.loc) # cell counts
    get_cell_counts(data = so_filt2, sample_name = capture[i], stage = "Retain only HTO singlets", save.loc) # cell counts
    seu <- so_filt2
    
    plot_HTODemux(
      seu, 
      sample_name = paste0(so[[capture[i]]]@project.name, ".onlySinglets"), 
      save.loc = here(save.loc, "3_HTODemux"),
      HTO_AssayName = "HTO", RNA_AssayName = "RNA", 
      hash_class_column = "HTO_classification.global", 
      hash.ident = "HTO_maxID", ncol = 4, nrow = 3,
      plot.width = 12, plot.height = 6
    )
    return(seu)
  }) %>% setNames(.,capture)
  
  # --- special exception for TIPTOE011 -----------------------------------------------------------
  if ("TIPTOE011" %in% capture) {
    message("TIPTOE011 resequenced HTO re-analysis...")
    seu <- so[["TIPTOE011"]]
    remove_cells <- colnames(seu)[which(
      seu$HTO_maxID %in% c("PBMC", "Primary-Tumour", "LN", "Colon-Normal-Tissue") &
        seu$HTO_secondID %in% c("PBMC", "Primary-Tumour", "LN", "Colon-Normal-Tissue")
    )]
    seu <- subset(seu, cells = remove_cells, invert = TRUE)
    seu$HTO_maxID[which(seu$HTO_maxID == "LN")] <- "idLN"
    so[["TIPTOE011"]] <- seu
    plot_HTODemux(
      seu, sample_name = "TIPTOE011_resequenced.onlySinglets", 
      save.loc = here(save.loc, "3_HTODemux"),
      HTO_AssayName = "HTO", RNA_AssayName = "RNA", 
      hash_class_column = "HTO_classification.global", 
      hash.ident = "HTO_maxID", ncol = 4, nrow = 3,
      plot.width = 12, plot.height = 6
    )
  }
  # fix metadata
  so <- lapply(so, fixMetadata, res = "RNA_snn_res.1")
  
  # save objects
  lapply(seq_along(capture), function(i) {
    message("[INFO] save filtered objects...")
    qs_save(so[[capture[i]]], here(save.loc, "3_HTODemux","int_obj", paste0(so[[capture[i]]]@project.name, ".filtered.qs")))
  })
}

# --- cluster_qc initial clustering --------------------------------------------------------------------------
if (!steps$cluster_qc) {
  message("[SKIP] Skipping clustering QC step...")
} else {
  
  if (!exists("so") || is.null(so) || length(so) == 0) {
    so <- load_so_from_stage(capture, "3_HTODemux", save.loc, pattern = "filtered")
  }
  
  # initial clustering
  so <- lapply(seq_along(capture), function(i) {
    message("[INFO] ", so[[capture[i]]]@project.name, " clustering analysis...")
    seu <- so[[capture[i]]]
    bk.list <- get_blacklist(seu, species = species) %>% unlist
    seu <- find_clusters_to_remove(
      seu, capture = capture[i], assay = "RNA",
      reduction = "pca", reduction.name = "umap", npcs = min.pc,
      seed.use = seed.use, prefix = "0", save.loc = here(save.loc, "5_cluster_qc_1"),
      resolutions = seq(0.5, 1, 0.1), features = NULL,
      vars.to.regress = vars.to.regress, bk.list = bk.list
    )
    return(seu)
    
    get_cell_counts(data = seu, sample_name = capture[i], stage = "Clustering analysis", save.loc) # cell counts
    
  }) %>% setNames(., capture)
  
  # add save.loc to make sure fixMetadata works
  # data_capture <- save.loc %>% str_split(pattern="/") %>% unlist %>% tail(1)
  so <- lapply(so, fixMetadata, res = "RNA_snn_res.1", fix_sampleClusters = TRUE)
  
  # save objects
  lapply(seq_along(capture), function(i) {
    message("[INFO] save clustered objects...")
    qs_save(so[[capture[i]]], here(save.loc, "5_cluster_qc_1", "int_obj", paste0(so[[capture[i]]]@project.name, ".clustered.qs")))
  })
}

return(so)
}
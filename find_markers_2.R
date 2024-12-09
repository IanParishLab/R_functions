find_markers <- function(SeuratObjectIN,
                         DefaultSeuratAssay,
                         SeuratIdents,
                         ident.1 = NULL,
                         ident.2 = NULL,
                         max.cells.per.ident = Inf, 
                         logfc.threshold = 0.125, 
                         min.pct = 0.05,
                         test.use = "LR",
                         latent.vars = paste0("nCount_",DefaultSeuratAssay),
                         save.loc="./full.markers/",
                         run_FindAllMarkers = TRUE,
                         run_FindMarkers = TRUE,
                         annotate_peaks = NULL,
                         annotations = NULL) {
  require(Seurat)
  require(Signac)
  # require(SeuratDisk)
  require(tidyverse)
  
  ################ tmp ################ 
  # SeuratObjectIN="/researchers/nicole.saw/projects/Sinead_Reading/R/TCF1plus.rds"
  # DefaultSeuratAssay = "RNA"
  # latent.vars = "nCount_RNA"
  # SeuratIdents = "RNA_snn_res.0.9"
  # save.loc="/researchers/nicole.saw/projects/Sinead_Reading/R/full.markers/"
  # run_FindAllMarkers = TRUE
  # run_FindMarkers = TRUE
  # test.use = "LR"
  # logfc.threshold = 0.125
  # min.pct = 0.05
  # annotate_peaks = TRUE
  # annotations = readRDS("/researchers/nicole.saw/projects/Sinead_Reading/R/EnsDb.Mmusculus.v79_annotations.rds")
  ################ tmp ################ 
  
  # setup #
  object.name <- deframe(strsplit(SeuratObjectIN,"/"))[length(deframe(strsplit(SeuratObjectIN,"/")))]
  object.name <- gsub(".rds","",object.name)
  saveRDSOutputName = paste0(c(object.name,DefaultSeuratAssay,test.use),collapse=".")
  
  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  # end setup #
  
  # check if saveDir exists
  print("[MSG] Checking save.loc, the output file path...")
  ifelse(!dir.exists(file.path(save.loc)), dir.create(file.path(save.loc)), FALSE)
  
  # read seurat object
  print(paste0("[MSG] reading in seurat object from ", SeuratObjectIN," ..."))
  so <- readRDS(SeuratObjectIN)
  
  # set seurat object default assay and identity to use
  print(paste0("[MSG] set default assay as ", DefaultSeuratAssay, " and Idents set as ", SeuratIdents))
  DefaultAssay(so) <- DefaultSeuratAssay
  
  # if SCT assay is used run this
  if(DefaultAssay(so) == "SCT"){
    so <- PrepSCTFindMarkers(so, assay = "SCT", verbose = FALSE)
  }
  
  # setting Idents and print clusters analysed
  Idents(so) <- SeuratIdents
  
  clusters = unique(so[[SeuratIdents]]) %>% deframe
  print((paste0("clusters analysed:", paste0(gtools::mixedsort(unlist(clusters)), collapse = ","))))
  
  # running FindAllMarkers
  if (run_FindAllMarkers == TRUE){  
    print(paste0("[MSG] run FindAllMarkers... with logfc.threshold = ", logfc.threshold, ", min.pct = ", min.pct, ", test.use = ", test.use, 
             ", latent.vars = ", paste0(latent.vars, collapse = ",")))

    FindAllMarkers.so <- FindAllMarkers(so, assay=DefaultSeuratAssay, 
                                        logfc.threshold = logfc.threshold, min.pct = min.pct, 
                                        test.use = test.use, latent.vars = latent.vars,
                                        max.cells.per.ident = max.cells.per.ident)
    # annotate peaks using ClosestFeature()
    if(!is.null(annotate_peaks))  {
      require(EnsDb.Mmusculus.v79)
      require(BSgenome.Mmusculus.UCSC.mm10)
      
      print(paste0("[MSG] Annotating peaks..."))
      annotations <- annotations
      colnames(FindAllMarkers.so)[which(colnames(FindAllMarkers.so) == "gene")] <- "peaks"
      dacr_ann <- ClosestFeature(so[[DefaultSeuratAssay]], annotation = annotations, regions = FindAllMarkers.so$peaks)
      colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
      da_cr <- full_join(FindAllMarkers.so, dacr_ann, by = "peaks") %>% distinct()
      
      print(paste0("[MSG] saving FindAllMarkers results in ","FindAllMarkers.",saveRDSOutputName,".annotated_dacr.rds"))
      saveRDS(da_cr, file.path(save.loc, paste0("FindAllMarkers.",saveRDSOutputName,".annotated_dacr.rds")))
    } else {
      
      print(paste0("[MSG] saving FindAllMarkers results in ",saveRDSOutputName,"..."))
      saveRDS(FindAllMarkers.so, file.path(save.loc, paste0("FindAllMarkers.",saveRDSOutputName,".rds")))
    }
  } 
  
  # running FindMarkers
  if (run_FindMarkers == TRUE){     
    # if all clusters
    # run findmarkers with logfc.threshold = 0, min.pct = 0
    print(paste0("[MSG] run FindMarkers... with logfc.threshold = 0, and min.pct = 0, and test.use = ", test.use,
                 ", and latent.vars = ", paste0(latent.vars, collapse = ",")))

    FindMarkers.so <- list()
    for (cluster in clusters){
      FindMarkers.so[[cluster]] <- FindMarkers(so, ident.1 = cluster,
                                               assay = DefaultSeuratAssay, 
                                               logfc.threshold = 0, 
                                               min.pct = 0, 
                                               test.use= test.use, 
                                               latent.vars = latent.vars
                                               )
      
      FindMarkers.so[[cluster]]$cluster <- cluster
      FindMarkers.so[[cluster]]$gene <- rownames(FindMarkers.so[[cluster]])
    }
    FindMarkers.so.rbind <- plyr::rbind.fill(FindMarkers.so)
    
    if(!is.null(annotate_peaks)) {
      print(paste0("[MSG] Annotating peaks..."))
      annotations <- annotations
      colnames(FindMarkers.so.rbind)[which(colnames(FindMarkers.so.rbind) == "gene")] <- "peaks"
      dacr_ann <- ClosestFeature(so[[DefaultSeuratAssay]], annotation = annotations, regions = FindMarkers.so.rbind$peaks)
      colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
      da_cr <- full_join(FindMarkers.so.rbind, dacr_ann, by = "peaks") %>% distinct()
      
      print(paste0("[MSG] saving FindMarkers. results in ","FindMarkers.",saveRDSOutputName,".annotated_dacr.rds"))
      saveRDS(da_cr, file.path(save.loc, paste0("FindMarkers.",saveRDSOutputName,".annotated_dacr.rds")))
    } else {
      
      print(paste0("[MSG] saving FindMarkers. results in ","FindMarkers.",saveRDSOutputName,".rds"))
      saveRDS(FindMarkers.so.rbind, file.path(save.loc, paste0("FindMarkers.",saveRDSOutputName,".rds")))
    }
  }
}
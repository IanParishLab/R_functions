find_markers <- function(seu_path,
                         assay,
                         group.by,
                         ident.1 = NULL,
                         ident.2 = NULL,
                         max.cells.per.ident = Inf, 
                         logfc.threshold = 0.125, 
                         min.pct = 0.05,
                         test.use = "LR",
                         latent.vars = paste0("nCount_",assay),
                         save.loc="find_markers",
                         run_FindAllMarkers = TRUE,
                         run_FindMarkers = TRUE,
                         annotate_peaks = NULL,
                         annotations = NULL,
                         ...) {
  suppressPackageStartupMessages({
    require(Seurat)
    require(tidyverse)
    # require(SeuratDisk)
  })
  
  ################ tmp ################ 
  # seu_path="/researchers/nicole.saw/projects/Sinead_Reading/R/TCF1plus.rds"
  # assay = "RNA"
  # latent.vars = "nCount_RNA"
  # group.by = "RNA_snn_res.0.9"
  # save.loc="/researchers/nicole.saw/projects/Sinead_Reading/R/full.markers/"
  # run_FindAllMarkers = TRUE
  # run_FindMarkers = TRUE
  # test.use = "LR"
  # logfc.threshold = 0.125
  # min.pct = 0.05
  # annotate_peaks = TRUE
  # annotations = readRDS("/researchers/nicole.saw/projects/Sinead_Reading/R/EnsDb.Mmusculus.v79_annotations.rds")
  ################ tmp ################ 
  
  # start setup #
  object_name <- deframe(strsplit(seu_path,"/"))[length(deframe(strsplit(seu_path,"/")))]
  object_name <- gsub(".rds","",object_name)
  output_name <- paste0(c(object_name,assay,test.use,group.by),collapse=".")
  
  # make output directory
  dir.create(save.loc, showWarnings = FALSE, recursive = TRUE)
  
  # end setup #
  
  # check if save.loc exists
  message("[MSG] Checking save.loc, the output file path...")
  ifelse(!dir.exists(file.path(save.loc)), dir.create(file.path(save.loc)), FALSE)
  
  # read Seurat object
  message("[MSG] reading in Seurat object from ", seu_path)
  so <- readRDS(seu_path)
  
  # set Seurat object default assay and identity to use
  message("[MSG] Default assay: ", assay, ", Idents: ", group.by)
  DefaultAssay(so) <- assay
  Idents(so) <- group.by
  # if SCT assay is used run this
  if(DefaultAssay(so) == "SCT"){
    so <- PrepSCTFindMarkers(so, assay = "SCT", verbose = FALSE)
  }
  
  # Print clusters analysed
  clusters = unique(so[[group.by]]) %>% deframe
  print((paste0("clusters analysed:", paste0(gtools::mixedsort(unlist(clusters)), collapse = ","))))
  
  # running FindAllMarkers
  if (run_FindAllMarkers == TRUE){  
    message("[MSG] run FindAllMarkers... with logfc.threshold = ", logfc.threshold, ", min.pct = ", min.pct, ", test.use = ", test.use, 
            ", latent.vars = ", paste0(latent.vars, collapse = ","))
    
    find_all_markers <- FindAllMarkers(so, assay=assay, 
                                       logfc.threshold = logfc.threshold, min.pct = min.pct, 
                                       test.use = test.use, latent.vars = latent.vars,
                                       max.cells.per.ident = max.cells.per.ident)
    # annotate peaks using ClosestFeature()
    if(!is.null(annotate_peaks))  {
      require(Signac)
      require(EnsDb.Mmusculus.v79)
      require(BSgenome.Mmusculus.UCSC.mm10)
      
      message("[MSG] Annotating peaks...")
      annotations <- annotations
      colnames(find_all_markers)[which(colnames(find_all_markers) == "gene")] <- "peaks"
      dacr_ann <- ClosestFeature(so[[assay]], annotation = annotations, regions = find_all_markers$peaks)
      colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
      da_cr <- full_join(find_all_markers, dacr_ann, by = "peaks") %>% distinct()
      
      message("[MSG] saving FindAllMarkers results in ","FindAllMarkers.",output_name,".annotated_dacr.rds")
      saveRDS(da_cr, file.path(save.loc, paste0("FindAllMarkers.",output_name,".annotated_dacr.rds")))
    } else {
      
      message("[MSG] saving FindAllMarkers results in ",output_name,"...")
      saveRDS(find_all_markers, file.path(save.loc, paste0("FindAllMarkers.",output_name,".rds")))
    }
  } 
  
  # running FindMarkers
  if (run_FindMarkers == TRUE){     
    # if all clusters
    # run findmarkers with logfc.threshold = 0, min.pct = 0
    message("[MSG] run FindMarkers... with logfc.threshold = 0, and min.pct = 0, and test.use = ", test.use,
            ", and latent.vars = ", paste0(latent.vars, collapse = ","))
    
    find_markers <- list()
    for (cluster in clusters){
      find_markers[[cluster]] <- FindMarkers(so, 
                                             ident.1 = cluster,
                                             assay = assay, 
                                             logfc.threshold = 0, 
                                             min.pct = 0, 
                                             test.use= test.use, 
                                             latent.vars = latent.vars,
                                             ...
      )
      
      find_markers[[cluster]]$cluster <- cluster
      find_markers[[cluster]]$gene <- rownames(find_markers[[cluster]])
    }
    find_markers.rbind <- plyr::rbind.fill(find_markers)
    
    if(!is.null(annotate_peaks)) {
      message("[MSG] Annotating peaks...")
      annotations <- annotations
      colnames(find_markers.rbind)[which(colnames(find_markers.rbind) == "gene")] <- "peaks"
      dacr_ann <- ClosestFeature(so[[assay]], annotation = annotations, regions = find_markers.rbind$peaks)
      colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
      da_cr <- full_join(find_markers.rbind, dacr_ann, by = "peaks") %>% distinct()
      
      message("[MSG] saving FindMarkers results in ","FindMarkers.",output_name,".annotated_dacr.rds")
      saveRDS(da_cr, file.path(save.loc, paste0("FindMarkers.",output_name,".annotated_dacr.rds")))
    } else {
      
      message("[MSG] saving FindMarkers results in ","FindMarkers.",output_name,".rds")
      saveRDS(find_markers.rbind, file.path(save.loc, paste0("FindMarkers.",output_name,".rds")))
    }
  }
}
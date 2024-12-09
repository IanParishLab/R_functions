find_markers_specific_clusters <- function(SeuratObjectIN,
                                           DefaultSeuratAssay,
                                           SeuratIdents,
                                           idents.1 = NULL,
                                           idents.2 = NULL,
                                           max.cells.per.ident = NULL, 
                                           test.use = "LR",
                                           latent.vars = paste0("nCount_",DefaultSeuratAssay),
                                           saveRDSPath="./full.markers/",
                                           annotate_peaks = NULL,
                                           annotations = NULL) {
  
  require(Seurat)
  require(Signac)
  require(SeuratDisk)
  require(tidyverse)
  
  ################ tmp ################ 
  # SeuratObjectIN="/researchers/nicole.saw/projects/Sinead_Reading/R/ikzf_filt1-customPeakList-chromVAR.rds"
  # DefaultSeuratAssay = "ATAC"
  # latent.vars = "nCount_ATAC"
  # SeuratIdents = "RNA_snn_res.0.9"
  # idents.1 = "1"
  # idents.2 = "4"
  # test.use = "LR"
  # saveRDSPath="/researchers/nicole.saw/projects/Sinead_Reading/R/full.markers/"
  # annotate_peaks = TRUE
  # annotations = readRDS("/researchers/nicole.saw/projects/Sinead_Reading/R/EnsDb.Mmusculus.v79_annotations.rds")
  ################ tmp ################ 
  
  # setup #
  object.name <- deframe(strsplit(SeuratObjectIN,"/"))[length(deframe(strsplit(SeuratObjectIN,"/")))]
  object.name <- gsub(".rds","",object.name)
  saveRDSOutputName = paste0(c(object.name,DefaultSeuratAssay,test.use),collapse=".")
  # end setup #
  
  # check if saveDir exists
  print("[MSG] Checking saveRDSPath, the output file path...")
  ifelse(!dir.exists(file.path(saveRDSPath)), dir.create(file.path(saveRDSPath)), FALSE)
  
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
  
  # if specific clusters only
  print(paste0("[MSG] run FindMarkers with clusters ",idents.1," and ",idents.2,", with logfc.threshold = 0, and min.pct = 0, and test.use = ", test.use))
  
  FindMarkers.so <- FindMarkers(so, ident.1 = idents.1, ident.2 = idents.2,
                                assay = DefaultSeuratAssay, 
                                logfc.threshold = 0, min.pct = 0, 
                                test.use = test.use, latent.vars = latent.vars)

  FindMarkers.so$cluster_1 <- ifelse(length(idents.1) > 1, paste0(idents.1, collapse = "_"),idents.1)
  FindMarkers.so$cluster_2 <- ifelse(length(idents.2) > 1, paste0(idents.2, collapse = "_"),idents.2)

  FindMarkers.so$gene <- rownames(FindMarkers.so)
  
  if(!is.null(annotate_peaks))  {
    print(paste0("[MSG] Annotating peaks..."))
    require(EnsDb.Mmusculus.v79)
    require(BSgenome.Mmusculus.UCSC.mm10)
    
    annotations <- annotations
    colnames(FindMarkers.so)[which(colnames(FindMarkers.so) == "gene")] <- "peaks"
    dacr_ann <- ClosestFeature(so[[DefaultSeuratAssay]], annotation = annotations, regions = FindMarkers.so$peaks)
    colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
    da_cr <- full_join(FindMarkers.so, dacr_ann, by = "peaks") %>% distinct()
    
    print(paste0("[MSG] saving FindMarkers. results in ","FindMarkers.",saveRDSOutputName,".annotated_dacr.rds"))
    saveRDS(da_cr, paste0(saveRDSPath,"FindMarkers.",saveRDSOutputName,".",idents.1,"_",paste0(idents.2, collapse = "_"),".annotated_dacr.rds"))
  } else {
    
    print(paste0("[MSG] saving FindMarkers results in ",saveRDSOutputName,"..."))
    saveRDS(FindMarkers.so ,paste0(saveRDSPath,"FindMarkers.",saveRDSOutputName,".",idents.1,"_",paste0(idents.2, collapse = "_"),".rds"))
    
  }
}
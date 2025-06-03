find_markers_specific_clusters <- function(seu_path,
                                           assay,
                                           group.by,
                                           idents.1 = NULL,
                                           idents.2 = NULL,
                                           max.cells.per.ident = NULL, 
                                           test.use = "LR",
                                           latent.vars = paste0("nCount_",assay),
                                           save.loc="full.markers/",
                                           annotate_peaks = NULL,
                                           annotations = NULL) {
  
  require(Seurat)
  require(tidyverse)
  # require(SeuratDisk)
  
  ################ tmp ################ 
  # seu_path="/researchers/nicole.saw/projects/Sinead_Reading/R/ikzf_filt1-customPeakList-chromVAR.rds"
  # assay = "ATAC"
  # latent.vars = "nCount_ATAC"
  # group.by = "RNA_snn_res.0.9"
  # idents.1 = "1"
  # idents.2 = "4"
  # test.use = "LR"
  # save.loc="/researchers/nicole.saw/projects/Sinead_Reading/R/full.markers/"
  # annotate_peaks = TRUE
  # annotations = readRDS("/researchers/nicole.saw/projects/Sinead_Reading/R/EnsDb.Mmusculus.v79_annotations.rds")
  ################ tmp ################c 
  
  # setup #
  object_name <- deframe(strsplit(seu_path,"/"))[length(deframe(strsplit(seu_path,"/")))]
  object_name <- gsub(".rds","",object_name)
  output_name <- paste0(c(object_name,assay,test.use),collapse=".")
  # end setup #
  
  # check if saveDir exists
  message("[MSG] Checking save.loc, the output file path...")
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  # read seurat object
  message("[MSG] reading in seurat object from ", seu_path)
  so <- readRDS(seu_path)
  
  # set seurat object default assay and identity to use
  message("[MSG] Default assay: ", assay, ", Idents: ", group.by)
  DefaultAssay(so) <- assay
  
  # if SCT assay is used run this
  if(DefaultAssay(so) == "SCT"){
    so <- PrepSCTFindMarkers(so, assay = "SCT", verbose = FALSE)
  }
  
  # setting Idents and print clusters analysed
  Idents(so) <- group.by
  
  # if specific clusters only
  message("[MSG] run FindMarkers with clusters ",idents.1," and ",idents.2,", with logfc.threshold = 0, and min.pct = 0, and test.use = ", test.use)
  
  find_markers <- FindMarkers(so, 
                              ident.1 = idents.1, ident.2 = idents.2,
                                assay = assay, 
                                logfc.threshold = 0, 
                                min.pct = 0, 
                                test.use = test.use, 
                                latent.vars = latent.vars)

  find_markers$cluster_1 <- ifelse(length(idents.1) > 1, paste0(idents.1, collapse = "_"),idents.1)
  find_markers$cluster_2 <- ifelse(length(idents.2) > 1, paste0(idents.2, collapse = "_"),idents.2)
  find_markers$gene <- rownames(find_markers)
  
  if(!is.null(annotate_peaks))  {
    message("[MSG] Annotating peaks...")
    require(Signac)
    require(EnsDb.Mmusculus.v79)
    require(BSgenome.Mmusculus.UCSC.mm10)
    
    annotations <- annotations
    colnames(find_markers)[which(colnames(find_markers) == "gene")] <- "peaks"
    dacr_ann <- ClosestFeature(so[[assay]], annotation = annotations, regions = find_markers$peaks)
    colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
    da_cr <- full_join(find_markers, dacr_ann, by = "peaks") %>% distinct()
    
    message("[MSG] saving FindMarkers. results in ","FindMarkers.",output_name,".annotated_dacr.rds")
    saveRDS(da_cr, paste0(save.loc,"FindMarkers.",output_name,".",idents.1,"_",paste0(idents.2, collapse = "_"),".annotated_dacr.rds"))
  } else {
    
    message("[MSG] saving FindMarkers results in ",output_name,"...")
    saveRDS(FindMarkers.so ,paste0(save.loc,"FindMarkers.",output_name,".",idents.1,"_",paste0(idents.2, collapse = "_"),".rds"))
    
  }
}
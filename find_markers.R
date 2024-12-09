find_markers <- function(SeuratObjectIN,
                         DefaultSeuratAssay,
                         SeuratIdents,
                         max.cells.per.ident = NULL, 
                         logfc.threshold = 0.125, 
                         min.pct = 0.05,
                         test.use = "LR",
                         latent.vars = "nCount_RNA",
                         saveRDSPath="./full.markers/",
                         run_FindAllMarkers = TRUE,
                         run_FindMarkers = TRUE) {
  require(Seurat)
  require(SeuratDisk)
  require(tidyverse)
  
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

  if(DefaultAssay(so) == "SCT"){
    so <- PrepSCTFindMarkers(so, assay = "SCT", verbose = FALSE)
  }

  Idents(so) <- SeuratIdents
  
  clusters = unique(so[[SeuratIdents]]) %>% deframe
  print(gtools::mixedsort(paste0("clusters analysed:", sort(unlist(clusters)))))
  
  if (run_FindAllMarkers == TRUE){  
    # run findallmarkers
    print(paste0("[MSG] run FindAllMarkers... with logfc.threshold = ", logfc.threshold, ", and min.pct = ", min.pct, ", and test.use = ", test.use, ", and latent.vars = ", latent.vars))
    FindAllMarkers.so <- FindAllMarkers(so, assay=DefaultSeuratAssay, logfc.threshold = logfc.threshold, min.pct = min.pct, test.use = test.use, latent.vars = latent.vars)
    print(paste0("[MSG] saving FindAllMarkers results in ",saveRDSOutputName,"..."))
    saveRDS(FindAllMarkers.so,paste0(saveRDSPath,"FindAllMarkers.",saveRDSOutputName,".rds")
    )}
  
  if (run_FindMarkers == TRUE){  
    # run findmarkers with logfc.threshold = 0, min.pct = 0
    print(paste0("[MSG] run FindMarkers... with logfc.threshold = 0, and min.pct = 0, and test.use = ", test.use))
    FindMarkers.so <- list()
    for (cluster in clusters){
      FindMarkers.so[[cluster]] <- FindMarkers(so, ident.1 = cluster, assay = DefaultSeuratAssay, logfc.threshold = 0, min.pct = 0, test.use= test.use, latent.vars = latent.vars)
      FindMarkers.so[[cluster]]$cluster <- cluster
      FindMarkers.so[[cluster]]$gene <- rownames(FindMarkers.so[[cluster]])
      # saveRDS(FindMarkers.so[[cluster]] ,paste0(saveRDSPath,"FindMarkers.",cluster,saveRDSOutputName,".rds"))
    }
    FindMarkers.so.rbind <- plyr::rbind.fill(FindMarkers.so)
    print(paste0("[MSG] saving FindMarkers results in ",saveRDSOutputName,"..."))
    saveRDS(FindMarkers.so.rbind ,paste0(saveRDSPath,"FindMarkers.",saveRDSOutputName,".rds"))
  }
}
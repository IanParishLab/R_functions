find_clusters_to_remove_multiome <- function(so, capture, assays, res = NULL, remove_clusters = NULL, bk.list, process_wsnn = TRUE, npcs, prefix, seed.use){
  require(Seurat)
  require(tidyverse)
  require(clustree)
  require(cowplot)
  require(Signac)
  
  if("RNA" %in% assays){
    assay = "RNA"
    print("processing RNA assay...")
    DefaultAssay(so) <- assay
    
    if(!is.null(remove_clusters)){
      remove.cells <- colnames(so)[so@meta.data[[res]] %in% remove_clusters]
      so <- subset(so, cells = remove.cells, invert = TRUE)
    }
    
    if(!("integrated" %in% assay)){
      so <- NormalizeData(object = so)
      so <- FindVariableFeatures(object = so)
      VariableFeatures(so) <- unique(setdiff(VariableFeatures(so), bk.list))
    } else {
      print(paste0("also processing `integrated` assay..."))
    }
    
    so <- ScaleData(object = so, vars.to.regress = c("percent.mito","CC.Difference"), features = rownames(so))
    so <- RunPCA(so, npcs = npcs, verbose = FALSE)
    so <- RunUMAP(so, dims = 1:npcs, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", return.model = TRUE, verbose = FALSE, seed.use = seed.use)
    print(paste0("seed.use = ",seed.use))
    
    resolutions<-seq(0.5,2,by=0.1)
    for (j in 1:length(resolutions)){
      so <- FindNeighbors(so, dims = 1:npcs, assay = assay, reduction = "pca", graph.name = paste0(assay,"_snn"), verbose = FALSE)
      so <- FindClusters(object = so, graph.name = paste0(assay,"_snn"), resolution = resolutions[j], verbose = FALSE)
    }
    
    pdf(paste0("qc_VlnPlot/", prefix ,"_", capture,"_",assay,"_snn_res_VlnPlot_res_0.5_2.pdf"), height = 9, width = 12)
    p<-list()
    for (j in 1:length(resolutions)) {
      so <- SetIdent(so, value = paste0(assay,"_snn_res.", resolutions[j]))
      p[[j]] <- VlnPlot(so, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0, ncol = 1)
      title <- ggdraw() + draw_label(paste0("resolution=", resolutions[j]), fontface = 'bold')
      print(cowplot::plot_grid(title, p[[j]], ncol = 1, rel_heights = c(0.1, 1)))
    }
    dev.off()
    
    pdf(paste0("clustree/", prefix ,"_",capture,"_",assay,"_clusterResTree.pdf"), height = 12, width = 10)
    print(clustree(so, prefix = paste0(assay,"_snn_res.")))
    dev.off()
  }
  
  # process ATAC assay
  if("ATAC" %in% assays){
    assay = "ATAC"
    print("processing ATAC assay...")
    DefaultAssay(so) <- assay
    if(!is.null(remove_clusters)){
      remove.cells <- colnames(so)[so@meta.data[[res]] %in% remove_clusters]
      so <- subset(so, cells = remove.cells, invert = TRUE)
    }
    
    so <- RunTFIDF(so)
    so <- FindTopFeatures(so, min.cutoff = 'q0')
    so <- RunSVD(so)
    so <- RunUMAP(so, reduction = 'lsi', dims = 2:npcs, reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE, seed.use = seed.use)
    print(paste0("seed.use = ",seed.use))
    
    for (j in 1:length(resolutions)){
      so <- FindNeighbors(so, dims = 2:npcs, assay = assay, reduction = "lsi", graph.name = paste0(assay,"_snn"), verbose = FALSE)
      so <- FindClusters(object = so, graph.name = paste0(assay,"_snn"), resolution = resolutions[j], algorithm = 3, verbose = FALSE)
    }
    
    resolutions<-seq(0.5,2,by=0.1)
    pdf(paste0("qc_VlnPlot/", prefix ,"_", capture,"_",assay,"_snn_res_VlnPlot_res_0.5_2.pdf"), height = 9, width = 12)
    p<-list()
    for (j in 1:length(resolutions)) {
      so <- SetIdent(so, value = paste0(assay,"_snn_res.", resolutions[j]))
      p[[j]] <- VlnPlot(so, features = c("nCount_ATAC", "nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0, ncol = 1)
      title <- ggdraw() + draw_label(paste0("resolution=", resolutions[j]), fontface = 'bold')
      print(cowplot::plot_grid(title, p[[j]], ncol = 1, rel_heights = c(0.1, 1)))
    }
    dev.off()
    
    pdf(paste0("clustree/", prefix ,"_",capture,"_",assay,"_clusterResTree.pdf"), height = 12, width = 10)
    print(clustree(so, prefix = paste0(assay,"_snn_res.")))
    dev.off()
  }
  
  if(process_wsnn == TRUE){
    print("processing MultiModal information and generating weighted nearest neighours (WNN) clustree...")
    print("reduction.name = wnn.umap...")
    
    resolutions<-seq(0.5,2,by=0.1)
    for (j in 1:length(resolutions)){
      so <- FindMultiModalNeighbors(so, reduction.list = list("pca", "lsi"), dims.list = list(1:npcs, 2:npcs))
      so <- RunUMAP(so, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
      so <- FindClusters(so, graph.name = "wsnn", algorithm = 3, resolution = resolutions[j], verbose = FALSE)
    }
    
    pdf(paste0("clustree/", prefix ,"_",capture,"_wsnn_clusterResTree.pdf"), height = 12, width = 10)
    print(clustree(so, prefix = "wsnn_res."))
    dev.off()
  }
  
  return(so)
}

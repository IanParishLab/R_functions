##todo 
# add features argument to customize features to include/exclude
# include scgate blacklist
find_clusters_to_remove <- function(so, capture, assay, reduction, reduction.name, res=NULL, vars.to.regress,
remove_clusters=NULL, npcs, prefix, seed.use, save.loc = "./cluster_qc", resolutions = seq(0.5,2,by = 0.1),
features, bk.list){
  
  require(Seurat)
  require(tidyverse)
  require(clustree)
  require(cowplot)
  require(scGate)
  
  # make output directory
  lapply(c("int_obj","plots"), function(d){
    ifelse(!dir.exists(file.path(save.loc,d)),
           dir.create(file.path(save.loc,d), recursive = TRUE), paste0(save.loc," directory exists"))
  })
  
  DefaultAssay(so) <- assay
  print(paste0("[MSG] assay = ", assay))
  
  if(!is.null(remove_clusters)){
    remove.cells <- colnames(so)[so@meta.data[[res]] %in% remove_clusters]
    so <- subset(so, cells = remove.cells, invert = TRUE)
  }
  
  # normalisation step
  if(assay == "RNA"){
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so)
    VariableFeatures(so) <- unique(setdiff(VariableFeatures(so), bk.list))
      if(!is.null(features)){
        VariableFeatures(so) <- unique(c(VariableFeatures(so), features))
      }
    so <- ScaleData(object = so, vars.to.regress = vars.to.regress)
  }
  
  if(assay == "integrated"){
    so <- ScaleData(object = so, vars.to.regress = vars.to.regress)
  }
  
  if(assay == "SCT"){
    SCTransform(so, vst.flavor = "v2", vars.to.regress = vars.to.regress, verbose = FALSE, seed.use = seed.use)
  }
  
  so <- RunPCA(so, npcs = 50, verbose = FALSE)
  so <- RunUMAP(so, dims = 1:npcs, reduction = reduction, reduction.name = reduction.name, verbose = FALSE, seed.use = seed.use)
  print(paste0("[MSG] reduction = ",reduction))
  print(paste0("[MSG] reduction.name = ",reduction.name))
  print(paste0("[MSG] seed.use = ",seed.use))
  
  for (j in 1:length(resolutions)){
    so <- FindNeighbors(so, graph.name = paste0(assay,"_snn"), dims = 1:npcs, verbose = FALSE)
    so <- FindClusters(object = so, graph.name = paste0(assay,"_snn"), resolution = resolutions[j])
  }
  
  pdf(file.path(save.loc, "plots", paste0(prefix ,"_", capture,"_",assay,"_snn_res_VlnPlot.pdf")), height = 9, width = 12)
  p<-list()
  for (j in 1:length(resolutions)) {
    so <- SetIdent(so, value = paste0(assay,"_snn_res.", resolutions[j]))
    p[[j]] <- VlnPlot(so, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0, ncol = 1)
    title <- ggdraw() + draw_label(paste0("resolution=", resolutions[j]), fontface = 'bold')
    print(cowplot::plot_grid(title, p[[j]], ncol = 1, rel_heights = c(0.1, 1)))
  }
  dev.off()
  
  library(clustree)
  pdf(file.path(save.loc, "plots", paste0(prefix ,"_",capture,"_clusterResTree.pdf")), height = 12, width = 10)
  print(clustree(so, prefix = paste0(assay,"_snn_res.")))
  dev.off()
  
  return(so)
}
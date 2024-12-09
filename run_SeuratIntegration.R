run_SeuratIntegration <- function(seu.obj.list, merged.seu.obj.name, min.pc, seed.use, 
                                  save.loc = "integration/", integration_method, integration_method_name,
                                  vars.to.regress = c("percent.mito","CC.Difference","RP.Score1"),
                                  conda_env_path = "/home/nsaw/miniconda3/envs/scvi-env/", 
                                  hash.ident = "HTO_maxID", bk.list,
                                  cols = list(HTO_col,capture_col)
){
  suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(reticulate)
  library(tidyverse)
  })
  
  # # test parameters
  # seu.obj.list = so[-2]
  # merged.seu.obj.name = "test"
  # min.pc = min.pc
  # seed.use = seed.use
  # save.loc = here("CD8", "8_integration","not_clean")
  # integration_method = scVIIntegration
  # integration_method_name = "scvi"
  # vars.to.regress = c("percent.mito","CC.Difference","RP.Score1")
  # conda_env_path = "/home/nsaw/miniconda3/envs/scvi-env/"
  # hash.ident = "HTO_maxID"
  # cols = c(HTO_col,capture_col)
  # bk.list = scGate::genes.blacklist.default$Hs$Ribo
  # # end test parameters

  # make output directory
  ifelse(!dir.exists(file.path(save.loc)),
         dir.create(file.path(save.loc), recursive = TRUE), paste0(save.loc," directory exists"))
  
  reduction <-  paste0("integrated.", integration_method_name)
  reduction.name <-  paste0("umap.", integration_method_name)
  
  # merge scRNAseq into 1 big dataset
  print("[MSG] Merging data & add.cell.ids with `seu.obj@project.name` ...")
  add.cell.ids = sapply(seu.obj.list, function(seu){ seu@project.name }) %>% unlist()
  merged.seu.obj <- merge(seu.obj.list[[1]], y = seu.obj.list[2:length(seu.obj.list)], add.cell.ids = add.cell.ids)
  
  print("[MSG] normalizeAndScaleData ...")
  merged.seu.obj <- normalizeAndScaleData(merged.seu.obj, 
                                          dry_run = FALSE,
                                          sampleName = merged.seu.obj.name, 
                                          min.pc = min.pc, 
                                          seed.use = seed.use, 
                                          save.loc = save.loc,
                                          vars.to.regress = vars.to.regress,
                                          bk.list = bk.list)
  
  merged.seu.obj <- RunUMAP(merged.seu.obj, dims = 1:min.pc, 
                            reduction = "pca", reduction.name = "umap.unintegrated")
  
  before_int_umap <- DimPlot(merged.seu.obj, reduction = "umap.unintegrated", 
                             group.by = c(hash.ident,"orig.ident"), 
                             cols = cols) & 
    theme(aspect.ratio = 1) & 
    NoAxes()
  
  print("[MSG] Saving unintegrated object ...")
  saveRDS(merged.seu.obj, file.path(save.loc, paste0(merged.seu.obj.name,".unintegrated_so.rds")))
  
  print("[MSG] IntegrateLayers ...")
  combined <- IntegrateLayers(object = merged.seu.obj,
                              method = integration_method,
                              orig.reduction = "pca",
                              new.reduction = reduction,
                              conda_env = conda_env_path,
                              verbose = TRUE)
  
  print("[MSG] IntegrateLayers ...")
  combined <- JoinLayers(combined)
  
  combined <- RunUMAP(combined, dims = 1:min.pc, seed.use = seed.use,
                      reduction = reduction, reduction.name = reduction.name, 
                      verbose =  FALSE)
  
  after_int_umap <- DimPlot(combined, 
                            reduction = reduction.name, 
                            group.by = c(hash.ident,"orig.ident"), 
                            cols = cols) &
    theme(aspect.ratio = 1) &
    NoAxes()
  
  print("[MSG] Saving integrated object ...")
  saveRDS(combined, file.path(save.loc, paste0(merged.seu.obj.name,".rds")))
  
  # save plot
  ggsave(before_int_umap/after_int_umap,
         filename = file.path(save.loc, paste0("before_after_integratelayers_umap.",merged.seu.obj.name,".pdf")),
         width = 10, height = 10, device = 'pdf')
  
}

run_ProjecTILs_make_reference <- function(seu, annot.by = "ClusterNames", filter_expr = NULL, name, save.loc, 
vars.to.regress = "percent.mito", 
                                          bk.list = unlist(scGate::genes.blacklist.default$Mm[c("Mito")]), 
                                          ndim = 30, seed = 100) {
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
    library(ProjecTILs)
  })

  # make output directory
  dir.create(file.path(save.loc), showWarnings = FALSE, recursive = TRUE)

  message("[MSG] Preparing subset: ", name)
  if (!is.null(filter_expr)) {
    seu <- subset(seu, subset = !!rlang::parse_expr(filter_expr))
  }
  
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  
  # Remove bk.list genes
  VariableFeatures(seu) <- setdiff(VariableFeatures(seu), bk.list)
  
  # Basic preprocessing
  seu <- ScaleData(seu, verbose = FALSE, vars.to.regress = vars.to.regress)
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = ndim)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:ndim, seed.use = seed)
  seu <- seu[-which(rowSums(seu) == 0), ]
  
  # Make reference
  ref <- make.reference(ref = seu, ndim = ndim, seed = seed, recalculate.umap = TRUE, annotation.column = annot.by)
  
  # Recalculate without bk.listed features in reference
  VariableFeatures(ref) <- setdiff(VariableFeatures(ref), bk.list)
  seu <- ScaleData(seu, verbose = FALSE, vars.to.regress = vars.to.regress)
  ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = ndim)
  ref <- RunUMAP(ref, reduction = "pca", dims = 1:ndim, seed.use = seed)
  
  # Save reference
  saveRDS(ref, file.path(save.loc, paste0("ref.", name, ".rds")))
  pdf(file.path(save.loc, paste0("ref.", name, ".umap.pdf")), width = 5, height = 5)
  print(DimPlot(ref, label = TRUE, repel = TRUE, label.size = 4, group.by = annot.by))
  dev.off()
  
  message("[MSG] Finished reference: ", name)
  return(invisible(ref))
}

# Example usage:
# joGiles_RNA <- readRDS("path/to/scRNAseq_longitudinal_LCMV_Arm_Cl13.Rds")
# save.loc <- "/researchers/nicole.saw/references/ProjecTILs/"
# 
# # Cl13_d30 only
# run_ProjecTILs_make_reference(joGiles_RNA, filter_expr = 'sampleID == "Cl13_d30"', name = "giles_Cl13_d30", save.loc = save.loc)
# 
# # All cells except Naive
# run_ProjecTILs_make_reference(joGiles_RNA, filter_expr = 'ClusterNames != "Naive"', name = "giles_no_NaiveCells", save.loc = save.loc)


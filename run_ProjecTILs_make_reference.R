run_ProjecTILs_make_reference <- function(seu, annot.by = "ClusterNames", filter_expr = NULL, name, outdir, 
                                          blacklist = unlist(scGate::genes.blacklist.default$Mm[c("Mito", "Ribo")]), 
                                          ndim = 30, seed = 100) {
  message("[MSG] Preparing subset: ", name)
  
  if (!is.null(filter_expr)) {
    seu <- subset(seu, subset = !!rlang::parse_expr(filter_expr))
  }
  
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, verbose = FALSE)
  
  # Remove blacklist genes
  VariableFeatures(seu) <- setdiff(VariableFeatures(seu), blacklist)
  
  # Basic preprocessing
  seu <- ScaleData(seu, verbose = TRUE)
  seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = ndim)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:ndim, seed.use = seed)
  seu <- seu[-which(rowSums(seu) == 0), ]
  
  # Make reference
  ref <- make.reference(ref = seu, ndim = ndim, seed = seed, recalculate.umap = TRUE, annotation.column = annot.by)
  
  # Recalculate without blacklisted features
  VariableFeatures(ref) <- setdiff(VariableFeatures(ref), blacklist)
  ref <- ScaleData(ref, verbose = TRUE)
  ref <- RunPCA(ref, features = VariableFeatures(ref), npcs = ndim)
  ref <- RunUMAP(ref, reduction = "pca", dims = 1:ndim, seed.use = seed)
  
  # Save reference
  saveRDS(ref, file.path(outdir, paste0("ref.", name, ".rds")))
  pdf(file.path(outdir, paste0("ref.", name, ".umap.pdf")), width = 5, height = 5)
  print(DimPlot(ref, label = TRUE, repel = TRUE, label.size = 4, group.by = annot.by))
  dev.off()
  
  message("[MSG] Finished reference: ", name)
  return(invisible(ref))
}

# Example usage:
# joGiles_RNA <- readRDS("path/to/scRNAseq_longitudinal_LCMV_Arm_Cl13.Rds")
# outdir <- "/researchers/nicole.saw/references/ProjecTILs/"
# 
# # Cl13_d30 only
# run_ProjecTILs_make_reference(joGiles_RNA, filter_expr = 'sampleID == "Cl13_d30"', name = "giles_Cl13_d30", outdir = outdir)
# 
# # All cells except Naive
# run_ProjecTILs_make_reference(joGiles_RNA, filter_expr = 'ClusterNames != "Naive"', name = "giles_no_NaiveCells", outdir = outdir)


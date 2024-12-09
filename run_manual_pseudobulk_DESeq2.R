run_manual_pseudobulk_DESeq2 <- function(so, idents, idents.1, idents.2 = NULL, assay = "customPeakList",
                                         fitType = 'mean',
                                         annotations_filepath = "EnsDb.Mmusculus.v79_annotations.rds",
                                         save_path = "full.markers/pseudobulk.notFindMarkers.annotated_dacr.rds"){
  
  require(Seurat)
  require(Signac)
  require(tidyverse)
  require(DESeq2)
  
  mtx <- as.matrix(so@assays[[assay]]$counts)
  ind <- which(apply(mtx, 1, function(x) all(x == 0)))
  if (length(ind) != 0){
    print("[MSG] removing peaks with rowSums() == 0...")
    rm_peaks <- row.names(so)[row.names(so) %in% names(ind)]
    tmp <- as.matrix(so@assays[[assay]]@counts) + 1
    tmp <- tmp[-which(rownames(tmp) %in% rm_peaks),]
    mtx <- tmp
  } else {
    print("[MSG] no peaks with rowSums() == 0...")
    mtx <- as.matrix(so@assays[[assay]]@counts) + 1
  }
  
  print(paste0("[MSG] using Idents(so) = ",idents,"..."))
  Idents(so) <- idents
  print(paste0("[MSG] using idents.1 = ",idents.1,"..."))
  cells.1 <- WhichCells(so, idents = idents.1)
  
  if(!is.null(idents.2)){
    print(paste0("[MSG] using idents.2 = ",idents.2,"..."))
    cells.2 <- WhichCells(so, idents = idents.2)
  } else {
    print(paste0("[MSG] comparing idents.1 = ",idents.1," vs all other idents..."))
    cells.2 <- WhichCells(so, idents = idents.1, invert = TRUE)
  }
  
  group.info <- data.frame(row.names = colnames(mtx))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  group.info <- group.info[!is.na(group.info$group),]
  
  print("[MSG] Run DESeq2...")
  dds1 <- DESeqDataSetFromMatrix(
    countData = mtx,
    colData = group.info,
    design = ~ group
  )
  
  dds1 <- estimateSizeFactors(object = dds1) #, type = "iterate")
  dds1 <- estimateDispersions(object = dds1, fitType = fitType, maxit = 100) #, fitType = "mean"/"parametric"
  
  # print plot
  print(plotDispEsts(dds1))
  
  dds1 <- nbinomWaldTest(object = dds1)
  
  resultsNames(dds1)
  res_pre <- results(dds1,
                     name = "group_Group2_vs_Group1",
                     alpha = 0.05)
  res <- lfcShrink(dds1,
                   coef = "group_Group2_vs_Group1",
                   res = res_pre,
                   type = "apeglm")
  res1 <- data.frame(res)
  res1$peaks <- rownames(res1)
  annotations <- readRDS(annotations_filepath)
  dacr_ann <- ClosestFeature(so[[assay]], annotation = annotations, regions = res1$peaks)
  colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region","query_region"))] <- c("closest_peaks","peaks")
  bulk.dacr <- full_join(res1, dacr_ann, by = "peaks") %>% distinct()
  
  print(paste0("[MSG] Saving file in ",save_path," ..."))
  saveRDS(bulk.dacr, save_path)
  return(bulk.dacr)
}

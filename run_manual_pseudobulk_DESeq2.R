run_manual_pseudobulk_DESeq2 <- function(so, idents, idents.1, idents.2 = NULL, assay,
                                         fitType = 'mean', test = "Wald",  # NEW PARAMETER
                                         annotations = readRDS("/researchers/nicole.saw/references/ensembl/EnsDb.Mmusculus.v79_annotations.rds"),
                                         save.loc = ".") {
  suppressPackageStartupMessages({
    require(Seurat)
    require(Signac)
    require(tidyverse)
    require(DESeq2)
  })
  
  # make output directory
  dir.create(file.path(save.loc), showWarnings = FALSE, recursive = TRUE)
  
  mtx <- so@assays[[assay]]$counts
  cell_bc <- colnames(so)[unlist(so[[idents]]) %in% c(idents.1,idents.2)]
  ind <- which(apply(mtx, 1, function(x) all(x == 0)))
  
  if (length(ind) != 0) {
    message("[MSG] removing peaks with rowSums() == 0...")
    mtx <- as.matrix(mtx[-ind,cell_bc]) + 1
  } else {
    message("[MSG] no peaks with rowSums() == 0...")
    mtx <- as.matrix(mtx[,cell_bc]) + 1
  }
  
  message("[MSG] using Idents(so) = ", idents, "...")
  Idents(so) <- idents
  message("[MSG] using idents.1 = ", idents.1, "...")
  cells.1 <- WhichCells(so, idents = idents.1)
  
  if (!is.null(idents.2)) {
    message("[MSG] using idents.2 = ", idents.2, "...")
    cells.2 <- WhichCells(so, idents = idents.2)
  } else {
    message("[MSG] comparing idents.1 = ", idents.1, " vs all other idents...")
    cells.2 <- WhichCells(so, idents = idents.1, invert = TRUE)
  }
  
  group.info <- data.frame(row.names = colnames(mtx))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)
  group.info <- group.info[!is.na(group.info$group), ]
  
  message("[MSG] Run DESeq2...")
  dds <- DESeqDataSetFromMatrix(
    countData = mtx,
    colData = group.info,
    design = ~ group
  )
  
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = fitType, maxit = 100)
  
  print(plotDispEsts(dds))
  
  if (tolower(test) == "lrt") {
    message("[MSG] Running LRT test...")
    dds <- DESeq(dds, test = "LRT", reduced = ~1)
    res <- results(dds)
  } else {
    message("[MSG] Running Wald test...")
    dds <- nbinomWaldTest(dds)
    resultsNames(dds)
    res_pre <- results(dds, name = "group_Group2_vs_Group1", alpha = 0.05)
    res <- lfcShrink(dds, coef = "group_Group2_vs_Group1", res = res_pre, type = "apeglm")
  }
  
  res1 <- data.frame(res)
  res1$peaks <- rownames(res1)
  
  dacr_ann <- ClosestFeature(so[[assay]], annotation = annotations, regions = res1$peaks)
  colnames(dacr_ann)[which(colnames(dacr_ann) %in% c("closest_region", "query_region"))] <- c("closest_peaks", "peaks")
  bulk.dacr <- full_join(res1, dacr_ann, by = "peaks") %>% distinct()
  
  message("[MSG] Saving file in ", save.loc, " ...")
  saveRDS(bulk.dacr, file.path(save.loc, paste0(idents.1," vs ",idents.2,".pseudobulk.DESeq2.annotated_dacr.rds")))
  
  return(bulk.dacr)
}

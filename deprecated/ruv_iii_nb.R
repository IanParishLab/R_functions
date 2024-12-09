# deprecated
# run_RUV_III_NB <- function(so, assay="RNA", genome, low_abundant_genes_limit = -5, libsize_drop_nmads = 2, mito_drop_nmads = 3){
# suppressWarnings(suppressPackageStartupMessages({
#   library(ruvIIInb) 
#   library(SingleCellExperiment)
#   library(scater)
#   library(scran)
#   library(scuttle)
#   library(edgeR)
#   library(SingleR)
#   library(celldex)
#   library(hrbrthemes)
#   library(tidyverse)
#   library(ggplot2)
#   library(uwot)
#   library(scMerge)
#   library(Seurat)
#   library(randomcoloR)
#   library(dittoSeq)
#   library(pheatmap)
#   library(gridExtra)
#   library(igraph)
#   library(DelayedArray)
# }))
# 
# dir.create(file.path(getwd(),"RUV-III-NB/"))
# setwd(file.path(getwd(),"RUV-III-NB/"))
# 
# print("Creating SingleCellExperiment...")
# # so is Seurat Object
# sce <- SingleCellExperiment(assays = list(counts = so@assays[[assay]]@counts), 
#                             colData = so@meta.data)
# 
# # 2.2.1 Filtering low quality cells
# # The cell-level quality control metrics such as the number of genes that have non-zero counts and the percentage of counts that comes from Mitochondrial genes for each cell can be computed and added to the SingleCellExperiment object as follows:
# if (genome == "mouse"){
#   sce <- addPerCellQCMetrics(x = sce,subsets=list(Mito=grep("Mt-",rownames(sce))))
# } else if (genome == "human"){
#   sce <- addPerCellQCMetrics(x = sce,subsets=list(Mito=grep("MT-",rownames(sce))))
# } else {
#   stop("genome not included - choose only human or mouse.")
# }
# 
# libsize_drop <- isOutlier(
#   metric = sce$total,
#   nmads = libsize_drop_nmads,
#   type = "lower",
#   log = TRUE)
# colData(sce)$libsize_drop<-libsize_drop
# 
# mito_drop <- isOutlier(
#   metric = colData(sce)$subsets_Mito_percent,
#   nmads = mito_drop_nmads,
#   type = "higher")
# colData(sce)$mito_drop<-mito_drop
# 
# # 2.2.2 Filtering low abundant genes
# # Similarly, we recommend computing gene-level quality control metrics such as the mean count across all cells for each gene, and the percentage of cells with non-zero counts for each gene. The code below adds these measures to the SingleCellExperiment object and flags genes with a low mean cell count (see Figure 2.3).
# sce <- addPerFeatureQCMetrics(x = sce)
# 
# #Remove genes with zero counts for each gene
# sce <- subset(sce, rowData(sce)$mean > 0 )
# 
# # detect low abundant genes
# lowcount_drop <- log(rowData(sce)$mean) < low_abundant_genes_limit
# 
# print("Subsetting SingleCellExperiment to filter out low quality cells...")
# # subset 
# print("Number of low counts genes:")
# print(table(lowcount_drop))
# print("Number of cell-level total counts with number of median absolute deviations (nMADs = 2) away from the median value:")
# print(table(libsize_drop))
# print("Number of cell-level mitochondrial percentage with number of median absolute deviations (nMADs = 3) away from the median value:")
# print(table(mito_drop))
# sce <- sce[!(lowcount_drop), !(libsize_drop | mito_drop)]
# 
# # 2.3 Normalising the data
# # 2.3.0.1 Using scHK as negative control genes
# # Reading in the control genes
# data(segList)
# if (genome == "mouse"){
#   ctl <- segList$mouse$mouse_scSEG
# } else if (genome == "human"){
#   ctl <- segList$human$human_scSEG
# } else {
#   stop("genome not included - choose only human or mouse")
# }
# 
# # Creating a logical vector to identify control genes
# rowData(sce)$ctlLogical<-rownames(assays(sce)$counts) %in% ctl
# 
# print("Perform initial clustering to identify pseudo-replicates...")
# # Perform initial clustering to identify pseudo-replicates
# sce <- computeSumFactors(sce,assay.type="counts")
# data_norm_pre <- sweep(assays(sce)$counts,2,sce$sizeFactor,'/')
# assays(sce, withDimnames=FALSE)$lognormcounts<- log(data_norm_pre+1)
# snn_gr_init <- buildSNNGraph(sce, assay.type = "lognormcounts")
# clusters_init <- igraph::cluster_louvain(snn_gr_init)
# sce$cluster_init <- factor(clusters_init$membership)
# 
# sce <- runUMAP(sce,exprs_values = "lognormcounts")
# sce <- runPCA(sce,exprs_values = "lognormcounts")
# 
# png("dittoDimPlot_cluster_init.png", width=10, height=10)
# print(dittoDimPlot(sce, "cluster_init"))
# dev.off()
# 
# 
# # Figure 2.5 indicates a potential association between biology (cell type) and an unwanted factor (in this case library size). This potential association will be taken into account when RUV-III-NB estimates and subsequently removes the variation due to the unwanted factors.
# sce$logLS <- log(colSums(assays(sce)$counts))
# 
# png("plotUMAP_cluster_init_logLS.png", width=10, height=10)
# plotUMAP(sce,colour_by='logLS')
# dev.off()
# 
# print("Running RUV-III-NB...")
# # 2.3.0.3 Running RUV-III-NB
# # The code shown below performs RUV-III-NB normalisation.
# # Construct the replicate matrix M using pseudo-replicates identified using initial clustering
# 
# M <- matrix(0,ncol(assays(sce)$counts),length(unique(sce$cluster_init)))
# cl<- sort(unique(as.numeric(unique(sce$cluster_init))))
# for(CL in cl){
#   M[which(as.numeric(sce$cluster_init)==CL),CL] <- 1
# }
# 
# #RUV-III-NB code
# ruv3nb_out <- fastruvIII.nb(Y=DelayedArray(assays(sce)$counts), # count matrix with genes as rows and cells as columns
#                             M=M, #Replicate matrix constructed as above
#                             ctl=rowData(sce)$ctlLogical, #A vector denoting control genes
#                             k=2, # dimension of unwanted variation factors
#                             ncores = 6
# )
# 
# sce_ruv3nb <- makeSCE(ruv3nb_out,cData=colData(sce))
# seurat_ruv3nb <- as.Seurat(sce_ruv3nb, counts = "counts", data ="logPAC")
# 
# print("Saving important objects...")
# saveRDS(ruv3nb_out, "ruv3nb_out.rds")
# saveRDS(sce_ruv3nb, "sce_ruv3nb.rds")
# saveRDS(seurat_ruv3nb, "seurat_ruv3nb.rds")
# 
# }

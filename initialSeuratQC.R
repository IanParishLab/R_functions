# initialSeuratQC(tmp, var_to_compare = "orig.ident", sampleName = "SHIELD004_vs_TIPTOE011", save.loc = "./")
initialSeuratQC <- function(seu.obj, var_to_compare, var_col, sampleName, save.loc){
  
  suppressPackageStartupMessages({
    library(Seurat)
    library(tidyverse)
  })
  
  seu.obj$var_to_compare <- seu.obj[[var_to_compare]]
  seu.obj$log10GenesPerUMI <- log10(seu.obj$nFeature_RNA) / log10(seu.obj$nCount_RNA)
  
  p <- list()
  # cell counts per sample
  p[[1]] <- seu.obj@meta.data %>% 
    ggplot(aes(x=var_to_compare, fill=var_to_compare)) + 
    geom_bar() +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    scale_fill_manual(values = var_col) +
    scale_color_manual(values = var_col) +
    ggtitle("Cell counts per sample")
  
  # Visualize the number UMIs/transcripts (nCount_RNA) per cell
  p[[2]] <- seu.obj@meta.data %>% 
    ggplot(aes(color=var_to_compare, x=nCount_RNA, fill= var_to_compare)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    ylab("Cell density") +
    # geom_vline(xintercept = 500) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    scale_fill_manual(values = var_col) +
    scale_color_manual(values = var_col) +
    ggtitle("nUMIs/transcripts\n(nCount_RNA) per cell")
  
  # Visualize the distribution of genes detected per cell via histogram
  p[[3]] <- seu.obj@meta.data %>% 
    ggplot(aes(color=var_to_compare, x=nFeature_RNA, fill= var_to_compare)) + 
    geom_density(alpha = 0.2) + 
    theme_classic() +
    scale_x_log10() + 
    # geom_vline(xintercept = 300) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    scale_fill_manual(values = var_col) +
    scale_color_manual(values = var_col) +
    ggtitle("Distribution of genes\n(nFeatures_RNA) detected per cell")
    
  
  # Visualize the distribution of genes detected per cell via boxplot
  p[[4]] <- seu.obj@meta.data %>% 
    ggplot(aes(x=var_to_compare, y=log10(nFeature_RNA), fill=var_to_compare)) + 
    geom_boxplot() + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    scale_fill_manual(values = var_col) +
    scale_color_manual(values = var_col) +
    ggtitle("Distribution of genes\ndetected per cell")
  
  # Visualize the correlation between genes detected and number of UMIs and determine whether 
  # strong presence of cells with low numbers of genes/UMIs
  p[[5]] <- seu.obj@meta.data %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mito)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~var_to_compare) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    ggtitle("Correlation between genes detected and nUMIs to detect strong presence of\ncells with low numbers of genes/UMIs")
  
  # Visualize the distribution of mitochondrial gene expression detected per cell
  p[[6]] <- seu.obj@meta.data %>% 
    ggplot(aes(color=var_to_compare, x=percent.mito, fill=var_to_compare)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    geom_vline(xintercept = 0.2) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    scale_fill_manual(values = var_col) +
    scale_color_manual(values = var_col) +
    ggtitle("Mitochondrial gene expression\ndetected per cell")
  
  # Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
  p[[7]] <- seu.obj@meta.data %>%
    ggplot(aes(x=log10GenesPerUMI, color=var_to_compare, fill=var_to_compare)) +
    geom_density(alpha = 0.2) +
    theme_classic() +
    geom_vline(xintercept = 0.8) +
    theme(plot.title = element_text(size = 10, face="bold")) +
    scale_fill_manual(values = var_col) +
    scale_color_manual(values = var_col) +
    ggtitle("Overall complexity of the gene expression\naccording to genes/UMI")
  
  pdf(file.path(save.loc, paste0(sampleName,".initialQC.pdf")), width = 12, height = 6)
  print(ggpubr::ggarrange(plotlist = p, ncol = 4, nrow = 2, common.legend = TRUE))
  dev.off()
  
  return(p)
}

# # example usage 1
# plotCombined <- readRDS("/researchers/nicole.saw/projects/Sara_Roth_snPATHOseq/R/plotCombined.rds")
# seqtype_col<- setNames(c("#007756", "#D99BBD"), c("snPATHOseq","scRNAseq"))
# p <- initialSeuratQC(
#   plotCombined, 
#   var_to_compare = "seqType",
#   sampleName = "TIPTOE003 scRNAseq vs snPATHOseq",
#   var_col = seqtype_col,
#   save.loc = "../../Sara_Roth_snPATHOseq/R/"
#   )
# 
# sampleName = "TIPTOE003 scRNAseq vs snPATHOseq"
# save.loc = "/researchers/nicole.saw/projects/Sara_Roth_snPATHOseq/R/"
# 
# a <- p[[1]]/p[[2]]/p[[3]] & scale_color_manual(values = seqtype_col)
# b <- p[[5]]
# pdf(file.path(save.loc, paste0(sampleName,".initialQC.subset.pdf")), width = 8, height = 5)
# # print(ggpubr::ggarrange(plotlist = p[c(1,5,2,3)], ncol = 2, nrow = 2, common.legend = TRUE))
# print((a | b + theme(aspect.ratio = 1)) + plot_layout(widths = c(0.15,0.25)) & theme(legend.position = 'none'))
# dev.off()

# # example usage 2
# so <- readRDS("/researchers/nicole.saw/projects/SHIELD/R/first_pass_merged/tiptoe_shield.rds")
# so$orig.ident[so$orig.ident == "TIPTOE"] <-  "NextGEM"
# so$orig.ident[so$orig.ident == "SHIELD"] <-  "GEM-X"
# seqtype_col<- setNames(c("#06A5FF","#00446B"), c("NextGEM","GEM-X"))
# p <- initialSeuratQC(
#   so, 
#   var_to_compare = "orig.ident",
#   sampleName = "NextGEM vs GEM-X",
#   var_col = seqtype_col,
#   save.loc = "/researchers/nicole.saw/projects/SHIELD/R/"
# )
# 
# sampleName = "NextGEM vs GEM-X"
# save.loc = "/researchers/nicole.saw/projects/SHIELD/R/"
# 
# a <- p[[1]]/p[[2]]/p[[3]] & scale_color_manual(values = seqtype_col)
# b <- p[[5]]
# pdf(file.path(save.loc, paste0(sampleName,".initialQC.subset.pdf")), width = 10, height = 5)
# # print(ggpubr::ggarrange(plotlist = p[c(1,5,2,3)], ncol = 2, nrow = 2, common.legend = TRUE))
# print((a | b + theme(aspect.ratio = 1)) + plot_layout(widths = c(0.15,0.25)) & theme(legend.position = 'none'))
# dev.off()
# 
# so1 <- merge(subset(plotCombined, seqType == "snPATHOseq"), so)
# so1$orig.ident[so1$seqType == "snPATHOseq"] <- "snPATHOseq"
# so1$orig.ident[so1$orig.ident == "TIPTOE"] <- "NextGEM"
# so1$orig.ident[so1$orig.ident == "SHIELD"] <- "GEM-X"
# # SHIELD vs TIPTOE
# seqtype_col<- setNames(c("#06A5FF","#00446B","#007756"), c("NextGEM","GEM-X","snPATHOseq"))
# p <- initialSeuratQC(
#   so1, 
#   var_to_compare = "orig.ident",
#   sampleName = "NextGEM vs GEM-X vs snPATHOseq",
#   var_col = seqtype_col,
#   save.loc = "/researchers/nicole.saw/projects/SHIELD/R/"
# )
# 
# sampleName = "NextGEM vs GEM-X vs snPATHOseq"
# save.loc = "/researchers/nicole.saw/projects/SHIELD/R/"
# 
# a <- p[[1]]/p[[2]]/p[[3]] & scale_color_manual(values = seqtype_col)
# b <- p[[5]]
# pdf(file.path(save.loc, paste0(sampleName,".initialQC.subset.pdf")), width = 10, height = 5)
# # print(ggpubr::ggarrange(plotlist = p[c(1,5,2,3)], ncol = 2, nrow = 2, common.legend = TRUE))
# print((a | b + theme(aspect.ratio = 1)) + plot_layout(widths = c(0.15,0.25)) & theme(legend.position = 'none'))
# dev.off()
# 

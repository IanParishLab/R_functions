# source: "/researchers/nicole.saw/projects/Sinead_Reading/R/7_plots.Rmd"

# plotChromVARHeatmap -------------------------------------------------------
plotChromVARHeatmap <- function(so, ident = "group", motif.assay = "MACS3", chromvar.assay = "chromvar", logfc_threshold = 0, abs = FALSE, ...) {
  #  identify differential motifs
  Idents(so) <- ident
  markers <- FindAllMarkers(so, assay = chromvar.assay, logfc.threshold = 0, min.pct = 0.1) %>%
    dplyr::filter(p_val_adj < 0.05)
  
  # abs?
  if(isTRUE(abs)){
    markers <- markers %>% dplyr::filter(abs(avg_log2FC) > logfc_threshold)
  } else {
    markers <- markers %>% dplyr::filter(avg_log2FC > logfc_threshold)
  }

  # Retrieve motif names
  Idents(so) <- ident
  DefaultAssay(so) <- motif.assay
  motif_obj <- Motifs(so)
  pfm <- GetMotifData(motif_obj, slot = "pwm")
  motif_obj <- SetMotifData(motif_obj, slot = "pwm", new.data = pfm)

  # Average chromVAR expression of differential motifs
  avg_expr <- AverageExpression(so, assays = chromvar.assay, group.by = ident, layer = "data",
                                return.seurat = TRUE, features = markers$gene)
  avg_expr[["hash.ID"]] <- rownames(avg_expr@meta.data)
  Idents(avg_expr) <- "hash.ID"

  chromvar_data <- GetAssayData(avg_expr, assay = chromvar.assay, layer = "data")
  motif_names <- ConvertMotifID(motif_obj, id = rownames(chromvar_data), assay = motif.assay)

  mat <- as.data.frame(cbind(motif.names = motif_names, chromvar_data))
  mat[,2:ncol(mat)] <-  mat[,2:ncol(mat)] %>% mutate_if(is.character,as.numeric)
  rownames(mat) <- motif_names

  # Heatmap
  mat_num <- as.matrix(as.data.frame(mat[,-1]))
  val_range <- seq(min(mat_num), max(mat_num), length.out = 11)
  col_fun <- colorRamp2(val_range, rev(RColorBrewer::brewer.pal(n = length(val_range), name = "RdYlBu")))
  
  ComplexHeatmap::Heatmap(
    mat_num,
    column_names_rot = 90, column_names_centered = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 10),
    col = col_fun,
    name = "TF motif\nchromVAR score",
    ...)
}

# example usage:
# plotChromVARHeatmap(so, ident = 'group', logfc_threshold = 0.75,
#                     cluster_rows = TRUE, cluster_columns = FALSE,
#                     show_row_names = TRUE, show_row_dend = FALSE, show_column_dend = FALSE)

# plotVolcanoPanels -------------------------------------------------------
plotVolcanoPanels <- function(res, label_list, label_colors) {
  plots <- vector("list", length(label_list))

  for (i in seq_along(label_list)) {
    group_genes <- label_list[[i]]
    keyvals <- ifelse(res$gene %in% group_genes, label_colors[[i]], "lightgrey")
    names(keyvals)[keyvals == label_colors[[i]]] <- names(label_colors)[i]

    plots[[i]] <- EnhancedVolcano(res, lab = res$gene,
                                  x = 'avg_log2FC', y = 'p_val_adj',
                                  title = "", subtitle = "", caption = "",
                                  ylab = bquote(~-Log[10] ~ italic(P) ~adjusted),
                                  pCutoff = 0.05, FCcutoff = 0,
                                  labSize = 3,
                                  pointSize = ifelse(res$gene %in% group_genes, 3, 1),
                                  colCustom = keyvals,
                                  selectLab = group_genes,
                                  drawConnectors = TRUE, widthConnectors = 0.75) +
      theme_minimal() + theme(legend.position = "none")
  }

  combined <- ggpubr::ggarrange(plotlist = plots, nrow = 1, ncol = 3)
  annotated <- ggpubr::annotate_figure(combined,
                                       top = text_grob("Pseudobulk DEGs in Exh-Prog/Pre I cells in Tcf1+ subset",
                                                       face = "bold", size = 14))
  return(annotated)
}

# deg_labels <- list(
#   Stem = c("Tcf4","Myb","Wls","Cux1","Tcf7","Cd28","Il7r"),
#   Signalling = c("Gzmk","Runx1","Klf13","Batf","Hif1a","Satb1"),
#   Effector = c("Cacna1d","Akt3","Itpkb","Pip4k2a","Prkcb","Prkch")
# )
# deg_colors <- c("Stem" = "#FCD116", "Signalling" = "#009E60", "Effector" = "#1eacbd")
# 
# volcano_panel <- plotVolcanoPanels(res, deg_labels, deg_colors)
# ggsave(file.path(save.loc, "volcano_pbdeg2.pdf"), volcano_panel, width = 12, height = 4)


# averaged_heatmap -------------------------------------------------------
averaged_heatmap <- function(markers=NULL, features=NULL, group.by, so, method = "aggregate", heatmap_title, assay = c("MACS3", "RNA"), genes_to_label = NULL, ...) {
  suppressPackageStartupMessages({
  library(tidyverse)
  library(circlize)
  library(colorRamp2)
  library(ComplexHeatmap)
  })

  DefaultAssay(so) <- assay
  
  # Choose feature column and heatmap name
  if (assay != "RNA") {
    if(is.null(features)){
      features <- if(!is.null(markers)) markers$peaks else rownames(so[[assay]])
    }
    feature_col <- "peaks"
    gene_name_col <- "gene_name"
    heatmap_label <- "Averaged\nchromatin\naccessibility"
  
    } else if (assay == "RNA"){
    if(is.null(features)){
      features <- if(!is.null(markers)) markers$gene else rownames(so[["RNA"]])
    }
    feature_col <- "gene"
    gene_name_col <- "gene"
    if(method == "aggregate"){
    heatmap_label <- "Aggregated\ngene\nexpression"
    } else if(method == "average"){
      heatmap_label <- "Averaged\ngene\nexpression"
    }
  }
  
  # Compute average expression or accessibility
  if(method == "aggregate"){
    averaged <- AggregateExpression(so,
                                    group.by = group.by,
                                    features = features,
                                    assays = assay)
  } else if(method == "average"){
    averaged <- AggregateExpression(so,
                                    group.by = group.by,
                                    features = features,
                                    layer = "data",
                                    assays = assay)
  }
  averaged <- as.data.frame(averaged[[assay]])
  averaged <- t(scale(t(averaged))) %>%
    as.data.frame() %>%
    na.omit()
  
  # Get gene/peak names for annotation
  if (assay == "MACS3") {
    if(is.null(markers)){
      stop("need `markers` to annotate peak with gene name")
    }
    gene_name_mapping <- setNames(markers[[feature_col]], markers[[gene_name_col]])
    allnames <- names(gene_name_mapping)
  } else if (assay == "RNA") {
    allnames <- rownames(averaged)
    gene_name_mapping <- setNames(allnames, allnames)
  }
  
  # Determine which labels to annotate
  if (!is.null(genes_to_label)) {
    genes_to_label <- names(gene_name_mapping)[names(gene_name_mapping) %in% genes_to_label]
    index_w_genes <- which(rownames(averaged) %in% gene_name_mapping[genes_to_label])
    ra <- rowAnnotation(foo = anno_mark(index_w_genes,
                                        labels =  genes_to_label,
                                        labels_gp = gpar(fontsize = 8), padding = unit(.1, "mm")))
    
  } else {
    genes_to_label <- names(gene_name_mapping)
    index_w_genes <- integer(0)
    ra <- NULL
  }
  
  hm <- ComplexHeatmap::Heatmap(
    averaged,
    right_annotation = ra,
    row_title = heatmap_title,
    row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    col = colorRamp2(c(min(averaged), 0, max(averaged)), c("#4575b4", "#ffffbf", "#d73027")),
    name = heatmap_label,
    ...)

  return(hm)
}

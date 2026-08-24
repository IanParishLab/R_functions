# Helper functions ---------------------------------------------------------------------------------------------------
# --- set colour palettes and define colour generating function --------------------------
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cluster_col <- setNames(dittoSeq::dittoColors(1),seq(from=0,to=39))

discrete_pal <- function(x) {
  cols <- x %>% unique() %>% length() %>% Seurat::DiscretePalette(n = .)
  return(setNames(cols, x))
}

# --- get_cell_counts function -------------------------------------------------------------
get_cell_counts <- function(data, sample_name, stage, save.loc, overwrite = FALSE) {
  # Ensure directory exists
  counts_dir <- here(save.loc)
  if (!dir.exists(counts_dir)) dir.create(counts_dir, recursive = TRUE)
  
  # Output file
  output_file <- file.path(counts_dir, paste0(sample_name, "_cell_counts.tsv"))
  
  # Load previous if exists
  if (file.exists(output_file)) {
    df <- read.delim(output_file)
    if (!overwrite && stage %in% df$stage) return(invisible(NULL))
  } else {
    df <- data.frame(stage = character(), n_cells = numeric(), stringsAsFactors = FALSE)
  }
  
  # Count cells
  n_cells <- if (inherits(data, "Seurat")) ncol(data) else if (inherits(data, "dgCMatrix")) ncol(data) else if (inherits(data, "SingleCellExperiment")) ncol(data) else NA
  
  # Add row
  df <- rbind(df, data.frame(stage = stage, n_cells = n_cells))
  
  # Save
  write.table(df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("[INFO] Recorded ", n_cells, " cells at stage: ", stage, " for sample: ", sample_name)
}

# plotChromVARHeatmap -------------------------------------------------------
# source: "/researchers/nicole.saw/projects/Sinead_Reading/R/7_plots.Rmd"
plotChromVARHeatmap <- function(so, ident = "group", averaged.by.ident = NULL, motif.assay = "MACS3", 
                                chromvar.assay = "chromvar", logfc_threshold = 0, abs = FALSE,
                                return_table = FALSE, motifs = NULL, ...) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(colorRamp2)
  })
  #  identify differential motifs
  Idents(so) <- ident
  markers <- FindAllMarkers(so, assay = chromvar.assay, logfc.threshold = logfc_threshold, min.pct = 0.1) %>%
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
  if(!is.null(averaged.by.ident)){
    group.by <- averaged.by.ident
  } else {
    group.by <- ident
  }
  avg_expr <- AverageExpression(so, assays = chromvar.assay, group.by = group.by, layer = "data",
                                return.seurat = TRUE, features = markers$gene)
  avg_expr[["hash.ID"]] <- rownames(avg_expr@meta.data)
  Idents(avg_expr) <- "hash.ID"
  
  chromvar_data <- GetAssayData(avg_expr, assay = chromvar.assay, layer = "data")
  motif_names <- ConvertMotifID(motif_obj, id = rownames(chromvar_data), assay = motif.assay)
  
  mat <- as.data.frame(cbind(motif.names = motif_names, chromvar_data))
  mat[,2:ncol(mat)] <-  mat[,2:ncol(mat)] %>% mutate_if(is.character,as.numeric)
  rownames(mat) <- motif_names
  
  if(isTRUE(return_table)){
    tbl <- as.data.frame(cbind(motif.names = motif_names, chromvar_data)) %>% 
      rownames_to_column(var = "gene") %>% 
      left_join(markers, ., by = "gene")
    return(tbl)
  }
  # Heatmap
  mat_num <- as.matrix(as.data.frame(mat[,-1]))
  if(!is.null(motifs)){
    mat_num <- mat_num[which(rownames(mat_num) %in% motifs),]
  }
  
  val_range <- seq(min(mat_num), max(mat_num), length.out = 11)
  col_fun <- colorRamp2(val_range, rev(RColorBrewer::brewer.pal(n = length(val_range), name = "RdYlBu")))
  
  ComplexHeatmap::Heatmap(
    mat_num,
    column_names_rot = 90, column_names_centered = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_names_gp = grid::gpar(fontsize = 10),
    col = col_fun,
    name = "Motif\ndeviation\nscore",
    ...)
}

# example usage:
# plotChromVARHeatmap(so, ident = 'group', logfc_threshold = 0.75,
#                     cluster_rows = TRUE, cluster_columns = FALSE,
#                     show_row_names = TRUE, show_row_dend = FALSE, show_column_dend = FALSE)

# plotVolcanoPanels -------------------------------------------------------
plotVolcanoPanels <- function(res, x = 'avg_log2FC', y = 'p_val_adj', label_column = "gene", label_list, label_colors, nrow = 1, ncol, plot_title,...) {
  library(EnhancedVolcano)
  plots <- vector("list", length(label_list))
  
  for (i in seq_along(label_list)) {
    group_genes <- label_list[[i]]
    keyvals <- ifelse(res[[label_column]] %in% group_genes, label_colors[[i]], "lightgrey")
    names(keyvals)[keyvals == label_colors[[i]]] <- names(label_colors)[i]
    
    plots[[i]] <- EnhancedVolcano(res, lab = res[[label_column]],
                                  x = x, y = y,
                                  title = "", subtitle = names(label_list)[i], caption = "",
                                  # pCutoff = 0.05, FCcutoff = 0,
                                  pointSize = ifelse(res[[label_column]] %in% group_genes, 3, 1),
                                  colCustom = keyvals,
                                  selectLab = group_genes,
                                  drawConnectors = TRUE, widthConnectors = 0.1,
                                  ...) +
      theme_minimal() + theme(legend.position = "none")
  }
  
  combined <- ggpubr::ggarrange(plotlist = plots, nrow = nrow, ncol = ncol)
  annotated <- ggpubr::annotate_figure(combined,
                                       top = text_grob(plot_title, face = "bold", size = 14))
  return(annotated)
}

# deg_labels <- list
#   Stem = c("Tcf4","Myb","Wls","Cux1","Tcf7","Cd28","Il7r"),
#   Signalling = c("Gzmk","Runx1","Klf13","Batf","Hif1a","Satb1"),
#   Effector = c("Cacna1d","Akt3","Itpkb","Pip4k2a","Prkcb","Prkch")
# )
# deg_colors <- c("Stem" = "#FCD116", "Signalling" = "#009E60", "Effector" = "#1eacbd")
# 
# volcano_panel <- plotVolcanoPanels(res, deg_labels, deg_colors)
# ggsave(file.path(save.loc, "volcano_pbdeg2.pdf"), volcano_panel, width = 12, height = 4)


# averaged_heatmap -------------------------------------------------------
averaged_heatmap <- function(
    markers = NULL, features = NULL, group.by, so,
    method = c("aggregate", "average"),
    heatmap_title, assay = c("MACS3", "RNA"),
    cols = c("#4575b4", "#ffffbf", "#d73027"),
    genes_to_label = NULL, return_table = FALSE, ...
) {
  suppressPackageStartupMessages({
    library(tidyverse)
    library(circlize)
    library(colorRamp2)
    library(ComplexHeatmap)
  })
  
  method <- match.arg(method)
  assay <- match.arg(assay)
  DefaultAssay(so) <- assay
  
  # Define method-specific metadata
  method_fns <- list(
    aggregate = AggregateExpression,
    average   = AverageExpression
  )
  label_templates <- list(
    RNA   = c(aggregate = "Aggregated\ngene\nexpression",
              average   = "Averaged\ngene\nexpression"),
    MACS3 = c(aggregate = "Aggregated\nchromatin\naccessibility",
              average   = "Averaged\nchromatin\naccessibility")
  )
  feature_cols <- list(
    RNA   = c(feature = "gene", gene_name = "gene"),
    MACS3 = c(feature = "peaks", gene_name = "gene_name")
  )
  
  # Select features
  if (is.null(features)) {
    if (assay == "RNA") {
      features <- if (!is.null(markers)) markers$gene else rownames(so[["RNA"]])
    } else {
      features <- if (!is.null(markers)) markers$peaks else rownames(so[[assay]])
    }
  }
  
  # Run averaging/aggregation
  averaged <- method_fns[[method]](
    so,
    group.by = group.by,
    features = features,
    assays = assay
  )[[assay]] %>%
    as.data.frame() %>%
    { t(scale(t(.))) } %>%
    as.data.frame() %>%
    na.omit()
  
  # Gene/peak name mapping
  feature_col <- feature_cols[[assay]]["feature"]
  gene_name_col <- feature_cols[[assay]]["gene_name"]
  
  if (assay == "MACS3") {
    if (is.null(markers)) {
      stop("need `markers` to annotate peak with gene name")
    }
    gene_name_mapping <- setNames(markers[[feature_col]], markers[[gene_name_col]])
    allnames <- names(gene_name_mapping)
  } else {
    allnames <- rownames(averaged)
    gene_name_mapping <- setNames(allnames, allnames)
  }
  
  # Annotate genes
  if (!is.null(genes_to_label)) {
    
    genes_to_label <- allnames[allnames %in% genes_to_label]
    index_w_genes <- which(allnames %in% genes_to_label)
    
    if (length(genes_to_label) != length(index_w_genes)){
      stop("`genes_to_label` & `index_w_genes` have different lengths")
    }
    
    ra <- rowAnnotation(
      foo = anno_mark(index_w_genes, labels = genes_to_label, 
                      labels_gp = gpar(fontsize = 8), padding = unit(.1, "mm") )
    )
  } else {
    ra <- NULL
  }
  
  # Build heatmap
  hm <- ComplexHeatmap::Heatmap(
    averaged,
    right_annotation = ra,
    row_title = heatmap_title,
    # row_names_gp = grid::gpar(fontsize = 10),
    column_names_gp = grid::gpar(fontsize = 10),
    col = colorRamp2(c(min(averaged), 0, max(averaged)), cols),
    name = label_templates[[assay]][[method]],
    ...
  )
  
  if (isTRUE(return_table)) {
    draw(hm)
    return(averaged)
  } else {
    return(hm)
  }
}

# --- get_max_clus function  ------------------------------------------------------------- 
get_max_clus <- function(seu.obj) { unique(seu.obj$seuratClusters) %>% levels() %>% as.numeric %>% max } 

# --- plot_perCluster function -------------------------------------------------------------
plot_perCluster <- function(seu.obj, clusters, ident = "seuratClusters", reduction = 'umap', cluster_col = NULL,
                            ncol, nrow, width = 12, height = 10, save.loc = NULL, print = FALSE, ...) {
  
  Idents(seu.obj) <- ident
  plot <- lapply(as.character(clusters), function(i) {
    DimPlot(seu.obj, 
            group.by = ident, 
            reduction = reduction,
            cells.highlight = WhichCells(seu.obj, idents = i), 
            cols.highlight = cluster_col[i], 
            order = TRUE, repel = TRUE, ...) + 
      NoLegend() + NoAxes() + ggtitle(i) + 
      theme(aspect.ratio = 1, title = element_text(size = 7), text = element_text(size = 7))
  }) %>% 
    ggpubr::ggarrange(plotlist = ., ncol = ncol, nrow = nrow)
  
  if(isTRUE(print)){
    return(plot)
  } else {
    png(file.path(save.loc, "plots", paste0(seu.obj@project.name,".", ident,".perCluster.png")), width = 12, height = 10, units = 'in', res = 300)
    print(plot)
    dev.off()
  }
  
  return(plot)
}

# --- plot_perFeature function -------------------------------------------------------------
plot_perFeature <- function(seu.obj, feature, ident = "PatientID", reduction = 'umap.harmony', 
                            cols = c("lightgrey", "darkred"), 
                            ncol = NULL, nrow = NULL, width = 10, height = 10, 
                            save.loc = NULL, print = FALSE, ...) {
  
  Idents(seu.obj) <- ident
  ident.list <- unique(Idents(seu.obj))
  
  # 1. Fetch ALL background coordinate limits so ggplot never clips any cells
  embed_coords <- as.data.frame(Embeddings(seu.obj, reduction = reduction))
  x_col <- colnames(embed_coords)[1]
  y_col <- colnames(embed_coords)[2]
  
  x_lims <- range(embed_coords[[x_col]], na.rm = TRUE)
  y_lims <- range(embed_coords[[y_col]], na.rm = TRUE)
  
  plot_list <- lapply(ident.list, function(i) {
    target_cells <- WhichCells(seu.obj, idents = i)
    
    # 2. Base FeaturePlot on target cells
    p <- FeaturePlot(seu.obj, 
                     features = feature,
                     cells = target_cells,
                     reduction = reduction,
                     cols = cols,
                     order = TRUE, 
                     ...)
    
    # 3. Add background layer FIRST & lock axis limits to the FULL dataset limits
    p$layers <- c(
      geom_point(data = embed_coords, 
                 aes(x = .data[[x_col]], y = .data[[y_col]]), 
                 color = "lightgrey", 
                 size = 0.5, 
                 inherit.aes = FALSE),
      p$layers
    )
    
    p <- p + 
      coord_cartesian(xlim = x_lims, ylim = y_lims) + # Forces full viewport
      NoLegend() + 
      NoAxes() + 
      ggtitle(paste0(feature, " - ", i)) + 
      theme(aspect.ratio = 1, 
            title = element_text(size = 8), 
            text = element_text(size = 8))
    
    return(p)
  })
  
  combined_plot <- ggpubr::ggarrange(plotlist = plot_list, ncol = ncol, nrow = nrow)
  
  if (!isTRUE(print) && !is.null(save.loc)) {
    dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
    png(file.path(save.loc, "plots", paste0(seu.obj@project.name, ".", feature, ".perFeature.png")), 
        width = width, height = height, units = 'in', res = 300)
    print(combined_plot)
    dev.off()
  }
  
  return(combined_plot)
}

# --- get_cell_counts function -------------------------------------------------------------
get_cell_counts <- function(data, sample_name, stage, save.loc, overwrite = FALSE) {
  # Ensure directory exists
  counts_dir <- here(save.loc)
  if (!dir.exists(counts_dir)) dir.create(counts_dir, recursive = TRUE)
  
  # Output file
  output_file <- file.path(counts_dir, paste0(sample_name, "_cell_counts.tsv"))
  
  # Load previous if exists
  if (file.exists(output_file)) {
    df <- read.delim(output_file)
    if (!overwrite && stage %in% df$stage) return(invisible(NULL))
  } else {
    df <- data.frame(stage = character(), n_cells = numeric(), stringsAsFactors = FALSE)
  }
  
  # Count cells
  if (any(inherits(data, c("Seurat", "dgCMatrix", "dgTMatrix", "SingleCellExperiment")))) {
    n_cells <- ncol(data)
  } else {
    n_cells <- NA
    stop("Did not record cell counts as data format is not `Seurat`, `dgCMatrix`, `dgTMatrix`, `SingleCellExperiment`")
  }
  
  # Add row
  df <- rbind(df, data.frame(stage = stage, n_cells = n_cells))
  
  # Save
  write.table(df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  message("[INFO] Recorded ", n_cells, " cells at stage: ", stage, " for sample: ", sample_name)
}

# --- get_blacklist function -------------------------------------------------------------
get_blacklist <- function(seu, species = "human", include_ribo = FALSE, unlist = TRUE) {
  
  # # include ribosomal genes and proteins in the blacklist?
  # include_ribo <- if(grepl("CD8", save.loc)) FALSE else TRUE
  
  # Validate species
  species <- tolower(species)
  if (!species %in% c("human", "mouse")) {
    stop("`species` must be either 'human' or 'mouse'")
  }
  
  # Define regex for IG genes and ribosomal genes by species
  ig_grep <- ifelse(species == "human", "^IG[HKL]", "^Ig[hlk]")
  ribo_grep <- ifelse(species == "human", "^RP[SL]", "^Rp[sl]")  # RPS/RPL vs Rps/Rpl
  annotations <- if (species == "human") qs_read("/researchers/nicole.saw/references/ensembl/EnsDb.Hsapiens.v86_annotations.qs") else qs_read("/researchers/nicole.saw/references/ensembl/EnsDb.Mmusculus.v79_annotations.qs") %>% data.frame
  
  # Build blacklist components
  blacklist <- list(
    Mito = grep("^[Mm][Tt]-", rownames(seu), value = TRUE),
    TCR = intersect(
      rownames(seu), 
      scGate::genes.blacklist.default[[ifelse(species == "human", "Hs", "Mm")]]$TCR
    ),
    Immunoglobulin = grep(ig_grep, rownames(seu), value = TRUE),
    Ribosomal = if (include_ribo) grep(ribo_grep, rownames(seu), value = TRUE) else character(0)
  )
  
  # further filtering by gene biotype
  ind <- if (species == "human") which(annotations$GENEBIOTYPE %in% c("protein_coding", "LRG_gene")) else which(annotations$gene_biotype %in% c("protein_coding"))
  genes <- if (species == "human") annotations$SYMBOL[!ind] %>% unique else annotations$gene_name[!ind] %>% unique
  genes2 <- if (species == "human") setdiff(rownames(seu), annotations$SYMBOL[ind] %>% unique)  else setdiff(rownames(seu), annotations$gene_name[ind] %>% unique)
  
  blacklist <- unique(c(blacklist, genes, genes2))
  
  # check if needs deduplicating
  if (isTRUE(unlist)){
    # # Combine and deduplicate
    unique(unlist(blacklist))
  } else {
    blacklist
  }
}

# --- plot alluvial function -------------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
  library(dittoSeq)
})

# make ditto barplot ----
make_ditto_barplot <- function(seu, var, group.by, split.by = NULL, color.panel, ggtitle = "", ncol = 1, data.out = NULL){
  theme_custom <- theme(text = element_text(size = 7), 
                        axis.text.x = element_text(angle = 90, size = 7, hjust = 1))
  
  p <- dittoBarPlot(seu, var = var, group.by = group.by, split.by = split.by, split.ncol = ncol,
                    xlab = "celltypist_label", retain.factor.levels = TRUE, color.panel = color.panel) +
    ggtitle(ggtitle) +
    theme_custom +
    guides(color = guide_legend(ncol = 2)) + coord_flip()
  
  print(p)
  
  if(isTRUE(data.out)){
    df <- dittoBarPlot(seu, var = var, group.by = group.by, split.by = split.by, split.ncol = ncol,
                       xlab = "celltypist_label", retain.factor.levels = TRUE, color.panel = color.panel, 
                       data.out = data.out)
    return(df)
  } else {
    return(p)
  }
}
# plot alluvial ----
plot_alluvial <- function(
    seu,
    var = "High.hierarchy.cell.types",
    group.by = "PatientID",
    split.by = "sampleSourceSuperset",
    color.panel = NULL,
    ncol = 5,
    x_angle = 60,
    alluvium_width = 1/6,
    alpha = 0.75
) {
  # 1. Extract data using dittoSeq
  plot_df <- make_ditto_barplot(
    seu = seu,
    var = var,
    group.by = group.by,
    split.by = split.by,
    color.panel = color.panel,
    data.out = TRUE
  )
  
  # 2. Construct ggplot using dynamic symbol evaluation (.data[[...]])
  p <- ggplot(
    plot_df$data,
    aes(
      x = .data[[split.by]],
      stratum = .data[["label"]],
      alluvium = .data[["label"]],
      y = .data[["percent"]],
      fill = .data[["label"]]
    )
  ) +
    geom_alluvium(width = alluvium_width, alpha = alpha) +
    geom_stratum(width = alluvium_width, color = "grey30") +
    facet_wrap(~ grouping, scales = "free_y", ncol = ncol) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.x = element_text(face = "bold", size = 11, angle = x_angle, hjust = 1),
      panel.grid.major.x = element_blank()
    )
  
  # 3. Add custom colors if provided
  if (!is.null(color.panel)) {
    p <- p + scale_fill_manual(values = color.panel)
  }
  
  return(p)
}
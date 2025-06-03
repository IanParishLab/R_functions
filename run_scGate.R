run_scGate <- function(seu.obj, models.list, bk.list, min.cells, seed.use, 
                       gsub_pattern = paste0(names(models.list),"_", collapse = "|"), 
                       name, save.loc, save.name, gate_cell_col = NULL, gate_model_col = NULL){
  
  suppressPackageStartupMessages({
    library(scGate)
    library(RColorBrewer)
    library(ggpubr)
  })
  
  # make output directory
   dir.create(file.path(save.loc, "int_obj"), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(save.loc, "plots"), showWarnings = FALSE, recursive = TRUE)
  
  # use scGate
  print("[MSG] run scGate...")
  for(i in 1:length(models.list)){
    print(names(models.list)[i])
    
    if("scGate_multi" %in% colnames(seu.obj@meta.data)){
      seu.obj@meta.data["scGate_multi"] <- NULL
    }
    
    seu.obj <- scGate(seu.obj, model = models.list[[i]], reduction = 'pca', keep.ranks = FALSE, save.levels = FALSE, 
                      genes.blacklist = bk.list, min.cells = min.cells, seed = seed.use)
    seu.obj[[names(models.list)[i]]] <- seu.obj$scGate_multi
  }
  
  #fix metadata
  message("[MSG] Fix metadata...")
  tmp <- seu.obj@meta.data[colnames(seu.obj@meta.data) %in% names(models.list)]
  for(i in 1:(length(colnames(tmp)))){
    tmp[,i][tmp[,i] == "Target"] <- colnames(tmp)[i]
  }
  tmp$scGate_model_def <- sapply(1:nrow(tmp), function(x){
    tmp[x,!is.na(tmp[x,])] %>% 
      paste0(.,collapse = ";") %>% 
      gsub(";not","",.)
  })
  tmp$contaminating_cell <- gsub(gsub_pattern, "", tmp$scGate_model_def)
  tmp$contaminating_cell <- lapply(1:nrow(tmp), function(x){
    ifelse(str_split(tmp$contaminating_cell[x], pattern = ";") %>% sapply(., unique) %>% length == 1,
           str_split(tmp$contaminating_cell[x], pattern = ";") %>% sapply(., unique),
           "Multi"
    )
  }) %>% unlist 
  tmp$contaminating_cell[which(nchar(tmp$contaminating_cell) == 0)] <- "not"
  tmp$scGate_model_def[which(nchar(tmp$scGate_model_def) == 0)] <- "not"
  
  # add contaminating_cell to seu.obj
  message("[MSG] Add `contaminating_cell` to seu.obj...")
  con_name_1 <- paste0('contaminating.', name)
  model_name_2 <- paste0('scGate_model_def.', name)
  if(all(rownames(tmp) == colnames(seu.obj))){
    seu.obj[[con_name_1]] <- tmp$contaminating_cell
    seu.obj[[model_name_2]] <- tmp$scGate_model_def
  }
  
  # rm old metadata
  seu.obj@meta.data[colnames(seu.obj@meta.data) %in% names(models.list)] <- NULL
  
  # save only the scGate labels
  message("[MSG] Save scGate cell labels...")
  predicted_scGate_labels <- seu.obj@meta.data[grep(paste0(name, collapse="|"), colnames(seu.obj@meta.data), value = TRUE)] 
  saveRDS(predicted_scGate_labels, file.path(save.loc, "int_obj", paste0(save.name,".scGate_models.list.result.rds")))
  
  # make plots
  message("[MSG] Make & save plots...")
  capture <- unique(seu.obj$orig.ident)  
  
  unique_gate_cell_type <- table(seu.obj[[con_name_1]]) %>% 
    unlist %>% 
    names %>% 
    gsub(paste0(capture,".", collapse="|"), "", .) %>%
    gsub("not",NA,.) %>%
    unique %>% 
    sort
  
  unique_gate_model <- table(seu.obj[[model_name_2]]) %>% 
    unlist %>% 
    names %>% 
    gsub(paste0(capture,".", collapse="|"), "", .)  %>%
    gsub("not",NA,.) %>%
    unique %>% 
    sort
  
  # set colours
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  if(is.null(gate_cell_col)){
    gate_cell_col <- setNames(gg_color_hue(length(unique_gate_cell_type)), unique_gate_cell_type)
  }
  if(is.null(gate_model_col)){
    gate_model_col <- setNames(gg_color_hue(length(unique_gate_model)), unique_gate_model)
  }
  theme = NoLegend() + NoAxes() + theme(aspect.ratio = 1, title = element_text(size = 7))
  
  # plot
  plot <- DimPlot(seu.obj, group.by = c(con_name_1, model_name_2), 
                  order = TRUE, label = FALSE, repel = TRUE, cols = c(gate_cell_col,gate_model_col)) & theme
  
  lgd1 = DimPlot(seu.obj, 
                 group.by = con_name_1,
                 cols = gate_cell_col) %>% get_legend %>% as_ggplot
  lgd2 = DimPlot(seu.obj, 
                 group.by = model_name_2,
                 cols = gate_model_col) %>% get_legend %>% as_ggplot
  
  # save plot & legends
  pdf(file.path(save.loc,"plots",paste0(save.name,".scGate_models.list.result.pdf")),height = 10, width = 6)
  print(plot)
  print(lgd1 | lgd2)
  dev.off()
  
  
  return(seu.obj)
}

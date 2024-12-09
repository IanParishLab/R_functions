run_enrichGO <- function(markers, split.by = "cluster", gene.col = "gene", direction = "both",
                         p.value.col = "p_val_adj", p.value.threshold = 0.05, logfc.col = "avg_log2FC", logfc.threshold = 0,
                         res.name = "DEMethod.projectName.assayName.DETest", organism = NULL, save.loc){
  require(clusterProfiler)
  require(enrichplot)
  
  if(organism == "mouse"){
    require(org.Mm.eg.db)
    OrgDb <-  org.Mm.eg.db
  } else if(organism == "human"){
    require(org.Hs.eg.db)
    OrgDb <-  org.Hs.eg.db
  }
  
  # markers = readRDS("full.markers.genotype/FindMarkers.TCF1plus.RNA.LR.Control_IkarosAiolos.rds")
  # split.by = NULL
  # direction = "both"
  # gene.col = "gene"
  # p.value.col = "p_val_adj"
  # p.value.threshold = 0.05
  # logfc.col = "avg_log2FC"
  # logfc.threshold = 0
  # res.name = "FindMarkers.TCF1plus.RNA.LR.Control_IkarosAiolos"
  
  ## read markers
  if (!is.null(split.by)){
    markers <- split(markers, markers[[split.by]])
  } else {
    markers <- list(markers)
  }
  
  # save df only
  ego_res_all_df <- lapply(markers, function(m){
    
    rm_gn <- grep("^Gm|Rik$", m[[gene.col]])
    if(length(rm_gn) > 1){
      m <- m[-rm_gn,]
    } else {
      m <- m
    }
    
    original_gene_list <- setNames(m[[logfc.col]], m[[gene.col]])
    gene_list <- na.omit(original_gene_list)
    gene_list <-  sort(gene_list, decreasing = TRUE)
    
    if (direction == "both" | direction == "up"){
      ## upregulated
      ego <- enrichGO(gene = m[[gene.col]][which(m[[p.value.col]] < 0.05 & m[[logfc.col]] > 0)],
                      universe = names(gene_list),
                      OrgDb = OrgDb,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05, 
                      pool = TRUE)
      ego <- data.frame(ego)
      if(nrow(ego)!=0){
        
        if (!is.null(split.by)){
          ego$cluster <- unique(m[[split.by]])
          ego$trend <- "upregulated"
        } else {
          ego$trend <- "upregulated"
        }
        
        terms <- ego$Description[ego$p.adjust < 0.05]
        print(length(terms))
        ego_up <- ego
      } else {
        ego_up <- data.frame()
      }
      
    }
    
    if (direction == "both" | direction == "dn"){
      ## downregulated
      ego <- enrichGO(gene = m[[gene.col]][which(m[[p.value.col]] < 0.05 & m[[logfc.col]] < 0)],
                      universe = names(gene_list),
                      OrgDb = OrgDb,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05, 
                      pool = TRUE)
      ego <- data.frame(ego)
      if(nrow(ego)!=0){
        
        if (!is.null(split.by)){
          ego$cluster <- unique(m[[split.by]])
          ego$trend <- "downregulated"
        } else {
          ego$trend <- "downregulated"
        }
        
        terms <- ego$Description[ego$p.adjust < 0.05]
        print(length(terms))
        ego_dn <- ego
      } else {
        ego_dn <- data.frame()
      }
      
    }
    
    ego <- plyr::rbind.fill(ego_up,ego_dn)
    return(ego)
  })
  ego_res_all_df <- plyr::rbind.fill(ego_res_all_df)
  saveRDS(ego_res_all_df, file.path(save.loc, paste0("enrichGOresDF.",res.name,".rds")))
  
  ## save result objects only
  if (direction == "both" | direction == "up"){
    # upregulated markers only
    ego_res_all_up <- lapply(markers, function(m){
      m <- m[grep("^Gm|Rik$", m[[gene.col]], invert = TRUE),]
      original_gene_list <- setNames(m[[logfc.col]], m[[gene.col]])
      gene_list <- na.omit(original_gene_list)
      gene_list <-  sort(gene_list, decreasing = TRUE)
      
      ## upregulated
      ego <- enrichGO(gene = m[[gene.col]][which(m[[p.value.col]] < 0.05 & m[[logfc.col]] > 0)],
                      universe = names(gene_list),
                      OrgDb = OrgDb,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05, 
                      pool = TRUE)
      return(ego)
    })
    saveRDS(ego_res_all_up, file.path(save.loc, paste0("enrichGOresObj.up.",res.name,".rds")))
  }
  
  # downregulated markers only
  if (direction == "both" | direction == "dn"){
    ego_res_all_dn <- lapply(markers, function(m){
      m <- m[grep("^Gm|Rik$", m[[gene.col]], invert = TRUE),]
      original_gene_list <- setNames(m[[logfc.col]], m[[gene.col]])
      gene_list <- na.omit(original_gene_list)
      gene_list <-  sort(gene_list, decreasing = TRUE)
      
      ## downregulated
      ego <- enrichGO(gene = m[[gene.col]][which(m[[p.value.col]] < 0.05 & m[[logfc.col]] < 0)],
                      universe = names(gene_list),
                      OrgDb = OrgDb,
                      keyType = "SYMBOL",
                      ont = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05, 
                      pool = TRUE)
      
      return(ego)
    })
    saveRDS(ego_res_all_dn, file.path(save.loc, paste0("enrichGOresObj.dn.",res.name,".rds")))
  }
}

# Cluster level comparison functions
library(Seurat)
library(tidyverse)
library(UpSetR)

read_metadata <- function(file) {
  seurat_obj <- readRDS(file)
  return(seurat_obj@meta.data)
}

get_jaccard <- function(set_a, set_b) {
  jaccardi <- length(intersect(set_a, set_b))/length(union(set_a, set_b))
  return(jaccardi)
}

get_cells_from_cluster <- function(seurat_obj_metadata, clust_name){
  clust_cells <- rownames(seurat_obj_metadata)[seurat_obj_metadata$seurat_clusters == clust_name] 
  return(clust_cells)
}

get_cluster_jaccard_df <- function(metadata1, metadata2){
  clust_by_clust_df <- expand.grid(unique(metadata1$seurat_clusters),
                                   unique(metadata2$seurat_clusters),stringsAsFactors = FALSE)
  
  clust_by_clust_df$ji <- apply(clust_by_clust_df,1,function(x,md1,md2){
    
    set_a <- get_cells_from_cluster(md1, x[1])
    set_b <- get_cells_from_cluster(md2, x[2])
    
    return(get_jaccard(set_a=set_a,set_b=set_b))
    
  },md1 = metadata1, md2 = metadata2)
  
  return(clust_by_clust_df)
}


plot_cluster_jaccard_matrix <- function(clust_by_clust_df){
  var1_order <- clust_by_clust_df[,c("Var1","ji")] %>%
    group_by(Var1) %>%
    summarise(max_ji = max(ji)) %>%
    arrange(desc(max_ji)) %>%
    pull(Var1)
  
  var2_order <- clust_by_clust_df[,c("Var2","ji")] %>%
    group_by(Var2) %>%
    summarise(max_ji = max(ji)) %>%
    arrange(desc(max_ji)) %>%
    pull(Var2)
  
  clust_by_clust_df %>%
    mutate(Var1 = ordered(Var1,levels = var1_order),
           Var2 = ordered(Var2,levels = var2_order)) %>%
    ggplot(aes(Var1, Var2,fill = ji)) +
    geom_tile() +
    scale_fill_viridis_c()
  
}

plot_cluster_jaccard_matrix_colorhigh <- function(clust_by_clust_df,cluster_high){
  var1_order <- clust_by_clust_df[,c("Var1","ji")] %>%
    group_by(Var1) %>%
    summarise(max_ji = max(ji)) %>%
    arrange(desc(max_ji)) %>%
    pull(Var1)
  
  var2_order <- clust_by_clust_df[,c("Var2","ji")] %>%
    group_by(Var2) %>%
    summarise(max_ji = max(ji)) %>%
    arrange(desc(max_ji)) %>%
    pull(Var2)
  
  cluster_high_lgl <- as.factor(var1_order %in% cluster_high)
  levels(cluster_high_lgl) <- c("black","red")
  clust_by_clust_df %>%
    mutate(Var1 = ordered(Var1,levels = var1_order),
           Var2 = ordered(Var2,levels = var2_order)) %>%
    ggplot(aes(Var1, Var2,fill = ji)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(colour = cluster_high_lgl))
  
}

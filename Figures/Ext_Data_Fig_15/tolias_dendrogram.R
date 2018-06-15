library(dplyr)
library(feather)
library(ggplot2)
library(dendextend)
options(stringsAsFactors = F)

source("prune_leaf_custom.R")

v1_alm_anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

v1_alm_anno <- v1_alm_anno %>%
  filter(cluster_id %in% 1:133)
# 
# dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")
# 
# alm_clusters <- unique(v1_alm_anno$cluster_id[grepl("ALM",v1_alm_anno$cluster_label)])
# #inh_clusters <- 1:60
# #nn_clusters <- 116:133
# 
# non_visp_clusters <- c(alm_clusters)
# visp_clusters <- setdiff(1:133, non_visp_clusters)
# 
# non_visp_cluster_labels <- unique(v1_alm_anno$cluster_label[v1_alm_anno$cluster_id %in% non_visp_clusters])
# 
# visp_dend <- dend %>%
#   prune.dendrogram(non_visp_cluster_labels)


load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/map.tree.df.rda")

v1_alm_cl_anno <- v1_alm_anno %>%
  select(cl, cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  mutate(cl = as.character(cl))

labels(V1.dend) <- v1_alm_cl_anno$cluster_label[match(labels(V1.dend), v1_alm_cl_anno$cl)]

node_n <- map.tree.df %>%
  mutate(cl = as.character(cl)) %>%
  group_by(cl) %>%
  summarise(n_cells = n()) %>%
  left_join(v1_alm_cl_anno) %>%
  mutate(cluster_label = ifelse(is.na(cluster_label), cl, cluster_label))

nodes <- as.ggdend(V1.dend)$nodes %>%
  mutate(cluster_label = get_nodes_attr(V1.dend,"label")) %>%
  left_join(node_n)

plot_nodes <- nodes %>%
  filter(!is.na(n_cells))

segments <- as.ggdend(V1.dend)$segments

node_labels <- plot_nodes %>%
  filter(y > 0)

leaf_labels <- nodes %>%
  filter(y == 0) %>%
  mutate(col = ifelse(is.na(n_cells),"#808080",col)) %>%
  left_join(v1_alm_cl_anno)

dendrogram_plot <- ggplot() +
  geom_segment(data = segments,
               aes(x = x, xend = xend,
                   y = y, yend = yend),
               lineend = "square",
               size = 0.5) +
  geom_point(data = plot_nodes,
             aes(x = x,
                 y = y,
                 size = n_cells),
             pch = 21,
             fill = "#00AEEF",
             color = "#2E3192",
             alpha = 1) +
  geom_text(data = node_labels,
            aes(x = x, y = y,
                label = cl),
            vjust = 1,
            size = 2*5/6) +
  geom_text(data = leaf_labels,
            aes(x = x, y = y - 0.01,
                label = cluster_label,
                color = col),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size_area(max_size = 2) +
  scale_y_continuous(limits = c(-0.1,1)) +
  theme_void()

ggsave("tolias_dendrogram_plot.pdf",dendrogram_plot, width = 8.75, height = 5, useDingbats = FALSE)

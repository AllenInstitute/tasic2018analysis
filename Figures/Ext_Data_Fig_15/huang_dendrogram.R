library(dplyr)
library(ggplot2)
library(dendextend)
options(stringsAsFactors = F)

source("prune_leaf_custom.R")

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

dend <- dend %>%
  prune(labels(dend)[c(1:55,116:133)])

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/external/mouse_GABA_Paul_2016_20171002/anno.feather")

node_n <- anno %>%
  group_by(cluster_label) %>%
  summarize(n_cells = n())

nodes <- as.ggdend(dend)$nodes %>%
  mutate(cluster_label = get_nodes_attr(dend,"label")) %>%
  left_join(node_n)

plot_nodes <- nodes %>%
  filter(!is.na(n_cells))

segments <- as.ggdend(dend)$segments

node_labels <- plot_nodes %>%
  filter(y > 0)

leaf_labels <- nodes %>%
  filter(y == 0) %>%
  mutate(col = ifelse(is.na(n_cells),"#808080",col))

dendrogram_plot <- ggplot() +
  geom_segment(data = segments,
               aes(x = x, xend = xend,
                   y = y, yend = yend)) +
  geom_point(data = plot_nodes,
             aes(x = x,
                 y = y,
                 size = n_cells),
             pch = 21,
             fill = "#00AEEF",
             color = "#2E3192",
             alpha = 0.7) +
  geom_text(data = node_labels,
            aes(x = x, y = y,
                label = cluster_label),
            vjust = 1) +
  geom_text(data = leaf_labels,
            aes(x = x, y = y - 0.01,
                label = cluster_label,
                color = col),
            angle = 90,
            hjust = 1,
            vjust = 0.3) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_size_area(breaks = c(10,25,50,100)) +
  scale_y_continuous(limits = c(-0.1,1)) +
  theme_void()

dendrogram_plot

ggsave("huang_dendrogram_plot.pdf",dendrogram_plot, width = 8, height = 6)

library(dplyr)
library(ggplot2)
library(feather)
library(scrattch)
library(Matrix)
options(stringsAsFactors = F)

# All cells
anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
anno <- anno %>% filter(cluster_id %in% 1:133)

downsample <- 100

anno_100 <- anno[0,]

cluster_ids <- unique(anno$cluster_id)
for(i in 1:length(cluster_ids)) {
  current_cluster_id <- cluster_ids[i]
  cluster_anno <- anno %>% filter(cluster_id == current_cluster_id)
  if(nrow(cluster_anno) > downsample) {
    cluster_anno <- cluster_anno %>%
      sample_n(downsample, replace = F)
  }
  anno_100 <- rbind(anno_100, cluster_anno)
}

anno_100 <- anno_100 %>%
  arrange(dendcluster_id)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/co.ratio.min.rda")

co.ratio.100 <- co.ratio.min[anno_100$sample_id,anno_100$sample_id]
co.ratio.100 <- as.matrix(co.ratio.100)

all.sample_ids <- rownames(co.ratio.100)
rownames(co.ratio.100) <- NULL
colnames(co.ratio.100) <- NULL

all_plot_data <- data.frame(x = numeric(),
                           y = numeric(),
                           val = numeric())

for(i in 1:ncol(co.ratio.100)) {
  
  col_data <- co.ratio.100[,i]
  
  cell_plot_data <- data.frame(x = i,
                               y = which(col_data > 0),
                               val = col_data[col_data != 0])
  
  all_plot_data <- rbind(all_plot_data,
                        cell_plot_data)
  
}

all_plot_data <- all_plot_data %>%
  mutate(fill = values_to_colors(val, minval = 0, colorset = c("white","black")))

all_tile_plot <- ggplot() +
  geom_tile(data = all_plot_data,
            aes(x = x, y = y,
                fill = fill)) +
  scale_fill_identity() +
  scale_y_reverse(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white"))

# margin of 0.0365 inches can't seem to be removed.
ggsave("all_coclustering_tile_100.png", all_tile_plot, width = nrow(anno_100)/1000 + 0.0365, height = nrow(anno_100)/1000 + 0.0365, dpi = 1000)


cocl_sample_id <- all.sample_ids

cocl_anno <- anno[match(cocl_sample_id, anno$sample_id),]

cocl_anno <- cocl_anno %>%
  mutate(cocl_pos = 1:n()) %>%
  group_by(cluster_id) %>%
  mutate(mid = mean(cocl_pos))

coarse_rects <- cocl_anno %>%
  group_by(subclass_label, subclass_color) %>%
  summarize(min_pos = min(cocl_pos),
            max_pos = max(cocl_pos))

all_tile_plot_clusters <- ggplot() +
  geom_tile(data = all_plot_data,
            aes(x = x, y = y,
                fill = fill)) +
  geom_rect(data = cocl_anno,
            aes(xmin = cocl_pos - 0.5,
                xmax = cocl_pos + 0.5,
                ymin = -100,
                ymax = 0.5,
                fill = cluster_color)) +
  geom_rect(data = cocl_anno,
            aes(ymin = cocl_pos - 0.5,
                ymax = cocl_pos + 0.5,
                xmin = max(cocl_anno$cocl_pos) + 0.5,
                xmax = max(cocl_anno$cocl_pos) + 100 + 0.5,
                fill = cluster_color)) + 
  geom_rect(data = coarse_rects,
            aes(xmin = min_pos, xmax = max_pos,
                ymin = min_pos, ymax = max_pos,
                color = subclass_color),
            fill = NA,
            size = 0.5) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_y_reverse(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white"))

ggsave("all_coclustering_tile_100_clusters_rects.png", all_tile_plot_clusters, width = nrow(anno_100)/1000 + 0.1 + 0.0365, height = nrow(anno_100)/1000 + 0.1 + 0.0365, dpi = 1000)

# Pvalb Zoom
pvalb_anno <- cocl_anno %>%
  filter(subclass_label == "Pvalb")

minpos <- min(pvalb_anno$cocl_pos)

pvalb_plot_data <- all_plot_data %>%
  filter(x %in% pvalb_anno$cocl_pos,
         y %in% pvalb_anno$cocl_pos) %>%
  mutate(x = x - minpos,
         y = y - minpos)

pvalb_anno <- pvalb_anno %>%
  mutate(cocl_pos = cocl_pos - minpos)

pvalb_tile_plot_clusters <- ggplot() +
  geom_tile(data = pvalb_plot_data,
            aes(x = x, y = y,
                fill = fill)) +
  geom_rect(data = pvalb_anno,
            aes(xmin = cocl_pos - 0.5,
                xmax = cocl_pos + 0.5,
                ymin = -25,
                ymax = 0.5,
                fill = cluster_color)) +
  geom_rect(data = pvalb_anno,
            aes(ymin = cocl_pos - 0.5,
                ymax = cocl_pos + 0.5,
                xmin = max(pvalb_anno$cocl_pos) + 0.5,
                xmax = max(pvalb_anno$cocl_pos) + 25 + 0.5,
                fill = cluster_color)) + 
  # geom_rect(data = coarse_rects,
  #           aes(xmin = min_pos, xmax = max_pos,
  #               ymin = min_pos, ymax = max_pos,
  #               color = coarse_color),
  #           fill = NA,
  #           size = 0.2) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_y_reverse(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_void() +
  theme(panel.background = element_rect(fill = "white"))

pvalb_tile_plot_clusters

ggsave("pvalb_coclustering_tile_100_clusters.png", pvalb_tile_plot_clusters, width = nrow(pvalb_anno)/1000 + 0.025 + 0.0365, height = nrow(pvalb_anno)/1000 + 0.025 + 0.0365, dpi = 1000)

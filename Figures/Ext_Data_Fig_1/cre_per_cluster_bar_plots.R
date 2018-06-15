library(dendextend)
library(ggplot2)
library(dplyr)
library(pvclust)
library(feather)
options(stringsAsFactors = F)

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

# Retrieve annotations feather to get n per cluster
anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

# cluster_id order doesn't match dendrogram order
anno <- anno %>%
  filter(cluster_id %in% 1:133)

n_clusters <- max(anno$cluster_id)

# convert to ggdend
dend_gg <- as.ggdend(dend)

dend_seg <- dend_gg$segments

dend_leaves <- dend_gg$labels %>%
  mutate(cluster_label = label)

cluster_anno <- anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  left_join(dend_leaves) %>%
  mutate(cluster_id = x)

anno <- anno %>%
  select(-cluster_id) %>%
  left_join(cluster_anno)

panel_pad <- 0.05

sections <- anno %>%
  group_by(coarse_id) %>%
  summarise(x = min(dendcluster_id - 0.5),
            xend = n_clusters + 0.5)

panel_segments <- data.frame(x = min(sections$x), xend = max(sections$xend),
                             y = c(-1*(0:3 + 1:4*panel_pad),
                                            -1*(1:3 + 1:3*panel_pad))) %>%
  mutate(yend = y)

wedge_lines <- data.frame(x = unique(c(sections$x, sections$xend)),
                          y = 0,
                          yend = -4.2) %>%
  mutate(xend = x)

# Region rectangles
av_rects <- anno %>%
  select(cluster_id, cluster_label, cluster_color, region_id, region_label, region_color) %>%
  group_by(cluster_id, region_id, region_label) %>%
  mutate(ly_n = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(region_id) %>%
  mutate(cluster_n = n(),
         ly_frac = ly_n/cluster_n) %>%
  unique() %>%
  arrange(region_id) %>%
  mutate(ly_cum_frac = cumsum(ly_frac)) %>%
  ungroup() %>%
  arrange(cluster_id, region_id) %>%
  group_by(cluster_id) %>%
  mutate(xmin = cluster_id - 0.5,
         xmax = cluster_id + 0.5,
         ymax = - lag(ly_cum_frac, default = 0) - panel_pad,
         ymin = - ly_cum_frac - panel_pad)

# Layer rectangles
layer_rects <- anno %>%
  select(cluster_id, cluster_label, cluster_color, layer_id, layer_label, layer_color) %>%
  group_by(cluster_id, layer_id, layer_label) %>%
  mutate(ly_n = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(layer_id) %>%
  mutate(cluster_n = n(),
         ly_frac = ly_n/cluster_n) %>%
  unique() %>%
  arrange(layer_id) %>%
  mutate(ly_cum_frac = cumsum(ly_frac)) %>%
  ungroup() %>%
  arrange(cluster_id, layer_id) %>%
  group_by(cluster_id) %>%
  mutate(xmin = cluster_id - 0.5,
         xmax = cluster_id + 0.5,
         ymax = -1 - lag(ly_cum_frac, default = 0) - panel_pad * 2,
         ymin = -1 - ly_cum_frac - panel_pad * 2)

# cre line rects
cre_rects <- anno %>%
  mutate(cre_id = ifelse(inj_type_label == "retrograde",0,cre_id),
         cre_label = ifelse(inj_type_label == "retrograde", "Retrograde", cre_label),
         cre_color = ifelse(inj_type_label == "retrograde", "#000000", cre_color)) %>%
  # mutate(cre_id = ifelse(cre_label == "0", cre_id + inj_roi_id/10, cre_id),
  #        cre_color = ifelse(cre_label == "0", inj_roi_color, cre_color),
  #        cre_label = ifelse(cre_label == "0", inj_roi_label, cre_label)) %>%
  select(cluster_id, cluster_label, cluster_color, cre_id, cre_label, cre_color) %>%
  group_by(cluster_id, cre_id, cre_label) %>%
  mutate(ly_n = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(cre_id) %>%
  mutate(cluster_n = n(),
         ly_frac = ly_n/cluster_n) %>%
  unique() %>%
  arrange(cre_id) %>%
  mutate(ly_cum_frac = cumsum(ly_frac)) %>%
  ungroup() %>%
  arrange(cluster_id, cre_id) %>%
  group_by(cluster_id) %>%
  mutate(xmin = cluster_id - 0.5,
         xmax = cluster_id + 0.5,
         ymax = -2 - lag(ly_cum_frac, default = 0) - panel_pad * 3,
         ymin = -2 - ly_cum_frac - panel_pad * 3)

n_rects <- anno  %>%
  group_by(cluster_id, cluster_color, cluster_label) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(adj_n = log10(n)) %>%
  mutate(f = adj_n/3) %>%
  mutate(xmin = cluster_id - 0.5,
         xmax = cluster_id + 0.5,
         ymin = -3 - f - panel_pad * 4,
         ymax = -3 - panel_pad * 4)

n_guides <- data.frame(y = seq(-4 - panel_pad * 4,-3 - panel_pad * 4,by = 1/6),
                       x = 0.5,
                       xend = n_clusters + 1,
                       label = seq(3, 0, by = -0.5)) %>%
  mutate(yend = y)
# 
# large_cluster_break <- data.frame(y = -4.05 - panel_pad * 4, yend = -4.05 - panel_pad * 4,
#                                   x = 0.5, xend = n_clusters + 1)
# 
# large_cluster_labels <- n_rects %>%
#   filter(n > 500) %>%
#   mutate(y = -4.2 - panel_pad * 4,
#          x = (xmin+xmax)/2)

region_frac <- anno %>%
  group_by(cluster_id) %>%
  summarise(v1_f = sum(region_label == "VISp")/n()) %>%
  mutate(color = "#FF0000",
         color = ifelse(v1_f == 1, "#000000",color),
         color = ifelse(v1_f == 0, "#FFFFFF",color),
         x = cluster_id,
         y = 0)

# Flat version

flat_plot <- ggplot() +
  geom_rect(data = cre_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = cre_color)) +
  # Leaf Labels
  geom_text(data = dend_leaves,
            aes(x = x,
                y = -3.2,
                label = label,
                color = col),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2) +
  #Vertical class separators
  # geom_segment(data = wedge_lines,
  #              aes(x = x,
  #                  xend = xend,
  #                  y = y,
  #                  yend = yend)) +
  # Borders between panels
  scale_size(range = c(0.5, 1)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0.25,0)) +
  scale_x_continuous(limits = c(-1,n_clusters + 1)) +
  theme_void()
flat_plot

ggsave("cre_per_cluster.pdf",
       flat_plot, 
       width = 16, height = 12, 
       useDingbats = F)

cre_anno <- anno %>%
  select(cre_id, cre_label, cre_color) %>%
  unique()

cre_legend <- ggplot() +
  geom_tile(data = cre_anno,
            aes(x = 1,
                y = cre_id,
                fill = cre_color)) +
  geom_text(data = cre_anno,
            aes(x = 1.51,
                y = cre_id,
                label = cre_label),
            hjust = 0,
            vjust = 0.3) +
  scale_fill_identity() +
  scale_x_continuous(limits = c(0.5, 3)) +
  scale_y_reverse() +
  theme_void()

cre_legend

ggsave("cre_legend.pdf",
       cre_legend,
       width = 3,
       height = 8)

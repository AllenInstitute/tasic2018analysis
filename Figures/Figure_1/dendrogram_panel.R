library(dendextend)
library(ggplot2)
library(dplyr)
library(pvclust)
library(feather)
options(stringsAsFactors = F)

# Load the dendrogram from the feather directory
dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.Rdata")

# 20170927: Use the dendrogram from Osnat with weights based on bootstrapping confidence
# dend <- dend %>%
#   pvclust_show_signif_gradient(res, signif_type = "bp", signif_col_fun=colorRampPalette(c("white","black"))) %>%
#   pvclust_show_signif(res, signif_type="bp", signif_value=c(2,1))
# labels(dend) <- labels(dend_feather)

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
  group_by(subclass_id) %>%
  summarise(x = min(dendcluster_id - 0.5),
            xend = n_clusters + 0.5)

panel_segments <- data.frame(x = min(sections$x), xend = max(sections$xend),
                             y = c(-1*panel_pad,     -1*panel_pad - 1,
                                   -2*panel_pad - 1, -2*panel_pad - 3,
                                   -3*panel_pad - 3)) %>%
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
         ymax = -1 - lag(ly_cum_frac, default = 0)*2 - panel_pad * 2,
         ymin = -1 - ly_cum_frac*2 - panel_pad * 2)

n_rects <- anno  %>%
  group_by(cluster_id, cluster_color, cluster_label) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  mutate(adj_n = ifelse(n <= 500, n, 500)) %>%
  mutate(f = adj_n/500) %>%
  mutate(xmin = cluster_id - 0.5,
         xmax = cluster_id + 0.5,
         ymin = -3 - f - panel_pad * 3,
         ymax = -3 - panel_pad * 3)

n_guides <- data.frame(y = seq(-4 - panel_pad * 3,-3 - panel_pad * 3,by = 1/5),
                       x = 0.5,
                       xend = n_clusters + 1,
                       label = seq(500, 0, by = -100)) %>%
  mutate(yend = y)

large_cluster_rects <- n_rects %>%
  filter(n > 500) %>%
  mutate(ymin = -4.2 - panel_pad*3,
         ymax = -4.05 - panel_pad*3)

large_cluster_labels <- n_rects %>%
  filter(n > 500) %>%
  mutate(y = -4.05 - panel_pad * 3,
         x = (xmin+xmax)/2)

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
  geom_segment(data = dend_seg,
               aes(x = x,
                   xend = xend,
                   y = y,
                   yend = yend,
                   size = lwd,
                   color = col),
               lineend = "square") +
  # Annotation panels
  geom_rect(data = av_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = region_color)) +
  geom_rect(data = layer_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = layer_color)) +
  geom_rect(data = n_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = cluster_color)) +
  # Large cluster tabs
  geom_rect(data = large_cluster_rects,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = cluster_color)) +
  # Large cluster number labels
  geom_text(data = large_cluster_labels,
            aes(x = x,
                y = y,
                label = n),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2) +
  # Dashed lines for N cells
  geom_segment(data = n_guides,
               aes(x = x,
                   xend = xend,
                   y = y,
                   yend = yend),
               linetype = "dashed") +
  # N Cells labels
  geom_text(data = n_guides,
            aes(x = 0,
                y = y,
                label = label),
            size = 2,
            hjust = 1) +
  # Leaf Labels
  geom_text(data = dend_leaves,
            aes(x = x,
                y = -4.4,
                label = label,
                color = col),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2) +
  # Borders between panels
  geom_segment(data = panel_segments,
               aes(x = x, xend = xend,
                   y = y, yend = yend)) +
  scale_size(range = c(0.5, 1)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0.25,0)) +
  scale_x_continuous(limits = c(-1,n_clusters + 1)) +
  theme_void()
flat_plot

ggsave("hierarchical_cluster_flat_plot_bars_t.pdf",
       flat_plot, 
       width = 16, height = 12, 
       useDingbats = F)


# Rotation to match page
flat_plot_r <- ggplot() +
  geom_segment(data = dend_seg,
               aes(x = -y,
                   xend = -yend,
                   y = -x,
                   yend = -xend,
               size = lwd,
               color = col),
               lineend = "square") +
  # Annotation panels
  geom_rect(data = av_rects,
            aes(xmin = -ymin,
                xmax = -ymax,
                ymin = -xmin,
                ymax = -xmax,
                fill = region_color)) +
  geom_rect(data = layer_rects,
            aes(xmin = -ymin,
                xmax = -ymax,
                ymin = -xmin,
                ymax = -xmax,
                fill = layer_color)) +
  geom_rect(data = n_rects,
            aes(xmin = -ymin,
                xmax = -ymax,
                ymin = -xmin,
                ymax = -xmax,
                fill = cluster_color)) +
  # Large cluster tabs
  geom_rect(data = large_cluster_rects,
            aes(xmin = -ymin,
                xmax = -ymax,
                ymin = -xmin,
                ymax = -xmax,
                fill = cluster_color)) +
  # Large cluster number labels
  geom_text(data = large_cluster_labels,
            aes(x = -y,
                y = -x,
                label = n),
            hjust = 0,
            vjust = 0.3,
            size = 2) +
  # Dashed lines for N cells
  geom_segment(data = n_guides,
               aes(x = -y,
                   xend = -yend,
                   y = -x,
                   yend = -xend),
               linetype = "dashed") +
  # N Cells labels
  geom_text(data = n_guides,
            aes(x = -y,
                y = 0,
                label = label),
            size = 2,
            hjust = 0.5,
            vjust = 0) +
  # Leaf Labels
  geom_text(data = dend_leaves,
            aes(x = 4.4,
                y = -x,
                label = label,
                color = col),
            angle = 0,
            hjust = 0,
            vjust = 0.3,
            size = 2) +
  #Vertical class separators
  geom_segment(data = wedge_lines,
               aes(x = -y,
                   xend = -yend,
                   y = -x,
                   yend = -xend)) +
  #Borders between panels
  geom_segment(data = panel_segments,
               aes(x = -y, xend = -yend,
                   y = -x, yend = -xend)) +
  scale_size(range = c(0.5, 1)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(expand = c(0.25,0)) +
  scale_y_continuous(limits = c(-n_clusters - 2,2)) +
  theme_void()
flat_plot_r

ggsave("hierarchical_cluster_flat_plot_bars.pdf",
       flat_plot_r, 
       width = 7.5, height = 10.5, 
       useDingbats = F)



dend_nodes <- dend_gg$nodes %>%
  mutate(label = get_nodes_attr(dend,"label")) %>%
  filter(y > 0)

nodes_plot <- ggplot() +
  geom_segment(data = dend_seg,
               aes(x = x,
                   xend = xend,
                   y = y,
                   yend = yend)) +
  geom_text(data = dend_leaves,
            aes(x = x,
                y = -0.01,
                label = paste(label, cluster_id),
                color = col),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2) +
  geom_text(data = dend_nodes,
            aes(x = x + 1,
                y = y + 0.02,
                label = label),
            size = 2) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0.25,0)) +
  scale_x_continuous(limits = c(-1,n_clusters + 1)) +
  theme_void()
nodes_plot

ggsave("hierarchical_cluster_node_labels.pdf",nodes_plot,height = 8, width = 12, useDingbats = F)




cre_anno <- anno %>%
  select(cre_id, cre_label, cre_color) %>%
  unique()

cre_legend <- ggplot() +
  geom_tile(data = cre_anno,
            aes(x = cre_id, y = 1,
                fill = cre_color)) +
  scale_fill_identity() +
  scale_x_continuous(breaks = cre_anno$cre_id, labels = cre_anno$cre_label) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust= 1, vjust = 0.3))
ggsave("cre_legend.pdf",cre_legend,height = 6, width = 12, useDingbats = F)


library(rbokeh)

b_av_rects <- av_rects
b_layer_rects <- layer_rects
names(b_layer_rects) <- names(b_av_rects)
b_cre_rects <- cre_rects
names(b_cre_rects) <- names(b_av_rects)

all_rects <- rbind(b_av_rects, b_layer_rects, b_cre_rects)

bokeh_plot <- figure(height = 800,
                     width = 1200,
                     padding_factor = 0,
                     xaxes = FALSE,
                     yaxes = FALSE,
                     xlab = "",
                     ylab = "") %>%
  ly_segments(data = dend_seg,
              x0 = x, x1 = xend,
              y0 = y, y1 = yend) %>%
  ly_rect(data = all_rects,
          xleft = xmin, xright = xmax,
          ybottom = ymin, ytop = ymax,
          color = region_color,
          fill_alpha = 1,
          hover = list("Cluster" = cluster_label,
                       "Annotation" = region_label,
                       "N Cells" = ly_n,
                       "Fraction" = round(ly_frac,2))) %>%
  ly_rect(data = n_rects,
          xleft = xmin, xright = xmax,
          ybottom = ymin, ytop = ymax,
          color = cluster_color,
          fill_alpha = 1,
          hover = list("Cluster" = cluster_label,
                       "N Cells" = n))
  
bokeh_plot

htmlwidgets::saveWidget(bokeh_plot, file = "all_cell_dendrogram.html")


# 
# flat_plot_av_ratio <- ggplot() +
#   geom_segment(data = dend_seg,
#                aes(x = x,
#                    xend = xend,
#                    y = y,
#                    yend = yend)) +
#   geom_segment(data = sections,
#                aes(x = x,
#                    xend = xend,
#                    y = y,
#                    yend = yend)) +
#   geom_rect(data = av_rects,
#             aes(xmin = xmin,
#                 xmax = xmax,
#                 ymin = ymin - 0.05,
#                 ymax = ymax - 0.05,
#                 fill = region_color)) +
#   geom_segment(data = wedge_lines,
#                aes(x = x,
#                    xend = xend,
#                    y = y,
#                    yend = yend - 0.3)) +
#   geom_text(data = dend_leaves,
#             aes(x = x,
#                 y = -1.07,
#                 label = label,
#                 color = col),
#             angle = 90,
#             hjust = 1,
#             vjust = 0.3,
#             size = 2) +
#   geom_point(data = region_frac,
#              aes(x = x, y = y, color = "#000000", fill = color),
#              size = 2,pch=21) +
#   scale_color_identity() +
#   scale_fill_identity() +
#   scale_y_continuous(expand = c(0.25,0)) +
#   scale_x_continuous(limits = c(-1,n_clusters + 1)) +
#   theme_void()
# flat_plot_av_ratio
# ggsave("hierarchical_cluster_flat_plot_ratio.pdf",flat_plot_av_ratio,height = 8, width = 12, useDingbats = F)

# 
# layer_heat <- anno %>%
#   left_join(layer_order) %>%
#   select(cluster_id, cluster_label, cluster_color, layer_order, layer_label, layer_color) %>%
#   group_by(cluster_id, layer_order, layer_label) %>%
#   mutate(ly_n = n()) %>%
#   ungroup() %>%
#   group_by(cluster_id) %>%
#   arrange(layer_order) %>%
#   mutate(cluster_n = n(),
#          ly_frac = ly_n/cluster_n) %>%
#   unique() %>%
#   arrange(layer_order) %>%
#   ungroup() %>%
#   arrange(cluster_id, layer_order) %>%
#   group_by(cluster_id) %>%
#   mutate(xmin = cluster_id - 0.5,
#          xmax = cluster_id + 0.5,
#          ymax = -1 - (layer_order - 1)/12,
#          ymin = -1 - layer_order/12,
#          ly_frac_color = values_to_colors(ly_frac, max = 1))
# 
# ggplot() +
#   geom_rect(data = av_rects,
#             aes(xmin = xmin, xmax = xmax,
#                 ymin = ymin, ymax = ymax,
#                 fill = region_color)) +
#   geom_rect(data = layer_heat,
#             aes(xmin = xmin, xmax = xmax,
#                 ymin = ymin, ymax = ymax,
#                 fill = ly_frac_color)) +
#   geom_text(data = av_rects,
#             aes(x = cluster_id,
#                 y = -2,
#                 label = cluster_label,
#                 color = cluster_color),
#             angle = 90, hjust = 1, vjust = 0.3) +
#   scale_fill_identity() +
#   scale_color_identity() +
#   scale_y_continuous(limits = c(-3,0)) +
#   theme_void()
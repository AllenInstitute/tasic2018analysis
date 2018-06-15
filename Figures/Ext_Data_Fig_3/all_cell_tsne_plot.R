library(dplyr)
library(ggplot2)
library(feather)
library(ggrepel)
library(purrr)
library(gridExtra)
options(stringsAsFactors = F)

# Load tSNE coordinates
tsne <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/tsne.feather")

# Load feather annotations for cluster colors
# Filter for ALM and VISp, and remove injections.
anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

point_colors <- anno %>%
  select(sample_id,cluster_color,cluster_label,class_color,region_color,genes_color,confusion_color,sex_color)

tsne_data <- tsne %>%
  filter(sample_id %in% anno$sample_id) %>%
  left_join(point_colors)

cluster_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = all_x,
                 y = all_y,
                 color = cluster_color)) +
  scale_color_identity() +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_void()

region_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = all_x,
                 y = all_y,
                 color = region_color)) +
  scale_color_identity() +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_void()

class_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = all_x,
                 y = all_y,
                 color = class_color)) +
  scale_color_identity() +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_void()

confusion_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = all_x,
                 y = all_y,
                 color = confusion_color)) +
  scale_color_identity() +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_void()

genes_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = all_x,
                 y = all_y,
                 color = genes_color)) +
  scale_color_identity() +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_void()

sex_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = all_x,
                 y = all_y,
                 color = sex_color)) +
  scale_color_identity() +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_void() 

ggsave("cluster_tsne.png",cluster_plot,width = 12, height = 12)
ggsave("region_tsne.png",region_plot,width = 12, height = 12)
ggsave("class_tsne.png",class_plot,width = 12, height = 12)
ggsave("confusion_tsne.png",confusion_plot,width = 12, height = 12)
ggsave("genes_tsne.png",genes_plot,width = 12, height = 12)
ggsave("sex_tsne.png",sex_plot,width = 12, height = 12)


outline_color <- function(x) {
  library(grDevices)
  
  rgb_x <- col2rgb(x)
  hsv_x <- rgb2hsv(rgb_x[1],
                   rgb_x[2],
                   rgb_x[3])
  new_hsv <- as.matrix(hsv_x)
  new_hsv[3] <- new_hsv[3]/2
  
  new_hex <- hsv(new_hsv[1],
                 new_hsv[2],
                 new_hsv[3])
  
  new_hex
}


# tSNE plots that highlight position of every cluster
tsne_highlight_plot <- function(select_id, highlight_group = "cluster", color_by = "cluster", tsne, anno) {
  print(select_id)
  highlight_id <- paste0(highlight_group,"_id")
  highlight_label <- paste0(highlight_group,"_label")
  color_column <- paste0(color_by,"_color")
  
  plot_label <- anno[[highlight_label]][anno[[highlight_id]] == select_id] %>% unique() %>% unlist()
  
  point_colors <- anno %>%
    select(one_of("sample_id", highlight_id, color_column))
  
  tsne_data <- tsne %>%
    filter(sample_id %in% anno$sample_id) %>%
    left_join(point_colors) %>%
    rename_("color" = color_column)
  
  fore <- tsne_data %>%
    filter_(paste0(highlight_id," == ",select_id)) %>%
    rowwise() %>%
    mutate(fill = color,
           color = outline_color(fill))
  
  back <- tsne_data  %>%
    filter_(paste0(highlight_id," != ",select_id)) %>%
    mutate(color = "#D1D3D4")
  
  p <- ggplot() +
    geom_point(data = back,
               aes(x = all_x,
                   y = all_y,
                   color = color),
               size = 0.5) +
    geom_point(data = fore,
               aes(x = all_x,
                   y = all_y,
                   fill = fill,
                   color = color),
               size = 0.5,
               pch = 21) +
    scale_color_identity() +
    scale_fill_identity() +
    scale_x_continuous("") +
    scale_y_continuous("") +
    theme_void() +
    ggtitle(plot_label) +
    theme(legend.position = "none")
  
  p
}

clusters <- anno %>%
  filter(cluster_id < 118) %>%
  select(cluster_id) %>%
  unique() %>%
  arrange(cluster_id) %>%
  unlist() %>%
  as.list()

cluster_highlight_plots <- map(clusters, tsne_highlight_plot, "cluster", "cluster" ,tsne, anno)

ggsave("tSNE_by_cluster_huge.png", plot = arrangeGrob(grobs = cluster_highlight_plots,
                                                 nrow = 8), width = 75, height = 40,
       limitsize = FALSE)

ggsave("tSNE_by_cluster.png", plot = arrangeGrob(grobs = cluster_highlight_plots,
                                                      nrow = 8), width = 37.5, height = 20,
       limitsize = FALSE)

cluster_confusion_plots <- map(clusters, tsne_highlight_plot, "cluster", "confusion" ,tsne, anno)

ggsave("tSNE_by_cluster_confusion.png", plot = arrangeGrob(grobs = cluster_confusion_plots,
                                                 nrow = 8), width = 37.5, height = 20)

anno_no_inj <- anno %>%
  filter(inj_type_label == "0") %>%
  arrange(cre_label)

cre_highlight_plots <- map(unique(anno_no_inj$cre_id), tsne_highlight_plot, "cre","cre", tsne, anno_no_inj)

ggsave("tsne_by_cre.png", plot = arrangeGrob(grobs = cre_highlight_plots,
                                             nrow = 6), width = 20, height = 20)

# tsne centroids
centroids <- tsne_data %>%
  group_by(label, color,id) %>%
  summarise(n = n(),
            x = mean(x),
            y = mean(y)) 

cent_plot <- ggplot() +
  geom_point(data = centroids,
             aes(x = x,
                 y = y,
                 color = "black",
                 fill = color,
                 size = n),
             pch = 21) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_area(max_size = 20) +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_bw()

ggsave("all_cell_centroid_disks.pdf", cent_plot, width = 12, height = 12, useDingbats = F)

tsne_cent_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = x,
                 y = y,
                 color = color), alpha = 0.5) +
  geom_point(data = centroids,
             aes(x = x,
                 y = y,
                 color = "black",
                 fill = color,
                 size = n),
             pch = 21) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_area(max_size = 20) +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_bw()

ggsave("all_cell_tsne_centroid_disks.pdf", tsne_cent_plot, width = 12, height = 12, useDingbats = F)

tsne_label_plot <- ggplot() +
  geom_point(data = tsne_data,
             aes(x = x,
                 y = y,
                 color = color)) +
  geom_text_repel(data = centroids,
             aes(x = x,
                 y = y,
                 label = label),
             fontface = "bold",
             size = 3) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_size_area(max_size = 20) +
  scale_x_continuous("") +
  scale_y_continuous("") +
  theme_bw()

ggsave("all_cell_tsne_labels.pdf", tsne_label_plot, width = 12, height = 12, useDingbats = F)


# tsne hulls
# not terribly useful.
hulls <- tsne_data[0]
groups <- unique(tsne_data$id)

for(i in 1:length(groups)) {
  group_tsne <- tsne_data %>%
    filter(id == groups[i])
  group_hull <- chull(group_tsne$x, group_tsne$y)
  hulls <- rbind(hulls, group_tsne[group_hull,])
}

tsne_hull_plot <- ggplot() +
  geom_polygon(data = hulls,
               aes(x = x, y = y, group = id, color = color),
               fill = NA) +
  scale_color_identity()

library(dplyr)
library(scrattch)
library(ggplot2)
library(feather)
library(reshape2)
options(stringsAsFactors = F)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/cell.marker.dat.rda")

anno_file <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather"

anno <- read_feather(anno_file)

cluster_anno <- anno %>%
  select(sample_id, dendcluster_id, dendcluster_color)

n_cells <- ncol(cell.marker.dat)
n_genes <- nrow(cell.marker.dat)

# remove combinatorial genes. ~ rows > 1110
#cell.marker.dat <- cell.marker.dat[c(966:1200,1:965),]
# shift pan markers. ~rows > 1040

gene_pos <- data.frame(gene = rownames(cell.marker.dat)) %>%
  mutate(ypos = n():1)

sample_pos <- data.frame(sample_id = colnames(cell.marker.dat)) %>%
  mutate(xpos = 1:n())

data_melt <- melt(cell.marker.dat)
names(data_melt)[1:2] <- c("gene","sample_id")

cluster_rects <- data.frame(sample_id = colnames(cell.marker.dat)) %>%
  mutate(xpos = 1:n()) %>%
  mutate(ymin = nrow(cell.marker.dat) + 5, ymax = nrow(cell.marker.dat) + 20) %>%
  left_join(cluster_anno)

plot_data <- data_melt %>%
  left_join(gene_pos) %>%
  left_join(sample_pos) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(fill = values_to_colors(value, colorset = c("white","black")))

all_tile_plot <- ggplot() +
  geom_tile(data = plot_data,
            aes(x = xpos, y = ypos,
                fill = fill)) +
  geom_rect(data = cluster_rects,
            aes(xmin = xpos - 0.5, xmax = xpos + 0.5,
                ymin = ymin, ymax = ymax,
                fill = dendcluster_color)) +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_void()

ggsave("zy_100sampled_truth_bw_rect.png", all_tile_plot, height = (n_genes*3)/600 + 20/600 + 0.0365, width = n_cells/600 + 0.0365, dpi = 600)

# Rearrange to match dendrogram
cell.marker.dat2 <- cell.marker.dat[c(1:240,
                                     436:780,
                                     241:435,
                                     781:1000),]


gene_pos <- data.frame(gene = rownames(cell.marker.dat2)) %>%
  mutate(ypos = n():1)

sample_pos <- data.frame(sample_id = colnames(cell.marker.dat2)) %>%
  left_join(cluster_anno) %>%
  arrange(dendcluster_id) %>%
  mutate(xpos = 1:n())

cell.marker.dat2 <- cell.marker.dat2[,sample_pos$sample_id]

data_melt <- melt(cell.marker.dat2)
names(data_melt)[1:2] <- c("gene","sample_id")

cluster_rects <- data.frame(sample_id = colnames(cell.marker.dat2)) %>%
  mutate(xpos = 1:n()) %>%
  mutate(ymin = nrow(cell.marker.dat2) + 5, ymax = nrow(cell.marker.dat2) + 20) %>%
  left_join(cluster_anno)

plot_data <- data_melt %>%
  left_join(gene_pos) %>%
  left_join(sample_pos) %>%
  mutate(value = ifelse(is.na(value), 0, value)) %>%
  mutate(fill = values_to_colors(value, colorset = c("white","black")))

all_tile_plot <- ggplot() +
  geom_tile(data = plot_data,
            aes(x = xpos, y = ypos,
                fill = fill)) +
  geom_rect(data = cluster_rects,
            aes(xmin = xpos - 0.5, xmax = xpos + 0.5,
                ymin = ymin, ymax = ymax,
                fill = dendcluster_color)) +
  scale_fill_identity() +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_void()

ggsave("dend_order_zy_100sampled_truth_bw_rect.png", all_tile_plot, height = (n_genes*3)/600 + 20/600 + 0.0365, width = n_cells/600 + 0.0365, dpi = 600)

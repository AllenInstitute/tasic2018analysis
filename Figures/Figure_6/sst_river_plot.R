library(ggplot2)
library(dplyr)
library(scrattch.io)
library(feather)
options(stringsAsFactors = F)

source("custom_annotate_cat.R")
source("sankey_functions.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

sst_cluster_anno <- anno %>%
  filter(subclass_label == "Sst") %>%
  select(cl, cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  arrange(cluster_id) %>%
  mutate(cluster_id = 1:n())

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.cell.cl.map.80.rda")
sst80 <- cell.cl.map.df %>%
  mutate(sample_id = rownames(.))

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.cell.cl.map.300.rda")
sst300 <- cell.cl.map.df %>%
  mutate(sample_id = rownames(.))

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/cell.cl.map.df.rda")
sst150 <- cell.cl.map.df[sst80$sample_id,] %>%
  mutate(sample_id = rownames(.)) %>%
  mutate(org.cl = as.numeric(as.character(org.cl))) %>%
  left_join(sst_cluster_anno, by = c("org.cl" = "cl"))

sst80_cl <- sst80 %>%
  select(sample_id, org.cl)
names(sst80_cl)[2] <- "sst80_cl"
sst150_cl <- sst150 %>%
  select(sample_id, org.cl, cluster_id, cluster_label, cluster_color)
names(sst150_cl)[2] <- "sst150_cl"
sst300_cl <- sst300 %>%
  select(sample_id, org.cl)
names(sst300_cl)[2] <- "sst300_cl"


river_data <- sst80_cl %>%
  left_join(sst150_cl) %>%
  left_join(sst300_cl) 

# Filter transitions of fewer than 3 cells
river_data2 <- river_data %>%
  group_by(sst80_cl, sst150_cl) %>%
  mutate(n_80_150 = n()) %>%
  group_by(sst150_cl, sst300_cl) %>%
  mutate(n_150_300 = n()) %>%
  ungroup() %>%
  filter(n_80_150 > 2 & n_150_300 > 2)

max_freq_name <- function(x) {
  names(table(x))[which(table(x) == max(table(x)))]
}

# Ordering:
# First, order 300 based on clusters
sst300_cl_anno <- river_data2 %>%
  group_by(sst300_cl) %>%
  summarise(sst300_cluster_weight = mean(cluster_id),
            sst300_cluster_color = color_mean(cluster_color),
            sst300_cluster_label = max_freq_name(cluster_label)) %>%
  arrange(sst300_cluster_weight) %>%
  mutate(sst300_cluster_id = 1:n())

river_data2 <- river_data2 %>%
  left_join(sst300_cl_anno)

# Then, order standard clustering based on 300
std_cl_anno <- river_data2 %>%
  group_by(cluster_id) %>%
  summarise(std_cluster_weight = mean(sst300_cluster_id),
            std_cluster_color = color_mean(sst300_cluster_color),
            std_cluster_label = cluster_label[1]) %>%
  arrange(std_cluster_weight) %>%
  mutate(std_cluster_id = 1:n())

river_data2 <- river_data2 %>%
  left_join(std_cl_anno)

sst80_cl_anno <- river_data2 %>%
  group_by(sst80_cl) %>%
  summarise(sst80_cluster_weight = mean(std_cluster_id),
            sst80_cluster_color = color_mean(std_cluster_color),
            sst80_cluster_label = max_freq_name(std_cluster_label)) %>%
  arrange(sst80_cluster_weight) %>%
  mutate(sst80_cluster_id = 1:n())

river_data2 <- river_data2 %>%
  left_join(sst80_cl_anno)

plot_nodes <- make_plot_nodes(make_group_nodes(river_data2, c("sst80_cluster","std_cluster","sst300_cluster")),pad = 0.2)

group_links <- make_group_links(river_data2, 
                                c("sst80_cluster","std_cluster","sst300_cluster"),
                                plot_nodes) %>%
  filter(n > 4)

plot_links <- make_plot_links(group_links,
                              fill = "std_cluster")

sst80_labels <- plot_nodes %>%
  filter(group == "sst80_cluster") %>%
  mutate(x = 1,
         y = (ymin + ymax) / 2)
std_labels <- plot_nodes %>%
  filter(group == "std_cluster") %>%
  mutate(x = 2,
         y = (ymin + ymax) / 2)
sst300_labels <- plot_nodes %>%
  filter(group == "sst300_cluster") %>%
  mutate(x = 3,
         y = (ymin + ymax) / 2)

river_plot <- build_river_plot_predefined(plot_nodes, plot_links) +
  geom_text(data = sst80_labels,
            aes(x = x, y = y,
                label = name),
            size = 2) +
  geom_text(data = std_labels,
            aes(x = x, y = y,
                label = name),
            size = 2) +
  geom_text(data = sst300_labels,
            aes(x = x, y = y,
                label = name),
            size = 2)

ggsave("sst_river_plot.pdf",
       river_plot,
       width = 6,
       height = 3.5)

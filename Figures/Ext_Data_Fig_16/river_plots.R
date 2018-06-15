library(dplyr)
library(ggplot2)
library(feather)
options(stringsAsFactors = F)
source("sankey_functions.R")


anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170405/anno.feather")

anno <- anno %>%
  filter(genotype_label == "Vipr2-IRES2-Cre/wt;Ai14(RCL-tdT)/wt")

# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
nodes <- make_plot_nodes(make_group_nodes(anno, c("genotype","cluster")))

left_nodes <- nodes %>%
  filter(group == "genotype") %>%
  mutate(name = paste0(name, " (n = ",n,")"))
right_nodes <- nodes %>%
  filter(group == "cluster") %>%
  mutate(name = paste0(name, " (n = ",n,")"))


river <- build_river_plot(anno, 
                          c("genotype","cluster"), 
                          fill_group = "cluster") +
  geom_text(data = right_nodes,
            aes(x = xmax + 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  geom_text(data = left_nodes,
            aes(x = xmin - 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  scale_color_identity()

river

ggsave("vipr2_ai14_river.pdf", river, height = 5, width = 5)

## Simulations for placeholders

anno <- anno %>%
  filter(cluster_id %in% c(1,29:47))

nodes <- make_plot_nodes(make_group_nodes(anno, c("genotype","cluster")))

left_nodes <- nodes %>%
  filter(group == "genotype") %>%
  mutate(name = paste0(name, " (n = ",n,")"))
right_nodes <- nodes %>%
  filter(group == "cluster") %>%
  mutate(name = paste0(name, " (n = ",n,")"))


river <- build_river_plot(anno, 
                          c("genotype","cluster"), 
                          fill_group = "cluster") +
  geom_text(data = right_nodes,
            aes(x = xmax + 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  geom_text(data = left_nodes,
            aes(x = xmin - 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  scale_color_identity()

river

ggsave("SIMULATION_vipr2_slc32a1_river.pdf", river, height = 5, width = 5)


anno <- anno %>%
  filter(cluster_id %in% c(47))

nodes <- make_plot_nodes(make_group_nodes(anno, c("genotype","cluster")))

left_nodes <- nodes %>%
  filter(group == "genotype") %>%
  mutate(name = paste0(name, " (n = ",n,")"))
right_nodes <- nodes %>%
  filter(group == "cluster") %>%
  mutate(name = paste0(name, " (n = ",n,")"))


river <- build_river_plot(anno, 
                          c("genotype","cluster"), 
                          fill_group = "cluster") +
  geom_text(data = right_nodes,
            aes(x = xmax + 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  geom_text(data = left_nodes,
            aes(x = xmin - 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  scale_color_identity()

river

ggsave("SIMULATION_vipr2_pvalb_river.pdf", river, height = 5, width = 5)

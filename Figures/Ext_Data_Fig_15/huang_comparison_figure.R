library(dplyr)
library(ggplot2)
library(feather)
options(stringsAsFactors = F)
source("sankey_functions.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/external/mouse_GABA_Paul_2016_20171002/anno.feather")

# rearrange the target clusters based on which query cluster has the highest freq
new_order <- anno %>%
  group_by(cluster_id, cell_class_id) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(-Freq) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(cell_class_id, cluster_id) %>%
  mutate(new_order = 1:n()) %>%
  select(cluster_id, new_order)

huang_to_new <- left_join(anno, new_order) %>%
  mutate(cluster_id = new_order)

# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
h2n_nodes <- make_plot_nodes(make_group_nodes(huang_to_new, c("cell_class","cluster")))

h2n_new_nodes <- h2n_nodes %>%
  filter(group == "cluster") %>%
  mutate(name = sub("^[0-9]+ ","",name),
         name = paste0(name, " (n = ",n,")"))
h2n_huang_nodes <- h2n_nodes %>%
  filter(group == "cell_class") %>%
  mutate(name = paste0(name, " (n = ",n,")"))

huang_to_new_river <- build_river_plot(huang_to_new, 
                                       c("cell_class","cluster"), 
                                       fill_group = "cell_class") +
  geom_text(data = h2n_new_nodes,
            aes(x = xmax + 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  geom_text(data = h2n_huang_nodes,
            aes(x = xmin - 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  scale_color_identity()

huang_to_new_river

ggsave("huang_to_new_river.pdf", huang_to_new_river, height = 10, width = 6)

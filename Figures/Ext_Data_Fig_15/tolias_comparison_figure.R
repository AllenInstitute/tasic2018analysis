library(dplyr)
library(ggplot2)
library(feather)
options(stringsAsFactors = F)
source("sankey_functions.R")

v1_alm_anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
v1_alm_cl_anno <- v1_alm_anno %>%
  select(cl, cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  mutate(cl = as.character(cl))

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/map.tree.df.rda")
map.tree.df <- map.tree.df %>%
  mutate(sample_id = rownames(.))

map.dat <- map.tree.df %>%
  mutate(cluster_label = as.character(cluster_label)) %>%
  mutate(cluster_label = ifelse(is.na(cluster_label), as.character(cl), cluster_label)) %>%
  mutate(cl = as.character(cl))

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/tolias_patchseq_20161102/anno.feather")

# Remove old mapping and add new mapping from map.tree.df
anno <- anno %>%
  select(-cluster_id, -cluster_label, -cluster_color) %>%
  left_join(map.dat) %>%
  rowwise() %>%
  mutate(cluster_color = ifelse(cl %in% v1_alm_cl_anno$cl,
                                v1_alm_cl_anno$cluster_color[v1_alm_cl_anno$cl == cl],
                                "#808080")) %>%
# %>%
#   mutate(cluster_id = ifelse(cl %in% v1_alm_cl_anno$cl,
#                                 v1_alm_cl_anno$cluster_id[v1_alm_cl_anno$cl == cl],
#                                 0)) %>%
  ungroup()

# fix cluster ids
anno_cluster_ids <- anno %>%
  select(cluster_label) %>%
  unique() %>%
  mutate(cluster_id = 1:n())

anno <- anno %>%
  left_join(anno_cluster_ids)

# rearrange the target clusters based on which query cluster has the highest freq
new_order <- anno %>%
  group_by(cluster_id, cluster_label, class_id, cluster_color) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(-Freq) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(class_id, cluster_id) %>%
  mutate(new_order = 1:n()) %>%
  select(cluster_id, cluster_label, cluster_color, Freq, new_order) %>%
  mutate(manual_order = c(2,1,7,4,3,8,5,9,11,6,10,12,13,14))

#edited_new_order <- edit(new_order)

tolias_to_new <- left_join(anno, new_order) %>%
  mutate(cluster_id = manual_order) %>%
  # select(-cluster_color) %>%
  # left_join(v1_alm_cl_anno) %>%
  mutate(cluster_color = ifelse(is.na(cluster_color), "#808080", cluster_color))

# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
t2n_nodes <- make_plot_nodes(make_group_nodes(tolias_to_new, c("class","cluster")),
                             pad = 0.2)

t2n_new_nodes <- t2n_nodes %>%
  filter(group == "cluster") %>%
  mutate(name = sub("^[0-9]+ ","",name),
         name = paste0(name, " (n = ",n,")"))
t2n_tolias_nodes <- t2n_nodes %>%
  filter(group == "class") %>%
  mutate(name = paste0(name, " (n = ",n,")"))

tolias_to_new_river <- build_river_plot(tolias_to_new, 
                                       c("class","cluster"), 
                                       fill_group = "class",
                                       pad = 0.2) +
  geom_text(data = t2n_new_nodes,
            aes(x = xmax + 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0,
            size = 2) +
  geom_text(data = t2n_tolias_nodes,
            aes(x = xmin - 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1,
            size = 2) +
  scale_color_identity()

tolias_to_new_river

ggsave("tolias_to_new_river.pdf", tolias_to_new_river, height = 3, width = 2)

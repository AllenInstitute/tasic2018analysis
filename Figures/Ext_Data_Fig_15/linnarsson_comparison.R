library(dplyr)
library(feather)
library(ggplot2)
options(stringsAsFactors = F)

source("sankey_functions.R")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process/WGCNA_result_anterograde/map.linnarsson.ss.df.rda")

ss.map.df <- ss.map.df %>%
  mutate(sample_id = rownames(.)) %>%
  rename(ll_label = cl)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170913/anno.feather")

cluster_anno <- anno %>%
  select(cl, cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  mutate(pred.cl = as.character(cl))



ll_anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/external/mouse_SS_CA1_Zeisel_2015_20170620/anno.feather")
ll_cluster_anno <- ll_anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique()
names(ll_cluster_anno) <- c("ll_id","ll_label","ll_color")


map_anno <- ss.map.df %>%
  left_join(cluster_anno) %>%
  left_join(ll_cluster_anno) %>%
  filter(complete.cases(.))

# rearrange the target clusters based on which query cluster has the highest freq
new_order <- map_anno %>%
  group_by(cluster_id, cluster_label, ll_id) %>%
  summarise(Freq = n()) %>%
  ungroup() %>%
  group_by(cluster_id) %>%
  arrange(-Freq) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(cluster_id,ll_id) %>%
  mutate(new_order = 1:n()) %>%
  select(cluster_id, cluster_label, Freq, new_order)

linnarsson_to_new <- left_join(map_anno, new_order) %>%
  mutate(cluster_id = new_order) %>%
  filter(!coarse_cl %in% c("interneurons","pyramidal SS"),
         ll_label != "(none)")

# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
l2n_nodes <- make_plot_nodes(make_group_nodes(linnarsson_to_new, c("ll","cluster")),
                             pad = 0.2)

l2n_new_nodes <- l2n_nodes %>%
  filter(group == "cluster") %>%
  mutate(name = sub("^[0-9]+ ","",name),
         name = paste0(name, " (n = ",n,")"))
l2n_linnarsson_nodes <- l2n_nodes %>%
  filter(group == "ll") %>%
  mutate(name = paste0(name, " (n = ",n,")"))

linnarsson_to_new_river <- build_river_plot(linnarsson_to_new, 
                                        c("ll","cluster"), 
                                        fill_group = "ll",
                                        pad = 0.2) +
  geom_text(data = l2n_new_nodes,
            aes(x = xmax + 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0,
            size = 2) +
  geom_text(data = l2n_linnarsson_nodes,
            aes(x = xmin - 0.01,
                y = (ymin + ymax)/2,
                label = name,
                color = "black"),
            hjust = 1,
            size = 2) +
  scale_color_identity()

linnarsson_to_new_river

ggsave("l2n_glia.pdf",width = 4, height = 12)


library(dplyr)
library(ggplot2)
library(feather)
options(stringsAsFactors = F)
source("sankey_functions.R")

# Load Zizhen's comparisons to the Tasic dataset
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/V1.ref.49.df.rda")
V1.ref.49.df <- map.df %>%
  mutate(pred.cl = as.character(pred.cl),
         org.cl = as.numeric(as.character(org.cl)),
         org.cl_label = as.character(org.cl_label))
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/V1.ref.101.df.rda")
V1.ref.101.df <- map.df %>%
  mutate(pred.cl = as.character(pred.cl),
         #org.cl = as.numeric(as.character(org.cl)),
         org.cl_label = as.character(org.cl))

# get annotations for Tasic 2016 and the current dataset
anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
anno_tasic <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_VISp_SMV1_1679/anno.feather")

cluster_anno_tasic <- anno_tasic %>%
  filter(final_id %% 1 == 0) %>%
  select(final_id, final_label, final_color) %>%
  unique()

cluster_anno <- anno %>%
  select(cl, dendcluster_id, dendcluster_label, dendcluster_color) %>%
  unique()

new_to_tasic <- V1.ref.49.df %>%
  rename(cl = org.cl,
         final_label = pred.cl) %>%
  mutate(cl = as.numeric(as.character(cl))) %>%
  left_join(cluster_anno, by = "cl") %>%
  left_join(cluster_anno_tasic, by = "final_label") %>%
  filter(prob > 0.9)


# rearrange the target clusters based on which query cluster has the highest freq
tasic_order <- new_to_tasic %>%
  group_by(dendcluster_id, final_id) %>%
  summarise(link_n = n()) %>%
  group_by(final_id) %>%
  arrange(-link_n) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(dendcluster_id, final_id) %>%
  mutate(tasic_order = 1:n()) %>%
  select(final_id, tasic_order) %>% 
  unique()

new_to_tasic <- left_join(new_to_tasic, tasic_order) %>%
  mutate(final_id = tasic_order)

# format modification for compatibility with anno-based scripts
repeat_rows <- function(df, reps) {
  out <- df[0,]
  for(i in 1:nrow(df)) {
    n <- df[i,reps]
    if(n > 0) {
      for(j in 1:n) {
        out <- rbind(out, df[i,])
      }
    }
  }
  out
} 


# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
n2t_nodes <- make_plot_nodes(make_group_nodes(new_to_tasic, c("dendcluster","final"))) %>%
  filter(n > 4)

n2t_new_nodes <- n2t_nodes %>%
  filter(group == "dendcluster")
n2t_tasic_nodes <- n2t_nodes %>%
  filter(group == "final")

n2t_links <- make_group_links(new_to_tasic, c("dendcluster","final"), n2t_nodes) %>%
  filter(n > 4)

n2t_plot_links <- make_plot_links(n2t_links, 
                                  fill = "dendcluster")

new_to_tasic_river <- build_river_plot_predefined(n2t_nodes,
                                                  n2t_plot_links) +
  geom_text(data = n2t_new_nodes,
            aes(x = xmin,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  geom_text(data = n2t_tasic_nodes,
            aes(x = xmax,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  scale_color_identity()

new_to_tasic_river

ggsave("randomforest_new_to_tasic_river_gt90.pdf", new_to_tasic_river, height = 24, width = 6)

# L4 Only
new_to_tasic_l4 <- new_to_tasic %>%
  filter(dendcluster_id == 7)

# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
n2t_nodes <- make_plot_nodes(make_group_nodes(new_to_tasic_l4, c("dendcluster","final")))

n2t_new_nodes <- n2t_nodes %>%
  filter(group == "dendcluster")
n2t_tasic_nodes <- n2t_nodes %>%
  filter(group == "final")

new_to_tasic_river <- build_river_plot(new_to_tasic_l4, 
                                       c("dendcluster","final"), 
                                       fill_group = "dendcluster") +
  geom_text(data = n2t_new_nodes,
            aes(x = xmin,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  geom_text(data = n2t_tasic_nodes,
            aes(x = xmax,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  scale_color_identity()

new_to_tasic_river

ggsave("randomforest_new_to_tasic_river_L4.pdf", new_to_tasic_river, height = 4, width = 6)




tasic_to_new <- V1.ref.101.df %>%
  rename(cl = pred.cl,
         final_label = org.cl) %>%
  mutate(cl = as.numeric(as.character(cl))) %>%
  left_join(cluster_anno, by = "cl") %>%
  left_join(cluster_anno_tasic, by = "final_label")

# rearrange the target clusters based on which query cluster has the highest freq
new_order <- tasic_to_new %>%
  group_by(dendcluster_id, final_id) %>%
  summarise(link_n = n()) %>%
  group_by(dendcluster_id) %>%
  arrange(-link_n) %>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  arrange(final_id, dendcluster_id) %>%
  mutate(new_order = 1:n()) %>%
  select(dendcluster_id, new_order) %>% 
  unique()

# new_order <- tasic_to_new %>%
#   group_by(dendcluster_id) %>%
#   arrange(-Rf.map.num) %>%
#   filter(row_number() == 1) %>%
#   ungroup() %>%
#   arrange(final_id, dendcluster_id) %>%
#   mutate(new_order = 1:n()) %>%
#   select(dendcluster_id, new_order)

tasic_to_new <- left_join(tasic_to_new, new_order) %>%
  mutate(dendcluster_id = new_order)

# 
# tasic_to_new_expanded <- tasic_to_new %>%
#   repeat_rows("Rf.map.num")


# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
t2n_nodes <- make_plot_nodes(make_group_nodes(tasic_to_new, c("final","dendcluster")),
                             pad = 0.2)

t2n_new_nodes <- t2n_nodes %>%
  filter(group == "dendcluster")
t2n_tasic_nodes <- t2n_nodes %>%
  filter(group == "final")

tasic_to_new_river <- build_river_plot(tasic_to_new, 
                                       c("final","dendcluster"), 
                                       fill_group = "final",
                                       pad = 0.2) +
  geom_text(data = t2n_new_nodes,
            aes(x = xmax,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  geom_text(data = t2n_tasic_nodes,
            aes(x = xmin,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  scale_color_identity()

tasic_to_new_river

ggsave("randomforest_tasic_to_new_river.pdf", tasic_to_new_river, height = 24, width = 6)


## Tasic L4 only

# make_plot_nodes converts group_nodes to rectangles for plotting with geom_rect()
tasic_to_new_l4 <- tasic_to_new %>%
  filter(dendcluster_id %in% 63:66)
t2n_nodes <- make_plot_nodes(make_group_nodes(tasic_to_new_l4, c("final","dendcluster")),
                             pad = 0.2)

t2n_new_nodes <- t2n_nodes %>%
  filter(group == "dendcluster")
t2n_tasic_nodes <- t2n_nodes %>%
  filter(group == "final")

tasic_to_new_river <- build_river_plot(tasic_to_new_l4, 
                                       c("final","dendcluster"), 
                                       fill_group = "final",
                                       pad = 0.2) +
  geom_text(data = t2n_new_nodes,
            aes(x = xmax,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 0) +
  geom_text(data = t2n_tasic_nodes,
            aes(x = xmin,
                y = (ymin + ymax)/2,
                label = name,
                color = color),
            hjust = 1) +
  scale_color_identity()

tasic_to_new_river

ggsave("randomforest_tasic_to_new_river_L4_L5IT.pdf", tasic_to_new_river, height = 4, width = 6)


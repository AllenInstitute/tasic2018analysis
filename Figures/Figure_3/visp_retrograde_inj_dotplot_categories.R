library(feather)
library(ggplot2)
library(dplyr)
library(dendextend)
library(scrattch)
options(stringsAsFactors = F)

source("prune_leaf_custom.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
#
anno <- anno %>%
  filter(cluster_id %in% 1:133)

# anno <- anno %>%
#   mutate(inj_mat_label = ifelse(donor_label == "309020","AAV2-retro-EF1a-Cre",inj_mat_label))

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

str <- "CP-c"
ctx <- c("MOp","VISp-c","RSP","ACA")
th <- c("LP","LD","LGd")
tec <- "SCs"
p <- "PG"
my <- ""

inj_cat_anno <- data.frame(inj_primary_label = paste("VISp :",c(str,ctx,th,tec,p,my)),
                           inj_cat_id = c(rep(2,length(str)),
                                          rep(1,length(ctx)),
                                          rep(3,length(th)),
                                          rep(4,length(tec)),
                                          rep(5,length(p)),
                                          rep(6,length(my))
                           ),
                           inj_cat_label = c(rep("STR",length(str)),
                                             rep("CTX",length(ctx)),
                                             rep("TH",length(th)),
                                             rep("TEC",length(tec)),
                                             rep("P",length(p)),
                                             rep("MY",length(my))),
                           inj_cat_color = c(rep("#97D6F9",length(str)),
                                             rep("#209F5C",length(ctx)),
                                             rep("#FF90A1",length(th)),
                                             rep("#FF7AFF",length(tec)),
                                             rep("#FFBA87",length(p)),
                                             rep("#FFB3D9",length(my))
                           ))

visp_retro_cells <- anno %>%
  filter(grepl("VISp : ",inj_primary_label)) %>%
  filter(inj_flag_label =="OK") %>%
  left_join(inj_cat_anno)

visp_clusters <- unique(anno$cluster_id[grepl("ALM",anno$cluster_label)])
inh_clusters <- 1:60
nn_clusters <- 116:133

non_visp_clusters <- c(visp_clusters, inh_clusters, nn_clusters)
visp_clusters <- setdiff(1:133, non_visp_clusters)

non_visp_cluster_labels <- unique(anno$cluster_label[anno$cluster_id %in% non_visp_clusters])

visp_dend <- dend %>%
  prune.dendrogram(non_visp_cluster_labels)

visp_cluster_anno <- anno %>%
  select(dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  filter(cluster_id %in% visp_clusters)

# Use the filtered dendrogram to assign cluster x-positions
visp_dend_leaves <- as.ggdend(visp_dend)$labels

visp_cluster_anno <- visp_cluster_anno %>%
  left_join(visp_dend_leaves, by = c("cluster_label"="label")) %>%
  ungroup()

visp_inj <- visp_retro_cells %>%
  left_join(visp_cluster_anno)

# visp_order <- data.frame(inj_primary_label = c("VISp","ACA","RSP","CP","PG","LD","LGd","LP","PRT")) %>%
#   mutate(inj_primary_y = 1:n())

visp_order <- visp_inj %>%
  group_by(inj_primary_label) %>%
  summarise(meanclust = mean(cluster_id)) %>%
  ungroup() %>%
  arrange(meanclust) %>%
  mutate(inj_primary_y = 1:n())

visp_inj <- left_join(visp_inj, visp_order)

# By class

# VISp retrograde dotplot
visp_inj_sum <- visp_inj %>%
  group_by(inj_cat_id, inj_cat_label, inj_cat_color) %>%
  summarise(n_cells_exc = sum(class_label == "Glutamatergic" & cluster_id %in% visp_clusters),
            n_cells_inh = sum(class_label == "GABAergic"),
            #n_cells_inh = ifelse(n_cells_inh == 0, NA, n_cells_inh),
            n_cells_non = sum(!class_label %in% c("Glutamatergic","GABAergic")),
            #n_cells_non = ifelse(n_cells_non == 0, NA, n_cells_non),
            n_expts = length(unique(donor_id))) %>%
  ungroup() %>%
  arrange(inj_cat_id) %>%
  mutate(ypos = n():1)

visp_inj_clust <- visp_inj %>%
  group_by(inj_cat_id, inj_cat_label, inj_cat_color, 
           cluster_id, cluster_label, cluster_color, 
           class_label, x) %>%
  filter(class_label == "Glutamatergic") %>%
  summarise(n_clust = n()) %>%
  ungroup()

visp_clust <- visp_inj_clust %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique()

visp_pd <- visp_inj_clust %>%
  left_join(visp_inj_sum) %>%
  left_join(visp_clust) %>%
  mutate(heat_color = values_to_colors(n_clust, 
                                       maxval = 300, 
                                       colorset = c("white","orangered","red"))) %>%
  mutate(heat_color = ifelse(n_clust == 1,"#CAF9F2",heat_color))

region_hlines <- visp_inj_sum %>%
  group_by(inj_cat_label) %>%
  summarise(y = min(ypos) - 0.5)

inh_nn_labels <- data.frame(x = max(visp_cluster_anno$x) + 2:3 + 0.5,
                            y = max(visp_inj_sum$ypos) + 1,
                            label = c("GABAergic","Non-neuronal"),
                            color = c("orangered","#808080"))

dend_seg <- as.ggdend(visp_dend)$segments %>%
  mutate(y = y/max(y)*5 + max(visp_inj_sum$ypos + 0.5),
         yend = yend/max(yend)*5 + max(visp_inj_sum$ypos) + 0.5)

heat_bg_rect <- data.frame(xmin = min(visp_cluster_anno$x, na.rm = T) - 0.5,
                           xmax = max(visp_cluster_anno$x, na.rm = T) + 0.5,
                           ymin = min(visp_pd$ypos) - 0.5,
                           ymax = max(visp_pd$ypos) + 0.5,
                           fill = "white")

visp_plot <- ggplot() +
  # Dendrogram segments
  geom_segment(data = dend_seg,
               aes(x = x, xend = xend,
                   y = y+0.4, yend = yend+0.4),
               lineend = "square") +
  geom_rect(data = heat_bg_rect,
          aes(xmin = xmin, xmax = xmax,
              ymin = ymin, ymax = ymax,
              fill = fill)) +
  geom_point(data = visp_pd,
             aes(x = x,
                 y = ypos,
                 size = n_clust,
                 fill = cluster_color),
             color = "#000000",
             pch = 21) +
  geom_text(data = visp_inj_sum,
            aes(x = 0,
                y = ypos,
                label = inj_cat_label,
                color = inj_cat_color),
            size = 3,
            hjust = 1, vjust = 0.3) +
  # Right side labels
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 1.5,
                y = ypos,
                label = n_cells_exc),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 3,
                y = ypos,
                label = n_cells_inh),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 4.5,
                y = ypos,
                label = n_cells_non),
            size = 3) +
  # cluster labels
  geom_text(data = visp_cluster_anno,
            aes(x = x,
                y = 0,
                label = cluster_label,
                color = cluster_color),
            size = 3,
            angle = 90, 
            hjust = 1,
            vjust = 0.3) +
  scale_size_area(max_size = 4, breaks = c(1,5,10,50,100,200)) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_void() +
  scale_y_continuous(expand = c(0,1), limits = c(-3, 15)) +
  scale_x_continuous(limits = c(-6, 40)) +
  ggtitle("VISp Retrograde Injections")

visp_plot

ggsave("visp_retro_by_class_dotplot.pdf",visp_plot, width = 11.21951, height = 6, useDingbats = F)

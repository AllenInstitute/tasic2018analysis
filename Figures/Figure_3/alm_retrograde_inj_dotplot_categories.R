library(feather)
library(ggplot2)
library(dplyr)
library(dendextend)
library(scrattch)
options(stringsAsFactors = F)

source("prune_leaf_custom.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

anno <- anno %>%
  filter(cluster_id %in% 1:133)

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

str <- "CP-c"
ctx <- c("ALM-c","MOp","SSp","PERI","ORBl-c","SSs","RSP")
th <- c("VAL","ZI","PF","MD")
tec <- c("SC","SCs")
p <- c("PG","PRNc")
my <- c("IRN","GRN","PARN")

inj_cat_anno <- data.frame(inj_primary_label = paste("ALM :",c(str,ctx,th,tec,p,my)),
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

alm_retro_cells <- anno %>%
  filter(grepl("ALM : ",inj_primary_label)) %>%
  filter(inj_flag_label =="OK") %>%
  left_join(inj_cat_anno) 

visp_clusters <- unique(anno$cluster_id[grepl("VISp",anno$cluster_label)])
inh_clusters <- 1:60
nn_clusters <- 116:133

non_alm_clusters <- c(visp_clusters, inh_clusters, nn_clusters)
alm_clusters <- setdiff(1:133, non_alm_clusters)

non_alm_cluster_labels <- unique(anno$cluster_label[anno$cluster_id %in% non_alm_clusters])

alm_dend <- dend %>%
  prune.dendrogram(non_alm_cluster_labels)

plot(alm_dend)

alm_cluster_anno <- anno %>%
  select(dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
  unique() %>%
  filter(cluster_id %in% alm_clusters)

# Use the filtered dendrogram to assign cluster x-positions
alm_dend_leaves <- as.ggdend(alm_dend)$labels

alm_cluster_anno <- alm_cluster_anno %>%
  left_join(alm_dend_leaves, by = c("cluster_label"="label")) %>%
  ungroup()

alm_inj <- alm_retro_cells %>%
  left_join(alm_cluster_anno)

# alm_order <- data.frame(inj_primary_label = c("ALM","ACA","RSP","CP","PG","LD","LGd","LP","PRT")) %>%
#   mutate(inj_primary_y = 1:n())

alm_order <- alm_inj %>%
  group_by(inj_primary_label) %>%
  summarise(meanclust = mean(cluster_id)) %>%
  ungroup() %>%
  arrange(meanclust) %>%
  mutate(inj_primary_y = 1:n())

alm_inj <- left_join(alm_inj, alm_order)

# By class

# ALM retrograde dotplot
alm_inj_sum <- alm_inj %>%
  group_by(inj_cat_id, inj_cat_label, inj_cat_color) %>%
  summarise(n_cells_exc = sum(class_label == "Glutamatergic" & cluster_id %in% alm_clusters),
            n_cells_inh = sum(class_label == "GABAergic"),
            #n_cells_inh = ifelse(n_cells_inh == 0, NA, n_cells_inh),
            n_cells_non = sum(!class_label %in% c("Glutamatergic","GABAergic")),
            #n_cells_non = ifelse(n_cells_non == 0, NA, n_cells_non),
            n_expts = length(unique(donor_id))) %>%
  ungroup() %>%
  arrange(inj_cat_id) %>%
  mutate(ypos = n():1)

alm_inj_clust <- alm_inj %>%
  group_by(inj_cat_id, inj_cat_label, inj_cat_color, 
           cluster_id, cluster_label, cluster_color, 
           class_label, x) %>%
  filter(class_label == "Glutamatergic") %>%
  summarise(n_clust = n()) %>%
  ungroup()

alm_clust <- alm_inj_clust %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique()

alm_pd <- alm_inj_clust %>%
  left_join(alm_inj_sum) %>%
  left_join(alm_clust) %>%
  mutate(heat_color = values_to_colors(n_clust, 
                                       maxval = 300, 
                                       colorset = c("white","orangered","red"))) %>%
  mutate(heat_color = ifelse(n_clust == 1,"#CAF9F2",heat_color))

region_hlines <- alm_inj_sum %>%
  group_by(inj_cat_label) %>%
  summarise(y = min(ypos) - 0.5)

inh_nn_labels <- data.frame(x = max(alm_cluster_anno$x) + 2:3 + 0.5,
                            y = max(alm_inj_sum$ypos) + 1,
                            label = c("GABAergic","Non-neuronal"),
                            color = c("orangered","#808080"))

dend_seg <- as.ggdend(alm_dend)$segments %>%
  mutate(y = y/max(y)*5 + max(alm_inj_sum$ypos + 0.5),
         yend = yend/max(yend)*5 + max(alm_inj_sum$ypos) + 0.5)

heat_bg_rect <- data.frame(xmin = min(alm_cluster_anno$x, na.rm = T) - 0.5,
                           xmax = max(alm_cluster_anno$x, na.rm = T) + 0.5,
                           ymin = min(alm_pd$ypos) - 0.5,
                           ymax = max(alm_pd$ypos) + 0.5,
                           fill = "white")

alm_plot <- ggplot() +
  # Dendrogram segments
  geom_segment(data = dend_seg,
               aes(x = x, xend = xend,
                   y = y+0.4, yend = yend+0.4),
               lineend = "square") +
  geom_rect(data = heat_bg_rect,
          aes(xmin = xmin, xmax = xmax,
              ymin = ymin, ymax = ymax,
              fill = fill)) +
  geom_point(data = alm_pd,
             aes(x = x,
                 y = ypos,
                 size = n_clust,
                 fill = cluster_color),
             color = "#000000",
             pch = 21) +
  # Left side labels
  geom_text(data = alm_inj_sum,
            aes(x = 0,
                y = ypos,
                label = inj_cat_label,
                color = inj_cat_color),
            size = 3,
            hjust = 1, vjust = 0.3) +
  # Right side labels
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 1.5,
                y = ypos,
                label = n_cells_exc),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 3,
                y = ypos,
                label = n_cells_inh),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 4.5,
                y = ypos,
                label = n_cells_non),
            size = 3) +
  # cluster labels
  geom_text(data = alm_cluster_anno,
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
  scale_x_continuous(limits = c(-6, 35)) +
  ggtitle("ALM Retrograde Injections")

alm_plot

ggsave("alm_retro_by_class_dotplot.pdf",alm_plot, width = 10, height = 6, useDingbats = F)

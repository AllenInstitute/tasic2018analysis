library(feather)
library(ggplot2)
library(dplyr)
library(dendextend)
options(stringsAsFactors = F)

source("prune_leaf_custom.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
#anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170905_exons_only/anno.feather")

anno <- anno %>%
  filter(cluster_id %in% 1:133)

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

alm_retro_cells <- anno %>%
  filter(grepl("ALM : ",inj_primary_label)) %>%
  filter(inj_flag_label == "OK")

visp_clusters <- unique(anno$cluster_id[grepl("VISp",anno$cluster_label)])
inh_clusters <- 1:60
nn_clusters <- 116:133

non_alm_clusters <- c(visp_clusters, inh_clusters, nn_clusters)
alm_clusters <- setdiff(1:133, non_alm_clusters)

non_alm_cluster_labels <- unique(anno$cluster_label[anno$cluster_id %in% non_alm_clusters])

alm_dend <- dend %>%
  prune.dendrogram(non_alm_cluster_labels)

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

alm_order <- data.frame(inj_primary_label = paste("ALM :",c("MOp","ALM-c","ORBl-c","SSp","PERI","SSs","RSP","CP-c",
                                              "GRN","PARN","IRN","PRNc","PG",
                                              "SC","MD","ZI","PF"))) %>%
  mutate(inj_primary_y = 1:n())

alm_inj <- alm_inj %>%
  mutate(inj_primary_label = ifelse(inj_primary_label == "ALM : VAL","ALM : MD",inj_primary_label)) %>%
  left_join(alm_order)

### By Target Region


# ALM retrograde dotplot
alm_inj_sum <- alm_inj %>%
  group_by(inj_primary_y, inj_primary_label) %>%
  summarise(n_cells_exc = sum(class_label == "Glutamatergic" & cluster_id %in% alm_clusters),
            n_cells_vis = sum(grepl("VISp", cluster_label)),
            n_cells_inh = sum(class_label == "GABAergic"),
            #n_cells_inh = ifelse(n_cells_inh == 0, NA, n_cells_inh),
            n_cells_non = sum(!class_label %in% c("Glutamatergic","GABAergic")),
            #n_cells_non = ifelse(n_cells_non == 0, NA, n_cells_non),
            n_expts = length(unique(donor_id)),
            n_aav_gfp = sum(inj_mat_label == "AAV-SL1-CAG-GFP"),
            n_aav_tdt = sum(inj_mat_label == "AAV-SL1-CAG-tdTomato"),
            n_aav2_cre = sum(inj_mat_label == "AAV2-retro-EF1a-Cre"),
            n_aav2_tdt = sum(inj_mat_label == "AAV2-retro-EF1a-tdTomato"),
            n_cav_cre = sum(inj_mat_label == "CAV2-Cre.V795"),
            n_rv_cre = sum(inj_mat_label == "RV-dGL-Cre")) %>%
  ungroup() %>%
  arrange(inj_primary_y) %>%
  mutate(ypos = n():1)

alm_inj_clust <- alm_inj %>%
  group_by(inj_primary_y, inj_primary_id, inj_primary_label, 
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
  left_join(alm_clust)

region_hlines <- alm_inj_sum %>%
  group_by(inj_primary_label) %>%
  summarise(y = min(ypos) - 0.5)

inh_nn_labels <- data.frame(x = max(alm_cluster_anno$x) + 2:3 + 0.5,
                            y = max(alm_inj_sum$ypos) + 1,
                            label = c("GABAergic","Non-neuronal"),
                            color = c("orangered","#808080"))

dend_seg <- as.ggdend(alm_dend)$segments %>%
  mutate(y = y/max(y)*5 + max(alm_inj_sum$ypos + 0.5),
         yend = yend/max(yend)*5 + max(alm_inj_sum$ypos) + 0.5)

alm_plot <- ggplot() +
  # Dendrogram segments
  geom_segment(data = dend_seg,
               aes(x = x, xend = xend,
                   y = y, yend = yend),
               lineend = "square") +
  # Dashed lines for clusters
  geom_segment(data = alm_cluster_anno,
               aes(x = x, xend = x,
                   y = 0.5, yend = max(alm_inj_sum$ypos)+0.5,
                   color = "#808080"),
               linetype = "dashed") +
  # Horizontal lines separating target regions
  geom_segment(data = region_hlines,
               aes(y = y, yend = y,
                   x = -6, xend = max(alm_cluster_anno$x) + 5),
               color = "#808080") +
  # Points for each injection-cluster result
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
                label = sub("^ALM : ","",inj_primary_label)),
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
                label = n_cells_vis),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 4.5,
                y = ypos,
                label = n_cells_inh),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 6,
                y = ypos,
                label = n_cells_non),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 7.5,
                y = ypos,
                label = n_aav_gfp),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 9,
                y = ypos,
                label = n_aav_tdt),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 10.5,
                y = ypos,
                label = n_aav2_cre),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 12,
                y = ypos,
                label = n_aav2_tdt),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 13.5,
                y = ypos,
                label = n_cav_cre),
            size = 3) +
  geom_text(data = alm_inj_sum,
            aes(x = max(alm_cluster_anno$x) + 15,
                y = ypos,
                label = n_rv_cre),
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
  scale_size_area(max_size = 5, breaks = c(5,10,20,40,80)) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_void() +
  scale_y_continuous(expand = c(0,1), limits = c(-5, 25)) +
  scale_x_continuous(limits = c(-6, 42)) +
  ggtitle("ALM Retrograde Injections")

alm_plot

ggsave("alm_retro_by_target.pdf",alm_plot, width = 7, height = 5, useDingbats = F)



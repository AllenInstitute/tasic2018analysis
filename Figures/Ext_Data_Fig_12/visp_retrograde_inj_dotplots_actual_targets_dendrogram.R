library(feather)
library(ggplot2)
library(dplyr)
library(dendextend)
options(stringsAsFactors = F)

source("prune_leaf_custom.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
#anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_visp_20170905_exons_only/anno.feather")
anno <- anno %>%
  filter(cluster_id %in% 1:133)

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

#inj_qc <- read.table("20171030_inj_qc.tsv", sep = "\t", header = T)

visp_retro_cells <- anno %>%
  filter(grepl("VISp : ",inj_primary_label)) %>%
  filter(inj_flag_label == "OK")
  # filter(inj_flag %in% c("OK","checking")) %>%
  # filter(inj_primary != "BLA")

alm_clusters <- unique(anno$cluster_id[grepl("ALM",anno$cluster_label)])
inh_clusters <- 1:60
nn_clusters <- 116:133

non_visp_clusters <- c(alm_clusters, inh_clusters, nn_clusters)
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

visp_order <- data.frame(inj_primary_label = paste("VISp :",c("ACA","RSP","VISp-c","CP-c","SCs","PG","LD","LP","LGd"))) %>%
  mutate(inj_primary_y = 1:n())

visp_inj <- left_join(visp_inj, visp_order) %>%
  mutate(inj_mat_label = ifelse(donor_label == "309020","AAV2-retro-EF1a-Cre",inj_mat_label))

### By Target Region


# visp retrograde dotplot
visp_inj_sum <- visp_inj %>%
  group_by(inj_primary_y, inj_primary_label) %>%
  summarise(n_cells_exc = sum(class_label == "Glutamatergic" & cluster_id %in% visp_clusters),
            n_cells_alm = sum(grepl("ALM",cluster_label)),
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

visp_inj_clust <- visp_inj %>%
  group_by(inj_primary_y, inj_primary_id, inj_primary_label, 
           cluster_id, cluster_label, cluster_color, 
           class_label, x) %>%
  filter(cluster_id %in% visp_clusters) %>%
  summarise(n_clust = n()) %>%
  ungroup()

visp_clust <- visp_inj_clust %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique()

visp_pd <- visp_inj_clust %>%
  left_join(visp_inj_sum) %>%
  left_join(visp_clust)

region_hlines <- visp_inj_sum %>%
  group_by(inj_primary_label) %>%
  summarise(y = min(ypos) - 0.5)

inh_nn_labels <- data.frame(x = max(visp_cluster_anno$x) + 2:3 + 0.5,
                            y = max(visp_inj_sum$ypos) + 1,
                            label = c("GABAergic","Non-neuronal"),
                            color = c("orangered","#808080"))

dend_seg <- as.ggdend(visp_dend)$segments %>%
  mutate(y = y/max(y)*5 + max(visp_inj_sum$ypos + 0.5),
         yend = yend/max(yend)*5 + max(visp_inj_sum$ypos) + 0.5)

visp_plot <- ggplot() +
  # Dendrogram segments
  geom_segment(data = dend_seg,
               aes(x = x, xend = xend,
                   y = y, yend = yend),
               lineend = "square") +
  # Dashed lines for clusters
  geom_segment(data = visp_cluster_anno,
               aes(x = x, xend = x,
                   y = 0.5, yend = max(visp_inj_sum$ypos)+0.5,
                   color = "#808080"),
               linetype = "dashed") +
  # Horizontal lines separating target regions
  geom_segment(data = region_hlines,
               aes(y = y, yend = y,
                   x = -6, xend = max(visp_cluster_anno$x) + 5),
               color = "#808080") +
  # Points for each injection-cluster result
  geom_point(data = visp_pd,
             aes(x = x,
                 y = ypos,
                 size = n_clust,
                 fill = cluster_color),
             color = "#000000",
             pch = 21) +
  # Left side labels
  geom_text(data = visp_inj_sum,
            aes(x = 0,
                y = ypos,
                label = sub("^VISp : ","",inj_primary_label)),
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
                label = n_cells_alm),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 4.5,
                y = ypos,
                label = n_cells_inh),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 6,
                y = ypos,
                label = n_cells_non),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 7.5,
                y = ypos,
                label = n_aav_gfp),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 9,
                y = ypos,
                label = n_aav_tdt),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 10.5,
                y = ypos,
                label = n_aav2_cre),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 12,
                y = ypos,
                label = n_aav2_tdt),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 13.5,
                y = ypos,
                label = n_cav_cre),
            size = 3) +
  geom_text(data = visp_inj_sum,
            aes(x = max(visp_cluster_anno$x) + 15,
                y = ypos,
                label = n_rv_cre),
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
  scale_size_area(max_size = 4, breaks = c(5,10,20,40,80)) +
  scale_color_identity() +
  scale_fill_identity() +
  theme_void() +
  scale_y_continuous(expand = c(0,1), limits = c(-5, 18)) +
  scale_x_continuous(limits = c(-6, 49)) +
  ggtitle("VISp Retrograde Injections")

visp_plot

ggsave("visp_retro_by_target.pdf",visp_plot, width = 8.5, height = 4.5, useDingbats = F)



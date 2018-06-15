library(dplyr)
library(feather)
library(ggplot2)
options(stringsAsFactors = F)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Lamp5.Vip.tsne.df.rda")

lv_tsne <- Lamp5.Vip.tsne.df

lv_samples <- rownames(lv_tsne)

lv_tsne <- lv_tsne %>%
  mutate(sample_id = lv_samples) %>%
  #mutate(cl = as.numeric(as.character(cl))) %>%
  left_join(anno, by = "sample_id")

lv_clust <- ggplot() +
  geom_point(data = lv_tsne,
             aes(x = Lim1,
                 y = Lim2,
                 color = cluster_color),
             size = 0.4) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("lv_clust_tsne.png", lv_clust, width = 3, height = 3, dpi = 1000)
ggsave("lv_clust_tsne.pdf", lv_clust, width = 3, height = 3, dpi = 1000, useDingbats = FALSE)

lv_clust_labels <- lv_tsne %>%
  group_by(cluster_id, cluster_label, cluster_color) %>%
  summarise(x = mean(Lim1, trim = 0.5),
         y = mean(Lim2, trim = 0.5))

lv_clust_label_plot <- ggplot() +
  geom_point(data = lv_tsne,
             aes(x = -Lim1,
                 y = Lim2,
                 color = cluster_color),
             size = 0.4) +
  geom_text(data = lv_clust_labels,
             aes(x = -x, y = y,
                 label = cluster_label)) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("lv_clust_tsne_labels.png", lv_clust_label_plot, width = 6, height = 6, dpi = 1000)
ggsave("lv_clust_tsne_labels.pdf", lv_clust_label_plot, width = 6, height = 6, dpi = 1000, useDingbats = FALSE)


lv_region <- ggplot() +
  geom_point(data = lv_tsne,
             aes(x = Lim1,
                 y = Lim2,
                 color = region_color),
             size = 0.4) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("lv_region_tsne.png", lv_region, width = 3, height = 3, dpi = 1000)
ggsave("lv_region_tsne.pdf", lv_region, width = 3, height = 3, dpi = 1000, useDingbats = FALSE)


load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.Pvalb.tsne.df.rda")

sp_tsne <- Sst.Pvalb.tsne.df

sp_samples <- rownames(sp_tsne)

sp_tsne <- sp_tsne %>%
  mutate(sample_id = sp_samples) %>%
  #mutate(cl = as.numeric(as.character(cl))) %>%
  left_join(anno, by = "sample_id")

sp_clust <- ggplot() +
  geom_point(data = sp_tsne,
             aes(x = Lim1,
                 y = Lim2,
                 color = cluster_color),
             size = 0.4) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("sp_clust_tsne.png", sp_clust, width = 3, height = 3, dpi = 1000)
ggsave("sp_clust_tsne.pdf", sp_clust, width = 3, height = 3, dpi = 1000, useDingbats = FALSE)

sp_clust_labels <- sp_tsne %>%
  group_by(cluster_id, cluster_label, cluster_color) %>%
  summarise(x = mean(Lim1, trim = 0.5),
            y = mean(Lim2, trim = 0.5))

sp_clust_label_plot <- ggplot() +
  geom_point(data = sp_tsne,
             aes(x = -Lim1,
                 y = Lim2,
                 color = cluster_color),
             size = 0.4) +
  geom_text(data = sp_clust_labels,
            aes(x = -x, y = y,
                label = cluster_label)) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("sp_clust_tsne_labels.png", sp_clust_label_plot, width = 6, height = 6, dpi = 1000)
ggsave("sp_clust_tsne_labels.pdf", sp_clust_label_plot, width = 6, height = 6, dpi = 1000, useDingbats = FALSE)



sp_region <- ggplot() +
  geom_point(data = sp_tsne,
             aes(x = Lim1,
                 y = Lim2,
                 color = region_color),
             size = 0.4) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("sp_region_tsne.png", sp_region, width = 3, height = 3, dpi = 1000)
ggsave("sp_region_tsne.pdf", sp_region, width = 3, height = 3, dpi = 1000, useDingbats = FALSE)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/ex.tsne.df.rda")

gluta_tsne <- ex.tsne.df

gluta_samples <- rownames(gluta_tsne)

gluta_tsne <- gluta_tsne %>%
  mutate(sample_id = gluta_samples) %>%
  mutate(cl = as.numeric(as.character(cl))) %>%
  left_join(anno, by = "sample_id")

gluta_clust <- ggplot() +
  geom_point(data = gluta_tsne,
             aes(x = Lim1,
                 y = Lim2,
                 color = cluster_color),
             size = 0.5) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("gluta_clust_tsne.png", gluta_clust, width = 4.5, height = 4.5, dpi = 1000)
ggsave("gluta_clust_tsne.pdf", gluta_clust, width = 3, height = 3, dpi = 1000, useDingbats = FALSE)

gluta_clust_labels <- gluta_tsne %>%
  group_by(cluster_id, cluster_label, cluster_color) %>%
  summarise(x = mean(Lim1, trim = 0.5),
            y = mean(Lim2, trim = 0.5))

gluta_clust_label_plot <- ggplot() +
  geom_point(data = gluta_tsne,
             aes(x = -Lim1,
                 y = Lim2,
                 color = cluster_color),
             size = 0.4) +
  geom_text(data = gluta_clust_labels,
            aes(x = -x, y = y,
                label = cluster_label)) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("gluta_clust_tsne_labels.png", gluta_clust_label_plot, width = 12, height = 12, dpi = 1000)


gluta_region <- ggplot() +
  geom_point(data = gluta_tsne,
             aes(x = Lim1,
                 y = Lim2,
                 color = region_color),
             size = 0.5) +
  scale_color_identity() +
  theme_void() +
  theme(panel.background = element_rect(fill = "#DAEFF1"))

ggsave("gluta_region_tsne.png", gluta_region, width = 4.5, height = 4.5, dpi = 1000)
ggsave("gluta_region_tsne.pdf", gluta_region, width = 3, height = 3, dpi = 1000, useDingbats = FALSE)



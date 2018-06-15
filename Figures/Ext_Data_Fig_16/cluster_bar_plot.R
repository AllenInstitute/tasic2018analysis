library(dplyr)
library(ggplot2)
library(feather)
options(stringsAsFactors = F)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")


vipr2_anno <- anno %>%
  filter(grepl("Vipr2",genotype_label)) %>%
  filter(!grepl("Pvalb",genotype_label)) %>%
  # Collapse two Vipr2;Ai14 genotypes
  mutate(genotype_label = ifelse(genotype_label == "Vipr2-IRES2-Cre/wt;PhiC31-neo/Ai14(RCL-tdT)",
                              "Vipr2-IRES2-Cre/wt;Ai14(RCL-tdT)/wt",
                              genotype_label)) 

cluster_xpos <- vipr2_anno %>%
  select(dendcluster_id) %>%
  unique() %>%
  arrange(dendcluster_id) %>%
  mutate(xpos = 1:n())

cre_ypos <- vipr2_anno %>%
  select(genotype_label) %>%
  unique() %>%
  mutate(ypos = c(0.5,0))

vipr2_clusters <- vipr2_anno %>%
  select(genotype_label, dendcluster_id, cluster_label, cluster_color, cre_color) %>%
  group_by(genotype_label) %>%
  mutate(n_genotype = n()) %>%
  ungroup() %>%
  group_by(genotype_label, n_genotype, dendcluster_id, cluster_label, cluster_color, cre_color) %>%
  summarise(n_cells = n(),
            frac_cells = n()/n_genotype[1]) %>%
  left_join(cluster_xpos) %>%
  left_join(cre_ypos)

ggplot() +
  geom_point(data = vipr2_clusters,
             aes(x = xpos, y = ypos,
                 size = n_cells,
                 color = cluster_color)) +
  scale_color_identity() +
  scale_size_area() +
  theme_bw()

y_guides <- data.frame(x = 0.5, 
                       xend = max(vipr2_clusters$xpos + 0.5),
                       y = c(0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9))

y_base <- data.frame(x = 0.5,
                     xend = max(vipr2_clusters$xpos + 0.5),
                     y = c(0, 0.5))

y_guide_labels <- data.frame(x = 0,
                             y = seq(0, .9, by = 0.1),
                             label = c(0,0.1,0.2,0.3,0.4,0,0.1,0.2,0.3,0.4))

p <- ggplot() +
  geom_segment(data = y_guides,
               aes(x = x, xend = xend, y = y, yend = y),
               size = 0.1,
               linetype = "dashed",
               color = "#808080") +
  geom_segment(data = y_base,
               aes(x = x, xend = xend, y = y, yend = y),
               size = 0.2) +
  geom_text(data = y_guide_labels,
            aes(x = x, y = y,
                label = label),
            hjust = 1,
            size = 1) +
  geom_rect(data = vipr2_clusters,
            aes(xmin = xpos - 0.4,
                xmax = xpos + 0.4,
                ymin = ypos,
                ymax = ypos + frac_cells,
                fill = cluster_color)) +
  geom_text(data = vipr2_clusters,
            aes(x = xpos,
                y = -0.02,
                label = cluster_label,
                color = cluster_color),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 1.5) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_continuous(limits = c(-0.5, 1)) +
  theme_void(6)

ggsave("vipr2_barplot.pdf", width = 3, height = 3, useDingbats = F)

## Plotted together
ggplot() +
  geom_bar(data = vipr2_clusters,
             aes(x = xpos, y = frac_cells, group = genotype_label,
                 fill = genotype_color),
           position = position_dodge(width = 0.5),
           stat = "identity") + 
  scale_fill_identity() +
  scale_size_area() +
  scale_x_continuous(breaks = vipr2_clusters$xpos, 
                     labels = vipr2_clusters$cluster_label,
                     limits = c(0.5, max(vipr2_clusters$xpos) + 0.5)) +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3,
                                   size = 1.5,
                                   color = vipr2_clusters$cluster_color))

ggsave("vipr2_barplot_dodge.pdf", width = 10, height = 3, useDingbats = F)

ggplot() +
  geom_rect(data = vipr2_clusters %>%
              filter(genotype_label == "Vipr2-IRES2-Cre/wt;Ai14(RCL-tdT)/wt"),
            aes(xmin = xpos - 0.9,
                xmax = xpos - 0.5,
                ymin = 0,
                ymax = frac_cells,
                fill = cre_color)) +
  geom_rect(data = vipr2_clusters %>%
             filter(genotype_label == "Vipr2-IRES2-Cre/wt;Slc32a1-T2A-FlpO/wt;Ai65(RCFL-tdT)/wt"),
           aes(xmin = xpos - 0.5,
               xmax = xpos - 0.1,
               ymin = 0,
               ymax = frac_cells,
               fill = cre_color)) +
  geom_rect(data = vipr2_clusters,
            aes(xmin = xpos - 1,
                xmax = xpos,
                ymin = -0.01,
                ymax = 0,
                fill = cluster_color)) +
  geom_vline(data = vipr2_clusters,
             aes(xintercept = xpos),
             linetype = "dashed") +
  scale_x_continuous(breaks = vipr2_clusters$xpos-0.5, 
                     labels = vipr2_clusters$cluster_label,
                     limits = c(0, max(vipr2_clusters$xpos) + 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3,
                                   size = 1.5,
                                   color = vipr2_clusters$cluster_color),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_fill_identity()
ggsave("vipr2_barplot_dodge.pdf", width = 10, height = 3, useDingbats = F)

library(dplyr)
library(feather)
library(ggplot2)
options(stringsAsFactors = F)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

anno <- anno %>%
  filter(cluster_id %in% 1:133)

cluster_anno <- anno %>%
  select(cluster_id, cluster_label, cluster_color) %>%
  unique()

cre_anno <- anno %>%
  select(cre_id, cre_label, cre_color) %>%
  unique()

plot_anno <- anno %>%
  filter(inj_type_label == "No Injection") %>%
  filter(facs_label != "RFP-negative")

plot_data <- plot_anno %>%
  group_by(cre_id, cre_label, cre_color,
           cluster_id, cluster_label, cluster_color) %>%
  summarise(n_cells = n())

max_cre <- max(anno$cre_id)
max_cluster <- max(anno$cluster_id)

bg_rects <- data.frame(xmin = 1:(max_cre/2)*2 - 0.5,
                       xmax = (1:(max_cre/2) + 1)*2 - 1.5,
                       ymin = 0.5,
                       ymax = max_cluster + 0.5)

class_breaks <- anno %>%
  group_by(subclass_id) %>%
  summarise(yintercept = max(cluster_id) + 0.5)

ggplot() +
  geom_rect(data = bg_rects,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "#CAD7D7") +
  geom_hline(data = class_breaks,
             aes(yintercept = yintercept),
             size = 0.2) +
  geom_point(data = plot_data,
             aes(x = cre_id, 
                 y = cluster_id,
                 size = n_cells,
                 fill = cluster_color,
                 color = "#FFFFFF"),
             pch = 21) +
  geom_text(data = cluster_anno,
            aes(x = 1,
                y = cluster_id,
                label = cluster_label,
                color = cluster_color),
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_reverse("",
                  limits = c(133.5, -0.5),
                  expand = c(0,0)) +
  scale_x_continuous("",
                     expand = c(0,0),
                     limits = c(-5, 57),
                     breaks = cre_anno$cre_id,
                     labels = cre_anno$cre_label,
                     position = "top") +
  scale_size_area(max_size = 5,
                  breaks = c(1,10,50,100,200,500)) +
  theme_classic(base_size = 7) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 0))

ggsave("cre_dotplot.pdf",width = 8.25, height = 10.5, useDingbats = F)


# Pan-broad_specific ordering
pan_cre <- c("Gad2-IRES-Cre","Snap25-IRES2-Cre","Slc32a1-IRES-Cre","Slc17a7-IRES2-Cre")
broad_cre <- c("Rbp4-Cre_KL100","Pvalb-IRES-Cre","Vip-IRES-Cre","Sst-IRES-Cre","Ctgf-T2A-dgCre","Htr3a-Cre_NO152","Trib2-F2A-CreERT2","Ntsr1-Cre_GN220","Tlx3-Cre_PL56","Ndnf-IRES2-dgCre","Rorb-IRES2-Cre")

plot_anno2 <- plot_anno %>%
  mutate(cre_class_id = ifelse(cre_label %in% pan_cre, 1, 
                               ifelse(cre_label %in% broad_cre,2,3)))


cre_anno2 <- plot_anno2 %>%
  select(cre_id, cre_label, cre_color, cre_class_id) %>%
  unique() %>%
  arrange(cre_class_id, cre_id) %>%
  mutate(new_cre_id = 2:(n()+1))

plot_anno2 <- plot_anno2 %>%
  left_join(cre_anno2)

plot_data2 <- plot_anno2 %>%
  group_by(new_cre_id, cre_label, cre_color,
           cluster_id, cluster_label, cluster_color) %>%
  summarise(n_cells = n())

ggplot() +
  geom_rect(data = bg_rects,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "#CAD7D7") +
  geom_point(data = plot_data2,
             aes(x = new_cre_id, 
                 y = cluster_id,
                 size = n_cells,
                 color = cluster_color)) +
  geom_text(data = cluster_anno,
            aes(x = 1,
                y = cluster_id,
                label = cluster_label,
                color = cluster_color),
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  geom_hline(data = class_breaks,
             aes(yintercept = yintercept),
             size = 0.2) +
  scale_color_identity() +
  scale_y_reverse("",
                  limits = c(116.5, -0.5),
                  expand = c(0,0)) +
  scale_x_continuous("",
                     expand = c(0,0),
                     limits = c(-5, 37),
                     breaks = cre_anno2$new_cre_id,
                     labels = cre_anno2$cre_label,
                     position = "top") +
  scale_size_area(max_size = 3,
                  breaks = c(10,50,100,200,500)) +
  theme_classic(base_size = 7) +
  theme(axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 0))

ggsave("cre_grouped_dotplot.pdf",width = 8.25, height = 10.5, useDingbats = F)

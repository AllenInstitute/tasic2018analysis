library(dplyr)
library(feather)
library(ggplot2)
options(stringsAsFactors = F)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/region.de.df.rda")

cluster_anno <- anno %>%
  select(cl, dendcluster_id, cluster_label, cluster_color, subclass_id, class_id) %>%
  unique()

region.de.df2 <- region.de.df %>%
  mutate(ALM.cl = as.character(ALM.cl),
         VISp.cl = as.character(VISp.cl))

# Custom class order
class_order <- data.frame(class_id = 1:4,
                          class_order = c(2,1,3,4))

de_data <- region.de.df2 %>%
  rowwise() %>%
  mutate(dend_order = sum(cluster_anno$dendcluster_id[match(ALM.cl,cluster_anno$cl)],
                          cluster_anno$dendcluster_id[match(VISp.cl,cluster_anno$cl)])) %>%
  ungroup() %>%
  mutate(alm_subclass = cluster_anno$subclass_id[match(ALM.cl, cluster_anno$cl)],
         visp_subclass = cluster_anno$subclass_id[match(VISp.cl, cluster_anno$cl)],
         class_order = class_order$class_order[match(cluster_anno$class_id[match(VISp.cl, cluster_anno$cl)], class_order$class_id)]) %>%
  arrange(class_order,visp_subclass,region.de.num,dend_order) %>%
  mutate(y = 1:n(),
         ymin = 1:n() - 0.4,
         ymax = 1:n() + 0.4) %>%
  #mutate(alm_xmin = ifelse(region.de.ALM.num > 0, -log10(region.de.ALM.num + 1), 0),
  mutate(alm_xmin = -region.de.ALM.num,
         alm_xmax = 0,
         v1_xmin = 0,
         v1_xmax = region.de.V1.num) %>%
         #v1_xmax = ifelse(region.de.V1.num > 0,log10(region.de.V1.num + 1),0)) %>%
  mutate(comp_label = paste(cluster_anno$cluster_label[match(ALM.cl, cluster_anno$cl)], 
                            "vs",
                            cluster_anno$cluster_label[match(VISp.cl, cluster_anno$cl)])) %>%
  mutate(alm_label = sub("ALM ","",cluster_anno$cluster_label[match(ALM.cl,cluster_anno$cl)]),
         visp_label = sub("VISp ","",cluster_anno$cluster_label[match(VISp.cl,cluster_anno$cl)]),
         alm_color = cluster_anno$cluster_color[match(ALM.cl,cluster_anno$cl)],
         visp_color = cluster_anno$cluster_color[match(VISp.cl,cluster_anno$cl)]) %>%
  mutate(trim_alm_xmin = ifelse(alm_xmin < -200, -200, alm_xmin),
         trim_alm_xmin_2 = ifelse(alm_xmin < -200, -230, NA),
         trim_alm_label = ifelse(alm_xmin < -200, region.de.ALM.num, NA))

same_cl <- de_data %>%
  filter(ALM.cl == VISp.cl)

diff_cl <- de_data %>%
  filter(ALM.cl != VISp.cl)

p_figure1 <- ggplot() +
  # Central boxes
  geom_rect(data = diff_cl,
            aes(xmin = -20, xmax = 0,
                ymin = ymin, ymax = ymax,
                fill = alm_color)) +
  geom_rect(data = diff_cl,
            aes(xmin = 0, xmax = 20,
                ymin = ymin, ymax = ymax,
                fill = visp_color)) +
  geom_vline(aes(xintercept = 0),
             size = 0.5) +
  geom_rect(data = same_cl,
            aes(xmin = -20, xmax = 20,
                ymin = ymin, ymax = ymax,
                fill = alm_color)) +
  # ALM Bars
  geom_rect(data = de_data,
            aes(xmin = trim_alm_xmin - 25, xmax = alm_xmax - 25,
                ymin = ymin, ymax = ymax,
                fill = "#212021")
            ) +
  geom_rect(data = de_data,
            aes(xmin = trim_alm_xmin_2 - 25, xmax = -210 - 25,
                ymin = ymin, ymax = ymax,
                fill = "#212021"))+
  geom_text(data = de_data,
            aes(x = trim_alm_xmin_2 - 35,
                y = y,
                label = trim_alm_label),
            size = 2,
            hjust = 1,
            vjust = 0.3) +
  # VISp Bars
  geom_rect(data = de_data,
            aes(xmin = v1_xmin + 25, xmax = v1_xmax + 25,
                ymin = ymin, ymax = ymax,
                fill = "#848EBC")) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_reverse(expand = c(0,0)) +
  scale_x_continuous("N DE Genes",
                     limits = c(-300,300),
                     breaks = c(-225,-125,-25,25,125,225),
                     labels = c(-200, -100, 0, 0, 100, 200)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.3))

ggsave("region_specific_degenes_barplot.pdf", p_figure1, width = 2, height = 8)

p <- ggplot() +
  # ALM Bars
  geom_rect(data = de_data,
            aes(xmin = alm_xmin, xmax = alm_xmax,
                ymin = ymin, ymax = ymax,
                fill = "#212021")
  ) +
  # ALM Labels
  geom_text(data = de_data,
            aes(y = y,
                x = alm_xmin - 5,
                #label = paste(alm_label,paste0("(",region.de.ALM.num,")")),
                label = alm_label,
                color = alm_color),
            hjust = 1,
            vjust = 0.3,
            size = 2) +
  # VISp Bars
  geom_rect(data = de_data,
            aes(xmin = v1_xmin, xmax = v1_xmax,
                ymin = ymin, ymax = ymax,
                fill = "#848EBC")) +
  # VISp Labels
  geom_text(data = de_data,
            aes(y = y,
                x = v1_xmax + 5,
                #label = paste(paste0("(",region.de.V1.num,")"),visp_label),
                label = visp_label,
                color = visp_color),
            hjust = 0,
            vjust = 0.3,
            size = 2) +
  geom_vline(aes(xintercept = 0)) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_y_reverse() +
  scale_x_continuous("N DE Genes",
                     limits = c(-500,500),
                     breaks = c(-400,-300,-200,-100,0,100,200,300,400)) +
  theme_bw()

ggsave("region_specific_degenes_barplot_labeled.pdf",width = 3, height = 10)


# Median calculations
region.de.df$VISp.cl <- as.numeric(as.character(region.de.df$VISp.cl))
region.de.df$ALM.cl <- as.numeric(as.character(region.de.df$VISp.cl))

gluta.de.df <- region.de.df[region.de.df$VISp.cl > 57,]
gluta.de.df <- gluta.de.df[gluta.de.df$VISp.cl < 126,]
gluta.de.df <- gluta.de.df[!grepl("Meis2",gluta.de.df$VISp.cl_label) & !grepl("CR",gluta.de.df$VISp.cl_label),]

median(gluta.de.df$region.de.num)

# Log10 version
# p <- ggplot() +
#   # ALM Bars
#   geom_rect(data = de_data,
#             aes(xmin = alm_xmin, xmax = alm_xmax,
#                 ymin = ymin, ymax = ymax,
#                 fill = "#212021")
#   ) +
#   # ALM Labels
#   geom_text(data = de_data,
#             aes(y = y,
#                 x = alm_xmin - 0.1,
#                 label = paste(alm_label,paste0("(",region.de.ALM.num,")")),
#                 color = alm_color),
#             hjust = 1,
#             vjust = 0.3,
#             size = 2) +
#   # VISp Bars
#   geom_rect(data = de_data,
#             aes(xmin = v1_xmin, xmax = v1_xmax,
#                 ymin = ymin, ymax = ymax,
#                 fill = "#848EBC")) +
#   # VISp Labels
#   geom_text(data = de_data,
#             aes(y = y,
#                 x = v1_xmax + 0.1,
#                 label = paste(paste0("(",region.de.V1.num,")"),visp_label),
#                 color = visp_color),
#             hjust = 0,
#             vjust = 0.3,
#             size = 2) +
#   geom_vline(aes(xintercept = 0)) +
#   scale_color_identity() +
#   scale_fill_identity() +
#   scale_y_reverse() +
#   scale_x_continuous("Log10(N DE Genes)",
#                      limits = c(-5,5),
#                      breaks = c(-3,-2,-1,0,1,2,3)) +
#   theme_bw()
# 
# ggsave("region_specific_degenes_barplot.pdf",width = 3, height = 10)
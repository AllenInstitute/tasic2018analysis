library(dplyr)
library(ggplot2)
library(feather)
library(purrr)
options(stringsAsFactors = F)

# Fix pair names for compatibility with Zizhen's matrix, which uses cl

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/de.summary.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/L4.de.summary.rda")

# Only keep the comparisons between the ends and the middle.
L4.de.summary <- L4.de.summary[c("1_3","1_5","3_5"),]

names(L4.de.summary)[1:3] <- c("de.num","de.lfc","de.q.diff")

l4.de.additional <- data.frame(cl1 = c(1001, 1001, 1002), cl2 = c(1002, 1003, 1003),
                               cl1_label = c("L4 VISp Rspo1 Ctxn3+ End",
                                             "L4 VISp Rspo1 Ctxn3+ End",
                                             "L4 VISp Rspo1 Central"),
                               cl2_label = c("L4 VISp Rspo1 Central",
                                             "L4 VISp Rspo1 Shisa3+ End",
                                             "L4 VISp Rspo1 Shisa3+ End"),
                               cl.cor = 0, cl.co = 0)

L4.de.summary <- cbind(L4.de.summary, l4.de.additional) %>%
  select(one_of(names(de.summary)))

de.summary$cl1 = as.numeric(as.character(de.summary$cl1))
de.summary$cl2 = as.numeric(as.character(de.summary$cl2))
de.summary <- rbind(de.summary, L4.de.summary)

l4.anno.additional <- data.frame(cl = c(1001, 1002, 1003),
                                 cluster_id = c(1001, 1002, 1003),
                                 cluster_label = c("L4 VISp Rspo1 Ctxn3+ End",
                                                   "L4 VISp Rspo1 Central",
                                                   "L4 VISp Rspo1 Shisa3+ End"),
                                 cluster_color = c("#1241FF", "#6808FF", "#C80052"),
                                 subclass_id = c(1001, 1002, 1003),
                                 subclass_label = c("custom1","custom2","custom3"),
                                 subclass_color = c("#1241FF", "#6808FF", "#C80052"))

cl_cluster1 <- anno %>%
  filter(cluster_id %in% 1:133) %>%
  select(cl, cluster_id, cluster_label, cluster_color, subclass_id, subclass_label, subclass_color) %>%
  unique() %>%
  rbind(l4.anno.additional)

cl_cluster2 <- cl_cluster1
names(cl_cluster1) <- c("cl1","cluster_id1","cluster_label1","cluster_color1","subclass_id1","subclass_label1","subclass_color1")
names(cl_cluster2) <- c("cl2","cluster_id2","cluster_label2","cluster_color2","subclass_id2","subclass_label2","subclass_color2")

plot_data <- de.summary %>%
  mutate(cl1 = as.numeric(as.character(cl1)),
         cl2 = as.numeric(as.character(cl2))) %>%
  left_join(cl_cluster1) %>%
  left_join(cl_cluster2)

arc_poly <- function(x0, 
                     y0, 
                     start = 0,
                     end = 2*pi,
                     r = 0.05, 
                     n = 60, 
                     aspect = 1) {
  theta <- seq(start, end, length = n)
  if(aspect > 1) {
    x <- x0 + r*cos(theta)/aspect
  } else {
    x <- x0 + r*cos(theta)
  }
  
  if(aspect < 1) {
    y <- y0 + r*sin(theta)/aspect
  } else {
    y <- y0 + r*sin(theta)
  }
  
  data.frame(x = x, y = y)
}

# plot_data_left <- map_dfr(1:nrow(plot_data),
#                           function(x) {
#                             in_data <- plot_data[x,]
#                             arc_data <- arc_poly(x0 = log10(in_data$de.num),
#                                                  y0 = in_data$de.lfc,
#                                                  start = pi/2,
#                                                  end = 3*pi/2,
#                                                  r = 0.01,
#                                                  n = 30,
#                                                  aspect = 4/8)
#                             out_df <- data.frame(group = paste(in_data$cl1, in_data$cl2, sep = "_"),
#                                                  x = arc_data$x,
#                                                  y = arc_data$y,
#                                                  fill = in_data$subclass_color1)
#                           })
# 
# plot_data_right <- map_dfr(1:nrow(plot_data),
#                           function(x) {
#                             in_data <- plot_data[x,]
#                             arc_data <- arc_poly(x0 = log10(in_data$de.num),
#                                                  y0 = in_data$de.lfc,
#                                                  start = 3*pi/2,
#                                                  end = 5*pi/2,
#                                                  r = 0.01,
#                                                  n = 30,
#                                                  aspect = 4/8)
#                             out_df <- data.frame(group = paste(in_data$cl1, in_data$cl2, sep = "_"),
#                                                  x = arc_data$x,
#                                                  y = arc_data$y,
#                                                  fill = in_data$subclass_color2)
#                           })
# 
# p <- ggplot() +
#   geom_polygon(data = plot_data_left,
#                aes(x = x, y = y, 
#                    group = group, 
#                    fill = fill)) +
#   geom_polygon(data = plot_data_right,
#                aes(x = x, y = y, 
#                    group = group, 
#                    fill = fill)) +
#   theme_bw() +
#   scale_fill_identity() +
#   scale_x_continuous("log10(N DE Genes)", limits = c(-.05, 4)) +
#   scale_y_continuous("Mean(Log2 Fold Change)",limits = c(2, 10))
# 
# ggsave("pairwise_de_genes_vs_lfc_subclass.pdf",
#        p,
#        width = 8, height = 8,
#        useDingbats = F)

# simplified version

# keep only a few colors
keep_cluster_ids <- c(67,68,70,1001, 1002, 1003)
plot_data2 <- plot_data %>%
  mutate(cluster_color1 = ifelse(cluster_id1 %in% keep_cluster_ids & cluster_id2 %in% keep_cluster_ids,
                                 cluster_color1,
                                 "#808080")) %>%
  mutate(cluster_color1 = ifelse(cluster_id1 == 68 & cluster_id2 == 70,
                                 "#808080",
                                 cluster_color1))

plot_data2 %>% filter(cluster_color1 != "#808080")

p2 <- ggplot() +
  geom_point(data = plot_data2,
             aes(x = log10(de.num), y = de.lfc,
                 color = cluster_color1)) +
  theme_bw() +
  scale_color_identity() +
  scale_x_continuous("log10(N DE Genes)", limits = c(-.05, 4)) +
  scale_y_continuous("Mean(Log2 Fold Change)",limits = c(0, 10))

p2
ggsave("pairwise_de_genes_vs_lfc_simplified_colors.pdf",
       p2,
       width = 8, height = 8,
       useDingbats = F
       )

# Inset plot
p3 <- ggplot() +
  geom_point(data = plot_data2,
             aes(x = log10(de.num), y = de.lfc,
                 color = cluster_color1)) +
  theme_bw() +
  scale_color_identity() +
  scale_x_continuous("log10(N DE Genes)", expand = c(0,0), limits = c(-.1, 2.5)) +
  scale_y_continuous("Mean(Log2 Fold Change)", expand = c(0,0), limits = c(1.8, 4.4))

ggsave("pairwise_de_genes_vs_lfc_simplified_colors_inset.pdf",
       p3,
       width = 4, height = 2.2,
       useDingbats = F
)

library(ggplot2)
library(dplyr)
library(feather)
library(scrattch.vis)
library(reshape2)
library(pals)
options(stringsAsFactors = F)

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

cl_anno <- anno %>%
  select(cl, dendcluster_id, cluster_id) %>%
  unique() %>%
  filter(cluster_id %in% 1:133) %>%
  arrange(dendcluster_id)


# Pairwise Pearson Correlation
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/cl.cor.rda")
cl.cor2 <- cl.cor[as.character(cl_anno$cl),as.character(cl_anno$cl)]
rownames(cl.cor2) <- cl_anno$dendcluster_id
colnames(cl.cor2) <- cl_anno$dendcluster_id

#colorset <- sub("FF$","",plasma(20))
colorset <- jet(20)

cl.cor3 <- melt(cl.cor2) %>%
  mutate(color = values_to_colors(value, colorset = colorset))

cl.cor_plot <- ggplot() +
  geom_tile(data = cl.cor3,
            aes(x = Var1,
                y = Var2,
                fill = color)) +
  scale_fill_identity() +
  scale_y_reverse() +
  theme_void()

ggsave("cl.cor_plot.pdf",
       cl.cor_plot,
       width = 2.25,
       height = 2.25)

# Pairwise DEGene expression
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/de.summary.rda")
de.summary <- de.summary %>%
  mutate(cl1 = as.numeric(as.character(cl1)),
         cl2 = as.numeric(as.character(cl2)))

cl_anno1 <- cl_anno
names(cl_anno1) <- paste0(names(cl_anno), "1")
cl_anno2 <- cl_anno
names(cl_anno2) <- paste0(names(cl_anno), "2")

de.summary2 <- de.summary %>%
  left_join(cl_anno1) %>%
  left_join(cl_anno2) %>%
  mutate(color = values_to_colors(log10(de.num+1), colorset = colorset))

diag_tiles <- data.frame(pos = 1:125,
                         color = colorset[1])

de.summary_plot <- ggplot() +
  geom_tile(data = de.summary2,
            aes(x = dendcluster_id1,
                y = dendcluster_id2,
                fill = color)) +
  geom_tile(data = de.summary2,
            aes(x = dendcluster_id2,
                y = dendcluster_id1,
                fill = color)) +
  geom_tile(data = diag_tiles,
            aes(x = pos,
                y = pos,
                fill = color))+
  scale_fill_identity() +
  scale_y_reverse() +
  theme_void()

ggsave("de.summary_plot.pdf",
       de.summary_plot,
       width = 2.25,
       height = 2.25)

# Pairwise coclustering
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/co.stats.rda")

cl.co.ratio <- co.stats$cl.co.ratio

cl.co.ratio2 <- cl.co.ratio[as.character(cl_anno$cl),as.character(cl_anno$cl)]
rownames(cl.co.ratio2) <- cl_anno$dendcluster_id
colnames(cl.co.ratio2) <- cl_anno$dendcluster_id

cl.co.ratio3 <- melt(cl.co.ratio2) %>%
  mutate(color = values_to_colors(value, colorset = colorset))

cl.co.ratio_plot <- ggplot() +
  geom_tile(data = cl.co.ratio3,
            aes(x = Var1,
                y = Var2,
                fill = color)) +
  scale_fill_identity() +
  scale_y_reverse() +
  theme_void()

ggsave("cl.co.ratio_plot.pdf",
       cl.co.ratio_plot,
       width = 2.25,
       height = 2.25)

legend_colors <- heatmap_legend_plot(0, 1, colorset = colorset) +
  theme_void()

ggsave("legend_colors.pdf",
       legend_colors,
       width = 1, height = 0.5)

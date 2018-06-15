library(dplyr)
library(ggplot2)
library(cowplot)
library(feather)
library(dendextend)
library(scrattch.vis)
library(scrattch.io)
options(stringsAsFactors = F)

source("color_functions.R")
source("prune_leaf_custom.R")
source("layer_and_violin_functions.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")


plot_anno <- data.frame(plot_id = c(1:4),
                        cluster_id = c(35,35,46,46),
                        region_id = c(1,2,1,2),
                        plot_color = rainbow(4),
                        plot_label = c("ALM Calb2 Necab1","VISp Calb2 Necab1",
                                       "ALM Esm1","VISp Esm1"))

sub_anno <- anno %>%
  filter(cluster_id %in% c(35,46)) %>%
  left_join(plot_anno) %>%
  mutate(dendcluster_id = plot_id) %>%
  mutate(cluster_color = plot_color)


layer_jitters <- build_layer_plot(sub_anno,
                              dendcluster_ids = 1:4)

## Violin plots

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

sst_pvalb_genes <- unique(split_cst("Sst Crh Calb2"))

data <- feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/data.feather")
gene_data <- data[,c("sample_id",sst_pvalb_genes)] %>%
  filter(sample_id %in% sub_anno$sample_id) %>%
  left_join(sub_anno) %>%
  mutate(xpos = plot_id)

sst_pvalb_markers <- group_violin_plot2(data = gene_data,
                                       group_by = "plot", 
                                       clusters = 1:4,
                                       genes = sst_pvalb_genes,
                                       logscale = TRUE,
                                       labelheight = 2,
                                       max_width = 10,
                                       fontsize = 5,
                                       showcounts = FALSE)

all_plots <- plot_grid(layer_jitters,
                       sst_pvalb_markers,
                       align = "v",
                       nrow = 2,
                       rel_widths = 1,
                       rel_heights = 1)

save_plot("layers_violins.pdf",
          all_plots,
          ncol = 2,
          nrow = 2,
          base_width = 1.5/2,
          base_height = 4/2)

count_summary <- sub_anno %>%
  group_by(plot_id) %>%
  summarise(n_cells = n())

max_summary <- data.frame(sst_max = max(gene_data$Sst),
                          crh_max = max(gene_data$Crh),
                          calb2_max = max(gene_data$Calb2))

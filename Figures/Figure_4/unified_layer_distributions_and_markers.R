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

cluster_anno <- anno %>%
  select(cl, dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
  unique()

cocl_in <- read.csv("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/co.stats.csv")
cocl_in$cluster_id.x <- cluster_anno$cluster_id[match(cocl_in$cl.x,cluster_anno$cl)]
cocl_in$cluster_id.y <- cluster_anno$cluster_id[match(cocl_in$cl.y,cluster_anno$cl)]

cocl <- cocl_in %>%
  select(cluster_id.x, cluster_id.y, co.ratio) %>%
  filter(cluster_id.x != cluster_id.y)

dend <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")

gluta_layers <- build_layer_plot(anno,
                                 dend,
                                 cocl,
                                 dendcluster_ids = 1:55)

## Violin plots

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

#gluta_genes <- split_cst("Slc30a3 Cux2 Stard8 Rorb Deptor Scnn1a Rspo1 Hsd11b1 Batf3 Colq Postn Wnt7b Lemd1 Rxfp2 Oprk1 Tunar Osr1 Car3 Fam84b Chrna6 Pvalb Stac Ctxn3 Pappa2 Foxp2 Slc17a8 Trhr Tshz2 Rapgef3 Trh Gpr139 Mup5 Nxph4 Efr3a Rprm Crym")
gluta_genes <- split_cst("Slc30a3 Cux2 Stard8 Rorb Deptor Scnn1a Rspo1 Hsd11b1 Batf3 Colq Postn Wnt7b Lemd1 Rxfp2 Oprk1 Tunar Osr1 Car3 Fam84b Chrna6 Pvalb Stac Ctxn3 Pappa2 Foxp2 Slc17a8 Trhr Tshz2 Rapgef3 Trh Gpr139 Mup5 Nxph4 Efr3a Rprm Crym")
remove_genes <- c("Efr3a","Stard8","Colq","Postn","Wnt7b","Lemd1","Rxfp2","Tunar","Stac","Ctxn3","Mup5")
gluta_genes <- setdiff(gluta_genes, remove_genes)
gluta_markers <- group_violin_plot2(data_source = fdir, 
                                       group_by = "dendcluster", 
                                       clusters = 1:55,
                                       genes = gluta_genes,
                                       logscale = TRUE,
                                       labelheight = 2,
                                       max_width = 10,
                                       fontsize = 5,
                                       showcounts = FALSE)

all_plots <- plot_grid(gluta_layers,
                       gluta_markers,
                       align = "v",
                       nrow = 2,
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = c("b","c"))

save_plot("panels_bc.pdf",
          all_plots,
          ncol = 2,
          nrow = 2,
          base_width = 7.5/2,
          base_height = 6/2)

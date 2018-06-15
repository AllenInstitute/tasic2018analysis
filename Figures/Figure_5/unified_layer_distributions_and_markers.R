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

sst_pvalb_layers <- build_layer_plot(anno,
                              dend,
                              cocl,
                              dendcluster_ids = 85:115)

lamp5_sncg_vip_layers <- build_layer_plot(anno,
                                     dend,
                                     cocl,
                                     dendcluster_ids = 56:84)

## Violin plots

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

sst_pvalb_genes <- c("Sst","Chodl","Nos1","Mme","Tac1","Tacr3","Calb2","Nr2f2","Myh8","Tac2","Hpse","Crhr2","Crh","Esm1","Rxfp1","Nts","Pvalb","Gabrg1","Th","Calb1","Akr1c18","Sema3e","Gpr149","Reln","Tpbg","Cpne5","Vipr2","Nkx2-1")
#sst_pvalb_genes <- unique(split_cst("Sst Chodl Nos1 Etv1 Il1rapl2 Myh8 Chrna2 Tac2 Crhr2 Etv1 Calb2 Hpse C1ql3 Crh Nts Pvalb Gabrg1 Th Prdm8 Calb1 Reln Gpr149 Cpne5 Vipr2 Nkx2-1"))
sst_pvalb_markers <- group_violin_plot2(data_source = fdir, 
                                       group_by = "dendcluster", 
                                       clusters = 85:115,
                                       genes = sst_pvalb_genes,
                                       logscale = TRUE,
                                       labelheight = 2,
                                       max_width = 10,
                                       fontsize = 5,
                                       showcounts = FALSE)

lamp5_sncg_vip_genes <- c("Lamp5","Ndnf","Krt73","Fam19a1","Pax6","Ntn1","Plch2","Lsp1","Lhx6","Nkx2-1","Vip","Sncg","Slc17a8","Nptx2","Gpr50","Itih5","Serpinf1","Igfbp6","Gpc3","Lmo1","Ptprt","Rspo4","Chat","Crispld2","Col15a1","Pde1a")
#lamp5_sncg_vip_genes <- split_cst("Lamp5 Pax6 Ndnf Egln3 Pdlim5 Slc35d3 Vax1 Lhx6 Nkx2-1 Vip Calb2 Serpinf1 Col14a1 Sncg Ptprk Crispld2 Slc17a8 Igfbp6 Reln Gpc3 Lmo1 Cck Rspo4 Cbln4 Htr1f C1ql1 Itih5 ")
lamp5_sncg_vip_markers <- group_violin_plot2(data_source = fdir, 
                                        group_by = "dendcluster", 
                                        clusters = 56:84,
                                        genes = lamp5_sncg_vip_genes,
                                        logscale = TRUE,
                                        labelheight = 2,
                                        max_width = 10,
                                        fontsize = 5,
                                        showcounts = FALSE)

all_plots <- plot_grid(sst_pvalb_layers,
                       lamp5_sncg_vip_layers,
                       sst_pvalb_markers,
                       lamp5_sncg_vip_markers,
                       align = "v",
                       nrow = 2,
                       rel_widths = 1,
                       rel_heights = 1,
                       labels = c("c","d","e","f"))

save_plot("panels_cdef.pdf",
          all_plots,
          ncol = 2,
          nrow = 2,
          base_width = 7.5/2,
          base_height = 6/2)

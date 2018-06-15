library(dplyr)
library(ggplot2)
library(feather)
library(scrattch.vis)
library(scrattch.io)
options(stringsAsFactors = F)
source("layer_and_violin_functions.R")

genes <- split_cst("Vipr2, Pvalb, Slc32a1")

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

p <- group_violin_plot2(data_source = fdir,
                        genes = genes,
                        group_by = "dendcluster",
                        clusters = 1:133,
                        logscale = FALSE,
                        labelheight = 2,
                        max_width = 10,
                        fontsize = 6,
                        showcounts = FALSE)

ggsave("marker_violins.pdf", width = 10, height = 1, useDingbats = F)


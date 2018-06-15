library(dplyr)
library(ggplot2)
library(feather)
library(scrattch.vis)
library(scrattch.io)
options(stringsAsFactors = F)
source("layer_and_violin_functions.R")

#genes <- unique(unlist(read.csv("markers.csv", header=F)))
genes <- unique(unlist(read.table("2018-05-25_markers.txt", header=F)))

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

p <- group_violin_plot2(data_source = fdir,
                        genes = genes,
                        group_by = "dendcluster",
                        clusters = 118:133,
                        logscale = TRUE,
                        labelheight = 10,
                        max_width = 10,
                        fontsize = 6,
                        showcounts = FALSE)

p

ggsave("marker_violins.pdf", p, width = 4.1, height = 8.2, useDingbats = F)


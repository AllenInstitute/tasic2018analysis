library(dplyr)
library(scrattch)
options(stringsAsFactors = F)

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170405/"

vipr2_genes <- c("Vipr2")
vipr2_markers <- group_violin_plot(data_source = fdir, 
                   group_by = "dendcluster", 
                   clusters = 1:112,
                   logscale = FALSE,
                   genes = vipr2_genes,
                   labelheight = 50,
                   showcounts = FALSE,
                   fontsize = 10)

vipr2_markers

ggsave("vipr2_marker_violins.pdf", vipr2_markers, width = 16, height = 2, useDingbats = F)

vipr2_boxplot <- group_box_plot(data_source = fdir, 
                                   group_by = "dendcluster", 
                                   clusters = 1:112,
                                   logscale = FALSE,
                                   genes = vipr2_genes,
                                   labelheight = 50,
                                   showcounts = FALSE,
                                   fontsize = 12)

vipr2_boxplot

ggsave("vipr2_boxplot.pdf", vipr2_boxplot, width = 16, height = 2, useDingbats = F)

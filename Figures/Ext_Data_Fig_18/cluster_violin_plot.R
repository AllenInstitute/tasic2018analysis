library(dplyr)
library(ggplot2)
library(feather)
library(scrattch)
options(stringsAsFactors = F)

genes <- split_cst("Rspo1, Scnn1a, Hsd11b1, Chrna6, Slc17a8, Batf3, Colq, Fam84b, Osr1, Foxp2, Slc17a7")

p <- group_violin_plot(data_source = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170913/",
                  genes = genes,
                  group_by = "dendcluster",
                  clusters = 1:128,
                  labelheight = 10)

ggsave("marker_violins.pdf", width = 10, height = 3, useDingbats = F)


library(dplyr)
library(ggplot2)
library(feather)
library(reshape2)

setwd("C:/Users/thucn/Dropbox/AIBS/Transcriptomics/Manuscripts/V1_alm/Figures/Figure_3_Vipr2/scripts")
source("coexpression_functions.R")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180312/anno.feather")
data <- feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180312/data.feather")

## VISp
anno <- anno %>%
  group_by(dendcluster_id) %>%
  # Filter annotations for only clusters with > 10% VISp cells
  mutate(group_visp_fraction = sum(region_label == "VISp")/n()) %>%
  filter(group_visp_fraction > 0.05) %>%
  # Only keep cells from those clusters that are from VISp
  filter(region_label == "VISp") %>%
  ungroup()

# colors and ordering for different combinations
frac_anno <- data.frame(frac_label = c("frac_1","frac_2","frac_3","frac_12","frac_23","frac_13","frac_123"),
                        frac_id = 1:7,
                        frac_color = c("#eb008b", "#006838", "cyan", "black", "magenta", "gray", "blue"))

keep_frac <- c("frac_1", "frac_12", "frac_2","frac_3","frac_12","frac_23","frac_13","frac_123")
keep_frac <- c("frac_12", "frac_3","frac_12","frac_23","frac_13","frac_123")

coexpression_barplot(condition_1 = c("Vipr2"),
                     condition_2 = c("Pvalb"),
                     condition_3 = NULL,
                     frac_anno = frac_anno,
                     fpkm_cutoff = 1,
                     group_by = "cluster",
                     groups = c(1:125),
                     anno = anno,
                     keep_frac = keep_frac,
                     data = data)

ggsave("coexpression_barplot_Crh_Sst.pdf", width = 8, height = 5, dpi = 120, useDingbats = FALSE)




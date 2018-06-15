library(dplyr)
library(feather)
library(ggplot2)
options(stringsAsFactors = F)

source("sankey_functions.R")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process/WGCNA_result_anterograde/map.linnarsson.ss.df.rda")

ss.map.df <- ss.map.df %>%
  mutate(sample_id = rownames(.)) %>%
  rename(ll_label = cl)

ss.stats <- ss.map.df %>%
  group_by(ll_label) %>%
  mutate(ll_n = n()) %>%
  ungroup() %>%
  group_by(pred.cl, pred_cluster_label, ll_label) %>%
  summarise(n_cells = n(),
            ll_frac = n()/ll_n[1]) %>%
  arrange(as.numeric(pred.cl),
          -n_cells) %>%
  filter(ll_frac > 0.1) %>%
  filter(n_cells > 1)

write.csv(ss.stats,"filtered_linnarsson_mapping.csv", quote = F, row.names = F)

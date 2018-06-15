library(dplyr)
library(ggplot2)
library(feather)
library(scrattch)
library(dendextend)
options(stringsAsFactors = F)

fdir <- "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"

anno <- read_feather(file.path(fdir,"anno.feather"))

markers <- read.csv("markers_2018-05-25.csv")

Global.markers <- unique(unlist(markers$Global[markers$Global != ""]))
Inh.markers <- unique(unlist(markers$Inh[markers$Inh != ""]))
Ex.markers <- unique(unlist(markers$Ex[markers$Ex != ""]))

all_clusters <- anno %>%
  filter(cluster_id %in% 1:133) %>%
  arrange(dendcluster_id) %>%
  select(dendcluster_id) %>%
  unique() %>% unlist()

inh_clusters <- anno %>%
  filter(class_label == "GABAergic") %>%
  arrange(dendcluster_id) %>%
  select(dendcluster_id) %>%
  unique() %>% unlist()

exc_clusters <- anno %>%
  filter(class_label == "Glutamatergic") %>%
  arrange(dendcluster_id) %>%
  select(dendcluster_id) %>%
  unique() %>% unlist()

all_plot <- group_heatmap_plot(data_source = fdir,
                   genes = Global.markers,
                   group_by = "dendcluster",
                   clusters = all_clusters,
                   calculation = "trimmed_mean",
                   labelheight = 13,
                   showcounts = F)

inh_plot <- group_heatmap_plot(data_source = fdir,
                               genes = Inh.markers,
                               group_by = "dendcluster",
                               clusters = inh_clusters,
                               calculation = "trimmed_mean",
                               labelheight = 13,
                               showcounts = F)

exc_plot <- group_heatmap_plot(data_source = fdir,
                               genes = Ex.markers,
                               group_by = "dendcluster",
                               clusters = exc_clusters,
                               calculation = "trimmed_mean",
                               labelheight = 13,
                               showcounts = F)

ggsave("broad_glia_markers.pdf",all_plot,height = 8, width = 15)
ggsave("inh_markers.pdf",inh_plot,height = 8, width = 7.5)
ggsave("exc_markers.pdf",exc_plot,height = 10, width = 7.5)

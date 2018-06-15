library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

anno <- anno %>%
  filter(cluster_id %in% 1:133) %>%
  mutate(core_int_label = ifelse(core_int_label == "core",
                                 "core",
                                 "int"))

plot_data <- anno %>%
  group_by(dendcluster_id, dendcluster_label, dendcluster_color) %>%
  summarise(n_core = sum(core_int_label == "core"),
            n_int = sum(core_int_label == "int"),
            n_cells = n()) %>%
  mutate(frac_core = n_core/n_cells,
         frac_int = n_int/n_cells)

hlines <- data.frame(yintercept = seq(0,1,by=0.2))

ggplot() +
  geom_rect(data = plot_data,
            aes(xmin = dendcluster_id - 0.5,
                xmax = dendcluster_id + 0.5,
                ymin = 0,
                ymax = frac_core,
                fill = dendcluster_color)) +
  geom_rect(data = plot_data,
            aes(xmin = dendcluster_id - 0.5,
                xmax = dendcluster_id + 0.5,
                ymin = frac_core,
                ymax = 1,
                fill = "#000000")) +
  geom_text(data = plot_data,
            aes(x = dendcluster_id,
                y = -0.02,
                color = dendcluster_color,
                label = dendcluster_label),
            angle = 90,
            vjust = 0.3,
            hjust = 1,
            size = 2) +
  geom_hline(data = hlines,
             aes(yintercept = yintercept),
             linetype = "dashed",
             size = 0.1) +
  scale_fill_identity() +
  scale_color_identity() +
  theme_void() +
  scale_y_continuous(limits = c(-1, 1))

ggsave("intermediate_fraction.pdf", width = 7.5, height = 2)

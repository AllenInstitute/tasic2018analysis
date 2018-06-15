library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(feather)
options(stringsAsFactors = FALSE)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

anno <- anno %>%
  mutate(dataset = "Tasic_2017") %>%
  filter(cluster_id %in% 1:133) %>%
  mutate(category = ifelse(class_label == "GABAergic", "GABA_Tasic_2017",
                           ifelse(class_label == "Glutamatergic", "Gluta_Tasic_2017",
                                  "Non_Tasic_2017")))

data_2017 <- anno %>%
  select(sample_id, dataset, category, genes_label)

anno_2016 <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_VISp_SMV1_1679/anno.feather")
# compute genes detected
data_2016 <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_VISp_SMV1_1679/data.feather")
mat_2016 <- t(data_2016[,names(data_2016) != "sample_id"])
genes_2016 <- data.frame(sample_id = data_2016$sample_id,
                         genes_label = apply(mat_2016,2,function(x) sum(x > 0)))
anno_2016 <- left_join(anno_2016, genes_2016)

# 1:23, inh; 24:42, exc; 43:49, nn
anno_2016 <- anno_2016 %>%
  mutate(dataset = "Tasic_2016") %>%
  mutate(category = ifelse(final_id < 24, "GABA_Tasic_2016",
                           ifelse(final_id >= 24 & final_id < 43, "Gluta_Tasic_2016",
                                  "Non_Tasic_2016")))
data_2016 <- anno_2016 %>%
  select(sample_id, dataset, category, genes_label)

plot_data <- rbind(data_2016, data_2017)

ggplot() +
  geom_violin(data = plot_data,
                   aes(x = dataset,
                       y = genes_label,
                       fill = dataset),
              size = 0.1) +
  stat_summary(data = plot_data,
               aes(x = dataset,
                   y = genes_label),
               fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               color = "black",
               size = 0.2) +
  scale_fill_manual(breaks = c("Tasic_2016","Tasic_2017"), values = c("dodgerblue","skyblue")) +
  scale_y_continuous(limits = c(0, 16100), breaks = c(0, 4000, 8000, 12000, 16000)) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank())


ggsave("tasic_2016_genes_detected_all.pdf", width = 2.5, height = 3, useDingbats = F)

plot_data %>%
  group_by(dataset) %>%
  summarise(med = median(genes_label),
            n_cells = n())

ggplot() +
  geom_violin(data = plot_data,
                   aes(x = category,
                       y = genes_label,
                       fill = dataset),
                   size = 0.1) + 
  stat_summary(data = plot_data,
               aes(x = category,
                   y = genes_label),
               fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               color = "black",
               size = 0.2) +
  scale_fill_manual(breaks = c("Tasic_2016","Tasic_2017"), values = c("dodgerblue","skyblue")) +
  scale_y_continuous(limits = c(0, 16100), breaks = c(0, 4000, 8000, 12000, 16000)) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank())

ggsave("tasic_2016_genes_detected.pdf", width = 4, height = 3, useDingbats = F)

data_2017 <- anno %>%
  select(sample_id, dataset, category, genes_label, core_int_label, region_label, total_reads_label)


plot_data %>%
  group_by(category) %>%
  summarise(med = median(genes_label),
            n_cells = n())

data_2017_ci <- data_2017 %>%
  mutate(core_int_label = ifelse(core_int_label == "core",
                                 "core",
                                 "int")) %>%
  mutate(category_ci = paste(category,core_int_label, sep = "_"))

ggplot() +
  geom_violin(data = data_2017_ci,
              aes(x = category_ci,
                  y = genes_label,
                  fill = core_int_label),
              size = 0.1) +
  stat_summary(data = data_2017_ci,
               aes(x = category_ci,
                   y = genes_label),
               fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               color = "black",
               size = 0.2) +
  scale_fill_manual(breaks = c("core","int"), values = c("skyblue","gray80")) +
  scale_y_continuous(limits = c(0, 16100), breaks = c(0, 4000, 8000, 12000, 16000)) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank())

ggsave("tasic_2018_core_int_genes_detected.pdf", width = 4, height = 3, useDingbats = F)

data_2017_ci %>%
  group_by(category_ci) %>%
  summarise(med = median(genes_label),
            n_cells = n())

data_2017_reg <- data_2017 %>%
  mutate(category_reg = paste(category, region_label, sep = "_"))

ggplot() +
  geom_violin(data = data_2017_reg,
              aes(x = category_reg,
                  y = genes_label,
                  fill = region_label),
              size = 0.1) +
  stat_summary(data = data_2017_reg,
               aes(x = category_reg,
                   y = genes_label),
               fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               color = "black",
               size = 0.2) +
  scale_fill_manual(breaks = c("ALM","VISp"), values = c("#212021","#848EBC")) +
  scale_y_continuous(limits = c(0, 16100), breaks = c(0, 4000, 8000, 12000, 16000)) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank())

ggsave("tasic_2018_region_genes_detected.pdf", width = 4, height = 3, useDingbats = F)

data_2017_reg %>%
  group_by(category_reg) %>%
  summarise(med = median(genes_label),
            n_cells = n())

library(ggExtra)
library(cowplot)

p1 <- ggplot() +
  geom_point(data = data_2017_reg %>% filter(category == "GABA_Tasic_2017"),
             aes(x = total_reads_label,
                 y = genes_label,
                 color = region_label),
             size = 0.2) +
  scale_color_manual(breaks = c("ALM","VISp"), values = c("#212021","#848EBC")) +
  scale_x_continuous(limits = c(0, 5e6)) +
  theme_bw(base_size = 5) +
  theme(legend.position = "none") +
  ggtitle("GABAergic cells")

p1 <- ggMarginal(p1, groupColour = TRUE, groupFill = TRUE)

p2 <- ggplot() +
  geom_point(data = data_2017_reg %>% filter(category == "Gluta_Tasic_2017"),
             aes(x = total_reads_label,
                 y = genes_label,
                 color = region_label),
             size = 0.2) +
  scale_color_manual(breaks = c("ALM","VISp"), values = c("#212021","#848EBC")) +
  scale_x_continuous(limits = c(0, 5e6)) +
  theme_bw(base_size = 5) +
  theme(legend.position = "none") +
  ggtitle("Glutamatergic cells")

p2 <- ggMarginal(p2, groupColour = TRUE, groupFill = TRUE)

p3 <- ggplot() +
  geom_point(data = data_2017_reg %>% filter(category == "Non_Tasic_2017"),
             aes(x = total_reads_label,
                 y = genes_label,
                 color = region_label),
             size = 0.2) +
  scale_color_manual(breaks = c("ALM","VISp"), values = c("#212021","#848EBC")) +
  scale_x_continuous(limits = c(0, 5e6)) +
  theme_bw(base_size = 5) +
  theme(legend.position = "none") +
  ggtitle("Non-neuronal cells")

p3 <- ggMarginal(p3, groupColour = TRUE, groupFill = TRUE)

all_plots <- plot_grid(p1, p2, p3,
                       nrow = 1,
                       rel_widths = 1,
                       rel_heights = 1)

ggsave("v1_alm_reads_vs_genes_panel.pdf",
       all_plots,
       width = 8.1,
       height = 2.7,
       useDingbats = F)

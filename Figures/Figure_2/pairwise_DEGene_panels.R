library(cowplot)
library(dplyr)
library(ggplot2)
library(feather)
library(purrr)
library(ggExtra)
options(stringsAsFactors = F)

source("color_functions.R")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/de.summary.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/region.de.summary.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/region.de.df.rda")

#region.de.summary <- cbind(region.de.summary, region.de.df)

region.de.summary <- region.de.df


anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
class_colors <- read.csv("class_comparison_colors_simplified.csv")

cl_cluster1 <- anno %>%
  filter(cluster_id %in% 1:133) %>%
  select(cl, cluster_id, cluster_label, cluster_color, subclass_id, subclass_label, subclass_color, class_label) %>%
  unique()

cl_cluster2 <- cl_cluster1
names(cl_cluster1) <- c("cl1","cluster_id1","cluster_label1","cluster_color1","subclass_id1","subclass_label1","subclass_color1","class_label1")
names(cl_cluster2) <- c("cl2","cluster_id2","cluster_label2","cluster_color2","subclass_id2","subclass_label2","subclass_color2","class_label2")

pairwise <- de.summary %>%
  mutate(comp_type = "pairwise") %>%
  select(comp_type, cl1, cl2, de.num, de.lfc) %>%
  mutate(cl1 = as.numeric(as.character(cl1)),
         cl2 = as.numeric(as.character(cl2))) %>%
  left_join(cl_cluster1) %>%
  left_join(cl_cluster2)  %>%
  mutate(cluster_region1 = ifelse(grepl("ALM", cluster_label1), "ALM", "none"),
         cluster_region1 = ifelse(grepl("VISp",cluster_label1), "VISp", cluster_region1)) %>%
  mutate(cluster_region2 = ifelse(grepl("ALM", cluster_label2), "ALM", "none"),
         cluster_region2 = ifelse(grepl("VISp",cluster_label2), "VISp", cluster_region2)) %>%
  left_join(class_colors) %>% 
  mutate(class_type = ifelse(subclass_label1 == subclass_label2,
                             paste0("insc_",class_type),
                             class_type))

region <- region.de.summary %>%
  mutate(cl1 = as.numeric(sub("[A-Z|p]+","",VISp.cl)),
         cl2 = as.numeric(sub("[A-Z|p]+","",ALM.cl))) %>%
  mutate(comp_type = "region") %>%
  select(comp_type, cl1, cl2, de.num, de.lfc) %>%
  left_join(cl_cluster1) %>%
  left_join(cl_cluster2)  %>%
  mutate(cluster_region1 = "VISp",
         cluster_region2 = "ALM") %>%
  left_join(class_colors) %>%
  filter(class_type != "non_non")

both <- rbind(pairwise, region) 


# P2: Best matches between regions
p1_data <- both %>%
  mutate(class_type = ifelse(comp_type == "region",
                             class_type,
                             "background")) %>%
  mutate(class_color = ifelse(class_type == "background",
                              "#A7A9AC",
                              class_color)) %>%
  mutate(class_color = ifelse(class_type == "both_gluta",
                              "#116F8C",
                              class_color))

p1_meds <- p1_data %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num)) %>%
  filter(!is.na(class_type))

p1 <- ggplot() +
  geom_point(data = p1_data,
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  scale_color_identity() +
  geom_hline(data = p1_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.25) +
  geom_vline(data = p1_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.25) +
  
  scale_x_continuous(limits = c(0, 4)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw(7) +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())

p1 <- ggMarginal(p1, groupColour = TRUE, groupFill = TRUE, 
                 yparams = list(bw = 0.1),
                 xparams = list(bw = 0.04))



# p2: within vs between subclasses

p2_data <- pairwise %>%
  mutate(class_type = ifelse(class_type %in% c("insc_both_gaba","insc_both_gluta"),
                             "in_subclass",
                             class_type)) %>%
  mutate(class_type = ifelse(class_type %in% c("both_gaba","both_gluta"),
                             "out_subclass",
                             class_type))%>%
  mutate(class_type = ifelse(class_type %in% c("non_non","insc_non_non"),
                             "non_non",
                             class_type))%>%
  mutate(class_type = ifelse(!class_type %in% c("in_subclass","out_subclass","non_non"),
                             "background", class_type)) %>%
  mutate(class_color = ifelse(class_type == "in_subclass", 
                              "#009444",
                              "#A7A9AC")) %>%
  mutate(class_color = ifelse(class_type == "out_subclass",
                              "#2E3192",
                              class_color)) %>%
  mutate(class_color = ifelse(class_type == "non_non",
                              "#5e5244",
                              class_color))

p2_meds <- p2_data %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num)) %>%
  filter(!is.na(class_type))

p2 <- ggplot() +
  geom_point(data = p2_data,
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  scale_color_identity() +
  geom_hline(data = p2_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.25) +
  geom_vline(data = p2_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.25) +
  scale_x_continuous(limits = c(0, 4)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw(7) +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())

p2 <- ggMarginal(p2, groupColour = TRUE, groupFill = TRUE, 
                 yparams = list(bw = 0.1),
                 xparams = list(bw = 0.04))


# P3: between classes

p3_data <- pairwise %>%
  mutate(class_type = gsub("insc_","",class_type)) %>%
  mutate(class_type = ifelse(!class_type %in% c("both_gluta","both_gaba","non_non"),
                             class_type,
                             "background")) %>%
  mutate(class_color = ifelse(class_type == "background",
                              "#A7A9AC",
                              class_color))

p3_meds <- p3_data %>%
  group_by(class_type, class_color) %>%
  summarise(lfc_med = median(de.lfc),
            den_med = median(de.num)) %>%
  filter(!is.na(class_type))

p3 <- ggplot() +
  geom_point(data = p3_data,
             aes(x = log10(de.num + 1),
                 y = de.lfc,
                 color = class_color),
             size = 0.1) +
  scale_color_identity() +
  geom_hline(data = p3_meds,
             aes(yintercept = lfc_med,
                 color = class_color),
             linetype = "dashed",
             size = 0.25) +
  geom_vline(data = p3_meds,
             aes(xintercept = log10(den_med + 1),
                 color = class_color),
             linetype = "dashed",
             size = 0.25) +
  scale_x_continuous(limits = c(0, 4)) +
  scale_y_continuous(limits = c(0, 10)) +
  theme_bw(7) +
  theme(panel.border = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank())

p3 <- ggMarginal(p3, groupColour = TRUE, groupFill = TRUE, 
                 yparams = list(bw = 0.1),
                 xparams = list(bw = 0.04))


all_plots <- plot_grid(p1, p2, p3,
                       nrow = 1,
                       rel_widths = 1,
                       rel_heights = 1)

ggsave("de_gene_panel_v5.pdf",
       width = 8.1,
       height = 2.7,
       useDingbats = F)

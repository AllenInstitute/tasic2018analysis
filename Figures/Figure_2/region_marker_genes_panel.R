library(dplyr)
library(feather)
library(ggplot2)
library(ggrepel)
options(stringsAsFactors = F)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/V1.ALM.diff.gene.rda")
exc.gene.df <- gene.df %>%
  mutate(class = "exc")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/V1.ALM.diff.inh.gene.rda")
inh.gene.df <- gene.df %>%
  mutate(class = "inh")

# Find DE genes for Inhibitory types
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/region.de.df.rda")
region.de.df[] <- lapply(region.de.df, function(x) ifelse(is.factor(x), as.character(x), x))

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

gene.df <- rbind(exc.gene.df, inh.gene.df) %>%
  filter(!grepl("^LOC",gene),
         !grepl("Rik$",gene))

alm_exc_label <- gene.df %>%
  filter(class == "exc") %>%
  arrange(-prop.diff) %>%
  head(25) %>%
  mutate(color = ifelse(class == "exc", "#212021", "#BE1E2D"),
         fill = "#E5E1E5") 

visp_exc_label <- gene.df %>%
  filter(class == "exc") %>%
  arrange(prop.diff) %>%
  head(25) %>%
  mutate(color = ifelse(class == "exc", "#848EBC", "#BE1E2D"),
         fill = "#D0D7EF") 

alm_inh_label <- gene.df %>%
  filter(class == "inh") %>%
  arrange(-prop.diff) %>%
  head(10) %>%
  mutate(color = ifelse(class == "exc", "#212021", "#DD6091"),
         fill = "#E5E1E5") 

visp_inh_label <- gene.df %>%
  filter(class == "inh") %>%
  arrange(prop.diff) %>%
  head(10) %>%
  mutate(color = ifelse(class == "exc", "#848EBC", "#DD6091"),
         fill = "#D0D7EF") 


all_labels <- rbind(alm_exc_label, visp_exc_label, alm_inh_label, visp_inh_label)
#all_labels <- rbind(alm_label, visp_label, inh_label)

no_label <- gene.df %>%
  filter(!gene %in% all_labels$gene) %>%
  mutate(color = ifelse(class == "exc", "#808080", "#EA9BA1"))


marker_plot <- ggplot() +
  geom_point(data = no_label,
             aes(x = -prop.diff,
                 y = log10(prop.max),
                 color = color),
             #alpha = 1,
             size = 0.2) +
  geom_point(data = all_labels,
             aes(x = -prop.diff,
                 y = log10(prop.max),
                 color = color),
             size = 1)+
  geom_label_repel(data = all_labels,
                   aes(x = -prop.diff,
                       y = log10(prop.max),
                       label = gene,
                       fill = fill,
                       color = color),
                   size = 2,
                   min.segment.length = 0) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-1.5, 1.2), 
                     breaks = seq(-1,1,by=0.25), 
                     labels = c(1,0.75,0.5,0.25,0,0.25,0.5,0.75,1)) +
  scale_y_continuous(limits = c(-3,0.5),
                     breaks = log10(c(0.001,0.01,0.1,0.25,0.5,0.75,1)), 
                     labels = c(0.001,0.01, 0.1, 0.25, 0.5, 0.75, 1)) +
  theme_bw(7) +
  theme(panel.grid.minor = element_blank())

marker_plot

ggsave("marker_plot.pdf",marker_plot, width = 7.5, height = 3.5, useDingbats = FALSE)

# One layer, with marginals
library(ggExtra)
all_points <- rbind(all_labels, no_label %>% mutate(fill = color))

margin_points <- all_points %>%
  mutate(color = ifelse(class == "exc", "#808080", "#EA9BA1"))

marker_plot2 <- ggplot() +
  geom_point(data = margin_points,
             aes(x = -prop.diff,
                 y = log10(prop.max),
                 color = color),
             size = 0.5) +
  geom_point(data = all_labels,
             aes(x = -prop.diff,
                 y = log10(prop.max),
                 color = color),
             size = 0.5) +
  geom_label_repel(data = all_labels,
                   aes(x = -prop.diff,
                       y = log10(prop.max),
                       label = gene,
                       fill = fill,
                       color = color),
                   size = 2) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(limits = c(-1.5, 1.2), 
                     breaks = seq(-1,1,by=0.25), 
                     labels = c(1,0.75,0.5,0.25,0,0.25,0.5,0.75,1)) +
  scale_y_continuous(limits = c(-3,0.5),
                     breaks = log10(c(0.001,0.01,0.1,0.25,0.5,0.75,1)), 
                     labels = c(0.001,0.01, 0.1, 0.25, 0.5, 0.75, 1)) +
  theme_bw(7) +
  theme(panel.grid.minor = element_blank())

marker_plot2 <- ggMarginal(marker_plot2, 
           margins = "x", 
           size = 10,
           groupColour = TRUE, 
           groupFill = TRUE, 
           bw = 0.05)

ggsave("marker_plot_marginals.pdf",marker_plot2, width = 7.5, height = 3.5, useDingbats = FALSE)

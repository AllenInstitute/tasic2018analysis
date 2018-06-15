library(dplyr)
library(ggplot2)
library(cowplot)
options(stringsAsFactors = F)

load("gradient.df.list.rda")

load("L4 IT VISp Rspo1 L5 IT VISp Batf3.rda")
L4_vs_Batf3 <- df %>%
  mutate(sample_id = rownames(.))
load("L4 IT VISp Rspo1 L5 IT VISp Hsd11b1 Endou.rda")
L4_vs_Hsd11b1 <- df %>%
  mutate(sample_id = rownames(.)) %>%
  mutate(prob = 1-prob) %>%
  mutate(eigen = -1 * eigen)
load("L4.df.rda")
L4_vs_L4 <- L4.df %>%
  mutate(sample_id = rownames(.)) %>%
  mutate(color = "#000000",
         color = ifelse(prob < 0.2, "#1241FF" , color),
         color = ifelse(prob > 0.4 & prob < 0.6, "#6808FF" , color),
         color = ifelse(prob > 0.8, "#C80052", color)) %>%
  mutate(eigen = -1 * eigen)

this_theme <- theme_bw(4) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#EBEBEB", size = 0.25),
        panel.border = element_blank(),
        axis.ticks = element_blank())

L4_vs_Batf3_points <- ggplot() +
  geom_point(data = L4_vs_Batf3,
             aes(x = prob,
                 y = eigen),
             size = 0.2) +
  ylab("") +
  scale_x_continuous("",
                     limits = c(-0.05, 1.05),
                     expand = c(0,0)) +
  this_theme

L4_vs_Batf3_ks_pval <- 0e0

L4_vs_Batf3_sampled <- L4_vs_Batf3 %>%
  filter(sample_id %in% rownames(df.list[["69_72"]]))

L4_vs_Batf3_bars <- ggplot(L4_vs_Batf3_sampled, aes(prob)) + 
  geom_histogram(aes(y=..count../sum(..count..)),
                 breaks = seq(0, 1 , by = 0.05),
                 binwidth = 0.5,
                 color = NA,
                 size = 0.25) +
  geom_hline(yintercept = 0.05,
             size = 0.5) +
  ylim(0, 0.3) +
  ylab("") +
  scale_x_continuous("Assignment probability",
                     limits = c(-0.05, 1.05),
                     expand = c(0,0)) +  
  this_theme

L4_vs_Hsd11b1_points <- ggplot() +
  geom_point(data = L4_vs_Hsd11b1,
             aes(x = prob,
                 y = eigen),
             size = 0.2) +
  ylab("") +
  scale_x_continuous("",
                     limits = c(-0.05, 1.05),
                     expand = c(0,0)) +  
  this_theme

L4_vs_Hsd11b1_ks_pval <- 4.84e-13

L4_vs_Hsd11b1_sampled <- L4_vs_Hsd11b1 %>%
  filter(sample_id %in% rownames(df.list[["69_70"]]))

L4_vs_Hsd11b1_bars <- ggplot(L4_vs_Hsd11b1_sampled, aes(prob)) + 
  geom_histogram(aes(y=..count../sum(..count..)),
                 breaks = seq(0, 1 , by = 0.05),
                 binwidth = 0.5,
                 color = NA,
                 size = 0.25) +
  geom_hline(yintercept = 0.05,
             size = 0.5) +
  ylim(0, 0.3) +
  ylab("") +
  scale_x_continuous("Assignment probability",
                     limits = c(-0.05, 1.05),
                     expand = c(0,0)) +  
  this_theme

L4_vs_L4_points <- ggplot() +
  geom_point(data = L4_vs_L4,
             aes(x = prob,
                 y = eigen,
                 color = color),
             size = 0.2) +
  scale_color_identity() +
  ylab("Eigengenvector position") +
  xlab("") +
  this_theme

L4_vs_L4_ks_pval <- 4.08e-3

L4_vs_L4_sampled <- L4_vs_L4 %>%
  filter(sample_id %in% rownames(df.list[["58"]]))

L4_vs_L4_bars <- ggplot(L4_vs_L4_sampled, aes(prob)) + 
  geom_histogram(aes(y=..count../sum(..count..)),
                 breaks = seq(0, 1 , by = 0.05),
                 binwidth = 0.5,
                 color = NA,
                 size = 0.25) +
  geom_hline(yintercept = 0.05,
             size = 0.5) +
  ylim(0, 0.3) +
  ylab("Fraction of Cells") +
  scale_x_continuous("Assignment probability",
                     limits = c(-0.05, 1.05),
                     expand = c(0,0)) +  
  this_theme

all_plots <- plot_grid(L4_vs_L4_points,
                       L4_vs_Hsd11b1_points,
                       L4_vs_Batf3_points,
                       L4_vs_L4_bars,
                       L4_vs_Hsd11b1_bars,
                       L4_vs_Batf3_bars,
                       align = "v",
                       nrow = 2,
                       rel_widths = 1,
                       rel_heights = c(1, 0.5))

save_plot("gradient_panels.pdf",
          all_plots,
          ncol = 3,
          nrow = 2,
          base_width = 3.5/3,
          base_height = 1,
          useDingbats = F)

library(feather)
library(ggplot2)
library(reshape2)
library(scrattch.vis)
library(viridisLite)
options(stringsAsFactors = F)

load("L4.df.rda")

genes <- c("Shisa3","Tmem215","Gpr88","Slc35g2","Rab3b","Scnn1a","Thsd7a","Nrsn2","Endou","Ccdc3",
           "Col8a1","Rassf2","Nog","Fmn1","Phactr2","Cnn3","Rreb1","Tshz2",
           "9530059O14Rik","Kitl","Cnr1","Pde1a","LOC105245781","Clec18a","St6galnac5","Galnt14","Doc2a","Xkr6","Lmo3","Stac2","Adam33","Chrm3","Ctxn3")

gene_data <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/data.feather",
                          columns = c("sample_id",genes))

l4_data <- L4.df %>%
  mutate(sample_id = rownames(.)) %>%
  left_join(gene_data) %>%
  arrange(prob,eigen) %>%
  mutate(xpos = 1:n()) %>%
  mutate(color = values_to_colors(prob, colorset = viridis(10)))

first_instance <- function(x, vals) {
  y <- numeric()
  for(val in vals) {
    y <- c(y, min(which(x >= val)))
  }
  y
}

x_breaks <- data.frame(val = seq(0, 1, by = 0.2),
                       xpos = first_instance(l4_data$prob, seq(0,1, by = 0.2)))

gene_ypos <- data.frame(gene = genes,
                        ypos = 1:length(genes))

l4_plot_data <- melt(l4_data[,c("xpos",genes)], "xpos") %>%
  mutate(color = values_to_colors(log10(value + 1))) %>%
  rename(gene = variable) %>%
  left_join(gene_ypos)

l4_heatmap <- ggplot() +
  geom_tile(data = l4_plot_data,
            aes(x = xpos,
                y = ypos,
                fill = color)) +
  geom_tile(data = l4_data,
            aes(x = xpos,
                y = -1,
                fill = color)) +
  scale_fill_identity() +
  scale_y_continuous("", expand = c(0,0),
                   breaks = 1:length(genes),
                   labels = genes) +
  scale_x_continuous("", expand = c(0,0),
                     breaks = x_breaks$xpos,
                     labels = x_breaks$val) +
  theme_bw(4) +
  theme(axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

ggsave("l4_heatmap.pdf",
       l4_heatmap,
       width = 2,
       height = 2)

library(dendextend)
library(dplyr)
library(feather)
library(ggplot2)
options(stringsAsFactors = F)

values_to_colors <- function(x, minval = NULL, maxval = NULL, colorset = c("darkblue","dodgerblue","gray80","orangered","red")) {
  
  heat_colors <- colorRampPalette(colorset)(1001)
  
  if(is.null(maxval)) {
    maxval <- max(x)
  }
  if (is.null(minval)) {
    minval <- min(x)
  }
  
  heat_positions <- unlist(round((x - minval) / (maxval - minval) * 1000 + 1, 0))
  
  colors <- heat_colors[heat_positions]
  
  colors
}

color_sum <- function(col1,col2) {
  
  rgbmat1 <- col2rgb(col1)/255
  rgbmat2 <- col2rgb(col2)/255
  
  mix <- rgbmat1 + rgbmat2
  
  rgb(mix[1],mix[2],mix[3])
  
}


# Load Zizhen's V1 to ALM comparisons
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/V1.ALM.compare.cor.rda")

# # Rotate the L5a types so that Ucma_2 is closer to L2/3
# rotate_labels <- labels(V1.dend)
# rotate_labels[4:12] <- rotate_labels[c(9:12,4:8)]
# V1.dend <- V1.dend %>% rotate(rotate_labels)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
cluster_anno <- anno %>% select(cl, cluster_id, cluster_label, cluster_color) %>% unique

labels(V1.dend) <- cluster_anno$cluster_label[match(labels(V1.dend), cluster_anno$cl)]
labels(ALM.dend) <- cluster_anno$cluster_label[match(labels(ALM.dend), cluster_anno$cl)]

comb.map$v1_label <- cluster_anno$cluster_label[match(comb.map$V1.cl, cluster_anno$cl)]
comb.map$alm_label <- cluster_anno$cluster_label[match(comb.map$ALM.cl, cluster_anno$cl)]

v1_ggdend <- as.ggdend(V1.dend)
alm_ggdend <- as.ggdend(ALM.dend)

# Invert and center dendrograms
# add 0.3 space between them
# 0.1 for leaves
# 0.1 on each side for labels

# V1: right-to-left
v1_seg <- v1_ggdend$segments
names(v1_seg)[1:4] <- c("y","x","yend","xend")
ymid <- max(v1_seg$y/2)
ymax <- max(v1_seg$y)
v1_seg <- v1_seg %>%
  mutate(y = (y - 2*(y - ymid))/ymax,
         yend = (yend - 2*(yend - ymid))/ymax,
         x = x + 0.15,
         xend = xend + 0.15)

v1_leaves <- v1_ggdend$labels
names(v1_leaves)[1:3] <- c("v1_y","v1_x","v1_label")
v1_leaves <- v1_leaves %>%
  mutate(v1_y = (v1_y - 2*(v1_y - ymid))/ymax,
         v1_x = v1_x + 0.05) %>%
  select(-col,-cex)

# ALM: left-to-right
alm_seg <- alm_ggdend$segments
names(alm_seg)[1:4] <- c("y","x","yend","xend")
ymid <- max(alm_seg$y/2)
ymax <- max(alm_seg$y)
alm_seg <- alm_seg %>%
  mutate(y = (y - 2*(y - ymid))/ymax,
         yend = (yend - 2*(yend - ymid))/ymax,
         x = -x - 0.15,
         xend = -xend - 0.15)

alm_leaves <- alm_ggdend$labels
names(alm_leaves)[1:3] <- c("alm_y","alm_x","alm_label")
alm_leaves <- alm_leaves %>%
  mutate(alm_y = (alm_y - 2*(alm_y - ymid))/ymax,
         alm_x = -alm_x - 0.05) %>%
  select(-col,-cex)

# build comparison segments
comb.map <- comb.map %>%
  mutate(Prob = as.numeric(Prob)) %>%
  left_join(alm_leaves) %>%
  left_join(v1_leaves)

v1_to_alm <- comb.map %>%
  filter(type == "V1.ref.ALM")
alm_to_v1 <- comb.map %>%
  filter(type == "ALM.ref.V1")

alm_v1_comp_plot <- ggplot() +
  geom_segment(data = v1_seg,
               aes(x = x, xend = xend,
                   y = y, yend = yend),
               lineend = "square") +
  geom_segment(data = alm_seg,
               aes(x = x, xend = xend,
                   y = y, yend = yend),
               lineend = "square") +
  geom_curve(data = alm_to_v1,
             aes(x = alm_x, xend = v1_x,
                 y = alm_y, yend = v1_y,
                 size = Prob),
             color = "#212021",
             curvature = -0.1,
             alpha = 0.5) + #arrow = arrow(length = unit(0.25,"cm"))) +
  geom_curve(data = v1_to_alm,
             aes(x = v1_x, xend = alm_x,
                 y = v1_y, yend = alm_y,
                 size = Prob),
             color = "#848EBC",
             curvature = -0.1,
             alpha = 0.5) +#, arrow = arrow(length = unit(0.25,"cm"))) +
  geom_text(data = alm_leaves,
            aes(x = alm_x - 0.01, y = alm_y,
                label = alm_label),
            hjust = 1,
            vjust = 0.3) +
  geom_text(data = v1_leaves,
            aes(x = v1_x + 0.01, y = v1_y,
                label = v1_label),
            hjust = 0,
            vjust = 0.3) +
  theme_void() +
  scale_size_continuous(range = c(0,3))

alm_v1_comp_plot

ggsave("v1_alm_crossmapping_cor.pdf",alm_v1_comp_plot,width = 12, height = 8)

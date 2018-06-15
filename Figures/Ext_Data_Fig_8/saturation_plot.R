library(dplyr)
library(ggplot2)
library(feather)
#library(scrattch)

color_mean <- function(x) {
  library(grDevices)
  
  rgb_x <- col2rgb(x)
  rgb_mean <- rowMeans(rgb_x)
  new_hex <- rgb(rgb_mean["red"]/255,
                 rgb_mean["green"]/255,
                 rgb_mean["blue"]/255)
  
  new_hex
  
}

ramp_colors <- function (x, minval = NULL, maxval = NULL, colorset = c("darkblue", 
                                                        "dodgerblue", "gray80", "orangered", "red")) 
{
  heat_colors <- colorRampPalette(colorset)(1001)
  if (is.null(maxval)) {
    maxval <- max(x)
  } else {
    x[x > maxval] <- maxval
  }
  if (is.null(minval)) {
    minval <- min(x)
  } else {
    x[x < minval] <- minval
  }
  heat_positions <- unlist(round((x - minval)/(maxval - minval) * 
                                   1000 + 1, 0))
  colors <- heat_colors[heat_positions]
  colors
}

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/saturation.df.rda")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

anno <- anno %>%
  filter(dendcluster_id < 134)

sat_xpos <- saturation.df %>%
  group_by(sample_size) %>%
  summarise(n_cl = sum(merge.cl == "nochange") + length(unique(merge.cl)) - 1) %>%
  mutate(sample_size = as.numeric(as.character(sample_size))) %>%
  arrange(desc(sample_size)) %>%
  mutate(xpos = 1:n() + 1) %>%
  mutate(xmin = xpos - 0.3,
         xmax = xpos + 0.3)

dendcluster_ypos <- anno %>%
  group_by(cl, dendcluster_id, dendcluster_label, dendcluster_color) %>%
  summarise(cl_size = n()) %>%
  ungroup() %>%
  arrange(-dendcluster_id) %>%
  mutate(ypos = 1:n()) %>%
  mutate(ymin = ypos - 0.5,
         ymax = ypos + 0.5)

full_pos <- dendcluster_ypos %>%
  mutate(xpos = 1,
         xmin = xpos - 0.3,
         xmax = xpos + 0.3) %>%
  mutate(cl_size = ifelse(is.na(cl_size), 0, cl_size)) %>%
  mutate(fill = ifelse(is.na(cl_size), "#FFFFFF", ramp_colors(cl_size, 0, 20, colorset = c("orange","white","skyblue"))))

sat_pos <- saturation.df %>%
  mutate(sample_size = as.numeric(as.character(sample_size))) %>%
  mutate(org.cl = as.numeric(as.character(org.cl))) %>%
  left_join(dendcluster_ypos, by = c("org.cl" = "cl")) %>%
  left_join(sat_xpos, by = "sample_size") %>%
  mutate(dendcluster_color = ifelse(merge.cl == "absent",NA,dendcluster_color)) %>%
  mutate(cl_size.x = ifelse(is.na(cl_size.x), 0, cl_size.x)) %>%
  mutate(outline = ifelse(merge.cl %in% c("absent","nochange"), NA, "#000000"),
         fill = ifelse(is.na(cl_size.x), "#FFFFFF", ramp_colors(cl_size.x, 0, 20, colorset = c("orange","white","skyblue"))))

sat_unmerged <- sat_pos %>%
  filter(merge.cl %in% c("absent","nochange"))

sat_merged <- sat_pos %>%
  filter(!merge.cl %in% c("absent","nochange")) %>%
  group_by(sample_size, merge.cl) %>%
  mutate(dendcluster_color = color_mean(dendcluster_color))

sat_curves <- sat_merged %>%
  group_by(sample_size, merge.cl) %>%
  arrange(ypos) %>%
  mutate(x = xmax, xend = xmax,
         y = ypos, yend = lead(ypos)) %>%
  
  ungroup() %>%
  filter(complete.cases(.))

merge_plot <- ggplot() +
  geom_rect(data = full_pos,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = fill,
                color = "white"),
            size = 0.1) +
  geom_text(data = full_pos,
            aes(x = xmax, y = ymax,
                label = cl_size),
            hjust = 0, vjust = 1,
            size = 2) +
  geom_text(data = dendcluster_ypos,
            aes(x = 0.5, y = ypos,
                label = dendcluster_label,
                color = dendcluster_color),
            hjust = 0, vjust = 0.3,
            size = 2*5/6) +
  geom_rect(data = sat_unmerged,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax,
                fill = fill,
                color = "white"),
            size = 0.1) +
  geom_text(data = sat_unmerged,
            aes(x = xmax, y = ymax,
                label = cl_size.x),
            hjust = 0, vjust = 1,
            size = 2*5/6) +
  geom_rect(data = sat_merged,
            aes(xmin = xmin + 0.2, xmax = xmax + 0.2,
                ymin = ymin, ymax = ymax,
                fill = fill,
                color = outline)) +
  geom_text(data = sat_merged,
            aes(x = xmax + 0.2, y = ymax,
                label = cl_size.x),
            hjust = 0, vjust = 1,
            size = 2*5/6) +
  geom_curve(data = sat_curves,
             aes(x = x + 0.2, xend = xend + 0.2,
                 y = y, yend = yend),
             curvature = -1) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_reverse(limits = c(12, 0)) +
  theme_void()

merge_plot

ggsave("merge_plot_numbers.pdf",merge_plot,width = 7.5, height = 10)

cluster_counts <- saturation.df %>%
  group_by(sample_size) %>%
  summarise(n_clusters = length(unique(as.character(merge.cl))) - 1 + sum(merge.cl == "nochange"))

cluster_counts

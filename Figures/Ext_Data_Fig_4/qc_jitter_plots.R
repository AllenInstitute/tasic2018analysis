library(dplyr)
library(ggplot2)
library(ggbeeswarm)
library(feather)

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE, roundall = F) {
  require(dplyr)
  # This does the summary. For each group return a vector with
  # N, mean, and sd
  
  names(data)[names(data) == measurevar] <- "measurevar"
  
  datac <- data %>%
    select(one_of(groupvars,"measurevar")) %>%
    filter(ifelse(na.rm == T, !is.na(measurevar), T)) %>%
    mutate(measurevar = as.numeric(measurevar)) %>%
    group_by_(c(groupvars)) %>%
    summarise(N = n(),
              median = median(measurevar),
              mean = mean(measurevar),
              max = max(measurevar),
              sd = ifelse(N == 1, 0, sd(measurevar)),
              q25 = as.numeric(quantile(measurevar, 0.25)),
              q75 = as.numeric(quantile(measurevar, 0.75))) %>%
    mutate(se = sd/sqrt(N))
  #%>%
  #  mutate(ci =  se * qt(conf.interval/2 + 0.5, N-1))
  
  
  if(roundall) {
    roundcols <- c("median","mean","max","sd","q25","q75","se","ci")
    datac[roundcols] <- round(datac[roundcols],3)
  }
  
  # datac <- datac %>%
  #   mutate(xpos = 1:n())
  
  return(datac)
}

#load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/V1_ALM_R_Objects/STAR2.0/V1_ALM_Star2.0_29002_samp.dat.Rdata")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

anno <- anno %>%
  filter(cluster_id %in% 1:133)

dendcluster_anno <- anno %>%
  select(dendcluster_id, dendcluster_label, dendcluster_color) %>%
  unique()

# dendcluster_id is the annotation for cluster ordering based on the current, bootstrapped dendrogram
genes_stats <- summarySE(data = anno,
                         measurevar = "genes_label",
                         groupvars = "dendcluster_id")


genes_plot <- ggplot() +
  # geom_quasirandom from the ggbeeswarm package
  # makes violin-shaped jittered point plots
  geom_quasirandom(data = anno,
                   aes(x = dendcluster_id,
                       y = genes_label),
                   color = "skyblue",
                   # Need to set position_jitter height = 0 to prevent
                   # jitter on the y-axis, which changes data representation
                   position = position_jitter(width = .3,height = 0),
                   size = 0.1) +
  # Errorbars built using genes_stats values
  geom_errorbar(data = genes_stats,
                aes(x = dendcluster_id,
                    ymin = q25,
                    ymax = q75),
                size = 0.2) +
  # Median points from genes_stats
  geom_point(data = genes_stats,
             aes(x = dendcluster_id,
                 y = median),
             color = "red",
             size = 0.5) +
  # Cluster labels as text objects
  geom_text(data = dendcluster_anno,
            aes(x = dendcluster_id,
                y = 0,
                label = dendcluster_label,
                color = dendcluster_color),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  # Median values next to cluster labels, since there's space there.
  geom_text(data = genes_stats,
            aes(x = dendcluster_id,
                y = 100,
                label = median),
            angle = 90,
            hjust = 0,
            vjust = 0.3,
            size = 2*5/6) +
  scale_color_identity() +
  # Expand the y scale so that the labels are visible
  scale_y_continuous("Genes Detected (thousands)",
                     limits = c(-3000, 16000), 
                     breaks = seq(0, 16000, 4000),
                     labels = seq(0,16,4)) +
  # Remove X-axis title
  scale_x_continuous("") +
  theme_bw() +
  # Theme tuning
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave("gene_detection_plot.pdf", genes_plot, width = 8, height = 3, useDingbats = F)


# dendcluster_id is the annotation for cluster ordering based on the current, bootstrapped dendrogram
reads_stats <- summarySE(data = anno,
                         measurevar = "total_reads_label",
                         groupvars = "dendcluster_id")


reads_plot <- ggplot() +
  # geom_quasirandom from the ggbeeswarm package
  # makes violin-shaped jittered point plots
  geom_quasirandom(data = anno,
                   aes(x = dendcluster_id,
                       y = total_reads_label/1e6),
                   color = "skyblue",
                   # Need to set position_jitter height = 0 to prevent
                   # jitter on the y-axis, which changes data representation
                   position = position_jitter(width = .3,height = 0),
                   size = 0.1) +
  # Errorbars built using genes_stats values
  geom_errorbar(data = reads_stats,
                aes(x = dendcluster_id,
                    ymin = q25/1e6,
                    ymax = q75/1e6),
                size = 0.2) +
  # Median points from genes_stats
  geom_point(data = reads_stats,
             aes(x = dendcluster_id,
                 y = median/1e6),
             color = "red",
             size = 0.5) +
  # Cluster labels as text objects
  geom_text(data = dendcluster_anno,
            aes(x = dendcluster_id,
                y = -0.5,
                label = dendcluster_label,
                color = dendcluster_color),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  # Median values next to cluster labels, since there's space there.
  geom_text(data = reads_stats,
            aes(x = dendcluster_id,
                y = -0.2,
                label = signif(median/1e6,4)),
            angle = 90,
            hjust = 0,
            vjust = 0.3,
            size = 2*5/6) +
  scale_color_identity() +
  # Expand the y scale so that the labels are visible
  scale_y_continuous("Total Sequenced Reads (millions)",
                     limits = c(-2,5),
                     breaks = c(0,1,2,3,4,5)) +
  # Remove X-axis title
  scale_x_continuous("") +
  theme_bw() +
  # Theme tuning
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

reads_plot

ggsave("total_reads_plot.pdf", reads_plot, width = 8, height = 3, useDingbats = F)

# Confusion score (for figure S3)

confusion_stats <- summarySE(data = anno,
                         measurevar = "confusion_label",
                         groupvars = "dendcluster_id")


confusion_plot <- ggplot() +
  # geom_quasirandom from the ggbeeswarm package
  # makes violin-shaped jittered point plots
  geom_quasirandom(data = anno,
                   aes(x = dendcluster_id,
                       y = confusion_label),
                   color = "skyblue",
                   # Need to set position_jitter height = 0 to prevent
                   # jitter on the y-axis, which changes data representation
                   position = position_jitter(width = .3,height = 0),
                   size = 0.1) +
  # Errorbars built using genes_stats values
  geom_errorbar(data = confusion_stats,
                aes(x = dendcluster_id,
                    ymin = q25,
                    ymax = q75),
                size = 0.2) +
  # Median points from genes_stats
  geom_point(data = confusion_stats,
             aes(x = dendcluster_id,
                 y = median),
             color = "red",
             size = 0.5) +
  # Cluster labels as text objects
  geom_text(data = dendcluster_anno,
            aes(x = dendcluster_id,
                y = -0.22,
                label = dendcluster_label,
                color = dendcluster_color),
            angle = 90,
            hjust = 1,
            vjust = 0.3,
            size = 2*5/6) +
  # Median values next to cluster labels, since there's space there.
  geom_text(data = confusion_stats,
            aes(x = dendcluster_id,
                y = -0.2,
                label = round(median,3)),
            angle = 90,
            hjust = 0,
            vjust = 0.3,
            size = 2*5/6) +
  scale_color_identity() +
  # Expand the y scale so that the labels are visible
  scale_y_continuous("Confusion Score",
                     limits = c(-0.5,1.25),
                     breaks = c(0,0.25,0.5,0.75,1,1.25)) +
  # Remove X-axis title
  scale_x_continuous("") +
  theme_bw() +
  # Theme tuning
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

confusion_plot

ggsave("confusion_plot.pdf", confusion_plot, width = 8, height = 2, useDingbats = F)

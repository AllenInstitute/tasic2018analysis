library(dendextend)
library(dplyr)
library(feather)
options(stringsAsFactors = F)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/dend.collapse.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/dend.bak.rda")

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

cl_anno <- anno %>%
  select(cl, cluster_id, cluster_label, cluster_color) %>%
  unique()

cl_to_cluster <- cl_anno %>%
  select(cl, cluster_id)

format_dend <- function(dend, cl_anno, cl_to_cluster) {
  # Rotate so that Glutamatergic samples are first:
  original_labels <- as.numeric(labels(dend))
  gaba_labels <- cl_to_cluster$cl[cl_to_cluster$cluster_id %in% 1:60]
  gluta_labels <- cl_to_cluster$cl[cl_to_cluster$cluster_id %in% 61:115]
  nn_labels <- cl_to_cluster$cl[cl_to_cluster$cluster_id %in% 116:133]
  
  dend <- dend %>%
    rotate(as.character(c(gluta_labels, gaba_labels, nn_labels)))
  
  labels_colors(dend) <- cl_anno$cluster_color[match(labels(dend), cl_anno$cl)]
  labels(dend) <- as.character(cl_anno$cluster_label[match(as.numeric(labels(dend)), cl_anno$cl)])
  
  labelDend <- function(dend,n=1)
  {
    if(is.null(attr(dend,"label"))){
      attr(dend, "label") =paste0("n",n)
      n= n +1
    }
    if(length(dend)>1){
      for(i in 1:length(dend)){
        tmp = labelDend(dend[[i]], n)
        dend[[i]] = tmp[[1]]
        n = tmp[[2]]
      }
    }
    return(list(dend, n))
  }
  
  
  dend <- labelDend(dend)[[1]]
  dend
}

# Original, uncollapsed dendrogram
dend_original <- format_dend(dend,
                             cl_anno, cl_to_cluster)

pdf("dend_0.pdf",
    width = 7.5, height = 4,
    useDingbats = F)
plot(dend_original)
dev.off()

# version used for figures
dend_0.4 <- readRDS("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/dend.RData")
pdf("dend_0.4.pdf",
    width = 7.5, height = 4,
    useDingbats = F)
plot(dend_0.4)
dev.off()

# Collapsed further
dend_0.8 <- format_dend(dend.collapse[[7]],
                        cl_anno, cl_to_cluster)
pdf("dend_0.8.pdf",
    width = 7.5, height = 4,
    useDingbats = F)
plot(dend_0.8)
dev.off()

# legend
dend_dfs <- as.ggdend(dend)

library(scrattch.vis)

colorset <- c("#FFFFFF","#C5C5C5","#921C1C","#000000")

legend <- heatmap_legend_plot(minval = 0, maxval = 1,
                    scale_name = "Bootstrapped Confidence",
                    colorset = colorset)

ggsave("legend.pdf",
       legend,
       width = 1,
       height = 1)

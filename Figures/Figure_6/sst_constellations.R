library(ggplot2)
library(dplyr)
source("plotting_functions.R")



cell.cl.map.df_to_transition.df.comb <- function(cell.cl.map.df, cl.df) {
  transition.df = with(cell.cl.map.df, as.data.frame(table(org.cl, transition.cl)))
  transition.df = transition.df[transition.df$Freq > 0,]
  transition.df$org.cl = as.character(transition.df$org.cl)
  transition.df$transition.cl = as.character(transition.df$transition.cl)
  
  ###combine transitions from both directions
  transition.df$cl.min = pmin(transition.df$org.cl, transition.df$transition.cl)
  transition.df$cl.max = pmax(transition.df$org.cl, transition.df$transition.cl)
  transition.df$cl.pair = paste(transition.df$cl.min, transition.df$cl.max)
  transition.df.comb= do.call("rbind",tapply(1:nrow(transition.df),transition.df$cl.pair, function(x){
    tmp = transition.df[x,][1,]
    tmp$Freq = sum(transition.df[x,"Freq"])
    tmp[,c(4,5,3)]
  }))
  cl.size = table(cell.cl.map.df$org.cl)
  transition.df.comb$cl.min.size = cl.size[transition.df.comb$cl.min]
  transition.df.comb$cl.max.size = cl.size[transition.df.comb$cl.max]
  transition.df.comb$ratio = with(transition.df.comb,Freq/pmin(cl.min.size,cl.max.size))
  transition.df.comb$cl1_label = cl.df[as.character(transition.df.comb$cl.min),"cluster_label"]
  transition.df.comb$cl2_label = cl.df[as.character(transition.df.comb$cl.max),"cluster_label"]
  colnames(transition.df.comb)[c(1:3,6)] = c("cl.x","cl.y","trans_n","trans_ratio")
  
  transition.df.comb
}
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.cl.80.rda")
cl.df$sst80_cl <- rownames(cl.df)

write.csv(cl.df,"sst80_cl_anno.csv")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.cl.300.rda")
cl.df$sst300_cl <- rownames(cl.df)

write.csv(cl.df,"sst300_cl_anno.csv")

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.cell.cl.map.80.rda")
samples <- rownames(cell.cl.map.df)
sst80_df <- cell.cl.map.df %>%
  mutate(org.cl = as.numeric(org.cl))

sst80_core_int <- sst80_df %>%
  group_by(org.cl) %>%
  summarise(n_core = sum(best.cl == org.cl))
names(sst80_core_int)[1] <- "cl"

sst80_cl_anno <- read.delim("sst80_coords.tsv", header = TRUE, sep = "\t")
rownames(sst80_cl_anno) <- sst80_cl_anno$sst80_cl
names(sst80_cl_anno)[1:4] <- c("cl","cluster_id","cluster_label","cluster_color")

sst80_cl_anno <- sst80_cl_anno %>%
  left_join(sst80_core_int)

sst80_transitions <- cell.cl.map.df_to_transition.df.comb(sst80_df, sst80_cl_anno) %>%
  mutate(cl.x = as.numeric(cl.x),
         cl.y = as.numeric(cl.y)) %>%
  filter(trans_ratio > 0.05) %>%
  select(cl.x, cl.y, trans_n)

sst80_interm <- net_coord_plot(sst80_cl_anno, sst80_transitions, sst80_cl_anno, cluster_ids = 1:30, split = F,
                               "trans_n", val_min = 0, link_color = "dodgerblue") +
  scale_y_reverse()

ggsave("sst80_constellation.pdf",
       sst80_interm,
       width = 4, height = 3,
       useDingbats = F)

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/Sst.cell.cl.map.300.rda")
sst300_df <- cell.cl.map.df %>%
  mutate(org.cl = as.numeric(org.cl))

sst300_core_int <- sst300_df %>%
  group_by(org.cl) %>%
  summarise(n_core = sum(best.cl == org.cl))
names(sst300_core_int)[1] <- "cl"

sst300_cl_anno <- read.delim("sst300_coords.tsv", header = TRUE, sep = "\t")
rownames(sst300_cl_anno) <- sst300_cl_anno$sst300_cl
names(sst300_cl_anno)[1:4] <- c("cl","cluster_id","cluster_label","cluster_color")

sst300_cl_anno <- sst300_cl_anno %>%
  left_join(sst300_core_int)

sst300_transitions <- cell.cl.map.df_to_transition.df.comb(sst300_df, sst300_cl_anno) %>%
  mutate(cl.x = as.numeric(cl.x),
         cl.y = as.numeric(cl.y)) %>%
  filter(trans_ratio > 0.05)

sst300_interm <- net_coord_plot(sst300_cl_anno, sst300_transitions, sst300_cl_anno, cluster_ids = 1:13, split = F,
                               "trans_n", val_min = 0, link_color = "dodgerblue") +
  scale_y_reverse()

ggsave("sst300_constellation.pdf",
       sst300_interm,
       width = 4, height = 3,
       useDingbats = F)


load("//allen/programs/celltypes/workgroups/rnaseqanalysis/V1_ALM/process_new/cell.cl.map.df.rda")
anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
sample_to_cl <- anno %>% select(sample_id, cl)

sst150_df <- cell.cl.map.df[samples,] %>%
  mutate(sample_id = rownames(.)) %>%
  mutate(org.cl = as.numeric(as.character(org.cl))) 

sst150_core_int <- sst150_df %>%
  group_by(org.cl) %>%
  summarise(n_core = sum(best.cl == org.cl))

sst150_cl_anno <- read.delim("sst150_coords.tsv", header = TRUE, sep = "\t")
rownames(sst150_cl_anno) <- sst150_cl_anno$sst150_cl
names(sst150_cl_anno)[1:5] <- c("cl","sst150_cl","cluster_id","cluster_label","cluster_color")

sst150_cl_anno <- sst150_cl_anno %>%
  left_join(sst150_core_int, by = c("cl"="org.cl")) 

sst150_transitions <- cell.cl.map.df_to_transition.df.comb(sst150_df, sst150_cl_anno) %>%
  mutate(cl.x = as.numeric(cl.x),
         cl.y = as.numeric(cl.y)) %>%
  filter(trans_ratio > 0.05)

sst150_interm <- net_coord_plot(sst150_cl_anno, sst150_transitions, sst150_cl_anno, cluster_ids = 30:50, split = F,
                                "trans_n", val_min = 0, link_color = "dodgerblue") +
  scale_y_reverse()

ggsave("sst150_constellation.pdf",
       sst150_interm,
       width = 4, height = 3,
       useDingbats = F)

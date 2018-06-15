library(feather)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = F)

anno <- read_feather("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")

# Fix missing cre data
# anno <- anno %>%
#   mutate(genotype_id = ifelse(donor_label == "296440", 6, genotype_id),
#          genotype_label = ifelse(donor_label == "296440", "Chat-IRES-Cre-neo/wt;Ai14(RCL-tdT)/wt", genotype_label),
#          genotype_color = ifelse(donor_label == "296440", "#FFAE00", genotype_color),
#          cre_id = ifelse(donor_label == "296440", 2,cre_id),
#          cre_label = ifelse(donor_label == "296440", "Chat-IRES-Cre-neo",cre_label),
#          reporter_id = ifelse(donor_label == "296440", 3 ,reporter_id),
#          reporter_label = ifelse(donor_label == "296440", "Ai14(RCL-tdT)" ,reporter_label))

sampling <- anno %>%
  filter(cluster_id %in% 1:133) %>%
  #filter(class_label != "Low Quality") %>%
  select(sample_id, cre_label, reporter_label, genotype_label, layer_label, layer_color, inj_type_label, region_label, facs_label)

genotype_sampling <- sampling %>%
  filter(inj_type_label == "No Injection" & facs_label != "RFP-negative")

genotype <- genotype_sampling %>%
  group_by(genotype_label, cre_label, reporter_label) %>%
  summarise(genotype_n = n()) %>%
  ungroup() %>%
  arrange(-genotype_n) %>%
  mutate(genotype_order = 1:n())

other_sampling <- sampling %>%
  filter(inj_type_label == "retrograde" | facs_label == "RFP-negative") %>%
  mutate(genotype_label = ifelse(inj_type_label != "No Injection","Retrograde Injection", genotype_label),
         cre_label = ifelse(inj_type_label != "No Injection","Retrograde", cre_label),
         reporter_label = ifelse(inj_type_label != "No Injection","Injection", reporter_label)) %>%
  mutate(genotype_label = ifelse(facs_label == "RFP-negative", paste0(genotype_label, "-neg"), genotype_label),
         cre_label = ifelse(facs_label == "RFP-negative", paste0(cre_label, "-neg"),cre_label),
         reporter_label = ifelse(facs_label == "RFP-negative","Negative",reporter_label))

other <- other_sampling %>%
  group_by(genotype_label, cre_label, reporter_label) %>%
  summarise(genotype_n = n()) %>%
  ungroup() %>%
  mutate(genotype_order = c(1,3,2) + max(genotype$genotype_order))

genotype <- rbind(genotype, other)
genotype_sampling <- rbind(genotype_sampling, other_sampling)

bar_width <- 0.8
bar_spacing <- 0.1

crl <- genotype_sampling %>%
  group_by(genotype_label, cre_label, reporter_label, region_label, layer_label, layer_color) %>%
  summarise(crl_n = n()) %>%
  left_join(genotype) %>%
  ungroup() %>%
  group_by(genotype_label, region_label) %>%
  arrange(desc(layer_label)) %>%
  mutate(layer_cum = cumsum(crl_n),
         layer_sum = sum(crl_n)) %>%
  mutate(genotype_xmax = 2*genotype_order,
         genotype_center = genotype_xmax - 1,
         alm_center = genotype_center - bar_width/2 - bar_spacing/2,
         v1_center = genotype_center + bar_width/2 + bar_spacing/2) %>%
  mutate(xmin = ifelse(region_label == "ALM", genotype_center - bar_width - bar_spacing/2, genotype_center + bar_spacing/2),
         xmax = ifelse(region_label == "ALM", genotype_center - bar_spacing/2, genotype_center + bar_width + bar_spacing/2),
         ymax = layer_cum,
         ymin = lag(layer_cum, default = 0)) %>%
  ungroup() %>%
  as.data.frame()

n_genotypes <- nrow(genotype)

region_labels <- data.frame(x = c((1:n_genotypes - 1)*2 + (1 - bar_width/2 - bar_spacing/2), 
                                  (1:n_genotypes - 1)*2 + (1 + bar_width/2 + bar_spacing/2)),
                            y = -1,
                            label = rep(c("A","V"), each = n_genotypes),
                            region_label = rep(c("ALM","VISp"), each = n_genotypes),
                            region_color = rep(c("#212021","#848EBC"), each = n_genotypes),
                            genotype_order = rep(c(1:n_genotypes), 2))

genotype_labels <- crl %>%
  select(genotype_label, cre_label, reporter_label, alm_center, v1_center, genotype_center) %>%
  unique() %>%
  mutate(label = paste(cre_label, reporter_label, sep = "\n"),
         x = genotype_center)

n_labels <- crl %>%
  select(region_label, genotype_order, alm_center, v1_center, layer_sum) %>%
  unique() %>%
  mutate(label = layer_sum,
         x = ifelse(region_label == "ALM", alm_center, v1_center),
         y = layer_sum + 20)

all_n_labels <- region_labels %>%
  left_join(n_labels, by = c("region_label","genotype_order")) %>%
  rename(y = y.y) %>%
  rename(x = x.x) %>%
  rename(label = label.y) %>%
  mutate(label = ifelse(is.na(layer_sum), 0, layer_sum),
         y = ifelse(is.na(layer_sum), 20, y))

# All genotype labelled cells
p <- ggplot() +
  geom_rect(data = crl,
            aes(xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = layer_color)) +
  geom_text(data = region_labels,
            aes(x = x,
                y = -1,
                label = label,
                color = region_color),
            vjust = 1,
            size = 2*5/6) +
  geom_text(data = genotype_labels,
            aes(x = alm_center,
                y = -125,
                label = cre_label),
            angle = 90, 
            hjust = 1, 
            vjust = 0.3,
            size = 2*5/6) +
  geom_text(data = genotype_labels,
            aes(x = v1_center,
                y = -125,
                label = reporter_label),
            angle = 90, 
            hjust = 1, 
            vjust = 0.3,
            size = 2*5/6) +
  geom_text(data = all_n_labels,
            aes(x = x,
                y = y,
                label = label),
            angle = 90, 
            hjust = 0, 
            vjust = 0.3,
            size = 2*5/6) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_y_continuous("N Cells",limits = c(-1000,3000), breaks = c(0, 500, 1000, 1500, 2000, 2500, 3000)) +
  scale_x_continuous("") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank())

p

ggsave("genotype_sampling_all.pdf", p, width = 8.5, height = 2.25)

# Pan lines
pan_p <- p +
  scale_y_continuous("",limits = c(-1000, 3000)) +
  scale_x_continuous("",limits = c(0,8))

ggsave("genotype_sampling_pan.pdf", pan_p, width = 6, height = 6)

# Broad lines
broad_p <- p +
  scale_y_continuous("",limits = c(-500, 1000)) +
  scale_x_continuous("",limits = c(8,16))

ggsave("genotype_sampling_broad.pdf", broad_p, width = 6, height = 6)


# Specific lines
spec_p <- p +
  scale_y_continuous("",limits = c(-150, 250)) +
  scale_x_continuous("",limits = c(16, 69))

ggsave("genotype_sampling_specific.pdf", spec_p, width = 12, height = 6)

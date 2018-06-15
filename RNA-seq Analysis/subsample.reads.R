load("exon.rda")
source("~/zizhen/My_R/subsample.R")
library(parallel)
frac= seq(0.1, 0.9, by=0.1)
tmp = mclapply(frac, function(f){
  print(f)
  dat = subsample_frac(exon, f)
  colSums(dat > 0)
},mc.cores=9)
tmp[[10]] = colSums(exon>0)
tmp.df = do.call("cbind",tmp)
colnames(tmp.df) = seq(0.1,1,by=0.1)
subsample.reads.df <- as.data.frame(as.table(tmp.df))
colnames(subsample.reads.df)= c("sample_id","Frac","gene.counts")

subsample.reads.df$class = cl.df[as.character(cl[as.character(subsample.reads.df$sample_id)]),"class_label"]
subsample.reads.df[subsample.reads.df$class=="Endothelial", "class"] = "Non-Neuronal"

subsample.reads.df= droplevels(subsample.reads.df)
subsample.reads.df$Frac = factor(subsample.reads.df$Frac)

g=ggplot(subsample.reads.df,
  aes(x = Frac,
      y = gene.counts)) +
  geom_violin(aes(fill = class),size = 0.1) +
  stat_summary(aes(x = Frac,
                   y = gene.counts),
               fun.ymin = function(z) { quantile(z,0.25) },
               fun.ymax = function(z) { quantile(z,0.75) },
               fun.y = median,
               color = "black",
               size = 0.1) +
  facet_wrap(~ class, nrow=3) + 
  scale_fill_manual(values=c("dodgerblue", "skyblue","gray80")) +
  scale_y_continuous(limits = c(0, 16100), breaks = c(0, 4000, 8000, 12000, 16000)) +
  theme_bw(base_size = 5) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        legend.position="none") +
  xlab("Sampling rate") +
  ylab("Gene counts")

ggsave("subsample.genes_detected.pdf", g, width=3, height=2.5, useDingbats=F)
save(subsample.reads.df, file="subsample.reads.df.rda")

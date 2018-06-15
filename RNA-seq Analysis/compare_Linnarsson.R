library(feather)
library(iterclust)

load("cl.clean.rda")
load("norm.dat.rda")
load("samp.dat.rda")
load("cl.med.rda")
load("V1.cl.rda")

###load mapping data and preprocessing
d="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/external/mouse_SS_CA1_Zeisel_2015_20170620"
linnarsson.samp.dat = as.data.frame(read_feather(file.path(d, "anno.feather")))
row.names(linnarsson.samp.dat)=linnarsson.samp.dat[["sample_id"]]
linnarsson.counts = as.data.frame(read_feather(file.path(d, "data_t.feather")))
row.names(linnarsson.counts)= linnarsson.counts[,1]
linnarsson.counts= as.matrix(linnarsson.counts[,-1])

###compute overlapping genes present in both dataset, and normalize data on overlapping genes
common.genes = intersect(row.names(linnarsson.counts), row.names(norm.dat))
cells= names(cl)
linnarsson.counts=linnarsson.counts[common.genes,]
linnarsson.dat = t(t(linnarsson.counts)*10^6/colSums(linnarsson.counts))
linnarsson.dat = log2(linnarsson.dat+1)

##Use our clusters as reference, compute markers
select.markers=intersect(V1.markers, row.names(linnarsson.dat))
select.markers = select.markers[rowMaxs(linnarsson.dat[select.markers,])>0]

map.result = map_sampling(norm.dat[,V1.cells], droplevels(cl[V1.cells]), linnarsson.dat, markers=select.markers)

map.df = map.result$map.df

map.df$pred_cluster_label = cl.df[as.character(map.df$pred.cl),"cluster_label"]
map.df$cl = linnarsson.samp.dat[match(row.names(map.df),linnarsson.samp.dat$sample_id),"cluster_label"]
map.df$coarse_cl= linnarsson.samp.dat[match(row.names(map.df),linnarsson.samp.dat$sample_id),"coarse_label"]
map.df$group = linnarsson.samp.dat[match(row.names(map.df),linnarsson.samp.dat$sample_id),"group_label"]
map.df$tissue = linnarsson.samp.dat[match(row.names(map.df),linnarsson.samp.dat$sample_id),"tissue_label"]

tmp.cor =  rowMaxs(cor(linnarsson.dat[select.markers,],cl.med[select.markers,]))
tmp.cor[is.na(tmp.cor)] = 0
map.df$cor = tmp.cor
save(map.df, file="map.linnarsson.df.rda")


ss.map.df = with(map.df,map.df[tissue == "sscortex" & cor > 0.4,])
save(ss.map.df, file="map.linnarsson.ss.df.rda")


tb = as.data.frame(with(ss.map.df%>%filter(prob > 0.9), table(coarse_cl, pred_cluster_label))) %>% filter(Freq > 1)

tb = as.data.frame(with(ss.map.df%>%filter(prob > 0.9), table(cl, pred_cluster_label))) %>% filter(Freq > 1)


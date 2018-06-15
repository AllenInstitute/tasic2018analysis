library(feather)
library(iterclust)
library(dplyr)
library(Matrix)


load("cl.final.rda")
load("cl.clean.rda")
load("norm.dat.rda")
load("samp.dat.rda")
load("V1.cl.rda")


d = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_VISp_SMV1_1679/"
V1.2016.samp.dat = as.data.frame(read_feather(file.path(d, "anno.feather")))
row.names(V1.2016.samp.dat)=V1.2016.samp.dat[["sample_id"]]

V1.2016.counts = as.matrix(round(read.csv("dat_1679/counts.csv",row.names=1,header=T)))
V1.2016.norm.dat = log2(t(t(V1.2016.counts)*10^6/Matrix::colSums(V1.2016.counts))+1)
V1.2016.cl = setNames(V1.2016.samp.dat[["final_label"]],V1.2016.samp.dat[["sample_id"]])
V1.2016.cl = V1.2016.cl[V1.2016.cl!=""]
tmp = unique(V1.2016.cl)
tmp = tmp[order(as.integer(gsub(" f.*$","",tmp)))]
V1.2016.cl = factor(V1.2016.cl,levels= tmp)
tmp.cl = V1.2016.cl
levels(tmp.cl) = 1:length(levels(tmp.cl))
tmp = select_markers(V1.2016.norm.dat, tmp.cl, n.markers=50, de.param = de_param(q1.th=0.4, q.diff.th=0.7))
V1.2016.de.genes  = tmp$de.genes
V1.2016.markers=tmp$markers
common.markers=intersect(V1.2016.markers, row.names(norm.dat))


map.result = map_sampling(V1.2016.norm.dat[,names(V1.2016.cl)], V1.2016.cl, norm.dat[, V1.cells], markers = common.markers)
 
map.df = map.result$map.df
map.df$org.cl = cl[row.names(map.df)]
map.df$org.cl_label = cl.df[as.character(cl[row.names(map.df)]), "cluster_label"]
tb = as.data.frame(with(map.df%>%filter(prob > 0.9), table(org.cl_label, pred.cl))) %>% filter(Freq > 1)
save(map.df, file="V1.ref.49.df.rda")
map.ref.2016.df = map.df


load("de.genes.rda")
tmp = select_markers(norm.dat, droplevels(cl[V1.cells]), de.genes=de.genes, n.markers=50)
common.markers=intersect(tmp$markers, row.names(V1.2016.norm.dat))
V1.markers= tmp$markers
save(V1.markers, file="V1.markers.rda")
map.result = map_sampling(as.matrix(norm.dat[common.markers,V1.cells]), droplevels(cl[V1.cells]), V1.2016.norm.dat, markers = common.markers)
map.df = map.result$map.df
map.df$org.cl = V1.2016.cl[row.names(map.df)]
map.df$pred.cl_label = cl.df[as.character(map.df$pred.cl), "cluster_label"]
map.df = map.df[!is.na(map.df$org.cl),]
tb = as.data.frame(with(map.df%>%filter(prob > 0.9), table(org.cl, pred.cl_label))) %>% filter(Freq > 1)

save(map.df, file="V1.ref.101.df.rda")








 

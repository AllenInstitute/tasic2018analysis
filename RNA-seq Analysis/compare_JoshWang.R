library(feather)
library(iterclust)

load("cl.clean.rda")
load("norm.dat.rda")
load("samp.dat.rda")
load("cl.med.rda")

###load mapping data and preprocessing
d = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/external/mouse_GABA_Paul_2016_20171002"
Wang.samp.dat = as.data.frame(read_feather(file.path(d, "anno.feather")))
row.names(Wang.samp.dat)=Wang.samp.dat[["sample_id"]]
Wang.counts = as.data.frame(read_feather(file.path(d, "data_t.feather")))
row.names(Wang.counts)= Wang.counts[,1]
Wang.counts= as.matrix(Wang.counts[,-1])

###compute overlapping genes present in both dataset, and normalize data on overlapping genes
common.genes = intersect(row.names(Wang.counts), row.names(norm.dat))

Wang.counts=Wang.counts[common.genes,]
Wang.dat = t(t(Wang.counts)*10^6/colSums(Wang.counts))
Wang.dat = log2(Wang.dat+1)

##Use our clusters as reference, compute markers
load("de.genes.rda")
inh.cl = droplevels(cl[cl %in% row.names(cl.df)[cl.df$class_label=="GABAergic"]])
inh.cells=names(inh.cl)
inh.markers = select_markers(norm.dat, inh.cl, de.genes = de.genes, n.markers=50)$markers
inh.markers=intersect(inh.markers, common.genes)
common.markers= inh.markers[rowMaxs(Wang.dat[inh.markers,])>0]

map.result = map_sampling(as.matrix(norm.dat[common.markers, inh.cells]), droplevels(cl[inh.cells]), Wang.dat, markers = common.markers)
map.df = map.result$map.df
map.freq = map.result$map.freq

tmp.cor = rowMaxs(cor(Wang.dat[common.markers,], cl.med[common.markers, droplevels(cl[inh.cells])]))
tmp.cor[is.na(tmp.cor)]=0
map.df$cor = tmp.cor

map.df$cluster_label = as.character(cl.df[as.character(map.df$pred.cl),"cluster_label"])
tmp = is.na(map.df$cluster_label)
map.df$cluster_label[tmp] = map.df$cl[tmp]
map.df$Wang.cl = Wang.samp.dat[match(row.names(map.df),Wang.samp.dat$sample_id),"cell_class_label"]
save(map.df, file="map.Wang.rda")

 


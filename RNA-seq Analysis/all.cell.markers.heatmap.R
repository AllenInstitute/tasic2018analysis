load("cl.final.rda")
library(feather)
anno <- read_feather("/data/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/anno.feather")
tmp.cl.df = as.data.frame(unique(anno[,c("cl","dendcluster_id")]))
dend.id = with(tmp.cl.df, setNames(dendcluster_id, cl))
cl.ord = names(dend.id)[order(dend.id)]
cl.ord = cl.ord[cl.ord %in% levels(cl.clean)]


load("de.genes.rda")
load("cl.med.rda")

load("cl.clean.rda")
cl.clean = setNames(factor(as.character(cl.clean), cl.ord), names(cl.clean))


####Compute heatmap at single cell level that cover all the clusters.
max.cl.size=50
all.markers= select_markers(norm.dat, cl.clean, de.genes=de.genes, n.markers=20)$markers
core.cells=names(cl)[is.na(cell.cl.map.df[names(cl.clean),"transition.cl"])]
tmp.cells = unlist(tapply(core.cells, cl.clean[core.cells], function(x){
  x= sample(x, min(length(x),max.cl.size))
},simplify=FALSE))
tmp.cl = setNames(factor(as.character(cl.clean[tmp.cells]), cl.ord), tmp.cells)
tmp.cl = sort(tmp.cl)

cl.med = cl.med[,levels(cl.clean)]
all.markers= all.markers[rowMaxs(cl.med[all.markers,]) > 2]
gene.min = apply(cl.med[all.markers,],1, function(x)min(which(x>2)))
gene.sd= apply(cl.med[all.markers,],1, function(x)sd(which(x>2)))
gene.sd[is.na(gene.sd)]=0

select.markers = all.markers[gene.sd < 20 | gene.sd > 20 & rowSums(cl.med[all.markers,] > 2) >  20]
select = gene.sd[select.markers] < 10 | rowSums(cl.med[select.markers,]>2) < 15
set1 = select.markers[select]
set2 = select.markers[!select]


tmp=setNames(rep(6, length(de.genes)),names(de.genes))
set3=unique(unlist(select_markers_pair(norm.dat[set2,], add.genes=tmp, de.genes=de.genes,rm.genes=setdiff(row.names(norm.dat),set2))))
set3=set3[hclust(dist(cl.med[set3,]),method="average")$order]

ord.genes = c(names(sort(gene.min[set1])),set3)
tmp.dat = as.matrix(norm.dat[ord.genes,names(tmp.cl)])
tmp.dat = tmp.dat/rowMaxs(tmp.dat)
cl.col = as.character(cl.df[as.character(cl[colnames(tmp.dat)]),"cluster_color"])

gray.red <-colorRampPalette(c("gray", "red"))

pdf("all.heatmap.pdf",height=15,width=12)
heatmap.3(tmp.dat,Rowv=NULL, Colv=NULL, ColSideColors=cl.col,trace="none",col=gray.red(100),labRow=FALSE, labCol=FALSE)
dev.off()

darkblue.orange <-colorRampPalette(c("darkblue", "orange"))
pdf("all.heatmap.2.pdf",height=15,width=12)
heatmap.3(tmp.dat,Rowv=NULL, Colv=NULL, ColSideColors=cl.col,trace="none",col=darkblue.orange(100),labRow=FALSE, labCol=FALSE)
dev.off()

cell.marker.dat = tmp.dat
save(cell.marker.dat, file="cell.marker.dat.rda")

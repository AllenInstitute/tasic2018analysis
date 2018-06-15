co.stats = get_cl_co_stats(cl, co.ratio)
save(co.stats, file="co.stats.rda")

select.cl= levels(cl.clean)
co.stats.df = as.data.frame(as.table(co.stats$cl.co.ratio[select.cl, select.cl]))
colnames(co.stats.df)=c("cl.x","cl.y","co.ratio")
co.stats.df = co.stats.df[co.stats.df[,3]>0.05,]
co.stats.df$cl.x.label = cl.df[as.character(co.stats.df$cl.x), "cluster_label"]
co.stats.df$cl.y.label = cl.df[as.character(co.stats.df$cl.y), "cluster_label"]
co.stats.df = co.stats.df[as.integer(co.stats.df[,1]) <= as.integer(co.stats.df[,2]),]
write.table(co.stats.df, file="co.stats.df.csv",sep=",",row.names=F,quote=F)

tmp.dat = co.stats$cl.co.ratio[select.cl, select.cl]
row.names(tmp.dat)=colnames(tmp.dat) = cl.df[colnames(tmp.dat),"cluster_label"]
pdf("cl.co.pdf",height=10,width=10)
heatmap.3(tmp.dat, Rowv=as.dendrogram(dend), Colv=as.dendrogram(dend), trace="none",col=blue.red(100),cexRow=0.5, cexCol=0.5)
dev.off()

select.cells = sample_cells(cl.clean, 50)

ord=order(match(as.character(cl.clean[select.cells]), labels(dend)), co.stats$cell.cl.confusion[select.cells])
select.cells= select.cells[ord]

tmp.dat = co.stats$cell.cl.co.ratio[select.cells,labels(dend)]
colnames(tmp.dat) = cl.df[colnames(tmp.dat),"cluster_label"]
library(gplots)
pdf("cell.cl.pdf")
heatmap.2(tmp.dat, Rowv=NULL, Colv=as.dendrogram(dend), trace="none",col=blue.red(100), labRow=FALSE, labCol=FALSE)
dev.off()

library(WGCNA)
tom = TOMsimilarity(as.matrix(co.ratio[select.cells, select.cells]))
colnames(tom)=row.names(tom)=select.cells
all.hc = hclust(as.dist(1-tom), method="average")
ord1 = all.hc$labels[all.hc$order]
ord = names(cl.clean)[order(cl.clean,match(names(cl.clean), ord1))]
sep = cl.clean[ord]
sep=which(sep[-1]!=sep[-length(sep)])
png("all.co.png",height=10,width=10)
heatmap.3(as.matrix(co.ratio[ord,ord]), col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
dev.off()


select.cells = names(cl.clean)[cl.clean %in% row.names(cl.df)[cl.df$top=="Inh"]]
select.cells = unlist(tapply(select.cells,cl.clean[select.cells],function(x)sample(x, min(length(x),100))))
ord1 = all.hc$labels[all.hc$order]
ord1 = ord1[ord1%in% select.cells]
ord = select.cells[order(cl.clean[select.cells],match(select.cells, ord1))]
Inh.co.ratio.100 = co.ratio[ord,ord]
sep = cl.clean[ord]
sep=which(sep[-1]!=sep[-length(sep)])
pdf("Inh.co.100.pdf",height=10,width=10)
#heatmap.3(Inh.co.ratio.100, col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
heatmap.3(Inh.co.ratio.100, col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL)
dev.off()
save(Inh.co.ratio.100, file="Inh.co.ratio.100.rda")

select.cells = names(cl.clean)[cl.clean %in% row.names(cl.df)[cl.df$top=="Ex"]]
select.cells = unlist(tapply(select.cells,cl.clean[select.cells],function(x)sample(x, min(length(x),100))))
ord1 = all.hc$labels[all.hc$order]
ord1 = ord1[ord1%in% select.cells]
ord = select.cells[order(cl.clean[select.cells],match(select.cells, ord1))]
sep = cl.clean[ord]
sep=which(sep[-1]!=sep[-length(sep)])
Ex.co.ratio.100 = co.ratio[ord,ord]
pdf("Ex.co.100.pdf",height=10,width=10)
heatmap.3(Ex.co.ratio.100, col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
dev.off()
save(Ex.co.ratio.100, file="Ex.co.ratio.100.rda")



select.cells = names(cl.clean)[cl.clean %in% row.names(cl.df)[cl.df$top%in%c("Glia","Endo")]]
select.cells = unlist(tapply(select.cells,cl.clean[select.cells],function(x)sample(x, min(length(x),100))))
ord1 = all.hc$labels[all.hc$order]
ord1 = ord1[ord1%in% select.cells]
ord = select.cells[order(cl.clean[select.cells],match(select.cells, ord1))]
sep = cl.clean[ord]
sep=which(sep[-1]!=sep[-length(sep)])
Noneuron.co.ratio.100 = co.ratio[ord,ord]
pdf("Non.neuronal.co.100.pdf",height=10,width=10)
heatmap.3(Noneuron.co.ratio.100, col = blue.red(100), trace="none", Rowv=NULL, Colv=NULL,colsep=sep,sepcolor="black")
dev.off()
save(Noneuron.co.ratio.100, file="Noneuron.co.ratio.100.rda")

select.cells = c(colnames(Inh.co.ratio.100),colnames(Ex.co.ratio.100),colnames(Noneuron.co.ratio.100))
tmp.cl = cl.clean[select.cells]
levels(tmp.cl) = cl.df[levels(tmp.cl),"cluster_label"]
tmp.cl = setNames(factor(as.character(tmp.cl), levels=labels(dend)),names(tmp.cl))
select.cells=names(tmp.cl)
ord1 = all.hc$labels[all.hc$order]
ord1 = ord1[ord1%in% select.cells]
ord = select.cells[order(tmp.cl,match(select.cells, ord1))]
co.ratio.100= co.ratio[ord,ord]
save(co.ratio.100, file="co.ratio.100.rda")


###cl, and cl.df are previous annotation of clusters
load("merge.result.rda")
compare.result = compare_annotate(merge.result$cl, ref.cl, ref.cl.df)
ggsave("map.pdf", compare.result$g,height=10,width=10)

cl = compare.result$cl
cl.df = compare.result$cl.df

cl.region =table(cl, samp.dat[names(cl), "Region"])
cl.region = as.data.frame.matrix(round(cl.region/rowSums(cl.region),digits=2))
cl.df[,colnames(cl.region)] = cl.region[row.names(cl.df),]

cl.df$region_label =  ""
cl.df$region_label[cl.df$VISp > 0.9] = "VISp"
cl.df$region_label[cl.df$ALM > 0.9] = "ALM"

cl.gene.counts=tapply(names(cl),cl, function(x)mean(samp.dat[x, "Genes.With.FPKM"]))
cl.df$gene.counts=round(cl.gene.counts[row.names(cl.df)])

write.csv(cl.df, "cl.df.csv")



de.genes = de_score(norm.dat, cl, de.param = de.param)
tmp = select_markers(norm.dat, cl, n.markers=50, de.genes=de.genes)
select.markers= tmp$markers
save(de.genes, file="de.genes.rda")


cl.cor= cor(cl.med[select.markers,])
de.genes.sym = de.genes
pairs= names(de.genes)
pairs = as.data.frame(do.call("rbind",strsplit(pairs,"_")))
row.names(pairs) = names(de.genes)
pairs$tmp = paste0(pairs[,2],"_",pairs[,1])
de.genes.sym[pairs$tmp] = de.genes[row.names(pairs)]
for(i in 1:nrow(pairs)){
  de.genes.sym[[pairs[i,"tmp"]]]$up.genes = de.genes[[row.names(pairs)[i]]]$down.genes
  de.genes.sym[[pairs[i,"tmp"]]]$down.genes = de.genes[[row.names(pairs)[i]]]$up.genes  
}

source("~/zizhen/My_R/hicat/R/doublet.R")
nn.df = find_doublet(cl.df, cl.sim = cl.cor, de.genes.sym=de.genes.sym)
doublet.df = nn.df[with(nn.df, which(max.olap.ratio1 > 0.80 & max.olap.ratio2 > 0.80 &  up.genes.score > 200 & max.olap.score2 > 100 & olap.cl.sim < 0.85)),]


cl.clean = row.names(cl.df)[cl.df$class_label!="Low Quality"]
nn.df = find_doublet(cl.df[cl.clean,], cl.sim = cl.cor[cl.clean,cl.clean], de.genes.sym=de.genes.sym)
  
low.df = nn.df[with(nn.df, which(up.genes < 3 & up.genes.score <50 &  nn.cl.gene.counts  - cl.gene.counts > 1000)),]

noise.cl = droplevels(cl[!cl %in% levels(cl.clean)])
###Drop clusters with confusion score > 0.6
refine.result=refine_cl(cl.clean, co.ratio=co.ratio.min, confusion.th=0.6, min.cells=100,verbose=1, tol.th=0.001)
cl.clean = droplevels(refine.result$cl)
cl = setNames(factor(c(as.character(cl.clean), as.character(noise.cl)),levels=row.names(cl.df)), c(names(cl.clean),names(noise.cl)))


cl.df = rbind(cl.df[levels(cl.clean),], cl.df[levels(noise.cl),])

cl.region =table(cl, samp.dat[names(cl), "Region"])
cl.region = as.data.frame.matrix(round(cl.region/rowSums(cl.region),digits=2))
cl.df[,colnames(cl.region)] = cl.region[row.names(cl.df),]


cl.gene.counts=tapply(names(cl),cl, function(x)mean(samp.dat[x, "Genes.With.FPKM"]))
cl.df$gene.counts=round(cl.gene.counts[row.names(cl.df)])
cl.size= table(cl)
cl.df$size = cl.size[row.names(cl.df)]
save(cl, cl.df, file="cl.final.rda")


de.genes = de_score(norm.dat, cl.clean, de.param = de.param)
select.markers = select_markers(norm.dat, cl.clean, n.markers=50, de.genes=de.genes)$markers
select.markers= tmp$markers
save(de.genes, file="de.genes.rda")
save(select.markers, file="select.markers.rda")


cl.med <- do.call("cbind",tapply(names(cl), cl, function(x){
  rowMedians(as.matrix(norm.dat[,x]))
}))
row.names(cl.med) = row.names(norm.dat)
save(cl.med, file="cl.med.rda")


l.rank = setNames(1:nrow(cl.df), row.names(cl.df))
l.color = setNames(as.character(cl.df$cluster_color),row.names(cl.df))
select.cl = levels(cl.clean)
dend.result <- build_dend(cl.med[select.markers,select.cl], l.rank=l.rank, l.color=l.color, nboot=100)
dend = dend.result$dend
cl.cor = dend.result$cl.cor
save(dend, file="dend.rda")

l.rank = setNames(1:nrow(cl.df), row.names(cl.df))

dend.collapse=sapply(seq(0.2, 0.9,by=0.1), function(conf){
  d=unbranch_by_conf(dend, conf)
  d=reorder_dend(d, l.rank)  
  d.labeled = d
  labels(d.labeled)=cl.df[labels(d.labeled),"cluster_label"]
  ggsave(paste0("dend.", conf, ".pdf"),plot(d.labeled), height=5, width=12)
  d
},simplify=F)
save(dend.collapse, file="dend.collapse.rda")
dend.p50   = dend.collapse[[1]]
save(dend.p50, file="dend.p50.rda")


clean.cl.df = cl.df[labels(dend.p50),]
noise.cl.df = cl.df[!row.names(cl.df) %in% row.names(clean.cl.df),]
cl.df = rbind(clean.cl.df, noise.cl.df)
cl.df$cluster_id = 1:nrow(cl.df)
save(cl, cl.df, file="cl.final.rda")

save(dend.p50, file="dend.p50.rda")

l.rank = setNames(1:nrow(cl.df), row.names(cl.df))
dend  = reorder_dend(dend, l.rank)
save(dend, file="dend.rda")
dend.labeled = dend
labels(dend.labeled) = cl.df[labels(dend.labeled),"cluster_label"]
ggsave("dend.pdf", plot(dend.labeled), height=7, width=14)



pdf("cl.cor.heatmap.pdf",height=12,width=12)
heatmap.3(cl.cor, Rowv=dend, Colv=dend, trace="none",col=jet.colors(100),cexCol=0.5,cexRow=0.5,breaks=c(-0.2, 0.2, seq(0.2, 1, length.out=99)))
dev.off()

set.seed(12345)

##Umap does not seem to work for this dataset. 
#umap.result <- smallvis(t(as.matrix(norm.dat[select.markers, names(cl.clean)])), method = "umap", perplexity = 25, eta = 0.01)
#umap.df = as.data.frame(umap.result)
#row.names(umap.df) = names(cl.clean)
#colnames(umap.df)=c("Lim1","Lim2")
#umap.result <- plot_tSNE_cl(norm.dat, select.markers, cl.clean, cl.df, theta=0.1, tsne.df=umap.df)
#ggsave("umap.pdf", umap.result$g, height=10, width=10)

fast_tsne <- function(dat,theta=0.05, nthreads=12, perplexity=20,...){
  fast.tsne.df <- fftRtsne(dat, theta=theta,nthreads = nthreads, perplexity=perplexity, fast_tsne_path="~/src/FIt-SNE/bin/fast_tsne",...)
  row.names(fast.tsne.df)= row.names(dat)
  colnames(fast.tsne.df)=c("Lim1","Lim2")
  fast.tsne.df = as.data.frame(fast.tsne.df)
}

source("~/src/FIt-SNE/fast_tsne.R")
tsne.result <- plot_tSNE_cl(norm.dat, select.markers, cl.clean, cl.df, theta=0.05)
ggsave("tsne.cl.pdf", tsne.result$g, height=10, width=10)
tsne.df = tsne.result$tsne.df
slow.tsne.df = tsne.df
save(slow.tsne.df, file="slow.tsne.df.rda")
                                
fast.tsne.df <- fast_tsne(t(as.matrix(norm.dat[select.markers, names(cl.clean)])))
fast.tsne.result <- plot_tSNE_cl(norm.dat, select.markers, cl.clean, cl.df, theta=0.05,tsne.df=fast.tsne.df)
ggsave("fast.tsne.cl.pdf", fast.tsne.result$g, height=10, width=10)
tsne.df = fast.tsne.df
save(fast.tsne.df, file="fast.tsne.df.rda")
save(tsne.df, file="tsne.df.rda")



load("co.ratio.min.rda")
co.stats = get_cl_co_stats(cl.clean, co.ratio = co.ratio.min)
save(co.stats, file="co.stats.rda")


select.cl=levels(cl.clean)
co.stats.df = as.data.frame(as.table(co.stats$cl.co.ratio[select.cl, select.cl]))
colnames(co.stats.df)=c("cl.x","cl.y","co.ratio")
co.stats.df = co.stats.df[co.stats.df[,3]>0.05,]
co.stats.df$cl.x.label = cl.df[as.character(co.stats.df$cl.x), "cluster_label"]
co.stats.df$cl.y.label = cl.df[as.character(co.stats.df$cl.y), "cluster_label"]
co.stats.df = co.stats.df[as.integer(co.stats.df[,1]) <= as.integer(co.stats.df[,2]),]
write.table(co.stats.df, file="co.stats.csv",sep=",",row.names=F,quote=F)

all.markers=select.markers
save(all.markers, file="all.markers.rda")


inh.cl = droplevels(cl[cl %in% row.names(cl.df)[cl.df$class_label=="GABAergic"]])
inh.markers = select_markers(norm.dat, inh.cl, n.markers=50, de.genes=de.genes)$markers
inh.tsne.df <- fast_tsne(t(as.matrix(norm.dat[inh.markers, names(inh.cl)])))
save(inh.tsne.df, file="inh.tsne.df.rda")
inh.tsne.result <- plot_tSNE_cl(norm.dat, inh.markers, inh.cl, cl.df, theta=0.1,tsne.df=inh.tsne.df)
ggsave("inh.tsne.pdf", inh.tsne.result$g, height=10, width=10)



ex.cl = droplevels(cl[cl %in% row.names(cl.df)[cl.df$class_label=="Glutamatergic"]])
ex.markers = select_markers(norm.dat, ex.cl, n.markers=50, de.genes=de.genes)$markers
ex.tsne.df <- fast_tsne(t(as.matrix(norm.dat[ex.markers, names(ex.cl)])))
save(ex.tsne.df, file="ex.tsne.df.rda")
ex.tsne.result <- plot_tSNE_cl(norm.dat, ex.markers, ex.cl, cl.df, theta=0.05,tsne.df=ex.tsne.df)
ggsave("ex.tsne.pdf", ex.tsne.result$g, height=10, width=10)


  
Lamp5.Vip.cl = droplevels(cl.clean[cl.clean %in% row.names(cl.df)[cl.df$subclass_label %in% c("Lamp5","Vip","Sncg")]])
Lamp5.Vip.markers = select_markers(norm.dat, Lamp5.Vip.cl, n.markers=50, de.genes=de.genes)$markers
Lamp5.Vip.tsne.df <- fast_tsne(t(as.matrix(norm.dat[Lamp5.Vip.markers, names(Lamp5.Vip.cl)])),theta=0.02)
save(Lamp5.Vip.tsne.df, file="Lamp5.Vip.tsne.df.rda")
Lamp5.Vip.tsne.result <- plot_tSNE_cl(norm.dat, Lamp5.Vip.markers,  Lamp5.Vip.cl, cl.df, theta=0.05,tsne.df=Lamp5.Vip.tsne.df)
ggsave("Lamp5.Vip.tsne.pdf", Lamp5.Vip.tsne.result$g, height=10, width=10)




Sst.Pvalb.cl = droplevels(cl[cl %in% row.names(cl.df)[cl.df$subclass_label %in% c("Sst","Pvalb")]])
Sst.Pvalb.markers = select_markers(norm.dat, Sst.Pvalb.cl, n.markers=50, de.genes=de.genes)$markers
Sst.Pvalb.tsne.df <- fast_tsne(t(as.matrix(norm.dat[Sst.Pvalb.markers, names(Sst.Pvalb.cl)])),theta=0.02)
save(Sst.Pvalb.tsne.df, file="Sst.Pvalb.tsne.df.rda")
Sst.Pvalb.tsne.result <- plot_tSNE_cl(norm.dat, Sst.Pvalb.markers,  Sst.Pvalb.cl, cl.df, tsne.df=Sst.Pvalb.tsne.df)
ggsave("Sst.Pvalb.tsne.pdf", Sst.Pvalb.tsne.result$g, height=10, width=10)



cl.df = read.csv("cl.df.csv",header=T)
row.names(cl.df) = cl.df$cl
l.rank = setNames(1:nrow(cl.df), row.names(cl.df))
dend = reorder_dend(dend, l.rank)
dend.labeled = dend
labels(dend.labeled) = cl.df[labels(dend), "cluster_label"]
ggsave(plot(dend.labeled), file="dend.pdf",height=7, width=13)
  

sum(sort(with(samp.dat[names(cl.clean),], table(cre_line[facs_population_plan!="RFP-negative"])),decreasing=TRUE)[1:4])

with(samp.dat[names(cl.clean),], table(facs_population_plan[cre_line=="Snap25"],Region[cre_line=="Snap25"]))




table(samp.dat[names(cl.clean),"Injection_type"])
with(cell.cl.cor.map, table(best.score + second.score  < 0.9 | org.cl!=best.score & org.cl!=second.cl & !is.na(second.cl)))
with(cell.cl.cor.map, head(cell.cl.cor.map[org.cl!=best.cl & org.cl!=second.cl & !is.na(second.cl),]))
int.df = cell.cl.cor.map[!cell.cl.cor.map$core,]
with(int.df, table(best.score + second.score > 0.9 & (org.cl==best.cl | is.na(second.cl) | org.cl==second.cl)))

cl.df.clean = droplevels(cl.df[levels(cl.clean),])
table(cl.df.clean$class_label)
table(cl.df[as.character(cl.clean),"class_label"])
    
tmp.cl = droplevels(cl[cl%in% c("118","120")])
tmp=display_cl(tmp.cl, norm.dat, prefix="Meis2-CR", col=all.col, de.genes=de.genes, n.markers=50)

table(samp.dat[names(cl.clean),"Injection_type"])
head(with(cl.df.clean, cl.df.clean[region_label!="",1:5]))

samp.dat.clean = samp.dat[names(cl.clean),]




    
###
d="/data/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520"
anno.df = read_feather(file.path(d, "anno.feather"))

clean.anno.df = as.data.frame(anno.df[anno.df$class_label!="Low Quality",])

non.retro = clean.anno.df %>% filter(inj_type_label != "retrograde")

cre.size= sort(table(anno.df$cre_label),decreasing=T)
pan.cre= names(cre.size)[c(1,3:5)]


table(non.retro$cre_label %in% pan.cre, non.retro$facs_label!="RFP-negative")

non.retro%>% filter(!cre_label %in% pan.cre, facs_label!="RFP-negative") %>% select(region_label) %>% table







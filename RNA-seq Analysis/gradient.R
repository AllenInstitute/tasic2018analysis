###
L4.cells = names(cl)[cl=="69"]
qc.index = setNames(samp.dat$percent_reads_aligned_to_introns/samp.dat$percent_reads_aligned_to_exons,row.names(samp.dat))

###normalize L4 cells first to remove QC bias. 
lm.result = lm_normalize(as.matrix(norm.dat[,L4.cells]), qc.index[L4.cells])
tmp.dat = lm.result[[1]]
fit= lm.result[[2]]
genes = row.names(fit)[fit$R_2 > 0.2]
qc.index= get_eigen(list(qc=genes), norm.dat, L4.cells)[[1]]
qc.index = setNames(qc.index[,1],row.names(qc.index))
lm.result = lm_normalize(as.matrix(norm.dat[,L4.cells]), qc.index[L4.cells])
tmp.dat = lm.result[[1]]
vg2 = findVG(tmp.dat[,L4.cells])
g=ggplot(vg2, aes(x=log10(g.means), y=dispersion)) + geom_point(size=0.5,alpha=0.5)
g=g+geom_smooth(data = vg2[vg2$g.means > 1,], method='loess')
pdf("vg2.pdf")
g
dev.off()
L4.norm.dat = tmp.dat

load("de.param.rda")
load("rm.eigen.rda")
###Run one  round of clustering on less stringent threshold.
de.param = de_param(q1.th=0.3, q.diff.th=0.7,de.score.th=100)
L4.result=onestep_clust(norm.dat=tmp.dat[select.genes,],select.cells=L4.cells, prefix="L4", de.param = de.param, vg.padj.th=0.9, rm.eigen=rm.eigen, rm.th=0.7)

de.param$q.diff.th=0.5
L4.cl= L4.result$cl
L4.result=display_cl(L4.cl, tmp.dat, prefix="L4", de.param=de.param, col=all.col[-5,L4.cells])
L4.markers= L4.result$markers

L4.eigen= get_eigen(list(L4=L4.markers),tmp.dat, L4.cells)[[1]]
L4.ord= row.names(L4.eigen)[order(-L4.eigen)]
tmp = cor(t(tmp.dat[L4.markers,L4.cells]), L4.eigen[L4.cells,])
tmp = tmp[abs(tmp[,1])>0.2,]
L4.markers=names(tmp)[order(tmp)]

###2016 datasets
select.cl = droplevels(V1.2016.cl[V1.2016.cl %in% levels(V1.2016.cl)[26:28]])
common.markers=L4.markers[L4.markers %in% row.names(V1.2016.norm.dat)]
eigen= get_eigen(list(L4=common.markers),V1.2016.norm.dat, names(select.cl))[[1]]
ord= order(select.cl, eigen[,1])
tmp.dat3 = V1.2016.norm.dat[common.markers, names(select.cl)[ord]]
tmp.dat3 = tmp.dat3/rowMaxs(tmp.dat3)
  
pdf("L4.2016.pdf",height=8,width=6)
heatmap.3(tmp.dat3, Colv=NULL,Rowv=NULL,col=blue.red(100), trace="none",cexRow=0.7)
dev.off()
L4.2016.dat = tmp.dat3
save(L4.2016.dat, file="L4.2016.dat.rda")

tmp.dat2=tmp.dat[common.markers, L4.ord]
tmp.dat2 = tmp.dat2/rowMaxs(tmp.dat2)


pdf("L4.pdf",height=8,width=6)
heatmap.3(tmp.dat2, Colv=NULL,Rowv=NULL,col=blue.red(100), trace="none",cexRow=0.7,ColSideColors=all.col["Injection_type",colnames(tmp.dat2)])
dev.off()
L4.dat = tmp.dat2
save(L4.dat, file="L4.dat.rda")






library(randomForest)
select.n=50
tmp.cl = setNames(rep(2:1, c(select.n,select.n)), c(head(L4.ord,select.n),tail(L4.ord, select.n)))
L4.rf = randomForest(t(tmp.dat[L4.markers,names(tmp.cl)]), as.factor(tmp.cl))
L4.prob = predict(L4.rf, t(tmp.dat[L4.markers, L4.ord]),type="prob")
L4.df = data.frame(prob=L4.prob[,1], eigen = L4.eigen[row.names(L4.prob),1])
L4.df$Injection_type = samp.dat[row.names(L4.df),"Injection_type"]
save(L4.df, file="L4.df.rda")

L4.test = sample(setdiff(L4.ord, names(tmp.cl)),400)
L4.sampled.df = L4.df[L4.test,]
ks.pval = format(ks.test(L4.sampled.df$prob,punif)$p.value,digits=3,scientific=TRUE)
ks.pval
save(L4.sampled.df, file="L4.sampled.df.rda")


pdf("L4.prob.pdf")
g=ggplot(L4.df, aes(x=prob, y= eigen)) + geom_point()
g = g + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g
g = ggplot(L4.df, aes(prob))+ geom_histogram(aes(y=..count../sum(..count..)),breaks=seq(0,1,by=0.05),binwidth=0.5,color="black")
g = g+geom_hline(yintercept=0.05)
g = g + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g = g +ylim(0,0.5)
g = g+ geom_text(aes(x=0.75, y=0.4, label=paste("KS test pvalue:",ks.pval)))
g
dev.off()


load("de.param.rda")
L4.cl = setNames(cut(L4.prob[,1],seq(0, 1, by=0.2),include.lowest=TRUE), row.names(L4.prob))
levels(L4.cl) = 1:length(levels(L4.cl))
L4.de.genes=display_cl(L4.cl, norm.dat, prefix="L4", de.param=de.param, col=all.col[-5,L4.cells])$de.genes
save(L4.de.genes, file="L4.de.genes.rda")
 
L4.de.num = sapply(L4.de.genes, function(x)length(x$genes))
L4.de.lfc = sapply(L4.de.genes, function(x){
  if(is.null(x$genes)){
    return(0)
  }
  top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])],100)
  mean(abs(x$de.df[top.genes, "lfc"]))
})
L4.de.q.diff = sapply(L4.de.genes, function(x){
  if(is.null(x$genes)){
    return(0)
  }
  top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])],100)
  mean(abs(x$de.df[top.genes, "q.diff"]))
})

L4.de.summary = data.frame(L4.de.num, L4.de.lfc, L4.de.q.diff)

cl1.genes= sapply(L4.de.genes, function(x)paste(head(x$up.genes,10), collapse=";"))
cl2.genes= sapply(L4.de.genes, function(x)paste(head(x$down.genes,10), collapse=";"))

L4.de.summary$cl1.genes = cl1.genes[row.names(L4.de.summary)]
L4.de.summary$cl2.genes = cl2.genes[row.names(L4.de.summary)]
save(L4.de.summary, file="L4.de.summary.rda")



###comparison between L4 VISp Rspo1 and L5 IT VISp Hsd11b1 Endou cluster, between L4 VISp Rspo1 and L5 IT VISp Batf3 
select.pairs = c("69_70", "69_72")
df.list=list()
for(p in select.pairs){
  select.cl = droplevels(cl.clean[cl.clean %in% c(as.character(pairs[p,1]),as.character(pairs[p,2]))])
  tmp.markers=select_markers(norm.dat, select.cl, de.genes=de.genes, n.markers=50)$markers
  eigen= get_eigen(list(markers=tmp.markers),norm.dat, names(select.cl))[[1]]
  ord= row.names(eigen)[order(eigen)]
  tmp.cl = setNames(rep(1:2, c(select.n,select.n)), c(head(ord,select.n),tail(ord, select.n)))
  
  rf = randomForest(t(as.matrix(norm.dat[tmp.markers,names(tmp.cl)])), as.factor(tmp.cl))
  prob = predict(rf, t(as.matrix(norm.dat[tmp.markers, ord])),type="prob")

  df = data.frame(prob=prob[,1], eigen = eigen[row.names(prob),1])
  df$cl = cl.clean[row.names(df)]
  test = setdiff(names(select.cl), names(tmp.cl))
  test = unlist(tapply(test, select.cl[test], function(x)sample(x, min(length(x), 200))))
  test.df = df[test,]

  prefix = paste(cl.df[levels(select.cl),"cluster_label"],collapse=":")
  save(df, file= paste0(prefix,".rda")) 
  df.list[[p]] = test.df
  
  ks.pval = format(ks.test(test.df$prob,punif)$p.value,digits=3,scientific=TRUE)
  g=ggplot(test.df, aes(x=prob, y= eigen)) + geom_point()
  g = g + theme_bw()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  pdf(paste0(prefix, ".prob.pdf"))
  plot(g)
  g = ggplot(test.df, aes(prob))+ geom_histogram(aes(y=..count../sum(..count..)),breaks=seq(0,1,by=0.05),binwidth=0.5,color="black")
  g = g+geom_hline(yintercept=0.05)
  g = g + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  g = g +ylim(0,0.5)
  g = g+ geom_text(aes(x=0.75, y=0.4, label=paste("KS test pvalue:",ks.pval)))
  plot(g)
  dev.off()
}
df.list[["69"]] = L4.sampled.df
save(df.list, file="gradient.df.list.rda")








load("cl.clean.rda")
load("de.genes.rda")
load("dend.rda")
pairs =  do.call("rbind",strsplit(names(de.genes),"_"))
row.names(pairs)= names(de.genes)
pairs = as.data.frame(pairs)
colnames(pairs)=c("P1","P2")
save(pairs, file="pairs.rda")

tmp=select_markers(norm.dat, cl.clean, de.genes=de.genes)
clean.de.genes=tmp$de.genes
up.de.num = sapply(clean.de.genes, function(x)length(x$up.genes))
up.score = sapply(clean.de.genes, function(x)x$up.score)
down.de.num = sapply(clean.de.genes, function(x)length(x$down.genes))
down.score = sapply(clean.de.genes, function(x)x$down.score)
head(sort(pmin(up.de.num, down.de.num)))
names(down.de.num) = with(pairs[names(down.de.num),],paste0(P2,"_",P1))

label = as.hclust(dend)$label 
de.num = c(up.de.num, down.de.num)
de.num.matrix <- convert_pair_matrix(de.num, directed=TRUE)
de.num.matrix <- de.num.matrix[label, label]
save(de.num.matrix, file="de.num.matrix.rda")
breaks=c(-1,seq(0.2,4,length.out=100))
colnames(de.num.matrix) = row.names(de.num.matrix) = cl.df[row.names(de.num.matrix),"cluster_label"]
tmp.dat=log10(de.num.matrix+1)
pdf("log10.de.num.pdf")
heatmap.3(tmp.dat, col=jet.colors(100), breaks=breaks,trace="none",Colv=dend, Rowv=dend,dendrogram="row",cexRow=0.3,cexCol=0.3)
dev.off()
full.de.num.matrix=de.num.matrix
full.up.num=up.de.num
full.down.num = down.de.num
full.de.num = de.num



de.num =  sapply(clean.de.genes, function(x)length(x$genes))
de.num.matrix <- convert_pair_matrix(de.num)
de.num.matrix <- de.num.matrix[label, label]
save(de.num.matrix, file="de.num.matrix.undirected.rda")
breaks=c(-1,seq(0.2,4,length.out=100))
colnames(de.num.matrix) = row.names(de.num.matrix) = cl.df[row.names(de.num.matrix),"cluster_label"]
tmp.dat=log10(de.num.matrix+1)
pdf("log10.de.num.undirected.pdf")
heatmap.3(tmp.dat, col=jet.colors(100), breaks=breaks,trace="none",Colv=dend, Rowv=dend,dendrogram="row",cexRow=0.3,cexCol=0.3)
dev.off()



load("cl.cor.rda")
de.num=sapply(clean.de.genes, function(x)length(x$genes))
cl.cor.df = data.frame(as.table(cl.cor))
row.names(cl.cor.df)=with(cl.cor.df,paste(Var1, Var2,sep="_"))
cl.cor.df = cl.cor.df[names(de.num),]
load("co.stats.rda")
cl.co.df = as.data.frame(as.table(co.stats$cl.co.ratio))
row.names(cl.co.df) = paste0(cl.co.df[,1], "_",cl.co.df[,2])

 
de.lfc = sapply(clean.de.genes, function(x){
  top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])],100)
  mean(abs(x$de.df[top.genes, "lfc"]))
})
de.q.diff = sapply(clean.de.genes, function(x){
  top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])],100)
  mean(abs(x$de.df[top.genes, "q.diff"]))
})

de.summary = data.frame(de.num, de.lfc, de.q.diff)
tmp=  pairs[row.names(de.summary),]
colnames(tmp)=c("cl1","cl2")
de.summary = cbind(de.summary, tmp)
de.summary$cl1_label = cl.df[as.character(de.summary$cl1),"cluster_label"]
de.summary$cl2_label = cl.df[as.character(de.summary$cl2),"cluster_label"]
de.summary$cl.cor = cl.cor.df[row.names(de.summary),"Freq"]
de.summary$cl.co = cl.co.df[row.names(de.summary),"Freq"]

save(de.summary, file="de.summary.rda")
cl1.genes= sapply(clean.de.genes, function(x)paste(head(x$up.genes,10), collapse=";"))
cl2.genes= sapply(clean.de.genes, function(x)paste(head(x$down.genes,10), collapse=";"))

de.summary$cl1.genes = cl1.genes[row.names(de.summary)]
de.summary$cl2.genes = cl2.genes[row.names(de.summary)]
save(de.summary, file="de.summary.rda")



load("gene.list.rda")
gene.list = sapply(gene.list, function(x)intersect(x, row.names(norm.dat)))
total.genes =  union(unlist(gene.list), unlist(lapply(clean.de.genes, function(x)x$genes)))
total.num = length(total.genes)
de.num.list=list()
de.pval.list=list()
pdf("log10.de.num.gene.list.pdf")
for(s in names(gene.list)){
  print(s)
  g = gene.list[[s]]
  up.de.num = sapply(clean.de.genes, function(x)length(intersect(x$up.genes,g)))
  down.de.num = sapply(clean.de.genes, function(x)length(intersect(x$down.genes,g)))
  names(down.de.num) = with(pairs[names(down.de.num),],paste0(P2,"_",P1))
  de.num = c(up.de.num, down.de.num)
  de.num.matrix <- convert_pair_matrix(de.num, directed=TRUE)
  de.num.matrix <- de.num.matrix[label, label]
  de.num.list[[s]]= de.num.matrix
  
  m = length(g)
  de.pval = sapply(names(de.num), function(p){
    x = de.num[[p]]
    n <- total.num - m
    k <- full.de.num[[p]]
    pval = phyper(x - 1, m, n, k, lower.tail = FALSE)
  })
  de.pval.matrix <- convert_pair_matrix(de.pval, directed=TRUE)
  de.pval.matrix <- de.pval.matrix[label, label]
  diag(de.pval.matrix)=1
  
   colnames(de.num.matrix) = row.names(de.num.matrix) = cl.df[row.names(de.num.matrix),"cluster_label"]
  colnames(de.pval.matrix) = row.names(de.pval.matrix) = cl.df[row.names(de.pval.matrix),"cluster_label"]
  m = max(de.num)
  main1=paste(s,"#DEG")
  main2=paste(s,"DEG Pval")
  if(m< 50){
    tmp.dat=de.num.matrix
    heatmap.3(tmp.dat, col=c("black",jet.colors(m-1)), breaks=seq(-1,m, length.out=m+1),trace="none",Colv=dend, Rowv=dend,cexRow=0.3,cexCol=0.3, main=main1)
  }
  else{
    breaks=c(-1,seq(0,log10(max(de.num)+1),length.out=100))
    tmp.dat=log10(de.num.matrix+1)
    heatmap.3(tmp.dat, col=c("black",jet.colors(99)), breaks=breaks,trace="none",Colv=dend, Rowv=dend,cexRow=0.3,cexCol=0.3, main=main1)
  }
  tmp.dat= -log10(de.pval.matrix)
  breaks=seq(0,max(tmp.dat)+1,length.out=101)
  heatmap.3(tmp.dat, col=c("black",jet.colors(99)), breaks=breaks,trace="none",Colv=dend, Rowv=dend,cexRow=0.3,cexCol=0.3, main=main2)
}
dev.off()
save(de.num.list, file="de.num.list.rda")


select.pair = c("40_41","41_48","18_41","41_53","41_65","1_41","146_41","126_41","127_41","133_41")
select.df = de.summary[select.pair,]
g=ggplot(de.summary, aes(de.num,de.lfc,color=de.q.diff)) + geom_point() + scale_color_gradient2(midpoint=0.85) + scale_x_log10()
g = g + geom_text(data=select.df, aes(de.num-0.02, de.lfc, label=paste(cl.df[as.character(cl1),"cluster_label"],cl.df[as.character(cl2),"cluster_label"],sep=":")),size=2,color="black")
g = g + geom_point(data=select.df, aes(de.num, de.lfc),color="red",pch=1)
g = g + xlab("Number of DE genes")
g = g + ylab("Mean log2(FC) of top 100 DE.genes")
ggsave("de.lfc.num.pdf",g)


 

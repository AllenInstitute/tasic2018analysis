####
load("samp.dat.rda")

library(dendextend)

pairs = do.call("rbind",strsplit(names(de.genes),"_"))
row.names(pairs) = names(de.genes)
Ex.types= row.names(cl.df)[cl.df$class_label=="Glutamatergic"]
Ex.pairs = pairs[pairs[,1] %in% Ex.types & pairs[,2] %in% Ex.types,]

Ex.V1.types = Ex.types[cl.df[Ex.types,"VISp"] >= 0.9]
Ex.ALM.types =  Ex.types[cl.df[Ex.types,"ALM"] >= 0.9]
Ex.mix.types=   Ex.types[cl.df[Ex.types,"ALM"] > 0.1 & cl.df[Ex.types, "VISp"]>0.1]

ALM.cl = droplevels(cl[cl %in% c(Ex.ALM.types, Ex.mix.types) & samp.dat[names(cl),"Region"]=="ALM"])
V1.cl = droplevels(cl[cl %in% c(Ex.V1.types,Ex.mix.types) & samp.dat[names(cl),"Region"]=="VISp"])

V1.markers = select_markers(norm.dat, V1.cl, de.genes=de.genes, n.markers=50)$markers
ALM.markers = select_markers(norm.dat, ALM.cl, de.genes=de.genes, n.markers=50)$markers


V1.ALM.comm.markers= intersect(V1.markers, ALM.markers)


ALM.ref.V1.result=map_cl_summary(norm.dat[V1.ALM.comm.markers,names(V1.cl)], V1.cl, norm.dat[V1.ALM.comm.markers,names(ALM.cl)], ALM.cl)
ALM.V1.map.df = ALM.ref.V1.result$cl.map.df


V1.ref.ALM.result=map_cl_summary(norm.dat[V1.ALM.comm.markers,names(ALM.cl)], ALM.cl, norm.dat[V1.ALM.comm.markers,names(V1.cl)], V1.cl)
V1.ALM.map.df = V1.ref.ALM.result$cl.map.df
 
V1.ALM.map.df$type = "V1.ref.ALM"
V1.ALM.map.df$V1.cl = V1.ALM.map.df$org.cl
V1.ALM.map.df$ALM.cl = V1.ALM.map.df$map.cl
ALM.V1.map.df$type = "ALM.ref.V1"
ALM.V1.map.df$V1.cl = ALM.V1.map.df$map.cl
ALM.V1.map.df$ALM.cl = ALM.V1.map.df$org.cl
 
comb.map = rbind(V1.ALM.map.df, ALM.V1.map.df)
comb.map$v1_label = cl.df[as.character(comb.map$V1.cl),"cluster_label"]
comb.map$alm_label = cl.df[as.character(comb.map$ALM.cl),"cluster_label"]

###Pruning
cell.map.df = rbind(ALM.ref.V1.result$map.df, V1.ref.ALM.result$map.df)
score.th = mean(cell.map.df$pred.score) - sd(cell.map.df$pred.score)* 1.64
comb.map.all=comb.map
comb.map=comb.map[comb.map$Prob>0.2 & comb.map$pred.score > score.th,]

load("dend.collapse.rda")
dend = dend.collapse[[3]]

V1.dend  = prune_dend(dend, setdiff(labels(tmp.dend),levels(V1.cl)))
ALM.dend  = prune_dend(dend, setdiff(labels(tmp.dend),levels(ALM.cl)))

V1.ALM.best.map=with(V1.ALM.map.df, tapply(1:nrow(V1.ALM.map.df), org.cl, function(x)as.character(map.cl)[x[which.max(Freq[x])]]))

ALM.V1.best.map=with(ALM.V1.map.df, tapply(1:nrow(ALM.V1.map.df), org.cl, function(x)as.character(map.cl)[x[which.max(Freq[x])]]))


##Reorder V1.dend according to the match to ALM
#match.score = with(comb.map.all,  tapply(Prob, list(V1.cl, ALM.cl), sum))
#match.score[is.na(match.score)]=0
#V1.ALM.best.map= setNames(colnames(match.score)[apply(match.score, 1, which.max)], row.names(match.score))
#ALM.V1.best.map= setNames(rownames(match.score)[apply(match.score, 2, which.max)], colnames(match.score))


l.rank = setNames(match(ALM.V1.best.map, labels(V1.dend)), names(ALM.V1.best.map))
ALM.dend = reorder_dend(ALM.dend, l.rank)

l.rank = setNames(match(V1.ALM.best.map, labels(ALM.dend)), names(V1.ALM.best.map))
V1.dend = reorder_dend(V1.dend, l.rank)



ALM.dend.labeled = ALM.dend
labels(ALM.dend.labeled) = as.character(cl.df[labels(ALM.dend),"cluster_label"])
ggsave("ALM.dend.pdf", plot(ALM.dend.labeled))

V1.dend.labeled = V1.dend
labels(V1.dend.labeled) = as.character(cl.df[labels(V1.dend),"cluster_label"])
ggsave("V1.dend.pdf", plot(V1.dend.labeled))

save(comb.map, V1.dend, ALM.dend, file = "V1.ALM.compare.cor.rda")
write.csv(comb.map, "comb.map.csv")  

###DEG between best map between V1 and ALM Ex types.
Ex.cells= c(names(V1.cl), names(ALM.cl))
Ex.region.cl =  setNames(paste0(cl[Ex.cells], samp.dat[Ex.cells, "Region"]), Ex.cells)
pairs = rbind(
  data.frame(VISp=paste0(names(V1.ALM.best.map), "VISp"), ALM=paste0(V1.ALM.best.map,"ALM")),
  data.frame(VISp=paste0(ALM.V1.best.map,"VISp"), ALM= paste0(names(ALM.V1.best.map),"ALM")))
pairs = pairs[!duplicated(pairs),]
pairs[,1] = as.character(pairs[,1])
pairs[,2] = as.character(pairs[,2])
Ex.region.pairs=pairs
Ex.region.de.result = de_score_pairs(norm.dat[,Ex.cells],cl=Ex.region.cl[Ex.cells], pairs=Ex.region.pairs, de.param = de.param)
Ex.region.de.genes = Ex.region.de.result$de.genes


median(sort(unlist(sapply(Ex.region.de.genes,function(x)length(x$genes)))))
sort(unlist(sapply(Ex.region.de.genes,function(x)x$score)))



tmp.cl= droplevels(cl[cl%in% c(Ex.V1.types, Ex.ALM.types)])
select.genes = select_markers(norm.dat,  tmp.cl, n.markers=50,de.genes=de.genes)$markers
V1.dat = as.matrix(norm.dat[select.genes, names(V1.cl)])
ALM.dat = as.matrix(norm.dat[select.genes, names(ALM.cl)])

load("low.th.rda")
V1.prop = rowSums(V1.dat > low.th[row.names(V1.dat)])/ncol(V1.dat)
ALM.prop = rowSums(ALM.dat > low.th[row.names(ALM.dat)])/ncol(ALM.dat)
prop.diff = (ALM.prop - V1.prop)/pmax(V1.prop,ALM.prop)

gene.df = data.frame(V1.prop, ALM.prop, prop.diff)
gene.df = gene.df[!is.na(gene.df$prop.diff),]
gene.df$prop.max = pmax(gene.df$V1.prop, gene.df$ALM.prop)
gene.df = gene.df[!is.na(gene.df$prop.diff),]
gene.df$gene = row.names(gene.df)

select.df1 = with(gene.df, gene.df[prop.diff > 0.75 & nchar(gene)< 10,])
select.df2 = with(gene.df, gene.df[prop.diff < -0.7 & nchar(gene)< 10,])

g= ggplot(gene.df, aes(sqrt(prop.max), prop.diff)) 
g = g+ geom_point() + scale_color_identity()
g = g + geom_text(data=select.df1, aes(x= prop.max, y=prop.diff, label=gene,size=0.5,hjust=0), position=position_jitter(height=0.02,width=0.01))
g = g + geom_text(data=select.df2, aes(x= prop.max, y=prop.diff, label=gene,size=0.5,hjust=0), position=position_jitter(height=0.02,width=0.01))
g = g + scale_x_log10()
ggsave("V1.ALM.gene.pdf",g,height=14,width=14)

save(gene.df, file="V1.ALM.diff.gene.rda")




#####Find regional genes among inhibitory clusters
load("low.th.rda")
Inh.types = row.names(cl.df)[cl.df$class_label=="GABAergic"]
Inh.V1.types= Inh.types[cl.df[Inh.types,"VISp"] > 0.9]
Inh.mix.types= setdiff(Inh.types, Inh.V1.types)
ALM.cl = droplevels(cl[cl %in% Inh.mix.types  & samp.dat[names(cl),"Region"]=="ALM"])
V1.cl = droplevels(cl[cl %in% c(Inh.V1.types,Inh.mix.types) & samp.dat[names(cl),"Region"]=="VISp"])
V1.markers = select_markers(norm.dat, V1.cl, de.genes=de.genes, n.markers=50)$markers
ALM.markers = select_markers(norm.dat, ALM.cl, de.genes=de.genes, n.markers=50)$markers
V1.ALM.comm.markers= intersect(V1.markers, ALM.markers)
ALM.ref.V1.result=map_cl_summary(norm.dat[V1.ALM.comm.markers,names(V1.cl)], V1.cl, norm.dat[V1.ALM.comm.markers,names(ALM.cl)], ALM.cl)
ALM.V1.map.df = ALM.ref.V1.result$cl.map.df
V1.ref.ALM.result=map_cl_summary(norm.dat[V1.ALM.comm.markers,names(ALM.cl)], ALM.cl, norm.dat[V1.ALM.comm.markers,names(V1.cl)], V1.cl)
V1.ALM.map.df = V1.ref.ALM.result$cl.map.df
 
V1.ALM.map.df$type = "V1.ref.ALM"
V1.ALM.map.df$V1.cl = V1.ALM.map.df$org.cl
V1.ALM.map.df$ALM.cl = V1.ALM.map.df$map.cl
ALM.V1.map.df$type = "ALM.ref.V1"
ALM.V1.map.df$V1.cl = ALM.V1.map.df$map.cl
ALM.V1.map.df$ALM.cl = ALM.V1.map.df$org.cl
 
comb.map = rbind(V1.ALM.map.df, ALM.V1.map.df)
comb.map$v1_label = cl.df[as.character(comb.map$V1.cl),"cluster_label"]
comb.map$alm_label = cl.df[as.character(comb.map$ALM.cl),"cluster_label"]


###Pruning
cell.map.df = rbind(ALM.ref.V1.result$map.df, V1.ref.ALM.result$map.df)
score.th = mean(cell.map.df$pred.score) - sd(cell.map.df$pred.score)* 1.64
comb.map.all=comb.map
comb.map=comb.map[comb.map$Prob>0.2 & comb.map$pred.score > score.th,]
V1.ALM.best.map=with(V1.ALM.map.df, tapply(1:nrow(V1.ALM.map.df), org.cl, function(x)as.character(map.cl)[x[which.max(Freq[x]+0.01*(as.character(org.cl[x])==as.character(map.cl[x])))]]))
ALM.V1.best.map=with(ALM.V1.map.df, tapply(1:nrow(ALM.V1.map.df), org.cl, function(x)as.character(map.cl)[x[which.max(Freq[x]+0.01*(as.character(org.cl[x])==as.character(map.cl[x])))]]))
save(comb.map, file="V1.ALM.compare.cor.Inh.rda")


Inh.cells= c(names(V1.cl), names(ALM.cl))
Inh.region.cl =  setNames(paste0(cl[Inh.cells], samp.dat[Inh.cells, "Region"]), Inh.cells)
pairs = rbind(
  data.frame(VISp=paste0(names(V1.ALM.best.map), "VISp"), ALM=paste0(V1.ALM.best.map,"ALM")),
  data.frame(VISp=paste0(ALM.V1.best.map,"VISp"), ALM= paste0(names(ALM.V1.best.map),"ALM")))
pairs = pairs[!duplicated(pairs),]
pairs[,1] = as.character(pairs[,1])
pairs[,2] = as.character(pairs[,2])
Inh.region.pairs=pairs
Inh.region.de.result = de_score_pairs(norm.dat[,Inh.cells],cl=Inh.region.cl[Inh.cells], pairs=Inh.region.pairs, de.param = de.param)
Inh.region.de.genes = Inh.region.de.result$de.genes

sort(unlist(sapply(Inh.region.de.genes,function(x)length(x$genes))))
sort(unlist(sapply(Inh.region.de.genes,function(x)x$score)))


select.genes = select_markers(norm.dat,  Inh.region.cl, n.markers=50,de.param = de.param)$markers
V1.dat = as.matrix(norm.dat[select.genes, names(V1.cl)])
ALM.dat = as.matrix(norm.dat[select.genes, names(ALM.cl)])

V1.prop = rowSums(V1.dat > low.th[row.names(V1.dat)])/ncol(V1.dat)
ALM.prop = rowSums(ALM.dat > low.th[row.names(ALM.dat)])/ncol(ALM.dat)
prop.diff = (ALM.prop - V1.prop)/pmax(V1.prop,ALM.prop)

gene.df = data.frame(V1.prop, ALM.prop, prop.diff)
gene.df = gene.df[!is.na(gene.df$prop.diff),]
gene.df$prop.max = pmax(gene.df$V1.prop, gene.df$ALM.prop)
gene.df = gene.df[!is.na(gene.df$prop.diff),]
gene.df$gene = row.names(gene.df)

select.df1 = with(gene.df, gene.df[prop.diff > 0.75 & nchar(gene)< 10,])
select.df2 = with(gene.df, gene.df[prop.diff < -0.75 & nchar(gene)< 10,])

g= ggplot(gene.df, aes(sqrt(prop.max), prop.diff)) 
g = g+ geom_point() + scale_color_identity()
g = g + geom_text(data=select.df1, aes(x= sqrt(prop.max), y=prop.diff, label=gene,size=0.5,hjust=0), position=position_jitter(height=0.02,width=0.01))
pdf("V1.ALM.Inh.gene.pdf",height=14,width=14)
g
dev.off()
save(gene.df, file="V1.ALM.diff.Inh.gene.rda")



###Differential genes for non-neuronal cluster
NN.types = row.names(cl.df)[cl.df$class_label=="Non-Neuronal"]
NN.V1.types= NN.types[cl.df[NN.types,"VISp"] > 0.9]
NN.mix.types= setdiff(NN.types, NN.V1.types)
ALM.cl = droplevels(cl[cl %in% NN.mix.types  & samp.dat[names(cl),"Region"]=="ALM"])
V1.cl = droplevels(cl[cl %in% c(NN.V1.types,NN.mix.types) & samp.dat[names(cl),"Region"]=="VISp"])
V1.markers = select_markers(norm.dat, V1.cl, de.genes=de.genes, n.markers=50)$markers
ALM.markers = select_markers(norm.dat, ALM.cl, de.genes=de.genes, n.markers=50)$markers
V1.ALM.comm.markers= intersect(V1.markers, ALM.markers)
ALM.ref.V1.result=map_cl_summary(norm.dat[V1.ALM.comm.markers,names(V1.cl)], V1.cl, norm.dat[V1.ALM.comm.markers,names(ALM.cl)], ALM.cl)
ALM.V1.map.df = ALM.ref.V1.result$cl.map.df
V1.ref.ALM.result=map_cl_summary(norm.dat[V1.ALM.comm.markers,names(ALM.cl)], ALM.cl, norm.dat[V1.ALM.comm.markers,names(V1.cl)], V1.cl)
V1.ALM.map.df = V1.ref.ALM.result$cl.map.df
 
V1.ALM.map.df$type = "V1.ref.ALM"
V1.ALM.map.df$V1.cl = V1.ALM.map.df$org.cl
V1.ALM.map.df$ALM.cl = V1.ALM.map.df$map.cl
ALM.V1.map.df$type = "ALM.ref.V1"
ALM.V1.map.df$V1.cl = ALM.V1.map.df$map.cl
ALM.V1.map.df$ALM.cl = ALM.V1.map.df$org.cl
 
comb.map = rbind(V1.ALM.map.df, ALM.V1.map.df)
comb.map$v1_label = cl.df[as.character(comb.map$V1.cl),"cluster_label"]
comb.map$alm_label = cl.df[as.character(comb.map$ALM.cl),"cluster_label"]

###Pruning
cell.map.df = rbind(ALM.ref.V1.result$map.df, V1.ref.ALM.result$map.df)
score.th = mean(cell.map.df$pred.score) - sd(cell.map.df$pred.score)* 1.64
comb.map.all=comb.map
comb.map=comb.map[comb.map$Prob>0.2 & comb.map$pred.score > score.th,]
V1.ALM.best.map=with(V1.ALM.map.df, tapply(1:nrow(V1.ALM.map.df), org.cl, function(x)as.character(map.cl)[x[which.max(Freq[x]+0.01*(as.character(org.cl[x])==as.character(map.cl[x])))]]))
ALM.V1.best.map=with(ALM.V1.map.df, tapply(1:nrow(ALM.V1.map.df), org.cl, function(x)as.character(map.cl)[x[which.max(Freq[x]+0.01*(as.character(org.cl[x])==as.character(map.cl[x])))]]))

NN.cells= c(names(V1.cl), names(ALM.cl))
NN.region.cl =  setNames(paste0(cl[NN.cells], samp.dat[NN.cells, "Region"]), NN.cells)
pairs = rbind(
  data.frame(VISp=paste0(names(V1.ALM.best.map), "VISp"), ALM=paste0(V1.ALM.best.map,"ALM")),
  data.frame(VISp=paste0(ALM.V1.best.map,"VISp"), ALM= paste0(names(ALM.V1.best.map),"ALM")))
pairs = pairs[!duplicated(pairs),]
pairs[,1] = as.character(pairs[,1])
pairs[,2] = as.character(pairs[,2])
NN.region.pairs=pairs
NN.region.de.result = de_score_pairs(norm.dat[,NN.cells],cl=NN.region.cl[NN.cells], pairs=NN.region.pairs, de.param = de.param)
NN.region.de.genes = NN.region.de.result$de.genes

sort(unlist(sapply(NN.region.de.genes,function(x)length(x$genes))))
sort(unlist(sapply(NN.region.de.genes,function(x)x$score)))





region.de.genes = c(Inh.region.de.genes, Ex.region.de.genes,nn.region.de.genes)
region.de.num=sapply(region.de.genes,function(x){
  if(length(x)==0) return(0)
  length(x$genes)
})
region.de.V1.num=sapply(region.de.genes,function(x){
  if(length(x)==0) return(0)
  length(x$up.genes)
})
region.de.ALM.num=sapply(region.de.genes,function(x){
  if(length(x)==0) return(0)
  length(x$down.genes)
})
region.de.score=sapply(region.de.genes,function(x){
  if(length(x)==0) return(0)
  x$score
})
region.de.V1.genes=  sapply(region.de.genes,function(x){
  if(is.null(x)) return("")
  return(paste(head(x$up.genes,8),collapse=" "))
})
region.de.ALM.genes=  sapply(region.de.genes,function(x){
  if(is.null(x)) return("")
  return(paste(head(x$down.genes,8),collapse=" "))
})
region.de.df = data.frame(region.de.num, region.de.V1.num, region.de.ALM.num,region.de.score, region.de.V1.genes, region.de.ALM.genes)
pairs = do.call("rbind",strsplit(names(region.de.num),"_"))
pairs[,1]=gsub("VISp","", pairs[,1])
pairs[,2]=gsub("ALM","", pairs[,2])
region.de.df = data.frame(VISp.cl=pairs[,1],ALM.cl=pairs[,2], region.de.df)
region.de.df$VISp_cl.label = cl.df[as.character(region.de.df$VISp.cl),"cluster_label"]
region.de.df$ALM_cl.label = cl.df[as.character(region.de.df$ALM.cl),"cluster_label"]


 
de.num=sapply(region.de.genes, function(x)length(x$genes))
de.lfc = sapply(region.de.genes, function(x){
  if(length(x$genes)==0){
    return(0)
  }
  top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])],100)
  mean(abs(x$de.df[top.genes, "lfc"]))
})
de.q.diff = sapply(region.de.genes, function(x){
  if(length(x$genes)==0){
    return(0)
  }
  top.genes = head(x$genes[order(x$de.df[x$genes, "pval"])],100)
  mean(abs(x$de.df[top.genes, "q.diff"]))
})

region.de.summary = data.frame(de.num, de.lfc, de.q.diff)
cl1.genes= sapply(region.de.genes, function(x)paste(head(x$up.genes,10), collapse=";"))
cl2.genes= sapply(region.de.genes, function(x)paste(head(x$down.genes,10), collapse=";"))

region.de.summary$cl1.genes = cl1.genes[row.names(region.de.summary)]
region.de.summary$cl2.genes = cl2.genes[row.names(region.de.summary)]

save(region.de.summary, file="region.de.summary.rda")

write.csv(region.de.summary, file="region.de.summary.csv")

region.de.df = cbind(region.de.df, region.de.summary[row.names(region.de.df),])
save(region.de.df, file="region.de.df.rda")
write.csv(region.de.df, file="region.de.csv")


Ex.shared.region.de.df = region.de.df[with(region.de.df,VISp.cl %in% Ex.types & as.character(VISp.cl)==as.character(ALM.cl)), ]
Ex.region.de.df = region.de.df[region.de.df$VISp.cl %in% Ex.types, ]
Inh.region.de.df = region.de.df[region.de.df$VISp.cl %in% Inh.types, ]
Inh.region.within.de.df = with(Inh.region.de.df, Inh.region.de.df[as.character(VISp.cl)==as.character(ALM.cl),])
summary(Inh.region.de.df)


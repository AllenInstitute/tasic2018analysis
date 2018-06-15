library(feather)
library(Matrix)
source("~/zizhen/My_R/map_river_plot.R")
source("~/zizhen/My_R/sankey_functions.R")

load("cl.clean.rda")
load("norm.dat.rda")
load("samp.dat.rda")
tmp=load("V1.cl.rda")

###load mapping data and preprocessing

hrvatin.samp.dat = read.csv(gzfile("GSE102827_cell_type_assignments.csv.gz"),row.names=1)
hrvatin.counts =  read.csv(gzfile("GSE102827_merged_all_raw.csv.gz"))
row.names(hrvatin.counts)= hrvatin.counts[,1]
hrvatin.counts= as.matrix(hrvatin.counts[,-1])
hrvatin.counts=Matrix(hrvatin.counts, sparse=TRUE)
hrvatin.cl = setNames(hrvatin.samp.dat$celltype,row.names(hrvatin.samp.dat))



###compute overlapping genes present in both dataset, and normalize data on overlapping genes
common.genes = intersect(row.names(hrvatin.counts), row.names(norm.dat))
cells= names(cl)
hrvatin.counts=hrvatin.counts[common.genes,]
hrvatin.dat = t(t(hrvatin.counts)*10^6/colSums(hrvatin.counts))
hrvatin.dat = log2(hrvatin.dat+1)

##Use our clusters as reference, compute markers
V1.cl = row.names(cl.df)[cl.df$ALM < 0.9]
V1.cl = droplevels(cl.clean[cl.clean %in% V1.cl])
V1.markers = select_markers(norm.dat, V1.cl, n.markers=50, de.genes=de.genes)$markers
V1.cells = intersect(names(V1.cl), row.names(samp.dat)[samp.dat$Region=="VISp"])
save(V1.cl, V1.cells, V1.markers, file="V1.cl.rda")

select.markers=intersect(V1.markers, row.names(hrvatin.dat))
select.markers = select.markers[rowMaxs(as.matrix(hrvatin.dat[select.markers,]))>0]

load("V1.cells.rda")
map.result = map_sampling(norm.dat[,V1.cells], droplevels(cl[V1.cells]), hrvatin.dat, markers=select.markers)

map.df = map.result$map.df

map.df$pred_cluster_label = cl.df[as.character(map.df$pred.cl),"cluster_label"]
map.df$cl = hrvatin.samp.dat[row.names(map.df),"celltype"]
map.df$sub.cl = hrvatin.samp.dat[row.names(map.df),"subtype"]
map.df$stim = hrvatin.samp.dat[row.names(map.df),"stim"]


tmp.cor =  rowMaxs(cor(as.matrix(hrvatin.dat[select.markers,]),cl.med[select.markers,]))
tmp.cor[is.na(tmp.cor)] = 0
map.df$cor = tmp.cor
colnames(map.df) = c("map_cluster_id","map_prob","map_cluster_label", "cluster_id","sub_cluster","stim","cor")
map.df$maintype= hrvatin.samp.dat[row.names(map.df),"maintype"]
tmp = compare_annotate(setNames(map.df$cluster_id, row.names(map.df)), setNames(map.df$map_cluster_id,row.names(map.df)), cl.df, reorder=FALSE)
map.df$cluster_label = map.df$cluster_id
map.df$cluster_color = tmp.cl.df[as.character(map.df$map_cluster_id),"cluster_color"]
map.df = droplevels(map.df)
save(map.df, file="map.hrvatin.df.rda")
g=plot_cl_meta_barplot(map.df$map_cluster_label, map.df$stim, c("yellow","green","blue"))
ggsave("hrvatin.stim.pdf",g,height=3,width=12)



ex.map.df = droplevels(map.df[which(map.df$maintype=="Excitatory"),])
g=river_plot(ex.map.df, min.cells=4, min.frac=0.1)
inh.map.df = droplevels(map.df[which(map.df$maintype=="Interneurons"),])
g=river_plot(inh.map.df, min.cells=20, min.frac=0.1)

nn.map.df = droplevels(map.df[map.df$maintype %in% c("Astrocytes","Macrophages","Microglia","Mural","Endothelial_SmoothMuscle","Oligodendrocytes"),])

hrvatin.color = rbind(rbind(unique(ex.map.df[,c("cluster_label","cluster_color")]),unique(inh.map.df[,c("cluster_label","cluster_color")])), unique(nn.map.df[,c("cluster_label","cluster_color")]))

hrvatin.color = setNames(hrvatin.color[,2], gsub("_","", hrvatin.color[,1]))
save(hrvatin.color, file="Hrvatin.color.rda")

#####map the our data to hrvatin.dat
select.cells=with(hrvatin.samp.dat, row.names(hrvatin.samp.dat))
                                        #[maintype %in% c("Interneurons","Excitatory") & !celltype %in% c("RSP","Hip")])

hrvatin.cl = droplevels(setNames(hrvatin.samp.dat[select.cells, "celltype"], select.cells))
new.name = setNames(gsub("_","", levels(hrvatin.cl)), levels(hrvatin.cl))
levels(hrvatin.cl)=new.name
de.param = de_param(q1.th=0.3)

tmp.cells = sample_cells(hrvatin.cl, 500)
select.markers = select_markers(hrvatin.dat[,tmp.cells],  hrvatin.cl[tmp.cells], n.markers=50, de.param = de.param)$markers
select.markers= intersect(select.markers, row.names(norm.dat))


map.result = map_sampling(hrvatin.dat, hrvatin.cl, norm.dat[,V1.cells], markers=select.markers,method="average")
map.df = map.result$map.df
colnames(map.df) = c("map_cluster_label","prob")
map.df$map_cluster_id = as.character(as.integer(map.df$map_cluster_label))
map.df$map_cluster_color= hrvatin.color[as.character(map.df$map_cluster_label)]
map.df$cluster_id = as.integer(as.character(cl[row.names(map.df)]))
map.df$cluster_label = cl.df[as.character(cl[row.names(map.df)]), "cluster_label"]
map.df$cluster_color = cl.df[as.character(cl[row.names(map.df)]), "cluster_color"]
map.df$class_label=  cl.df[as.character(cl[row.names(map.df)]), "class_label"]
save(map.df, file="map.ref.hrvatin.df.rda")
map.df = map.df[map.df$prob > 0.95,]
ex.map.df = droplevels(map.df[which(map.df$class_label=="Glutamatergic"),])
ggsave("ref.hrvatin.ex.riverplot.pdf",river_plot(ex.map.df, min.cells=4, min.frac=0.1))
inh.map.df = droplevels(map.df[which(map.df$class_label=="GABAergic" & map.df$cluster_label!="Meis2 Adamts19"),])
ggsave("ref.hrvatin.inh.riverplot.pdf",river_plot(inh.map.df, min.cells=4, min.frac=0.1))
nn.map.df = droplevels(map.df[which(map.df$class_label%in%c("Non-Neuronal", "Endothelial")),])
save(nn.map.df, file="Hrvatin.nn.map.df.rda")
ggsave("ref.hrvatin.nn.riverplot.pdf",river_plot(nn.map.df, min.cells=2, min.frac=0.05))




map.hrvatin.cl = setNames(map.df$map_cluster_id, row.names(map.df))
 
hrvatin.stim.cl = setNames(paste0(hrvatin.cl,hrvatin.samp.dat[names(hrvatin.cl),"stim"]), names(hrvatin.cl))

stim.early.pairs = data.frame(stim=paste0(levels(hrvatin.cl),"1h"),nostim=paste0(levels(hrvatin.cl),"0h"),stringsAsFactors=FALSE)
stim.early.pairs$cl = levels(hrvatin.cl)
stim.early.pairs$type = "early"
stim.late.pairs = data.frame(stim=paste0(levels(hrvatin.cl),"4h"),nostim=paste0(levels(hrvatin.cl),"0h"),stringsAsFactors=FALSE)
stim.late.pairs$cl = levels(hrvatin.cl)
stim.late.pairs$type = "late"
stim.pairs=rbind(stim.early.pairs, stim.late.pairs)
row.names(stim.pairs) = paste0(stim.pairs$stim,"_", stim.pairs$nostim)

select.cl=with(hrvatin.samp.dat, gsub("_", "",unique(celltype[maintype %in% c("Interneurons","Excitatory") & !celltype %in% c("RSP","Hip","Sub")])))
                                        
select.pairs = stim.pairs[stim.pairs$cl  %in% select.cl,]


de.param = de_param(q1.th=0.05, q.diff.th=0.5, lfc.th=1, padj.th=0.05)
hrvatin.stim.de.df = DE_genes_pairs(as.matrix(hrvatin.dat[common.genes,names(hrvatin.stim.cl)]), cl=hrvatin.stim.cl, pairs=stim.pairs, use.voom = TRUE, counts= as.matrix(hrvatin.counts[common.genes, names(hrvatin.stim.cl)]))
save(hrvatin.stim.de.df, file="hrvatin.stim.de.df.voom.rda")


 
stim.genes <- sapply(c("early","late"),function(cat){
  sapply(select.cl,function(x){
    print(x)
    pairs = row.names(select.pairs)[select.pairs$cl==x & select.pairs$type==cat]
    stim.markers = sapply(pairs, function(s){
      tmp = hrvatin.stim.de.df[[s]]
      tmp = tmp[order(tmp$padj),]
      tmp = tmp[with(tmp, which(abs(lfc)>1 & padj < 0.05)),]
      up.de.df= tmp[tmp$lfc > 0,]
      head(row.names(up.de.df),20)
    },simplify=F)
    markers = intersect(unlist(stim.markers),row.names(norm.dat))  
  },simplify=F)
},simplify=F)


stim.gene.counts=sapply(stim.genes, function(x){
  sapply(x, length)
})
select.cl= row.names(stim.gene.counts)[rowMins(stim.gene.counts)>=4]

tb = table(map.df$cluster_id, map.df$map_cluster_label)
tb = tb[,select.cl]
tb = tb[rowMaxs(tb)>5,]
tb = tb/rowSums(tb)
tb = tb[rowMaxs(tb)>0.9,]
match.cl = setNames(colnames(tb)[apply(tb, 1, which.max)], row.names(tb))


####compute average centered expression value for LRG and EGR 
avg = list(early=c(),late=c())
for(cat in c("early","late")){
  for(x in select.cl){
    select.genes= stim.genes[[cat]][[x]]
    tmp.cl = names(match.cl)[match.cl==x]
    select.cells = with(map.df, row.names(map.df)[which(map_cluster_label==x & cluster_id %in% tmp.cl)])
    tmp.dat = norm.dat[select.genes, select.cells]
    avg[[cat]][select.cells] <-  colMeans(tmp.dat - rowMeans(tmp.dat))
  }
}




map.df[names(avg$early),"ERG.avg"] = avg$early
map.df[names(avg$late),"LRG.avg"] = avg$late
select.map.df = map.df[names(avg$early),]
select.map.df = select.map.df[select.map.df$cluster_id %in% names(match.cl),]
tmp = match.cl[as.character(select.map.df$cluster_id)] == select.map.df$map_cluster_label
select.map.df = droplevels(select.map.df[tmp,])
levels= row.names(cl.df)[row.names(cl.df) %in% names(match.cl)]
select.map.df$cluster_label = factor(as.character(select.map.df$cluster_label), cl.df[levels, "cluster_label"])

match.cl = match.cl[match.cl %in% select.map.df$map_cluster_label]


tmp.df1 = select.map.df[,c("cluster_label", "map_cluster_label", "ERG.avg")]
tmp.df2 = select.map.df[,c("cluster_label", "map_cluster_label", "LRG.avg")]
colnames(tmp.df1)[3] = colnames(tmp.df2)[3] = "avg.exp"
tmp.df1$cat = "ERG"
tmp.df2$cat = "LRG"
tmp.df = rbind(tmp.df1, tmp.df2)
hrvatin.avg.df = tmp.df

match.cl= match.cl[names(match.cl)%in% as.character(select.map.df$cluster_id)]
stats=list(early=list(),late=list())
for(cat in c("early","late")){
  for(x in select.cl){
    tmp.cl = names(match.cl)[match.cl==x]
    if(length(tmp.cl)>1){
      all.cells = intersect(with(map.df, row.names(map.df)[which(map_cluster_label==x & cluster_id %in% tmp.cl)]),names(avg[[cat]]))
    
      for(y in tmp.cl){
        print(y)
        select.cells = intersect(with(map.df, row.names(map.df)[which(map_cluster_label==x & cluster_id ==y )]),names(avg[[cat]]))
        stats[[cat]][[y]] =t.test(avg[[cat]][select.cells], avg[[cat]][setdiff(all.cells, select.cells)] )
        
      }
    }
  }
}
t.stats = sapply(stats, function(x)sapply(x, function(y)y$statistic))
pval = sapply(stats, function(x)sapply(x, function(y)y$p.value))
diff = sapply(stats, function(x)sapply(x, function(y)(y$estimate[1]  - y$estimate[2])))
stats.df=as.data.frame(as.table(pval))
colnames(stats.df)=c("cl","cat","pval")
stats.df$t.stats = as.vector(t.stats)
stats.df$diff = as.vector(diff)
stats.df$cluster_label = cl.df[as.character(stats.df$cl),"cluster_label"]
save(stats.df, file="hrvatin.genes.stats.df.rda")


select.stats.df =   stats.df %>% filter(pval < 10^-5 & abs(diff) >1)
select.stats.df$sig_label = "*"
select.stats.df$sig_label[select.stats.df$pval < 10^-10] = "**"
select.stats.df$sig_label[select.stats.df$pval < 10^-20] = "***"
select.stats.df$avg.exp = 5
select.stats.df[select.stats.df$diff < -1,"avg.exp"] =-5
levels(select.stats.df$cat)=c("ERG","LRG")

tmp.level= c("IntSst2","IntPv","ExcL23","ExcL4","ExcL51","ExcL53","ExcL6")

hrvatin.avg.df$map_cluster_label = factor(hrvatin.avg.df$map_cluster_label, levels= tmp.level)
g = ggplot(hrvatin.avg.df, aes(x=cluster_label, y=avg.exp)) + geom_violin(aes(color=map_cluster_label))
g = g + geom_text(aes(x= cluster_label, y=avg.exp,label=sig_label),data=select.stats.df)
g = g + facet_wrap(~cat, nrow=2)
g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("hrvatin.avg.pdf",g, height=5, width=12)
save(hrvatin.avg.df, stats.df, file="hrvatin.avg.df.rda")
 
library(viridis)
for(x in c("IntSst2","ExcL23","ExcL51")){
  ERG = stim.genes[["early"]][[x]]
  LRG = stim.genes[["late"]][[x]]
  common = intersect(ERG, LRG)
  ERG = setdiff(ERG, common)
  LRG = setdiff(LRG, common)
  select.cells = with(select.map.df, row.names(select.map.df)[which(map_cluster_label==x)])
  tmp.dat= norm.dat[c(ERG,common, LRG),select.cells]
  cl.means= get_cl_means(tmp.dat - rowMeans(tmp.dat), droplevels(cl[select.cells]))
  colnames(cl.means)= cl.df[colnames(cl.means),"cluster_label"]
  select.cells2 = names(hrvatin.cl)[which(hrvatin.cl==x)]
  tmp.dat2= hrvatin.dat[c(ERG,common, LRG),select.cells2]
  
  hrvatin.cl.means = get_cl_means(tmp.dat2 - rowMeans(tmp.dat2), setNames(hrvatin.samp.dat[select.cells2, "stim"], select.cells2))

  df1 = as.data.frame(as.table(cl.means))
  df1$cat = "AIBS"
  df2 = as.data.frame(as.table(hrvatin.cl.means))
  df2$cat = "Hrvatin"

  df = rbind(df1, df2)
  colnames(df)=c("gene","cl","Exp","cat")
  df$Exp[df$Exp > 4]=4
  g = ggplot(df, aes(cl, gene)) + geom_tile(aes(fill=Exp)) + facet_wrap(~cat,scales="free_x")
  #g = g + scale_fill_gradient2(low = "darkgreen", mid ="",high = "lightyellow")
  g = g + scale_fill_viridis()
  g = g + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(x, ".heatmap.pdf"),g)
}
    









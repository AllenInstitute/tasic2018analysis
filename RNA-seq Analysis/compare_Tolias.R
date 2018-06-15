library(iterclust)
load("norm.dat.rda")
load("cl.final.rda")
tmp=load("V1.cl.rda")
load("cl.med.rda")

d = "/allen/programs/celltypes/workgroups/rnaseqanalysis/osnat/patchseq_simple_mapping/results/Tolias"
tmp=load(file.path(d, "quick_map.Tolias.rda"))

shiny.d="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/tolias_patchseq_20161102"
Tolias.df = as.data.frame(read_feather(file.path(shiny.d, "anno.feather")))
row.names(Tolias.df) = Tolias.df$sample_id


Tolias.counts = as.matrix(read.csv(file.path(d,"counts.csv"),row.names=1))

genes=row.names(Tolias.counts)
genes[genes=="6330527O06Rik"]="Lamp5"
genes[genes=="A930038C07Rik"]="Ndnf"
genes[genes=="Tmem90b"]="Syndig1"
genes[genes=="Cxcr7"]="Ackr3"
row.names(Tolias.counts)=genes


common.genes= intersect(genes, row.names(norm.dat))
Tolias.counts <- as.matrix(Tolias.counts[common.genes,])
Tolias.dat =log2(t(t(Tolias.counts)*10^6/colSums(Tolias.counts))+1)



common.markers= intersect(V1.markers, common.genes)
map.result = map_sampling(as.matrix(norm.dat[common.markers,V1.cells]), droplevels(cl[V1.cells]), Tolias.dat, markers = common.markers)
map.df = map.result$map.df
map.freq = map.result$map.freq

tmp.cor = rowMaxs(cor(Tolias.dat[common.markers,], cl.med[common.markers, droplevels(cl[V1.cells])]))
tmp.cor[is.na(tmp.cor)]=0
map.df$cor = tmp.cor

map.df$cluster_label = as.character(cl.df[as.character(map.df$pred.cl),"cluster_label"])
tmp = is.na(map.df$cluster_label)
map.df$cluster_label[tmp] = map.df$cl[tmp]
map.df$class_label = Tolias.df[row.names(map.df), "class_label"]

dend=readRDS("dend.RData")
labels(dend) = row.names(cl.df)[match(labels(dend), cl.df$cluster_label)]
V1.dend = prune_dend(dend, setdiff(labels(dend), levels(V1.cl)))


summarize_cl <- function(dend, map.freq,conf.th=0.7)
  {
    node.height=setNames(get_nodes_attr(dend, "height"),get_nodes_attr(dend, "label"))
    dend.list = dend_list(dend)
    dend.list = dend.list[sapply(dend.list, length)>1]
    memb = sapply(names(dend.list),function(x){
      rowSums(map.freq[,intersect(labels(dend.list[[x]]),colnames(map.freq)),drop=F])
    })
    memb = cbind(memb,map.freq)
    memb.th= lapply(row.names(memb),function(cell){
      ###Check all the node with confidence > conf.th
      x = memb[cell,]
      mapped.node = colnames(memb)[which(x>conf.th)]
               
      ###mapped nodes not met the minimal gene number/ratio constraints 
      ###Choose the deepest nodes that pass all the criteria. 
      mapped.node=mapped.node[order(node.height[mapped.node])]
      i=mapped.node[1]
      ###Get the markers on every mapped nodes. 
      c(cl=i, score=x[i])
    })
    memb.th = do.call("rbind",memb.th)
    row.names(memb.th) = row.names(memb)
    colnames(memb.th)=c("cl","score")
    memb.df = as.data.frame(memb.th)
    memb.df$resolution.index = 1- (node.height[memb.df$cl]/attr(dend,"height"))
    return(memb.df) 
  }

 
map.tree.df = summarize_cl(V1.dend, map.freq/100, 0.7)
map.tree.df$cluster_label = as.character(cl.df[as.character(map.tree.df$cl),"cluster_label"])
tmp = is.na(map.tree.df$cluster_label)
map.tree.df$cluster_label[tmp] = as.character(map.tree.df$cl[tmp])
map.tree.df$class_label = Tolias.df[row.names(map.tree.df), "class_label"]
with(map.tree.df, table(class_label, cl))

save(map.tree.df, V1.dend, file="map.tree.df.rda")


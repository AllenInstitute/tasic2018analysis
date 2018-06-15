cl.df$subclass_region_label = paste(cl.df$subclass_label, cl.df$region_label)


#####Find markers to annotate cluster
cl.present.counts= get_cl_sums(norm.dat > de.param$low.th, cl)
#####Find tau markers among sibling
low.th = setNames(pmax(rowMaxs(cl.med) - 10, 1), row.names(cl.med))
de.param$low.th = low.th
de.param$q1.th=0.4
specific.markers <- within_group_specific_markers(levels(cl.clean), norm.dat, cl.clean, de.param = de.param, cl.present.counts=cl.present.counts, n.markers=5)
cl.df.clean = cl.df[levels(cl.clean),]
within.class.specific.markers <- do.call("rbind",sapply(unique(cl.df.clean$subclass_region_label), function(g){
  cl.g = row.names(cl.df.clean)[cl.df.clean$subclass_region_label==g]
  print(g)
  if(length(cl.g)>1){
    df=within_group_specific_markers(cl.g, norm.dat, cl.clean, de.param = de.param, cl.present.counts=cl.present.counts, n.markers=5)
    df$class = g
    df
  }
},simplify=F))

specific.markers = with(specific.markers, specific.markers[specificity > 0.95,])
write.csv(specific.markers, file="specific.markers.csv")


within.class.specific.markers =  with(within.class.specific.markers, within.class.specific.markers[specificity > 0.75,])
write.csv(within.class.specific.markers, file="within.class.specific.markers.csv")


marker.list1 = with(specific.markers, split(g, cl))
marker.list1.str = sapply(marker.list1, function(x){paste(head(x,4), collapse=",")})

marker.list2 = with(within.class.specific.markers, split(g, cl))
marker.list2.str = sapply(marker.list2, function(x){paste(head(x,4), collapse=",")})


dend.list = dend_list(dend)
class.region.list <- split(row.names(cl.df.clean),cl.df.clean$subclass_region_label)
within.class.nodes = sapply(names(dend.list), function(x){
  class.region.label= unique(cl.df[labels(dend.list[[x]]), "subclass_region_label"])
  if(length(class.region.label)==1){
    if(length(labels(dend.list[[x]]))  < length(class.region.list[[class.region.label]])){
      return(dend.list[[x]])
    }
  }
  NULL
})
library(dendextend)
within.class.nodes = within.class.nodes[sapply(within.class.nodes, length)>0]
within.class.nodes = within.class.nodes[sapply(within.class.nodes, function(x) nleaves(x)>1 & nleaves(x)<5)]

within.class.node.specific.markers <- do.call("rbind",sapply(within.class.nodes, function(x){
  cl.g = labels(x)
  tmp.cl = droplevels(cl.clean[cl.clean %in% row.names(cl.df.clean)[cl.df.clean$subclass_region_label %in% cl.df.clean[cl.g,"subclass_region_label"]]])
  df=group_specific_markers(cl.g, norm.dat, cl=tmp.cl, de.param = de.param, cl.present.counts=cl.present.counts, n.markers=5)
  if(!is.null(df)){
    df$g = row.names(df)
    df$node = attr(x, "label")
    df$groups = paste(cl.g, collapse=",")
  }
  df
},simplify=F))
within.class.node.specific.markers = within.class.node.specific.markers[within.class.node.specific.markers$specificity > 0.75,]
write.csv(within.class.node.specific.markers,"within.class.node.specific.markers.csv")

cl.size=table(cl)
cl.present = t(t(cl.present.counts) / as.vector(cl.size[colnames(cl.present.counts)]))
tmp = do.call("rbind", sapply(names(within.class.nodes), function(n){
  cl.g = labels(dend.list[[n]])
  print(cl.g)
  if(length(cl.g) <= 3) {
    tmp.df = do.call("rbind",sapply(cl.g, function(l){
      tmp.df = within.class.node.specific.markers[within.class.node.specific.markers$node==n,]
      g = tmp.df$g[cl.present[tmp.df$g, l] >= de.param$q1.th]
      tmp.df = tmp.df[tmp.df$g %in%g,]
      if(nrow(tmp.df)>0){
        tmp.df$cl = l
        tmp.df$group.size = length(cl.g)
      }
      tmp.df
    },simplify=F))
  }
  else{
    NULL
  }
}))
tmp = with(tmp, tmp[order(as.integer(cl), -group.size, -specificity),])
tmp = tmp[!duplicated(paste(tmp$cl, tmp$g)),]
marker.list3 = split(tmp$g, tmp$cl)
marker.list3.str = sapply(marker.list3, function(x){paste(head(x,4), collapse=",")})


cl.df$specific.markers = marker.list1.str[row.names(cl.df)]
cl.df$within.class.markers = marker.list2.str[row.names(cl.df)]
cl.df$within.class.node.markers = marker.list3.str[row.names(cl.df)]




tmp=get_gene_score(de.genes)
up.gene.score=tmp$up.gene.score
down.gene.score=tmp$down.gene.score
node.markers <- findNodeSpecificMarkers(dend[[1]], norm.dat, cl.clean, cl.df, de.genes=de.genes, up.gene.score=up.gene.score, down.gene.score=down.gene.score)


node.vs.sibling.markers <- node_vs_sibling_markers(dend.list, norm.dat, cl,cl.df, de.param=de.param, cl.present.counts=cl.present.counts,n.markers=10)
node.vs.sibling.markers = with(node.vs.sibling.markers, node.vs.sibling.markers[specificity > 0.75 & node %in% names(within.class.nodes),])
tmp = do.call("rbind", tapply(1:nrow(node.vs.sibling.markers), node.vs.sibling.markers$node, function(x){
  n = unique(node.vs.sibling.markers[x, "node"])
  cl.g = labels(dend.list[[n]])
  if(length(cl.g) <= 3) {
    tmp.df = do.call("rbind",sapply(cl.g, function(l){
      tmp.df = node.vs.sibling.markers[x,]
      tmp.df$cl = l
      tmp.df
    },simplify=F))
  }
  else{
    NULL
  }
}))
present = get_pair_matrix(cl.present, as.character(tmp$g), as.character(tmp$cl))
tmp$cl.freq  = present
tmp = tmp[present > de.param$q1.th,]
marker.list3 = tapply(1:nrow(tmp), tmp$cl, function(x){
  tmp=split(tmp[x,"g"], tmp[x, "node"])
  unlist(lapply(tmp, function(x)head(x, 2)))
})


cl.sibling.markers=sapply(marker.list3, function(x)paste(x, collapse=","))
cl.df$sibling.markers = cl.sibling.markers[row.names(cl.df)]

write.csv(cl.df, file="cl.df.csv")

###find pairs that are fine splits
pairs =split(row.names(cl.df.clean), cl.df.clean$ref.cl)
pairs = pairs[sapply(pairs, length)==2]
pairs[["51"]] = c("58","59")
p = sapply(pairs, paste, collapse="_")
p.name = sapply(pairs, function(x)paste(cl.df[x,"cluster_label"],collapse="_"))
tmp.de.genes = de.genes.sym[p]
up.genes= sapply(tmp.de.genes, function(x)paste(head(x$up.genes, 4),collapse=","))
down.genes= sapply(tmp.de.genes, function(x)paste(head(x$down.genes, 4),collapse=","))
split.pair = data.frame(up.genes, down.genes)
row.names(split.pair)=p.name
write.csv(split.pair, file="split.pair.csv")




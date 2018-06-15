library(iterclust)
load("cl.final.rda")
load("norm.dat.rda")
load("top.cl.rda")
load("rm.eigen.rda")
ref.cl = cl
ref.cl.df = cl.df



select.cells=names(cl)[cl  %in% row.names(cl.df)[cl.df$subclass_label=="Sst"]]
save(select.cells, file="Sst.select.cells.rda")
all.cells=select.cells

top.result = list(cl=as.factor(top.cl[select.cells]),markers=NULL)
de.param = de_param(q1.th=0.5, q.diff.th=0.7, min.cells=4)

library(parallel)
nodes <- makeCluster(12, type="FORK")
clusterExport(nodes,c("norm.dat", "all.cells","rm.eigen","top.result"))

for(de.score.th in c(80,300)){
  de.param$de.score.th = de.score.th 
  clusterExport(nodes,c("de.param"))
  d = paste0("subsample_pca_",de.score.th)
  dir.create(d)
  clusterExport(nodes,c("d"))
  tmp= parSapply(nodes, 1:100, function(i){
    library(iterclust)
    prefix = paste("iter",i,sep=".")
    outfile= file.path(d, paste0("result.",i,".rda"))
      if(file.exists(outfile)){
      return(NULL)
    }
    select.cells=sample(all.cells, round(length(all.cells)*0.8))
    save(select.cells, file=file.path(d, paste0("cells.",i,".rda")))
    result <- iter_clust(norm.dat=norm.dat, select.cells=select.cells,prefix=prefix, split.size = 10, de.param = de.param, dim.method="pca",result= top.result,rm.eigen=rm.eigen, rm.th=0.7)
    save(result, file=outfile)
  })
  result.files=file.path(d, dir(d, "result.*.rda"))
  PCA.co.result <- collect_co_matrix_sparseM(norm.dat, result.files, all.cells)
  
  d = paste0("subsample_WGCNA_", de.score.th)
  dir.create(d)
  clusterExport(nodes,c("d"))
  tmp= parSapply(nodes, 1:100, function(i){
    library(iterclust)
    prefix = paste("iter",i,sep=".")
    outfile= file.path(d, paste0("result.",i,".rda"))
    if(file.exists(outfile)){
      return(NULL)
    }
    select.cells=sample(all.cells, round(length(all.cells)*0.8))
    save(select.cells, file=file.path(d, paste0("cells.",i,".rda")))
    result <- iter_clust(norm.dat=norm.dat, select.cells=select.cells,prefix=prefix, split.size = 10, de.param = de.param, dim.method="WGCNA",result= top.result,rm.eigen=rm.eigen, rm.th=0.7)
    save(result, file=outfile)
  })
  result.files=file.path(d, dir(d, "result.*.rda"))
  WGCNA.co.result <- collect_co_matrix_sparseM(norm.dat, result.files, all.cells)

  cl.list = c(PCA.co.result$cl.list, WGCNA.co.result$cl.list)
  co.ratio.min = pmin(PCA.co.result$co.ratio, WGCNA.co.result$co.ratio)
  cl.mat = cbind(PCA.co.result$cl.mat, WGCNA.co.result$cl.mat)
  
  consensus.result = iter_consensus_clust(co.ratio.min, cl.list, norm.dat, select.cells=all.cells, de.param = de.param)

  refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio.min, tol.th=0.005, confusion.th=0.8)
  merge.result= merge_cl(norm.dat=norm.dat, cl=refine.result$cl, rd.dat=t(norm.dat[consensus.result$markers,]), de.param = de.param,return.markers=FALSE)
  compare.result = compare_annotate(merge.result$cl, ref.cl, ref.cl.df)
  cl = compare.result$cl
  cl.df = compare.result$cl.df
  cl = setNames(factor(as.character(cl),levels=row.names(cl.df)), names(cl))
  save(cl, cl.df, file=paste0("Sst.cl.", de.score.th, ".rda"))
}

ref.cl = droplevels(ref.cl[all.cells])
ref.cl.df = ref.cl.df[levels(ref.cl),]




source("~/zizhen/My_R/map_river_plot.R")
source("~/zizhen/My_R/sankey_functions.R")
map.df = ref.cl.df[as.character(ref.cl),c(2,3,5:11)]
row.names(map.df) = names(ref.cl)
colnames(map.df) = paste0("map_",colnames(map.df))

load("Sst.cl.80.rda")
cl.df$cluster_id = as.integer(row.names(cl.df))
compare.result = compare_annotate(cl, ref.cl, ref.cl.df)

ggsave("Sst.cl.80.map.pdf", compare.result$g)
map.80.df = cbind(map.df, cl.df[as.character(cl[row.names(map.80.df)]),])
save(map.80.df, file="map.80.df.rda")


g80 <- river_plot(map.80.df, min.cells=4, min.frac=0.15)
ggsave("Sst.80.river.pdf", g80)

load("Sst.cl.300.rda")
map.300.df = cbind(map.df, cl.df[as.character(cl[row.names(map.80.df)]),c(3:4,6:12)])


g300 <- river_plot(map.300.df, min.cells=4, min.frac=0.1)
ggsave("Sst.300.river.pdf", g300)
save(map.300.df, file="map.300.df.rda")


source("~/zizhen/My_R/hicat/R/cl.transition.R")
for(de.score.th in c(80,300)){
  load(paste0("Sst.cl.", de.score.th, ".rda"))
  select.markers= select_markers(norm.dat, cl, n.markers=50, de.param = de.param)$markers
  print(de.score.th)
  print(length(cl))
  print(length(select.markers))
  cell.cl.map.df = get_core_transition(norm.dat, cl, select.markers, n.bin=5, n.iter=100, mc.cores=10)
  save(cell.cl.map.df, file=paste0("Sst.cell.cl.map.", de.score.th, ".rda"))
}

  

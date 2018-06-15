library(iterclust)
library(matrixStats)

tmp=load("norm.dat.rda")
load("select.cells.rda")
load("top.cl.rda")
load("rm.eigen.rda")

library(feather)
anno= read_feather("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170913/anno.feather")
ref.cl = setNames(anno$cluster_id, anno$sample_id)
ref.cl.df = unique(as.data.frame(anno[,c("cluster_id","cluster_color","cluster_label","class_id","class_label","class_color","category_id","category_label","category_color")]))
ref.cl.df = ref.cl.df[order(ref.cl.df$cluster_id),]
row.names(ref.cl.df) = ref.cl.df$cluster_id
load("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170913/select.markers.rda")
ref.markers= select.markers

all.cells=select.cells
top.result = list(cl=as.factor(top.cl[select.cells]),markers=NULL)
de.param = de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=4)

d = paste0("subsample_pca")
dir.create(d)


library(parallel)
cl <- makeCluster(12, type="FORK")
clusterExport(cl,c("norm.dat", "de.param", "all.cells","rm.eigen","top.result"))

tmp= parSapply(cl, 1:100, function(i){
  library(iterclust)
  d = paste0("subsample_pca")
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

d = paste0("subsample_WGCNA")
dir.create(d)
tmp= parSapply(cl, 1:100, function(i){
  library(iterclust)
  d = paste0("subsample_WGCNA")
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

####Computing PCA-based consensus clustering result
d = "subsample_pca/"
result.files=file.path(d, dir(d, "result.*.rda"))
co.result <- collect_co_matrix_sparseM(norm.dat, result.files, all.cells)
save(co.result, file="pca.co.result.rda")
co.ratio = co.result[[1]]
   
consensus.result = iter_consensus_clust(co.ratio, co.result$cl.list, norm.dat, select.cells=all.cells, de.param = de.param)

save(consensus.result, file="pca.consensus.result.rda")
cl = consensus.result$cl


refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio, tol.th=0.005, confusion.th=0.8)
cl = refine.result$cl
merge.result= merge_cl(norm.dat=norm.dat, cl=cl, rd.dat=t(norm.dat[consensus.result$markers,]), de.param = de.param,return.markers=FALSE)
compare.result = predict_annotate_cor(merge.result$cl, ref.markers, ref.cl, ref.cl.df,norm.dat)
ggsave("map.PCA.pdf", compare.result$g)
cl = compare.result$cl
cl.df = compare.result$cl.df

#tmp.cl = droplevels(cl[cl %in% as.character(c(29:50))])
#tmp= display_cl(cl=tmp.cl, norm.dat=norm.dat, prefix="Sst",de.param=de.param, max.cl.size=300, col=all.col)

cl.df = read.csv("cl.df.csv",row.names=1)
cl = setNames(factor(as.character(cl),levels=row.names(cl.df)), names(cl))
save(cl, cl.df, file="pca.cl.final.rda")


####Computing WGCNA-based consensus clustering result
d = "subsample_WGCNA/"
result.files=file.path(d, dir(d, "result.*.rda"))
co.result <- iterclust::collect_co_matrix_sparseM(norm.dat, result.files, all.cells)
save(co.result, file="WGCNA.co.result.rda")
co.ratio = co.result[[1]]

consensus.result = iter_consensus_clust(co.ratio, co.result$cl.list, norm.dat, select.cells=all.cells, de.param = de.param)
save(consensus.result, file="WGCNA.consensus.result.rda")
cl = consensus.result$cl
refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio, tol.th=0.005, confusion.th=0.8)
cl = refine.result$cl

merge.result= merge_cl(norm.dat=norm.dat, cl=cl, rd.dat=t(norm.dat[consensus.result$markers,]), de.param = de.param,return.markers=FALSE,verbose=1)
compare.result = predict_annotate_cor(merge.result$cl, ref.markers, ref.cl, ref.cl.df,norm.dat)
ggsave("map.WGCNA.pdf", compare.result$g,height=10, width=10)
cl = compare.result$cl
cl.df = compare.result$cl.df
save(cl, cl.df, file="WGCNA.cl.final.rda")


WGCNA.cl = cl
WGCNA.cl.df = cl.df

load("pca.cl.final.rda")
write.csv(cl.df, "PCA.cl.df.csv")
PCA.cl = cl
PCA.cl.df = cl.df

compare.result = compare_annotate(WGCNA.cl, PCA.cl, PCA.cl.df)
ggsave("WGCNA.PCA.compare.pdf", compare.result$g, height=10, width=10)



####Reach consensus clusters
load("pca.co.result.rda")
PCA.co.result = co.result
load("WGCNA.co.result.rda")
WGCNA.co.result = co.result

cl.list = c(PCA.co.result$cl.list, WGCNA.co.result$cl.list)
co.ratio.min = pmin(PCA.co.result$co.ratio, WGCNA.co.result$co.ratio)

cl.mat = cbind(PCA.co.result$cl.mat, WGCNA.co.result$cl.mat)
#co.ratio = crossprod(t(cl.mat))
#co.ratio@x = co.ratio@x/length(cl.list)


 
de.param = de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=4)
consensus.result = iter_consensus_clust(co.ratio.min, cl.list, norm.dat, select.cells=all.cells, de.param = de.param)

load("bak/cl.final.rda")
compare.result = compare_annotate(consensus.result$cl, cl, cl.df)
ggsave("map.new3.pdf", compare.result$g,height=10,width=10)


refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio, tol.th=0.005, confusion.th=0.8)
merge.result= merge_cl(norm.dat=norm.dat, cl=refine.result$cl, rd.dat=t(norm.dat[consensus.result$markers,]), de.param = de.param,return.markers=FALSE)

compare.result = compare_annotate(merge.result$cl, cl, cl.df)
ggsave("map.new4.pdf", compare.result$g,height=10,width=10)

cl = compare.result$cl
cl.df = compare.result$cl.df


cl.region =table(cl, samp.dat[names(cl), "Region"])
cl.region = as.data.frame.matrix(round(cl.region/rowSums(cl.region),digits=2))
cl.df[,colnames(cl.region)] = cl.region[row.names(cl.df),]

cl.gene.counts=tapply(names(cl),cl, function(x)mean(samp.dat[x, "Genes.Detected"]))
cl.df$gene.counts=round(cl.gene.counts[row.names(cl.df)])
write.csv(cl.df, "cl.df.new.csv")
save(cl, cl.df, file="cl.rda")

all.col[is.na(all.col)]="black"
#####Examine the clusters
tmp.cl = droplevels(cl[cl %in% as.character(80:88)])
tmp= display_cl(cl=tmp.cl, norm.dat=norm.dat, prefix="tmp",de.param=de.param, max.cl.size=300, col=all.col)

de.genes = de_score(norm.dat, cl, de.param = de.param)
cl.clean = droplevels(cl[cl %in% row.names(cl.df)[cl.df$category_label!="Noise"]])

tmp = select_markers(norm.dat, cl.clean, n.markers=50, de.genes=de.genes)
select.markers= tmp$markers
de.genes.clean = tmp$de.genes

cl.med <- do.call("cbind",tapply(names(cl), cl, function(x){
  rowMedians(as.matrix(norm.dat[,x]))
}))
row.names(cl.med) = row.names(norm.dat)
save(cl.med, file="cl.med.rda")

load("cl.med.rda")
l.label = setNames(cl.df.cluster_labels, row.names(cl.df))
l.rank = setNames(1:nrow(cl.df), row.names(cl.df))
l.color = setNames(as.character(cl.df$cluster_color),row.names(cl.df))
select.cl = levels(cl.clean)
cl.med = cl.med[,select.cl]
colnames(cl.med)= l.label[colnames(cl.med)]
dend.result <- build_dend(cl.med[select.markers,], l.rank, l.color, nboot=0)
dend = dend.result$dend
cl.cor = dend.result$cl.cor
 
save(cl.cor, file="cl.cor.rda")
pdf("cl.cor.pdf",height=5,width=12)
plot(dend)
dev.off()
pdf("cl.cor.heatmap.pdf",height=12,width=12)
heatmap.3(cl.cor, Rowv=dend, Colv=dend, trace="none",col=jet.colors(100),cexCol=0.5,cexRow=0.5,breaks=c(-0.2, 0.2, seq(0.2, 1, length.out=99)))
dev.off()
save(dend, file="dend.rda")



common.cl = cl
common.cl.df = cl.df
common.cl = setNames(factor(as.character(common.cl), levels=row.names(common.cl.df)),names(common.cl))
tmp=load("pca.cl.final.rda")
compare.result = compare_annotate(cl, common.cl, common.cl.df)
ggsave("PCA.compare.pdf", compare.result$g, height=10, width=10)
pca.cl = compare.result$cl
pca.cl.df = compare.result$cl.df







tmp.cl = droplevels(common.cl[common.cl %in% as.character(126:149)])
tmp = de_score(norm.dat, cl=tmp.cl, de.param = de.param)

tmp.cells = names(tmp.cl)
consensus.result = iter_consensus_clust(co.ratio.min, cl.list, norm.dat, select.cells= tmp.cells, de.param = de.param, all.col=all.col, verbose=1)
compare.result = compare_annotate(consensus.result$cl, common.cl, common.cl.df)
ggsave("tmp.pdf", compare.result$g)

compare.result = compare_annotate(result$cl, common.cl, common.cl.df)
ggsave("tmp.pdf",compare.result$g)


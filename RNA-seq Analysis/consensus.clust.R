library(iterclust)
library(matrixStats)

load("norm.dat.rda")
load("select.cells.rda")
load("top.cl.rda")
load("rm.eigen.rda")

###Use previous version of clustering for annotation. 
library(feather)
anno= read_feather("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180312/anno.feather")
ref.cl = setNames(anno$cluster_id, anno$sample_id)
ref.cl.df = unique(as.data.frame(anno[,c("cluster_id","cluster_color","cluster_label","class_id","class_label","class_color","subclass_id","subclass_label","subclass_color")]))
ref.cl.df = ref.cl.df[order(ref.cl.df$cluster_id),]
row.names(ref.cl.df) = ref.cl.df$cluster_id
load("../process_22679/select.markers.rda")
ref.markers= select.markers


####Prepare parameters for clustering. 
all.cells=select.cells
top.result = list(cl=as.factor(top.cl[select.cells]),markers=NULL)
de.param = de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=4)



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
compare.result = compare_annotate(merge.result$cl, ref.cl, ref.cl.df)
ggsave("map.PCA.pdf", compare.result$g)
cl = compare.result$cl
cl.df = compare.result$cl.df
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
compare.result = compare_annotate(merge.result$cl, ref.cl, ref.cl.df)
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


co.ratio.min = pmin(PCA.co.result$co.ratio, WGCNA.co.result$co.ratio)
save(co.ratio.min, file="co.ratio.min.rda")

cl.list = c(PCA.co.result$cl.list, WGCNA.co.result$cl.list)
cl.mat = cbind(PCA.co.result$cl.mat, WGCNA.co.result$cl.mat)
co.ratio.comb = crossprod(t(cl.mat))
co.ratio.comb@x = co.ratio.comb@x/length(cl.list)
save(co.ratio.comb, file="co.ratio.comb.rda")

de.param = de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=4)
consensus.result = iter_consensus_clust(co.ratio.min, cl.list, norm.dat, select.cells=all.cells, de.param = de.param)

refine.result = refine_cl(consensus.result$cl, co.ratio=co.ratio.min, tol.th=0.005, confusion.th=0.6)
merge.result= merge_cl(norm.dat=norm.dat, cl=refine.result$cl, rd.dat=t(norm.dat[consensus.result$markers,]), de.param = de.param,return.markers=FALSE)
save(merge.result, file="merge.result.rda")


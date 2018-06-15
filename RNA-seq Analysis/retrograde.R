tmp=with(samp.dat, is.na(Injection_type) | Injection_type == "" )
samp.dat[tmp,"Injection_type"]="0"
tb = table(ex.cl, samp.dat[names(ex.cl),"Injection_type"])
select.cl = row.names(tb)[rowMins(tb) > 20]
select.cells= names(cl)[cl %in% select.cl]

de.param = de_param(q1.th=0.5, q.diff.th=0.7, de.score.th=150, min.cells=4)
###DEG between V1 and ALM within a mixed Inhibitory types.

retro.cl = setNames(paste0(cl[select.cells], samp.dat[select.cells, "Injection_type"]), select.cells)

tmp= levels(droplevels(cl[select.cells]))
retro.pairs = data.frame(null=paste0(tmp,"0"), retro=paste0(tmp,"retrograde"), stringsAsFactors=FALSE)

row.names(retro.pairs) = paste(retro.pairs[,1], retro.pairs[,2], sep="_")

load("select.genes.rda")
retro.de.genes = de_score_pairs(norm.dat[select.genes,select.cells], cl=retro.cl, pairs=retro.pairs, de.param = de.param)$de.genes
unlist(sapply(retro.de.genes, function(x)x$num))

tmp=sapply(retro.de.genes, function(x)x$genes)
tmp[sapply(tmp,length)>0]
tmp=sapply(retro.de.genes, function(x)x$up.genes)
tmp[sapply(tmp,length)>0]
tmp=sapply(retro.de.genes, function(x)x$down.genes)
tmp[sapply(tmp,length)>0]


retro.cells = intersect(row.names(samp.dat)[samp.dat$Injection_type=="retrograde"], names(cl.clean))
tb=table(droplevels(cl.clean[retro.cells]), samp.dat[retro.cells, "injection_roi"])


  
tb = tb[rowSums(tb > 20) > 1,]
select.cl= row.names(tb)
de.df=sapply(select.cl, function(x){
  select.roi = colnames(tb)[tb[x, ] > 20]
  tmp.cells= retro.cells[as.character(cl.clean[retro.cells])==x & samp.dat[retro.cells,"injection_roi"] %in% select.roi]
  select.roi = setNames(samp.dat[tmp.cells, "injection_roi"], tmp.cells)
  de.df = DE_genes_pw(norm.dat, select.roi, counts=as.matrix(exon[select.genes, names(select.roi)], use.voom=TRUE))
  de.df=sapply(de.df, function(x){
    x=x[order(x$padj,-abs(x$lfc)),]
    x=x[x$padj < 0.05 & abs(x$lfc)>1,]
  },simplify=F)
},simplify=F)
save(de.df, file="retro.roi.de.df.rda")
tmp=sapply(de.df, function(x)sapply(x, nrow))

cl.low = droplevels(cl[!names(cl) %in% names(cl.clean)])
save(cl.low, file="cl.low.rda")

retro.cells = intersect(names(cl.clean), row.names(samp.dat)[samp.dat$Injection_type=="retrograde"]),
non.retro.cells = setdiff(names(cl.clean),retro.cells)
  
clean.cells=split(non.retro.cells, paste(samp.dat[non.retro.cells,"Region"],cl.df[as.character(cl.clean[non.retro.cells]),"class_label"],sep="_"))
clean.cells=split(non.retro.cells, paste(samp.dat[non.retro.cells,"Region"],samp.dat[non.retro.cells, "full_genotype"]))

cells.list = list(control=control.cells, fail.qc = fail.qc.cells, doublet = rm.cells, cl.low = names(cl.low), retro = retro.cells, non.retro.cells = non.retro.cells)
save(cells.list, file="cells.list.rda")


old.cl = as.data.frame(read_feather("/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170913/anno.feather"))
row.names(old.cl) = old.cl[,1]
old.cl[[1]]=NULL
old.cl= old.cl[,c("class_label","cluster_label","region_label","inj_type_label")]

absent.cells= setdiff(row.names(old.cl), names(cl))
tb= as.data.frame(table(old.cl[absent.cells, "cluster_label"]))
tb2 = as.data.frame(table(old.cl[,"cluster_label"]))

absent.df = inner_join(tb, tb2, by="Var1")
colnames(absent.df)=c("cluster","absent","total")
write.csv(absent.df, "absent.cells.csv")




  







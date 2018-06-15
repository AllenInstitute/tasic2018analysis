#####internal heterogenity#######
pca_hetero <- function(norm.dat,select.cells, perm.num=10,vg.padj=0.5,rm.gene.mod= NULL, rm.eigen=NULL,rm.th=0.6)
  {
    vg= findVG(2^norm.dat[,select.cells]-1)
    select.genes = row.names(vg)[vg$loess.padj < vg.padj]
    if(length(select.genes)==0){
      return(0)
    }
    if(length(select.genes)==1){
      return(1)
    }
    pca = prcomp(t(as.matrix(norm.dat[select.genes, select.cells])),tol=0.01)
    if(!is.null(rm.gene.mod)|!is.null(rm.eigen)){
      if(is.null(rm.eigen)){
        rm.eigen = get_eigen(rm.gene.mod, norm.dat, select.cells)[[1]]
      }
      else{
        rm.eigen = rm.eigen[select.cells,]
      }
      rm.cor=cor(pca$x, rm.eigen)
      rm.score =rowMaxs(abs(rm.cor))
      select.pca = which(rm.score < rm.th)[1]
    }
    pca.var = summary(pca)$importance[2,select.pca]
    ###shuffle expression of every gene.
    perm.pca.var = sapply(1:10, function(x){
      perm.dat = norm.dat[select.genes, select.cells]
      for(i in 1:nrow(perm.dat)){
        perm.dat[i, ] = sample(perm.dat[i, ])
      }
      perm.pca = prcomp(t(as.matrix(perm.dat)),tol=0.5)
      summary(perm.pca)$importance[2,1]
    })
    vg = vg[select.genes,]
    vg = vg[order(vg$loess.padj),]
    list(pca.ratio=pca.var/mean(perm.pca.var), pca.var = pca.var, vg.num=length(select.genes), vg=vg)
  }


WGCNA_hetero <- function(norm.dat,select.cells, perm.num=10,vg.padj=0.5,rm.gene.mod=NULL)
  {
    vg= findVG(2^norm.dat[,select.cells]-1)
    select.genes = row.names(vg)[vg$loess.padj < vg.padj]
    dat = norm.dat[select.genes,select.cells]
    adj=adjacency(t(dat), power = 4,type="unsigned")
    adj[is.na(adj)]=0
    TOM = TOMsimilarity(adj,TOMType="unsigned")
    dissTOM = as.matrix(1-TOM)
    row.names(dissTOM)= colnames(dissTOM) = row.names(dat)
    geneTree = flashClust(as.dist(dissTOM), method =  "average")
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, cutHeight=0.99, deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = 5)
    gene.mod = split(row.names(dissTOM), dynamicMods)
    gene.mod = gene.mod[setdiff(names(gene.mod),"0")]
    rm.eigen=NULL
    if(!is.null(rm.gene.mod)){
      rm.eigen = getEigen(rm.gene.mod, norm.dat, select.cells)
    }
    gm= filterGeneMod(norm.dat, select.cells, gene.mod, minModuleSize=5, rm.eigen=rm.eigen,min.padj = 40, padj.th=0.05,min.cells=3, ...)
    
    if(length(select.genes)==0){
      return(0)
    }
    if(length(select.genes)==1){
      return(1)
    }
    rd.WGCNA(norm.dat, select.genes, select.cells, minModuleSize=5,...)

    pca.var = summary(pca)$importance[2,1]
 
    ###shuffle expression of every gene.
    perm.pca.var = sapply(1:10, function(x){
      perm.dat = norm.dat[select.genes, select.cells]
      for(i in 1:nrow(perm.dat)){
        perm.dat[i, ] = sample(perm.dat[i, ])
      }
      perm.pca = prcomp(t(as.matrix(perm.dat)),tol=0.5)
      summary(perm.pca)$importance[2,1]
    })
    vg = vg[select.genes,]
    vg = vg[order(vg$loess.padj),]
    list(pca.ratio=pca.var/mean(perm.pca.var), pca.var = pca.var, vg.num=length(select.genes), vg=vg)
  }


qc.index = setNames(samp.dat$percent_reads_aligned_to_genome_only/samp.dat$percent_reads_aligned_to_mrna,row.names(samp.dat))

###Add oligodendrocyte contamination genes
save(rm.gene.mod, file="rm.gene.mod.hetero.rda")
pca.hetero= sapply(levels(clean.cl),function(x){
  print(x)
  select.cells = names(clean.cl)[clean.cl==x]
  #tmp.norm = lm_normalize(norm.dat[,select.cells], qc.index[select.cells])[[1]]
  tmp=pca_hetero(norm.dat,select.cells, perm.num=10,vg.padj=0.5, rm.eigen=rm.eigen)
},simplify=F)
pca.ratio = unlist(sapply(pca.hetero, function(x)x[1]))

names(pca.ratio) = paste(names(pca.hetero), cl.df[names(pca.hetero),"cluster_label"])




library(iterclust)

subsample.de.genes=list()
load("norm.dat.rda")
load("cl.final.rda")
load("select.genes.rda")
load("cl.med.rda")
load("cl.clean.rda")

dir.create("subsample_merge")
load("de.param.rda")
n.bin = 10
bins = sample(1:n.bin, length(cl.clean),replace=TRUE)
for(i in 1:n.bin){
  ####subsample cells and recompute de.genes 
  select.cells = names(cl.clean)[bins<= i]
  print(length(select.cells))
  rd.dat = t(as.matrix(norm.dat[select.genes, select.cells]))
  merge.result= merge_cl(norm.dat, cl.clean[select.cells], rd.dat=rd.dat, de.param=de.param)
  save(merge.result, file=file.path("subsample_merge/",paste0(i,".result.rda")))
}



merge.files <- dir("subsample_merge/",pattern="result.rda")
merge.stats <- do.call("rbind",sapply(merge.files, function(f){
  load(file.path("subsample_merge/",f))
  tmp.cl = merge.result$cl
  cl.size = table(cl[names(tmp.cl)])
  tb=table(tmp.cl, cl.clean[names(tmp.cl)])
  tb.df = as.data.frame(tb)
  tb.df = tb.df[tb.df$Freq>0,]
  tmp=split(tb.df[,2],tb.df[,1])
  merge.size=sapply(tmp, length)
  merge.cl =names(merge.size)[merge.size>1]
  if(length(merge.cl)>0){
    tb.df = tb.df[tb.df[,1]%in% merge.cl,]
    tb.df =tb.df[order(tb.df[,1]),]
    #tb.df[,"mergeCl"] = cl.df[as.character(tb.df[,1]),"cluster_label"]
    #tb.df[,"orgCl"] = cl.df[as.character(tb.df[,2]),"cluster_label"]
    colnames(tb.df) = c("mergeCl", "orgCl")
    absent.cl = colnames(tb)[colSums(tb)==0]
    merge.df = data.frame(org.cl = as.character(tb.df$orgCl), merge.cl = as.character(tb.df$mergeCl),stringsAsFactors=FALSE)
    merge.df = rbind(merge.df, data.frame(org.cl = as.character(absent.cl), merge.cl=rep("absent", length(absent.cl)),stringsAsFactors=FALSE))

    nochange.cl = setdiff(colnames(tb), c(absent.cl, as.character(tb.df[,2])))
    merge.df = rbind(merge.df, data.frame(org.cl = as.character(nochange.cl,"cluster_label"), merge.cl=rep("nochange", length(nochange.cl)),stringsAsFactors=FALSE))
    merge.df$sample_size = length(tmp.cl)
    merge.df$cl_size = cl.size[as.character(merge.df$org.cl)]
    return(merge.df)
  }
  return(NULL)
},simplify=F))


merge.stats$org.cl = factor(as.character(merge.stats$org.cl), labels(dend))
merge.stats$merge.cl = droplevels(factor(merge.stats$merge.cl, levels=c("absent","nochange",labels(dend))))

merge.stats$org.cl_label = merge.stats$org.cl
levels(merge.stats$org.cl_label) = cl.df[levels(merge.stats$org.cl_label),"cluster_label"]
merge.stats$merge.cl_label = merge.stats$merge.cl
levels(merge.stats$merge.cl_label) = c("absent", "nochange", as.character(cl.df[levels(merge.stats$merge.cl_label)[-(1:2)],"cluster_label"]))

merge.stats$sample_size = factor(merge.stats$sample_size)
saturation.df = merge.stats
save(saturation.df, file="saturation.df.rda")

cl.color = setNames(cl.df$cluster_color, cl.df$cluster_label)
cl.color = c("black","white",as.character(cl.color[levels(merge.stats$merge.cl)[-(1:2)]]))

g=ggplot(saturation.df, aes(org.cl, sample_size, fill = merge.cl)) + geom_raster() + scale_fill_manual(values=cl.color)
g = g + theme(axis.text.x = element_text(angle=90, hjust=1))
pdf("saturation.pdf",height=6,width=14)
g
dev.off()



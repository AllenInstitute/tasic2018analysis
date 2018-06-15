library(limma)
library(iterclust)
library(Matrix)
library(matrixStats)
library(RColorBrewer)
library(dplyr)

jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
blue.red <-colorRampPalette(c("blue", "white", "red"))
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/R_Object/20180515_RSC-11-153_mouse_star_samp.dat.Rdata")

exclude_genotypes <- c("Ai90","Ai94",
                       "Ai139","Ai140",
                       "Fezf2-CreER/wt;Ai148",
                       "Sst-IRES-Cre/wt;Ai14\\(RCL-tdT\\)/wt;Ai148",
                       "Rbp4-Cre_KL100/wt;Ai14\\(RCL-tdT\\)/wt;Ai148",
                       "Ai162","Ai163","Ai166",
                       "Ai173","Ai174","Ai175",
                       "Snap25-T2A-GCaMP6s")

genotype_filter <- paste0(exclude_genotypes, collapse = "|")
visp_ep_donors <- c("362849","362845")
mouse_genetics_donors <- c("366232","369426")

samp.dat_keep <- samp.dat %>%
    filter(Region %in% c("VISp","ALM")) %>%
    filter(Type == "Cells") %>%
    filter(!grepl(genotype_filter, full_genotype)) %>%
    filter(Injection_type != "anterograde" | is.na(Injection_type)) %>%
    filter(!external_donor_name %in% visp_ep_donors) %>%
    filter(!external_donor_name %in% mouse_genetics_donors) %>%
    # These are reagent validation runs, and should be excluded.
    filter(!(full_genotype == "Rbp4-Cre_KL100/wt;Ai14(RCL-tdT)/wt" & batch %in% c("R8S4-180314","R8S4-180328"))) %>%
    # This is the bad 2018 lot, and should be excluded.
    filter(S4.Lot != "1709695A")

nrow(samp.dat_keep)


save(samp.dat_keep, file="samp.dat_keep.rda")
samp.dat = samp.dat_keep
row.names(samp.dat) = samp.dat[,1]
samp.dat = samp.dat[,-1]
all.cells= row.names(samp.dat)
save(all.cells, file="all.cells.rda")

d = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/facs/V1_ALM_R_Objects/STAR2.0"
load(file.path(d, "V1_ALM_Star2.0_29002_samp.dat.Rdata"))
row.names(samp.dat) = samp.dat[,1]
samp.dat = samp.dat[all.cells,]
save(samp.dat, file="samp.dat.rda")

load(file.path(d,"V1_ALM_Star2.0_29002_exon.Rdata"))
genechr<-read.csv("~/zizhen/mm10/gene_chr.csv",as.is=T,header=F)
genes = row.names(exon)
select.genes = setdiff(genes, c(grep("^Gm",genes,value=T),genechr[genechr[,1]=="NC_005089.1",2]))
select.genes = setdiff(select.genes, c(grep("^Rps",select.genes,value=T),grep("^Rpl",select.genes,value=T)))
load("../process_21729/rm.gene.mod.rda")
sex.genes= rm.gene.mod[["6"]]
select.genes = setdiff(select.genes, sex.genes)
save(select.genes, file="select.genes.rda")



exclude.cells=with(samp.dat, row.names(samp.dat)[Genes.With.CPM < 1000 | percent_reads_aligned_total < 75 | complexity_cg > 0.5 | total_reads < 100000])
qc.cells= setdiff(all.cells, exclude.cells)

save(qc.cells, file="qc.cells.rda")




cpm = t(t(exon)*10^6/colSums(exon))
save(cpm, file="cpm.rda")

exon = exon[select.genes, qc.cells]
norm.dat = t(t(exon)*10^6/colSums(exon))
norm.dat = Matrix(norm.dat, sparse=TRUE)
norm.dat@x = log2(norm.dat@x +1)
save(norm.dat, file="norm.dat.rda")


###Remove duplets
load("top.gene.mod.rda")
select.cells= qc.cells
top.eigen=get_eigen(top.gene.mod, norm.dat[,select.cells], select.cells)[[1]]

eigen.min = colMins(top.eigen)
eigen.max = colMaxs(top.eigen)
eigen.norm = t((t(top.eigen)-eigen.min)/(eigen.max - eigen.min))
top.cl = apply(eigen.norm, 1, which.max)


eigen.mean = do.call("rbind",tapply(names(top.cl), top.cl, function(x) colMeans(eigen.norm[x,,drop=F])))
eigen.sd =  do.call("rbind",tapply(names(top.cl), top.cl, function(x) colSds(eigen.norm[x,,drop=F])))
th = sapply(1:ncol(eigen.mean[-1,]), function(i){
  match.cl=row.names(eigen.mean)[which.max(eigen.mean[,i])]
  bg.cells= names(top.cl)[!top.cl %in% match.cl]
  bg.mean = mean(eigen.norm[bg.cells,i])
  bg.sd = sd(eigen.norm[bg.cells,i])
  return(c(bg.mean, bg.sd))
})
th = th[1,] + 3 * th[2,]

eigen.offset = t(t(eigen.norm) - as.vector(th))
###Merge non-neuronal types
eigen.offset = cbind(eigen.offset[,1:2], non.neuronal=rowMaxs(eigen.offset[,-(1:2)]))

rm.cells = row.names(eigen.offset)[rowSums(eigen.offset>0)>1]
save(rm.cells, file="rm.cells.rda")
save(eigen.offset, file="eigen.offset.rda")

select.cells = setdiff(qc.cells, rm.cells)
save(select.cells, file="select.cells.rda")
top.cl = top.cl[select.cells]
###Merge non-neuronal cells
top.cl[top.cl %in% c(3:7)]=3
save(top.cl, file="top.cl.rda")
####



###Compute color bars for heatmap
cols = c("Region","roi","batch_vendor_name","Injection_type","injection_roi")
meta.df = samp.dat[,cols]
all.col.map=sapply(c(cols), function(x){
  tmp = droplevels(as.factor(samp.dat[[x]]))
  setNames(jet.colors(length(levels(tmp))), levels(tmp))
},simplify=FALSE)


d="/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20170913"
anno = read_feather(file.path(d, "anno.feather"))
tmp=as.data.frame(unique(anno[,c("layer_label","layer_color")]))
all.col.map$layer = setNames(tmp$layer_color, tmp$layer_label)

layer = gsub("^.*_","", samp.dat$roi)
layer[layer %in% c("L1-L2-3-L4-L5-L6","L1-L2-3-L5-L6","L2-3-L5")] = "L1-L6"
layer[layer %in% c("L2-3-L4")]= "L2/3-L4"
layer[layer %in% c("L4-L5-L6")]= "L4-L6"
layer[layer %in% c("L1-L2-L3","L1-L2-3")]= "L1-L2/3"
layer[layer %in% c("L1-L2-3-L4")]= "L1-L4"
layer=gsub("-3","/3",layer)
layer[layer %in% c("L4-6")]= "L4-L6"
layer[layer %in% c("L5-6")]= "L5-L6"
layer[layer %in% c("L2")]="L2/3"
layer[!layer %in% tmp[,1]]=""
meta.df$layer = layer


continuous.cols = c("percent_reads_aligned_to_introns", "percent_reads_aligned_to_exons", "Genes.With.FPKM")
meta.df.continuous = data.frame(sapply(continuous.cols, function(x){
  val = samp.dat[[x]]
  r = quantile(val, c(0.05, 0.95))
  breaks = c(min(val), seq(r[1],r[2],length.out=99), max(val))
  bin = cut(val, breaks)
},simplify=F))

continuous.col.map = sapply(continuous.cols, function(x){
  print(x)
  bin = meta.df.continuous[[x]]
  setNames(blue.red(length(levels(bin))),levels(bin))
},simplify=F)

meta.df = cbind(meta.df, meta.df.continuous)
all.col.map = c(all.col.map, continuous.col.map)
all.col = sapply(names(all.col.map), function(x){
  all.col.map[[x]][as.character(meta.df[[x]])]
})
all.col = t(all.col)
colnames(all.col)=row.names(samp.dat)
all.col[is.na(all.col)]="black"
save(meta.df, all.col, all.col.map, file="all.col.rda")


####compute rm.eigen for masking
log2Genes=log2(samp.dat$Genes.Detected.CPM)
log2Genes  = setNames(log2Genes - mean(log2Genes), row.names(samp.dat))


load("batch.gene.rda")
batch.effect = get_eigen(batch.gene, norm.dat, colnames(norm.dat))[[1]]
rm.eigen = cbind(log2Gene=log2Genes[colnames(norm.dat)], batch.effect)
save(rm.eigen, file="rm.eigen.rda")




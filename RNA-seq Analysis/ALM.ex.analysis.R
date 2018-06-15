ALM.ex.cl = droplevels(ex.cl[samp.dat[names(ex.cl),"Region"]=="ALM"])
ALM.ex.cl = ALM.ex.cl[cl.df[as.character(ALM.ex.cl),"region_label"]!="VISp"]
ALM.ex.markers = select_markers(norm.dat, ALM.ex.cl, n.markers=50, de.genes=de.genes)$markers
ALM.ex.tsne.result <- plot_tSNE_cl(norm.dat, ALM.ex.markers, ALM.ex.cl, cl.df, theta=0.05)
ALM.ex.tsne.df = ALM.ex.tsne.result$tsne.df
save(ALM.ex.tsne.df, file="ALM.ex.tsne.df.rda")
ggsave("ALM.ex.tsne.pdf", ALM.ex.tsne.result$g, height=10, width=10)
save(ALM.ex.cl, file="ALM.ex.cl.rda")
save(ALM.ex.markers, file="ALM.ex.markers.rda")
ALM.ex.cell.cl.map.df =  get_core_transition(norm.dat, ALM.ex.cl, ALM.ex.markers, n.bin=5, n.iter=100, mc.cores=10)
save(ALM.ex.cell.cl.map.df, file="ALM.ex.cell.cl.map.df.rda")

de.param$q1.th=0.4
ALM.ex.PT.cl = droplevels(ALM.ex.cl[ALM.ex.cl %in% row.names(cl.df)[cl.df$subclass_label=="L5 PT"]])




ALM.ex.cl.present.counts = get_cl_sums(norm.dat>0, droplevels(ALM.ex.cl))
cl.present.counts = get_cl_sums(norm.dat>0, cl.clean)

df = within_group_specific_markers(levels(ALM.ex.PT.cl), norm.dat, ALM.ex.PT.cl, de.param = de.param, cl.present.counts=ALM.ex.cl.present.counts)


PT = levels(ALM.ex.PT.cl)
all = levels(ALM.ex.cl)
pairs <- list(c1=list(95, PT), c2=list(96:97, PT), c3=list(PT, all), c4=list(95, all), c5=list(96:97,all), c6=list(96,96:97), c7=list(96, all), c8=list(97, 96:97),c9=list(97, all))
 
de.param$q1.th=0.4
ALM.PT.gene.list= sapply(pairs, function(x){
  print(x)
  df=group_specific_markers(as.character(x[[1]]), norm.dat, droplevels(cl.clean[cl.clean%in% as.character(x[[2]])]), de.param=de.param, n.markers=10, cl.present.counts = cl.present.counts)
  df=df[order(df$pval),]
  head(df,20)
},simplify=F)
save(ALM.PT.gene.list, file="ALM.PT.gene.list.top.20.rda")
PT.pairs = read.csv("PT_DE_pairs",row.names=1,header=FALSE)

names(ALM.PT.gene.list) = PT.pairs[[1]]

ALM.PT.gene.df = do.call("rbind", ALM.PT.gene.list)
ALM.PT.gene.df$pair = gsub("\\..*","",row.names(ALM.PT.gene.df))
ALM.PT.gene.df = ALM.PT.gene.df[,c(8,1:7)]
row.names(ALM.PT.gene.df) = NULL
write.csv(ALM.PT.gene.df, file="ALM.PT.gene.df.csv")

de.param = de_param(q1.th=0.4)
tmp=display_cl(ALM.ex.PT.cl,norm.dat, prefix="ALM.L5.PT", de.param = de.param, n.markers=30)
ALM.ex.PT.markers= tmp$markers
top.genes= unique(unlist(lapply(ALM.PT.gene.list, row.names)))

setdiff(unique(unlist(lapply(ALM.PT.gene.list, row.names))), ALM.ex.PT.markers)
tmp=display_cl(ALM.ex.PT.cl,norm.dat, prefix="ALM.L5.PT.tmp", de.param = de.param, markers=top.genes)
save(ALM.ex.PT.markers, file="ALM.ex.PT.markers.rda")


source("~/zizhen/My_R/map_river_plot.R")
source("~/zizhen/My_R/sankey_functions.R")

colnames(map.df) = c("map_cluster_id","map_prob","map_cluster_label", "cluster_id","sub_cluster","stim","cor")
map.df$map_cluster_color = cl.df[as.character(map.df$map_cluster_id),"cluster_color"]

tmp = compare_annotate(setNames(map.df$cluster_id, row.names(map.df)), setNames(map.df$map_cluster_id,row.names(map.df)), cl.df, reorder=FALSE)
tmp.cl.df = tmp$cl.df
map.df$cluster_label = map.df$cluster_id
map.df$cluster_id = as.integer(as.factor(map.df$cluster_label))
map.df$cluster_color = tmp.cl.df[as.character(map.df$cluster_label),"cluster_color"]
map.df$cluster_color=as.character(map.df$cluster_color)
map.df$map_cluster_color=as.character(map.df$map_cluster_color)

map.df$maintype = hrvatin.samp.dat$maintype
inh.map.df = droplevels(map.df[which(map.df$maintype=="Interneurons"),])
g=river_plot(inh.map.df, min.cells=4, min.frac=0.1)
ggsave(g, file="Hrvatin.inh.map.pdf")

ex.map.df = droplevels(map.df[which(map.df$maintype=="Excitatory" & !map.df$cluster_label %in% c("RSP","Hip","Sub")),])
g=river_plot(ex.map.df, min.cells=10, min.frac=0.1)
ggsave(g, file="Hrvatin.ex.map.pdf")

tb=with(ex.map.df, table(map_cluster_label, stim))
tb.df = as.data.frame(tb)
table(map.df$map_cluster_label)
stim.pval = sapply(names(de.num), function(p){
  x = de.num[[p]]
  n <- total.num - m
  k <- full.de.num[[p]]
  pval = phyper(x - 1, m, n, k, lower.tail = FALSE)
})









                








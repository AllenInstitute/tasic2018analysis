VISp.ex.cl = droplevels(ex.cl[samp.dat[names(ex.cl),"Region"]=="VISp"])
VISp.ex.cl = VISp.ex.cl[cl.df[as.character(VISp.ex.cl),"region_label"]!="ALM"]
VISp.ex.markers = select_markers(norm.dat, VISp.ex.cl, n.markers=50, de.genes=de.genes)$markers
VISp.ex.tsne.result <- plot_tSNE_cl(norm.dat, VISp.ex.markers, VISp.ex.cl, cl.df, theta=0.05)
VISp.ex.tsne.df = VISp.ex.tsne.result$tsne.df
save(VISp.ex.tsne.df, file="VISp.ex.tsne.df.rda")
ggsave("VISp.ex.tsne.pdf", VISp.ex.tsne.result$g, height=10, width=10)
save(VISp.ex.cl, file="VISp.ex.cl.rda")
save(VISp.ex.markers, file="VISp.ex.markers.rda")
VISp.ex.cell.cl.map.df =  get_core_transition(norm.dat, VISp.ex.cl, VISp.ex.markers, n.bin=5, n.iter=100, mc.cores=10)
save(VISp.ex.cell.cl.map.df, file="VISp.ex.cell.cl.map.df.rda")





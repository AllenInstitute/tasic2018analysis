load(file.path(d, "20180228_RSC-11-142_mouse_star_exon.Rdata"))
load(file.path(d, "20180228_RSC-11-142_mouse_star_intron.Rdata"))
load(file.path(d, "20180228_RSC-11-142_mouse_star_samp.dat.Rdata"))
row.names(samp.dat) = samp.dat[[1]]

cell.cl = cl
cell.cl.df = cl.df


load("~/zizhen/projects/Mouse/Nuclei/M1_V1/cl.final.rda")
nuclei.cl = cl
nuclei.cl.df = cl.df

L45.cells = names(cell.cl)[cell.cl %in% as.character(c(58:61,130,137))]
L45.nuclei = names(nuclei.cl)[nuclei.cl %in% as.character(c(37:39,42))]

gene.cell.intron.ratio <- rowSums(intron[,L45.cells])/rowSums(exon[,L45.cells]+intron[,L45.cells])

gene.nuclei.intron.ratio <- rowSums(intron[,L45.nuclei])/rowSums(exon[,L45.nuclei]+intron[,L45.nuclei])

nuclei.ratio = gene.nuclei.intron.ratio/gene.cell.intron.ratio
nuclei.ratio=nuclei.ratio[!is.na(nuclei.ratio)]


exon.cpm = t(t(exon)*10^6/Matrix::colSums(exon))
cell.exon = rowMeans(exon.cpm[,L45.cells])
nuclei.exon = rowMeans(exon.cpm[,L45.nuclei])
cell.nuclei.exon = data.frame(cell.exon, nuclei.exon)
cell.nuclei.exon$cell.nuclei.ratio = cell.nuclei.exon[,1]/cell.nuclei.exon[,2]
  
exon.ratio = rowMeans(exon[,L45.cells])/

cell.norm.dat = norm.dat

tmp=load("~/zizhen/projects/Mouse/Nuclei/M1_V1/norm.dat.rda")
tmp=load("~/zizhen/projects/Mouse/Nuclei/M1_V1/all.col.rda")
all.col[is.na(all.col)]="black"
nuclei.norm.dat= norm.dat
nuclei.all.col = all.col[,L45.nuclei]
tmp=display_cl(droplevels(nuclei.cl[L45.nuclei]), nuclei.norm.dat[,L45.nuclei],prefix="L45.nuclei.new",col=nuclei.all.col)

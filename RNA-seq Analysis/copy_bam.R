select.cl = c(121, 90, 85, 57, 42, 31, 18, 6)
tmp.cl = droplevels(cl[cl %in% as.character(select.cl)])

d = "select_bam"
dir.create(d)





select.cells = sample_cells(tmp.cl, 10)
bam_dir1 = setNames(file.path("/allen/programs/celltypes/workgroups/rnaseqanalysis/STARforLIMS/Mouse/star_out/", select.cells, paste0(select.cells,"_Aligned.sortedByCoord.out.bam")), select.cells)
tmp1=file.exists(bam_dir1)

bam_dir2 = setNames(with(samp.dat[select.cells,], file.path(fpkm_dir,paste0("ar_", ar_id, "_STAR_Aligned.sortedByCoord.out.bam"))), select.cells)
tmp2= file.exists(bam_dir2)

bam_dir = bam_dir1
bam_dir[!tmp1] = bam_dir2[!tmp1]
table(file.exists(bam_dir))


tmp=file.copy(bam_dir, file.path(d, paste0(select.cells,".bam")))
write.csv(cl.df, file="cl.df.csv")
write.csv(cl, file="cl.csv")


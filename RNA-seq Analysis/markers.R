Inh.markers <- c("Gad2","Slc32a1","Prox1","Adarb2","Nfix","Nfib","Cacna2d1",
                                     "Cxcl14","Tnfaip8l3","Cplx3","Lamp5","Cd34","Pax6","Krt73",
                                     "Scrg1","Egln3","Ndnf","Tmem182","Ntn1","Pde11a","Pdlim5",
                                     "Lsp1","Slc35d3","Nkx2-1","Serpinf1","Col14a1","Vip","Sncg",
                                     "Crabp1","Slc10a4","Cldn10","Bhlhe22","Crispld2","Slc17a8",
                                     "Cyb5r2","Nr1h4","Wnt7b","Prss12","Igfbp6","Calb2","Grpr",
                                     "Pthlh","Elfn1","Rspo1","Slc18a3","Lmo1","Rspo4","Sostdc1",
                                     "Chat","Cbln4","Gsx2","Gpc3","Mab21l1","C1ql1","Itih5","Mybpc1",
                                     "Myl1","Lhx6","Sox6","Sst","Chodl","Calb1","Cbln4","Etv1","Edn1",
                                     "Kl","Il1rapl2","Myh8","Ptprk","Chrna2","Myh13","Ptgdr","Crhr2",
                                     "Hpse","Igsf9","C1ql3","Tacstd2","Th","Col6a1","Nts","Tac1","Pvalb",
                                     "Gabrg1","Il7","Bche","Prdm8","Syt2","Ostn","Pdlim3","C1ql1",
                                     "Gpr149","Vipr2","Meis2","Adamts19","Cpa6","Lgr6")
Inh.comb.markers <- c("Reln","Cnr1","Nr2f2","Cck","Npy","Crh","Tac2")
                                     
                 

Ex.markers <- unique(c("Slc17a7","Rtn4rl2","Slc30a3","Cux2","Stard8","Otof","Rrad","Penk","Agmat",
                                "Emx2","S100a3","Macc1","Rorb","Scnn1a","Whrn","Endou","Col26a1",
                                "Rspo1","Fezf2","Hsd11b1","Batf3","Arhgap25","Colq","Pld5","Olfr78",
                                "Tcap","Fgf17","Wfdc18","Wfdc17","Aldh1a7","Tgfb1","Ctsc","Rxfp2",
                                "Prss35","Rgs12","Osr1","Oprk1","Cd52","Col23a1","Col18a1","Car1",
                                "Car3","Fam84b","Chrna6","Chrnb3","Fn1","Tac1","Lce3c","Erg",
                                "Cdc42ep5","Bmp5","Pvalb","Depdc7","Stac","C1ql2","Ptgfr","Slco2a1",
                                "Pappa2","Dppa1","Npsr1","Htr2c","Hpgd","Nxph3","Sla2","Tshz2",
                                "Rapgef3","Slc17a8","Trh","Nxph2","Foxp2","Col12a1","Syt6","Col5a1",
                                "Gpr139","Ly6d","Sla","Cpa6","Ppp1r18","Faim3","Ctxn3","Nxph4",
                                "Cplx3","Ctgf","Col8a1","Mup5","Ngf","Fam150a","F2r","Serpinb11","Fbxl7","P2ry12",
                                "Crh","Kynu","Hsd17b2","Mup3","Tlcd1","Lhx5","Trp73","Cpa6","Gkn1","Col18a1","Lce3c","Erg","Bmp5","Stac","C1ql2","Slco2a1","Lrrc9","Trhr","Myzap","Krt80","H60b","Fam150a","Clic5","Kcnj5","Olfr110","Olfr111"))
Ex.comb.markers <- c("Reln","Cdh13","Cpne7","Alcam","Rprm","Marcksl1")
                                


Global.markers <- c("Fez1","Phyhipl","Aplp1","Gnao1","Caly","Snap25","Atp1a3","Camk2b",
                                        "Syt1","Gabrg2","Fabp3","Stmn2","Kif5c","Slc32a1","Gad2","Dlx1","Dlx5","Dlx2","Dlx6os1",
                                        "Slc6a1","Sox2","Slc17a7","Nrn1","Neurod2","Sv2b","Satb2","Tbr1","Vsig2","Cmtm5","Kcnj10",
                                        "S100a16","S100a13","S1pr1","Gja1","Gjb6","Aqp4","Lcat","Acsbg1","Olig1","Sox10","Neu4",
                                        "Sapcd2","Gpr17","Plp1","Cldn11","Mag","Mog","Nkx6-2","Enpp6","9630013A20Rik","Brca1",
                                        "Mog","Opalin","Gjb1","Hapln2","Cyba","Ctsh","Ifitm3","Sparc","S100a11","Dcn","Col1a1",
                                        "Pltp","Vtn","Slc6a13","Spp1","Slc13a3","Col15a1","Slc47a1","Tgtp2","Ifi47","Esam",
                                        "Slco1a4","Slc38a5","Cldn5","H2-Q7","Slc38a11","Art3","Ace2","Acta2","Myh11","Pln",
                                        "Gja5","Kcnj8","Atp13a5","Aoc3","Ctss","C1qb","C1qc","C1qa","Cbr2","F13a1","Pf4",
                                        "Mrc1","Siglech","Selplg")



cl.present = get_cl_means(norm.dat > 1, cl.clean)
Ex.markers=Ex.markers[order(apply(cl.present[Ex.markers, row.names(cl.df)[c(53:109,111)]]>0.4 , 1, function(x)which(x)[1]))]
write.csv(Ex.markers, file="markers_plot2.csv")

Inh.markers=Inh.markers[order(apply(cl.present[Inh.markers, row.names(cl.df)[c(1:52,110)]]>0.4 , 1, function(x)which(x)[1]))]
write.csv(Inh.markers, file="markers_plot3.csv")


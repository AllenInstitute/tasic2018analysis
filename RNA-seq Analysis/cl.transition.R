library(iterclust)
library(matrixStats)

cl.clean = droplevels(cl[cl %in% row.names(cl.df)[cl.df$class_label!="Low Quality"]])  
 
test_cv_cor <- function(norm.dat, cl, markers, n.bin=5,g.perc=1){
  bins=unlist(tapply(names(cl), cl, function(x){
    if(length(x) > n.bin){
      tmp=rep_len(1:n.bin, length(x))
    }else{
      tmp = sample(1:n.bin, length(x))
    }
    setNames(tmp[sample(length(tmp))], x)
  }))
  names(bins) = gsub(".*\\.", "", names(bins))
  bins= bins[names(cl)]
  pred.cl = setNames(rep(NA, length(cl)), names(cl))
  for(i in 1:n.bin){
    print(i)
    train.cells = names(cl)[bins!=i]
    test.cells =names(cl)[bins==i]
    select.markers=sample(markers, round(length(markers)*g.perc))
    map.result <- map_by_cor(norm.dat[select.markers,], cl[train.cells], norm.dat[select.markers, test.cells])$pred.df
    pred.cl[test.cells] = as.character(map.result[test.cells, "pred.cl"])
  }
  return(pred.cl)
}


get_core_transition <- function(norm.dat, cl, select.markers, n.bin=5, n.iter=100, mc.cores=10)
  {
    cl.cv <- mclapply(1:n.iter, function(i){
      tmp=test_cv_cor(norm.dat, cl, select.markers, n.bin=n.bin)
    }, mc.cores=mc.cores)
    
    cell.cl.cor.map = do.call("rbind",sapply(cl.cv, function(x){
      df = data.frame(cell=names(x),cl=x)
    },simplify=F))
    cell.cl.cor.map = table(cell.cl.cor.map[,1],cell.cl.cor.map[,2])
    cell.cl.cor.map = cell.cl.cor.map / rowSums(cell.cl.cor.map)

    cell.cl.map.df = data.frame(org.cl= as.character(cl[row.names(cell.cl.cor.map)]),best.score=rowMaxs(cell.cl.cor.map), best.cl = colnames(cell.cl.cor.map)[apply(cell.cl.cor.map, 1, which.max)], stringsAsFactors=FALSE)
    row.names(cell.cl.map.df) = row.names(cell.cl.cor.map)
    tmp=get_pair_matrix_coor(cell.cl.cor.map, row.names(cell.cl.map.df), as.character(cell.cl.map.df$best.cl))
    tmp1 = cell.cl.cor.map
    tmp1[tmp]= 0
    cell.cl.map.df$second.score = rowMaxs(tmp1)
    cell.cl.map.df$second.cl =colnames(tmp1)[apply(tmp1,1, which.max)]
    cell.cl.map.df$second.cl[cell.cl.map.df$second.score ==0] = NA
    
    cell.cl.map.df$transition.cl = NA
    tmp = with(cell.cl.map.df, org.cl!=best.cl | best.score < 0.9)
    cell.cl.map.df[tmp,"transition.cl"] = as.character(cell.cl.map.df[tmp,"best.cl"])
    tmp = with(cell.cl.map.df, which(org.cl==transition.cl))
    cell.cl.map.df$transition.cl[tmp] = as.character(cell.cl.map.df[tmp,"second.cl"])
    
    cl.med <- do.call("cbind",tapply(names(cl), cl, function(x){
      rowMedians(as.matrix(norm.dat[select.markers,x]))
    }))
    row.names(cl.med) = select.markers
    
    cell.cl.cor=cor(as.matrix(norm.dat[select.markers, row.names(cell.cl.map.df)]), cl.med[select.markers,])
    cell.cl.map.df$cor = with(cell.cl.map.df, get_pair_matrix(cell.cl.cor, row.names(cell.cl.map.df),as.character(org.cl)))
    cell.cl.map.df$core = is.na(cell.cl.map.df$transition.cl)
    return(cell.cl.map.df)
  }

load("select.markers.rda")
cell.cl.map.df = get_core_transition(norm.dat, cl.clean, select.markers, n.bin=5, n.iter=100, mc.cores=10)

corrected.cells = with(cell.cl.map.df, org.cl!=best.cl & best.score > 0.9)
cl.clean.correct = cl.clean
cl.clean.correct[corrected.cells] = cell.cl.map.df[corrected.cells, "best.cl"]

tmp.cells = intersect(names(cl.clean)[cl.clean=="104"], row.names(samp.dat)[samp.dat$Region=="VISp"])

tmp1.cells = intersect(names(cl.clean)[cl.clean=="104"], row.names(samp.dat)[samp.dat$Region=="ALM"])



save(cell.cl.map.df, file="cell.cl.map.df.rda")
 
transition.df = with(cell.cl.map.df, as.data.frame(table(org.cl, transition.cl)))
transition.df = transition.df[transition.df$Freq > 0,]
transition.df$org.cl = as.character(transition.df$org.cl)
transition.df$transition.cl = as.character(transition.df$transition.cl)
save(transition.df, file="transition.df.rda")

###combine transitions from both directions
transition.df$cl.min = pmin(transition.df$org.cl, transition.df$transition.cl)
transition.df$cl.max = pmax(transition.df$org.cl, transition.df$transition.cl)
transition.df$cl.pair = paste(transition.df$cl.min, transition.df$cl.max)
transition.df.comb= do.call("rbind",tapply(1:nrow(transition.df),transition.df$cl.pair, function(x){
  tmp = transition.df[x,][1,]
  tmp$Freq = sum(transition.df[x,"Freq"])
  tmp[,c(4,5,3)]
}))
cl.size = table(cl.clean)
transition.df.comb$cl.min.size = cl.size[transition.df.comb$cl.min]
transition.df.comb$cl.max.size = cl.size[transition.df.comb$cl.max]
transition.df.comb$ratio = with(transition.df.comb,Freq/pmin(cl.min.size,cl.max.size))
save(transition.df.comb, file="transition.df.comb.rda")
transition.df.comb$cl1_label = cl.df[as.character(transition.df.comb$cl.min),"cluster_label"]
transition.df.comb$cl2_label = cl.df[as.character(transition.df.comb$cl.max),"cluster_label"]
transition.df.comb[transition.df.comb$ratio > 0.1 & transition.df.comb$Freq > 1,]
colnames(transition.df.comb)[1:2] = c("cl1","cl2")
save(transition.df.comb, file="transition.df.comb.rda")







fast_tsne <- function(dat, ...)
  {
    result <- fftRtsne(dat,fast_tsne_path = "~/src/FIt-SNE/bin/fast_tsne", ...)
    df = as.data.frame(result)
    row.names(df) = row.names(dat)
    colnames(df)=c("Lim1","Lim2")
    df
  }


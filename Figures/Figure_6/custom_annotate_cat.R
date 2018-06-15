annotate_cat <- function(df,
                         col = NULL, base = NULL,
                         sort_label = T, na_val = "ZZ_Missing",
                         colorset = "varibow", color_order = "sort") {
  
  library(dplyr)
  library(viridis)
  
  if(class(try(is.character(col), silent = T)) == "try-error") {
    col <- lazyeval::expr_text(col)
  } else if(class(col) == "NULL") {
    stop("Specify a column (col) to annotate.")
  }
  
  if(class(try(is.character(base), silent = T)) == "try-error") {
    base <- lazyeval::expr_text(base)
  } else if(class(base) == "NULL") {
    base <- col
  }
  
  if(!is.character(df[[col]])) {
    df[[col]] <- as.character(df[[col]])
  }
  
  df[[col]][is.na(df[[col]])] <- na_val
  
  x <- df[[col]]
  
  annotations <- data.frame(label = unique(x), stringsAsFactors = F)
  
  if(sort_label == T) {
    annotations <- annotations %>% arrange(label)
  }
  
  annotations <- annotations %>%
    mutate(id = 1:n())
  
  if(colorset == "varibow") {
    colors <- varibow(nrow(annotations))
  } else if(colorset == "rainbow") {
    colors <- sub("FF$","",rainbow(nrow(annotations)))
  } else if(colorset == "viridis") {
    colors <- sub("FF$","",viridis(nrow(annotations)))
  } else if(colorset == "magma") {
    colors <- sub("FF$","",magma(nrow(annotations)))
  } else if(colorset == "inferno") {
    colors <- sub("FF$","",inferno(nrow(annotations)))
  } else if(colorset == "plasma") {
    colors <- sub("FF$","",plasma(nrow(annotations)))
  } else if(colorset == "terrain") {
    colors <- sub("FF$","",terrain.colors(nrow(annotations)))
  } else if(is.character(colorset)) {
    colors <- colorRampPalette(colorset)(nrow(annotations))
  }
  
  if(color_order == "random") {
    
    colors <- sample(colors, length(colors))
    
  }
  
  annotations <- mutate(annotations, color = colors)
  
  names(annotations) <- paste0(base, c("_label","_id","_color"))
  
  names(df)[names(df) == col] <- paste0(base,"_label")
  
  df <- left_join(df, annotations, by = paste0(base, "_label"))
  
  df
}

color_mean <- function(x) {
  library(grDevices)
  
  rgb_x <- col2rgb(x)
  rgb_mean <- rowMeans(rgb_x)
  new_hex <- rgb(rgb_mean["red"]/255,
                 rgb_mean["green"]/255,
                 rgb_mean["blue"]/255)
  
  new_hex
  
}
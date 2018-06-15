hsv_palette <- function(n_colors, 
                        hue_start = 0, 
                        hue_end = max(1, n-1)/n,
                        sat_start = 0.55,
                        sat_end = 1,
                        sat_steps = 4,
                        val_start = 1,
                        val_end = 0.8,
                        val_steps = 3) {
  sats <- rep_len(seq(sat_start, sat_end, length.out = sat_steps),length.out = n_colors)
  vals <- rep_len(seq(val_start, val_end, length.out = val_steps),length.out = n_colors)
  if(hue_end < hue_start) {
    rev(sub("FF$","",rainbow(n_colors, s = sats, v = vals, start = hue_end, end = hue_start)))
  } else {
    sub("FF$","",rainbow(n_colors, s = sats, v = vals, start = hue_start, end = hue_end))
  }
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

check_s <- function(x, min_sat = 0) {
  library(grDevices)
  
  hsv_x <- rgb2hsv(col2rgb(x))
  
  hsv_x["s",] < min_sat
  
}

check_v <- function(x, min_val = 0) {
  library(grDevices)
  
  hsv_x <- rgb2hsv(col2rgb(x))
  
  hsv_x["v",] < min_val
  
}

adjust_s <- function(x, shift_sat = 0.3) {
  library(grDevices)
  
  hsv_x <- rgb2hsv(col2rgb(x))
  
  hsv_x["s", ] <- hsv_x["s", ] + shift_sat
  if(hsv_x["s", ] > 1) {
    hsv_x["s", ] <- 1
  } else if (hsv_x["s",] <- 0) {
    hsv_x["s", ] <- 0
  }
  
  hsv(hsv_x[1,], hsv_x[2,], hsv_x[3,])
}

adjust_v <- function(x, shift_val = 0.3) {
  library(grDevices)
  
  hsv_x <- rgb2hsv(col2rgb(x))
  
  hsv_x["v", ] <- hsv_x["v", ] + shift_val
  if(hsv_x["v", ] > 1) {
    hsv_x["v", ] <- 1
  } else if (hsv_x["v",] <- 0) {
    hsv_x["v", ] <- 0
  }
  
  hsv(hsv_x[1,], hsv_x[2,], hsv_x[3,])
}

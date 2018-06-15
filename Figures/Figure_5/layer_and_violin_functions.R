get_node_dend <- function(x, match_attr, match_value) {
  
  output <- NULL
  
  for(i in seq_len(length(x))) {
    
    if(attr(x[[i]], match_attr) == match_value) {
      
      output <- x[[i]]
      
    } else {
      if(is.list(x[[i]])) {
        nest <- get_node_dend(x[[i]], match_attr, match_value)
        if(!is.null(nest)) {
          output <- nest
        }
      }
    }
    
  }
  return(output)
  
}

build_layer_plot <- function(anno,
                             dend,
                             cocl,
                             dendcluster_ids,
                             seed_val = 42,
                             right_pad = 10) {
  
  cluster_anno <- anno %>%
    select(cl, dendcluster_id, cluster_id, cluster_label, cluster_color) %>%
    unique()
  
  keep_layers <- c("L1","L2/3","L4","L5","L6")
  #keep_layers <- c("L2/3","L4","L5","L6","L6b")
  
  xpad <- 0.1
  ypad <- 0.05
  
  layer_ranges <- data.frame(layer_label = rev(keep_layers),
                             ymin = (1:5 - 1) + ypad,
                             ymax = 1:5 - ypad)
  
  filtered_anno <- anno %>%
    filter(dendcluster_id %in% dendcluster_ids) %>%
    filter(layer_label %in% keep_layers)
  
  cluster_ranges <- filtered_anno %>%
    select(cluster_id, cluster_color, cluster_label, dendcluster_id) %>%
    unique() %>%
    arrange(dendcluster_id) %>%
    mutate(xmin = 1:n() - 1 + xpad,
           xmax = 1:n()     - xpad,
           xmid = 1:n() - 0.5)
  
  set.seed(seed_val)
  
  plot_data <- filtered_anno %>%
    select(sample_id,
           dendcluster_id, cluster_color, cluster_label,
           layer_id, layer_color, layer_label) %>%
    left_join(layer_ranges) %>%
    left_join(cluster_ranges) %>%
    group_by(dendcluster_id, layer_id) %>%
    mutate(x = runif(n(),xmin + xpad, xmax - xpad),
           y = runif(n(),ymin + ypad, ymax - ypad),
           fill_color = cluster_color)
  
  # Layer background rectangles
  layer_rects <- layer_ranges %>%
    mutate(xmin = min(cluster_ranges$xmin) - xpad, xmax = max(cluster_ranges$xmax) + xpad) %>%
    mutate(fill = c("#ECE09C","#F7F2DA","#A7D7DF","#C1E5E7","#D4EDED"))
    #mutate(fill = c("#FDE4DF","#ECE09C","#F7F2DA","#A7D7DF","#C1E5E7"))
    
  # Cluster color highlights at bottom of the plot
  cluster_rects <- cluster_ranges %>%
    mutate(ymin = -ypad, ymax = ypad)
  
  # Filter the dendrogram
  prune_dend_labels <- labels(dend)[!labels(dend) %in% filtered_anno$cluster_label]
  filtered_dend <- dend %>%
    prune.dendrogram(prune_dend_labels)
  dend_seg <- as.ggdend(filtered_dend)$segments %>%
    mutate(y = (y/max(y))*3 + max(layer_rects$ymax) + ypad,
           yend = (yend/max(yend))*3 + max(layer_rects$ymax) + ypad,
           x = x - 0.5,
           xend = xend - 0.5)
  
  pad_rect <- data.frame(ymin = min(layer_rects$ymin),
                         ymax = max(layer_rects$ymax),
                         xmin = max(layer_rects$xmax),
                         xmax = max(layer_rects$xmax) + max(layer_rects$xmax)*(right_pad/100)/(1 - right_pad/100))
  
  p <- ggplot() +
    # right side padding for alignment
    geom_rect(data = pad_rect,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = "#FFFFFF",
                  color = "#FFFFFF")) +
    # dendrogram segments
    geom_segment(data = dend_seg,
                 aes(x = x, xend = xend,
                     y = y, yend = yend,
                     size = lwd,
                     color = col),
                 lineend = "square") +
    # layer background rectangles
    geom_rect(data = layer_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill)) +
    # cluster label rectangles
    geom_rect(data = cluster_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = cluster_color)) +
    geom_rect(data = cluster_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = -2 - ypad, ymax = -2,
                  fill = cluster_color)) +
    # jittered cell points
    geom_point(data = plot_data,
               aes(x = x,
                   y = y,
                   color = cluster_color),
               size = 0.1) +
    # cluster name labels
    geom_rect(data = cluster_ranges,
              aes(xmin = xmid - 0.5 + xpad/2,
                  xmax = xmid + 0.5 - xpad/2,
                  ymax = 0 - ypad,
                  ymin = -2),
              fill = "#CAD7D7")+
    geom_text(data = cluster_ranges,
              aes(x = xmid,
                  y = -2 + ypad,
                  label = cluster_label),
              angle = 90,
              vjust = 0.3,
              hjust = 0,
              size = 1.5) +
    # Plot settings
    scale_color_identity() +
    scale_size(range = c(0.5,1), guide = FALSE) +
    scale_fill_identity() +
    scale_y_continuous(limits = c(-2.1, 8)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_void(base_size = 7) +
    theme(text = element_text(size = 6),
          legend.box.spacing = unit(0,"pt"))
  
  p
}

group_violin_plot2 <- function (genes = c("Hspa8", "Snap25", "Gad2", "Vip"), group_by = "final", 
                                clusters = 1:10, data_source = "internal", sort = F, logscale = F, 
                                showcounts = T, rotatecounts = F, fontsize = 7, labelheight = 25, 
                                max_width = 10) 
{
  library(dplyr)
  library(ggplot2)
  genes <- rev(genes)
  if (data_source == "internal") {
    data <- get_internal_data(genes, group_by, clusters)
  }
  else if (is.list(data_source)) {
    data <- get_list_data(data_source, genes, group_by, clusters)
  }
  else if (grepl("\\.db$", data_source)) {
    data <- get_db_data(data_source, genes, group_by, clusters)
  }
  else if (file.exists(paste0(data_source, "/anno.feather"))) {
    data <- get_feather_data(data_source, genes, group_by, 
                             clusters)
  }
  else {
    stop("Cannot identify data_source.")
  }
  genes <- sub("-", ".", genes)
  genes <- genes[genes %in% names(data)]
  data <- data %>% select(-xpos) %>% mutate(xpos = plot_id)
  genes[grepl("^[0-9]", genes)] <- paste0("X", genes[grepl("^[0-9]", 
                                                           genes)])
  names(data)[grepl("^[0-9]", genes)] <- paste0("X", names(data)[grepl("^[0-9]", 
                                                                       genes)])
  
  ngenes <- length(genes)
  nclust <- length(unique(data$plot_id))
  nsamples <- nrow(data)
  max_vals <- data %>% select(one_of(genes)) %>% summarise_each(funs(max)) %>% 
    unlist()
  data[genes] <- data[genes] + runif(nrow(data), 0, 1e-05)
  for (i in 1:length(genes)) {
    gene <- genes[i]
    gene_max <- max_vals[i]
    if (logscale) {
      data[gene] <- log10(data[gene] + 1)/log10(gene_max + 
                                                  1) * 0.9 + i
    }
    else {
      data[gene] <- data[gene]/gene_max * 0.9 + i
    }
  }
  print(names(data))
  header_labels <- build_header_labels(data = data, grouping = "plot", 
                                       ymin = ngenes + 1,
                                       labelheight = labelheight, 
                                       label_type = "simple")
  max_labels <- data.frame(x = (nclust + 0.5) * 1.01, y = 1:ngenes + 
                             0.5, label = sci_label(max_vals))
  max_header <- data.frame(x = (nclust + 0.5) * 1.01, y = ngenes + 
                             1, label = "Max value")
  max_width <- nclust * (max_width/100)/(1 - max_width/100)
  label_y_size <- max(header_labels$ymax) - min(header_labels$ymin)
  cluster_data <- data %>% group_by(plot_label, plot_color, 
                                    plot_id) %>% summarise(cn = n()) %>% as.data.frame(stringsAsFactors = F) %>% 
    arrange(plot_id) %>% mutate(labely = ngenes + label_y_size * 
                                  0.05, cny = max(header_labels$ymax) - 0.1 * label_y_size, 
                                xpos = plot_id)
  background_rects <- data.frame(ymin = 1:length(genes),
                                 ymax = 1:length(genes) + 1,
                                 xmin = 0.5,
                                 xmax = nclust + 0.5) %>%
    mutate(fill = ifelse(row_number() %% 2 == 0, "#CAD7D7", "#FFFFFF"))
  
  p <- ggplot() + 
    scale_fill_identity() + 
    scale_y_continuous("", 
                       breaks = 1:length(genes) + 0.45, 
                       labels = genes, 
                       expand = c(0,0)) + 
    scale_x_continuous("", expand = c(0, 0)) + 
    theme_classic(fontsize) + 
    theme(axis.text = element_text(size = rel(1)), 
          axis.text.x = element_blank(), 
          axis.ticks = element_blank(), 
          axis.text.y = element_text(face = "italic", color = "#000000"),
          legend.position = "none",
          axis.line = element_line(size = 0.1)) +
    geom_rect(data = background_rects,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = fill))
  
  for (i in 1:length(genes)) {
    
    p <- p + geom_violin(data = data, 
                         aes_string(x = "xpos", 
                                    y = genes[i], 
                                    fill = "plot_color"), 
                         scale = "width", 
                         size = 0.1,
                         adjust = 2) + 
      stat_summary(data = data, aes_string(x = "xpos", 
                                           y = genes[i]), 
                   fun.y = "median", 
                   fun.ymin = "median", 
                   fun.ymax = "median", 
                   geom = "point", 
                   size = 0.1)
    
    
  }
  p <- p + geom_rect(data = header_labels, 
                     aes(xmin = xmin, 
                         ymin = ymin, 
                         xmax = xmax, 
                         ymax = ymax, 
                         fill = color)) + 
    geom_rect(aes(xmin = nclust +  0.5, 
                  xmax = nclust + 0.5 + max_width, 
                  ymin = 1, 
                  ymax = max(header_labels$ymax)), 
              fill = "white") + 
    geom_text(data = max_labels, 
              aes(x = x, 
                  y = y, 
                  label = label), 
              hjust = 0, 
              vjust = 0.35, 
              size = pt2mm(fontsize), 
              parse = TRUE)
  
  if (showcounts) {
    if (rotatecounts) {
      p <- p + geom_text(data = cluster_data, aes(y = cny, 
                                                  x = xpos, label = cn), angle = 90, vjust = 0.35, 
                         hjust = 1, size = pt2mm(fontsize))
    }
    else {
      p <- p + geom_text(data = cluster_data, aes(y = cny, 
                                                  x = xpos, label = cn), size = pt2mm(fontsize))
    }
  }
  p
}
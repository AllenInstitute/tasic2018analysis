net_graph_plot <- function(net_graph,
                           net_df,
                           anno,
                           val_col, val_min = 0, val_scale = 1,
                           link_color = "gray40",
                           split = FALSE) {
  
  node_pos <- data.frame(cl = as.numeric(names(V(g))),
                         x = l[,1],
                         y = l[,2])
  
  net_pos <- net_df %>%
    left_join(node_pos, by = c("cl.x" = "cl")) %>%
    rename(xend = x, yend = y) %>%
    left_join(node_pos, by = c("cl.y" = "cl")) %>%
    filter(cl.x != cl.y)
  
  net_pos[[val_col]] <- net_pos[[val_col]] * val_scale
  
  net_pos <- net_pos[net_pos[[val_col]] >= val_min,]
  
  net_pos[[val_col]] <- net_pos[[val_col]]**2
  
  node_pos <- left_join(node_pos, anno)
  
  x_range <- range(node_pos$x) + c(-0.5,0.5)
  y_range <- range(node_pos$y) + c(-0.5,0.5)
  
  if(split) {
    p_nodes <- ggplot() +
      geom_point(data = node_pos,
                 aes(x = x, 
                     y = y, 
                     size = n_core, 
                     color = cluster_color)) +
      geom_text(data = node_pos,
                aes(x = x, 
                    y = y, 
                    label = cluster_label),
                size = 2) +
      scale_size_area(max_size = 15, breaks = c(50,100,250,500,1000)) +
      scale_color_identity() +
      scale_fill_identity() +
      theme_void() +
      xlim(x_range) +
      ylim(y_range)
    
    p_links <- ggplot() +
      geom_segment(data = net_pos,
                   aes_string(x = "x", xend = "xend",
                              y = "y", yend = "yend",
                              size = val_col),
                   color = link_color,
                   alpha = 0.4,
                   lineend = "round") +
      scale_color_identity() +
      scale_fill_identity() +
      theme_void() +
      xlim(x_range) +
      ylim(y_range)
    
    if(val_col == "cor") {
      
      p_links <- p_links + scale_size_area(max_size = 4, 
                 breaks = c(0.2,0.4,0.6,0.8,1.0)**2,
                 labels = c(0.2,0.4,0.6,0.8,1.0))
      
    } else if(val_col == "cocl") {
      
      p_links <- p_links + scale_size_area(max_size = 5, 
                          breaks = c(0.052,0.1,0.2,0.3,0.4,0.5)**2,
                          labels = c(0.052,0.1,0.2,0.3,0.4,0.5))
      
    } else if(val_col == "trans_n") {
      p_links <- p_links + scale_size_area(max_size = 5,
                                      breaks = c(1,3, 10, 50, 100, 200)**2,
                                      labels = c(1,3, 10, 50, 100, 200))
    }
        
    list(nodes = p_nodes,
         links = p_links)
  } else {
    p <- ggplot() +
      geom_segment(data = net_pos,
                   aes_string(x = "x", xend = "xend",
                              y = "y", yend = "yend",
                              size = val_col),
                   color = link_color,
                   alpha = 0.4,
                   lineend = "round") +
      geom_point(data = node_pos,
                 aes(x = x, 
                     y = y, 
                     size = n_core, 
                     color = cluster_color)) +
      geom_text(data = node_pos,
                aes(x = x, 
                    y = y, 
                    label = cluster_label),
                size = 2) +
      scale_size_area(max_size = 15, breaks = c(50,100,250,500,1000)) +
      scale_color_identity() +
      scale_fill_identity() +
      theme_void()
    
    p
  }

  
}

net_coord_plot <- function(coords,
                           net_df,
                           anno,
                           cluster_ids = 1:25,
                           val_col, val_min = 0, val_scale = 1,
                           link_color = "gray40",
                           split = FALSE) {
  
  keep_cl <- unique(anno$cl[anno$cluster_id %in% cluster_ids])
  
  node_pos <- data.frame(cl = coords$cl,
                         x = coords$const_x,
                         y = coords$const_y) %>%
    filter(cl %in% keep_cl)
  
  net_pos <- net_df %>%
    left_join(node_pos, by = c("cl.x" = "cl")) %>%
    rename(xend = x, yend = y) %>%
    left_join(node_pos, by = c("cl.y" = "cl")) %>%
    filter(cl.x %in% keep_cl & cl.y %in% keep_cl) %>%
    filter(cl.x != cl.y)
  
  net_pos[[val_col]] <- net_pos[[val_col]] * val_scale
  
  net_pos <- net_pos[net_pos[[val_col]] >= val_min,]
  
  node_pos <- left_join(node_pos, anno)
  
  x_range <- range(node_pos$x) + c(-0.5,0.5)
  y_range <- range(node_pos$y) + c(-0.5,0.5)
  
  if(split) {
    net_pos[[val_col]] <- net_pos[[val_col]]**2
    
    p_nodes <- ggplot() +
      geom_point(data = node_pos,
                 aes(x = x, 
                     y = y, 
                     size = n_core, 
                     color = cluster_color)) +
      geom_text(data = node_pos,
                aes(x = x, 
                    y = y, 
                    label = cluster_label),
                size = 2) +
      scale_size_area(max_size = 15, range = c(0.1, 15), breaks = c(50,100,250,500,1000)) +
      scale_color_identity() +
      scale_fill_identity() +
      theme_void() +
      xlim(x_range) +
      ylim(y_range)
    
    p_links <- ggplot() +
      geom_segment(data = net_pos,
                   aes_string(x = "x", xend = "xend",
                              y = "y", yend = "yend",
                              size = val_col),
                   color = link_color,
                   alpha = 0.4,
                   lineend = "round") +
      scale_color_identity() +
      scale_fill_identity() +
      theme_void() +
      scale_y_reverse() +
      xlim(x_range) +
      ylim(y_range)
    
    if(val_col == "cor") {
      
      p_links <- p_links + scale_size_area(max_size = 4, 
                                           breaks = c(0.2,0.4,0.6,0.8,1.0)**2,
                                           labels = c(0.2,0.4,0.6,0.8,1.0))
      
    } else if(val_col == "cocl") {
      
      p_links <- p_links + scale_size_area(max_size = 10, 
                                           breaks = c(0.052,0.1,0.2,0.3,0.4,0.5)**2,
                                           labels = c(0.052,0.1,0.2,0.3,0.4,0.5))
      
    } else if(val_col == "interm") {
      p_links <- p_links + scale_size_area(max_size = 10,
                                           breaks = c(3, 10, 50, 100, 200)**2,
                                           labels = c(3, 10, 50, 100, 200))
    }
    
    list(nodes = p_nodes,
         links = p_links)
  } else {
    net_pos[[val_col]] <- net_pos[[val_col]]/2
    
    scaler <- data.frame(x = 2, y = 2, size = 500, color = "white")
    
    p <- ggplot() +
      geom_point(data = scaler,
                 aes(x = x, y = y, size = size, color = color)) +
      geom_segment(data = net_pos,
                   aes_string(x = "x", xend = "xend",
                              y = "y", yend = "yend",
                              size = val_col),
                   color = link_color,
                   alpha = 0.4,
                   lineend = "round") +
      geom_point(data = node_pos,
                 aes(x = x, 
                     y = y, 
                     size = n_core, 
                     color = cluster_color)) +
      geom_text(data = node_pos,
                aes(x = x, 
                    y = y, 
                    label = cluster_label),
                size = 2) +
      scale_size_area(max_size = 11, breaks = c(1.5,3,5,10,25,50,100,125,250,500,1000)) +
      scale_color_identity() +
      scale_fill_identity() +
      theme_void()
    
    p
  }
  
  
}

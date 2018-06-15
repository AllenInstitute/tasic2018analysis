condition_to_summarise_call <- function(pos_condition, neg_condition, fpkm_cutoff) {
  call <- paste0("sum(")
  pos_condition_comparisons <- paste(paste0(pos_condition," > ",fpkm_cutoff), collapse = " & ")
  if(!is.null(neg_condition)) {
    neg_condition_comparisons <- paste(paste0(neg_condition," > ",fpkm_cutoff), collapse = " & ")
    call <- paste0(call, pos_condition_comparisons, " & !(", neg_condition_comparisons,")")
  } else {
    call <- paste0(call, pos_condition_comparisons)
  }
  call <- paste0(call,")/n()")
  call
}


coexpression_barplot <- function(condition_1, condition_2, condition_3,
                                 frac_anno, fpkm_cutoff = 1, group_by = "dendcluster",
                                 groups = 1:128,
                                 anno, keep_frac = frac_anno$frac_label, data) {
  
  
  group_id <- paste0(group_by,"_id")
  group_label <- paste0(group_by,"_label")
  
  group_filter <- paste0(group_id," %in% c(",paste0(groups,collapse = ","),")")
  
  anno <- anno %>%
    filter_(group_filter)
  
  # recompute filtered positions
  group_order <- data.frame(group = groups) %>%
    mutate(group_order = 1:n())
  names(group_order)[1] <- group_id
  
  group_anno <- anno %>%
    ungroup() %>%
    select(one_of(paste0(group_by,c("_id","_label","_color")))) %>%
    unique() %>%
    left_join(group_order) %>%
    arrange(group_order) %>%
    mutate(group_xpos = 1:n())
  
  anno <- left_join(anno,group_anno)
  
  # pull data for all genes used in conditions
  
  cond_data <- data[,c("sample_id",condition_1,condition_2,condition_3)]
  
  # join gene values to the annotations
  anno_data <- anno %>%
    left_join(cond_data, by = "sample_id")
  
  # for each group, calculate the fraction of cells that match the conditions
  
  if(!is.null(condition_3)) {
    frac_data <- anno_data %>%
      group_by(group_xpos) %>%
      summarise_(frac_1 = condition_to_summarise_call(condition_1, c(condition_2, condition_3), fpkm_cutoff),
                 frac_2 = condition_to_summarise_call(condition_2, c(condition_1, condition_3), fpkm_cutoff),
                 frac_3 = condition_to_summarise_call(condition_3, c(condition_1, condition_2), fpkm_cutoff),
                 frac_12 = condition_to_summarise_call(c(condition_1, condition_2), condition_3, fpkm_cutoff),
                 frac_23 = condition_to_summarise_call(c(condition_2, condition_3), condition_1, fpkm_cutoff),
                 frac_13 = condition_to_summarise_call(c(condition_1, condition_3), condition_2, fpkm_cutoff),
                 frac_123 = condition_to_summarise_call(c(condition_1, condition_2, condition_3), NULL, fpkm_cutoff))
  } else {
    frac_data <- anno_data %>%
      group_by(group_xpos) %>%
      summarise_(frac_1 = condition_to_summarise_call(condition_1, condition_2, fpkm_cutoff),
                 frac_2 = condition_to_summarise_call(condition_2, condition_1, fpkm_cutoff),
                 frac_12 = condition_to_summarise_call(c(condition_1, condition_2), NULL, fpkm_cutoff))
    
    frac_anno <- frac_anno %>%
      filter(!grepl("3",frac_label))
  }
  
  frac_anno <- frac_anno %>%
    filter(grepl(paste(paste0(keep_frac, "$"), collapse="|"), frac_label))
  
  plot_data <- melt(frac_data, "group_xpos") %>%
    filter(value > 0) %>%
    rename_("frac_label" = "variable") %>%
    left_join(frac_anno) %>%
    arrange(group_xpos, frac_id) %>%
    group_by(group_xpos) %>%
    mutate(xmin = group_xpos - 0.4,
           xmax = group_xpos + 0.4) %>%
    mutate(csum = cumsum(value),
           ymin = lag(csum, default = 0),
           ymax = csum)
  
  ggplot() + 
    geom_rect(data = plot_data,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = frac_color)) +
    scale_fill_identity(breaks = frac_anno$frac_color,
                        labels = frac_anno$frac_label,
                        guide = "legend") +
    scale_x_continuous(breaks = group_anno$group_xpos, 
                       labels = group_anno[[group_label]],
                       limits = c(0.5, max(group_anno$group_xpos) + 0.5)) +
    scale_y_continuous("Fraction of Cells in Cluster", expand = c(0,0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3, colour=group_anno$cluster_color))
}

prune.dendrogram <- function(dend, leaves, reindex_dend = TRUE, ...) {
  leaves <- as.character(leaves)
  
  for(i in seq_along(leaves))
  {
    # this function is probably not the fastest - but it works...
    dend <- prune_leaf(dend, leaves[i])	# move step by stem to remove all of these leaves...
  }
  
  if(reindex_dend) dend <- reindex_dend(dend)
  
  return(dend)
}

stats_midcache.dendrogram <- function (x, type = "hclust", quiet = FALSE) 
{
  type <- match.arg(type)
  stopifnot(inherits(x, "dendrogram"))
  setmid <- function(d, type) {
    if (is.leaf(d)) 
      return(d)
    k <- length(d)
    if (k < 1) 
      stop("dendrogram node with non-positive #{branches}")
    r <- d
    midS <- 0
    for (j in 1L:k) {
      r[[j]] <- unclass(setmid(d[[j]], type))
      midS <- midS + .midDend(r[[j]])
    }
    if (!quiet && type == "hclust" && k != 2) 
      warning("midcache() of non-binary dendrograms only partly implemented")
    attr(r, "midpoint") <- (.memberDend(d[[1L]]) + midS)/2
    r
  }
  setmid(x, type = type)
}

stats_.midDend <- function (x) {
  if (is.null(mp <- attr(x, "midpoint"))) 0 else mp
}
.midDend <- stats_.midDend # copied so that they would work inside the various functions here...

stats_.memberDend <- function (x) 
{
  r <- attr(x, "x.member")
  if (is.null(r)) {
    r <- attr(x, "members")
    if (is.null(r)) 
      r <- 1L
  }
  r
}
.memberDend <- stats_.memberDend


prune_leaf <- function(dend, leaf_name,...)
{
  labels_dend <- labels(dend)
  
  if(length(labels_dend) != length(unique(labels_dend)))	warning("Found dubplicate labels in the tree (this might indicate a problem in the tree you supplied)")
  
  if(!(leaf_name %in% labels_dend)) {	# what to do if there is no such leaf inside the tree
    warning(paste("There is no leaf with the label", leaf_name , "in the tree you supplied", "\n" , "Returning original tree", "\n" ))
    return(dend)
  }
  
  if(sum(labels_dend %in% leaf_name) > 1) {	# what to do if there is no such leaf inside the tree      
    warning(paste("There are multiple leaves by the name of '", leaf_name , "' in the tree you supplied.  Their locations is:",
                  paste(which(labels_dend %in% leaf_name), collapse = ","),"\n" , "Returning original tree", "\n" ))
    return(dend)
  }
  
  is.father.of.leaf.to.remove <- function(dend, leaf_name)
  {
    # this function checks if the leaf we wish to remove is the direct child of the current branch (dend) we entered the function
    is.father <- FALSE
    for(i in seq_len(length(dend)))
    {
      if(is.leaf(dend[[i]]) == TRUE  &&  labels(dend[[i]]) == leaf_name) is.father <- TRUE
    }
    return(is.father)
  }
  
  
  remove_leaf_if_child <- function(dend, leaf_name)
  {
    # print(labels(dend))
    if(all(labels(dend) != leaf_name))
    {	# if the leaf we want to remove is not in this branch, simply return the branch without going deeper intoit.
      return(dend)
    } else {	# but if the leaf we want to remove is here somewhere, go on searching
      attr(dend, "members") <- attr(dend, "members") - 1 
      
      if(!is.father.of.leaf.to.remove(dend, leaf_name))	# if you are not the father, then go on and make this function work on each child
      {
        for(i in seq_len(length(dend)))
        {
          dend[[i]] <- remove_leaf_if_child(dend[[i]], leaf_name)
        }
      } else { # we'll merge 
        if(length(dend) == 2) {
          leaf_location <- 1 
          # if leaf location is 1, then move branch in leaf 2 to be the new x
          if(is.leaf(dend[[leaf_location]]) == T  &&  labels(dend[[leaf_location]]) == leaf_name) {

            branch_to_bumpup <- 2
            dend <- dend[[branch_to_bumpup]]
          } else { # else - the leaf location must be located in position "2"

            branch_to_bumpup <- 1
            dend <- dend[[branch_to_bumpup]]
          }
        } else if(length(dend) > 2) {
          # If more than 2 branches, check if any are leaves
          dend_leaves <- unlist(lapply(dend, is.leaf))
          if(sum(dend_leaves) > 0) {
            # If so, check for matching labels to the leaf to prune
            dend_labels <- unlist(lapply(dend, function(x) attr(x, "label")))
            dend_matches <- dend_labels == leaf_name
            # Return a list containing the non-matching branches
            dend[dend_leaves & dend_matches] <- NULL
            # Note that in some cases, the following DOES NOT yield a correct result:
            # dend <- dend[!(dend_leaves & dend_matches)]
            
            # If the length is now 1, it can be bumped up
            if(length(dend) == 1) {
              dend <- dend[[1]]
            }
            
          }
        }
      }
    }
    return(dend)
  }
  
  
  new_dend <- remove_leaf_if_child(dend, leaf_name)
  new_dend <- suppressWarnings(stats_midcache.dendrogram(new_dend)) # fixes the attributes
  #   new_x <- fix_members_attr.dendrogram(new_x) # fix the number of memebers attr for each node
  return(new_dend)
}
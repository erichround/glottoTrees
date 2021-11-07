# Manipulating tree topoplogy


#' Glottolog trees by version
#'
#' Returns a multiPhylo object containing all, or a requested subset, of the
#' glottolog trees.
#'
#' By default, trees are returned from the most recent version of glottolog.
#' Alternatively, an older version of glottolog can be specified.
#'
#' @param family A character vector. Elements are names of glottolog families
#'   whose trees are to be returned. If \code{family} is left unspecified, all
#'   trees are returned.
#' @inheritParams get_glottolog_languages
#' @return A \code{phylo} object containing one glottolog tree, or a
#'   \code{multiPhylo} object containing multiple glottolog trees.
#' @examples 
#' 
#' library(ape)
#' tree_totonacan <- get_glottolog_trees(family = "Totonacan")
#' tree_totonacan_v4.3 <- get_glottolog_trees("Totonacan", "4.3")
#' plot(tree_totonacan)
#' plot(tree_totonacan_v4.3)
#' trees <- get_glottolog_trees(family = c("Caddoan", "Tangkic"))
#' plot(trees[[1]])
#' plot(trees[[2]])
get_glottolog_trees = function(
  family,
  glottolog_version
) {
  
  # Check glottolog_version
  if (missing(glottolog_version)) {
    glottolog_version <- .get_newest_version()
  } else {
    error_msg <- .check_glottolog_version(glottolog_version)
    if (!is.na(error_msg)) { stop(error_msg) }
    glottolog_version <- as.character(glottolog_version)
  }
  
  # Choose appropriate dataset
  if (glottolog_version == "4.0") {
    phy <- glottolog_trees_v4.0
  } else if (glottolog_version == "4.1") {
    phy <- glottolog_trees_v4.1
  } else if (glottolog_version == "4.2") {
    phy <- glottolog_trees_v4.2
  } else if (glottolog_version == "4.3") {
    phy <- glottolog_trees_v4.3
  } else if (glottolog_version == "4.4") {
    phy <- glottolog_trees_v4.4
  }
  
  # If family is missing, return the whole dataset
  if (missing(family)) { return(phy) }
  
  # Else, select the requested subset of the families
  
  # Check family
  check_result <- .check_family(family, glottolog_version)
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }

  phy <- phy[which_tree(family, glottolog_version)]
  
  if (length(phy) == 1) { 
    # If only one family, return as phylo object.
    phy[[1]]
  } else {
    # Otherwise, return a multiPhylo.
    phy
  }
}


#' Bind trees as a high-level rake
#'
#' Takes a multiPhylo object containing multiple trees and combines them into a
#' single tree with a rake structure at its root, below which each tree appears
#' on its own branch.
#'
#' @param phy A multiphylo object containing the trees to be combined.
#' @return A phylo object, a single tree.
#' @examples 
#' 
#' library(ape)
#' 
#' arnhem_hypothesis <- 
#'   c("Gunwinyguan", "Mangarrayi-Maran", "Maningrida",
#'     "Kungarakany", "Gaagudju")
#' trees <- get_glottolog_trees(arnhem_hypothesis)
#' simple_rake <- assemble_rake(trees)
#' rake_in_rake <- assemble_rake(c(assemble_rake(trees[1:3]), trees[4:5]))
#' plot(simple_rake)
#' plot(rake_in_rake)
#' 
#' # If `phy` contains only one tree, a warning is issued.
#' mono_rake <- assemble_rake(trees[3])
#' plot(mono_rake)
assemble_rake = function(phy) {
  
  if (class(phy) != "multiPhylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class multiPhylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  n_trees <- length(phy)
  if (n_trees == 1) {
    warning("`phy` contained only one tree. No changes were made to it.")
    return(phy[[1]])
  }
  
  n_tip_vec <- Ntip(phy)
  n_tips <- sum(n_tip_vec)
  root <- n_tips + 1
  tip_offset_vec <- cumsum(c(0, n_tip_vec[-n_trees]))
  
  n_node_vec <- Nnode(phy) - (n_tip_vec == 1) # Remove extraneous node above isolate
  n_node_cum <- cumsum(c(1, n_node_vec[-n_trees]))
  node_offset_vec <- n_node_cum  + n_tips - n_tip_vec
  
  edge <-
    lapply(1:n_trees, function(i) {
      # Adjust node numbers
      e <- eo <- phy[[i]]$edge
      nt <- n_tip_vec[i]
      e[eo <= nt] <- eo[eo <= nt] + tip_offset_vec[i]
      e[eo > nt] <- eo[eo > nt] + node_offset_vec[i]
      # Add edge connecting to new root
      edge_to_root <- c(root, root + n_node_cum[i])
      e <- rbind(edge_to_root, e, deparse.level = 0)
      # Remove extraneous node above isolates
      if (n_tip_vec[i] == 1) {
        e <- matrix(c(root, e[2,2]), nrow = 1)
      }
      e
    }) %>%
    do.call(rbind, .)
    
  edge_lengths <- 
    lapply(phy, function(x) {
      if (!("edge.length" %in% names(x))) {
        x$edge.length <- rep(1, Nedge(x))
      }
      e <- c(1, x$edge.length)
      # Remove extraneous node above isolates
      if (Ntip(x) == 1) { e <- 1 }
      e
    } ) %>% 
    unlist()
  node_labels <- 
    lapply(phy, function(x) {
      if (Ntip(x) == 1) { 
        # Remove extraneous node above isolates
        NULL 
        } else { x$node.label }
    }) %>% 
    unlist() %>% c("", .)
  tip_labels <- unlist(lapply(phy, function(x) x$tip.label))
  
  bound_tree <- list(
    edge = edge,
    Nnode = length(node_labels),
    node.label = node_labels,
    tip.label = tip_labels,
    edge.length = edge_lengths
  )
  class(bound_tree) <- "phylo"
  bound_tree
}


#' Create a glottolog super-tree
#'
#' Combining glottolog family trees into one large tree. Families can be
#' assembled directly below a rake structure at the root, or can be grouped, so
#' that the root first branches into groups, and the families then branch out
#' below the group nodes.
#'
#' Grouping is controlled by the \code{macro_groups} parameter. Groups can
#' comprise a single glottolog macroarea, or multiple macroareas. Current
#' macroareas are \code{Africa}, \code{Australia}, \code{Eurasia}, \code{North
#' America}, \code{Papunesia} and \code{South America}. Setting
#' \code{macro_groups} to \code{NULL} causes the tree to be assembled without
#' groups.
#'
#' @param macro_groups A list of character vectors, in which each vector
#'   contains the names of one or more macroareas which define a group.
#'   Alternatively, setting \code{macro_groups} to \code{NULL} causes the tree
#'   to be assembled without groups.
#' @inheritParams get_glottolog_languages
#' @examples 
#' 
#' # Supertree whose first order branches are the glottolog macroareas
#' supertree <- assemble_supertree()
#' supertree_v.4.3 <- assemble_supertree(glottolog_version = "4.3")
#' 
#' # Supertree whose first order branches are glottolog families
#' supertree <- assemble_supertree(macro_groups = NULL)
#' 
#' # Supertree whose first order branches are the African & Eurasian macroareas
#' supertree <- assemble_supertree(macro_groups = list("Africa", "Eurasia"))
#' 
#' # Supertree whose first order branches are the glottolog macroareas, but
#' # with the Americas combined:
#' supertree <- assemble_supertree(
#'   macro_groups = list("Africa", "Australia", "Eurasia", "Papunesia", 
#'                       c("South America", "North America"))
#'   )
assemble_supertree = function(
  macro_groups,
  glottolog_version
) {
  
  # Check glottolog_version
  if (missing(glottolog_version)) {
    glottolog_version <- .get_newest_version()
  } else {
    error_msg <- .check_glottolog_version(glottolog_version)
    if (!is.na(error_msg)) { stop(error_msg) }
    glottolog_version <- as.character(glottolog_version)
  }
  
  phy <- get_glottolog_trees(glottolog_version = glottolog_version)
  
  # Check macro_groups
  available_macros <-
    unique(get_glottolog_families(glottolog_version)$main_macroarea)
  if (missing(macro_groups)) { 
    # No argument given; use the available macro groups without
    # further grouping
    macro_groups <- as.list(available_macros) 
  }
  if (!is.null(macro_groups)){
    if (!is.list(macro_groups)) {
      if (is.character(macro_groups) & 
          all(macro_groups %in% available_macros)) {
        stop(str_c("`macro_groups` must be a list.\n",
                     "You have supplied character vector.\n",
                     "Did you mean `list(", 
                     str_c(str_c("'", macro_groups, "'"), 
                           collapse = ","),
                     ")`?"))
      }
      mc <- class(macro_groups)
      stop(str_c("`macro_groups` must be a list.\n",
                   "You have supplied an object of class ", mc, "."))
    }
  }

  # Get the predominant macro_area for each tree
  main_macro <- get_glottolog_families(glottolog_version)$main_macroarea
  
  # Group family trees by macroarea
  if (is.null(macro_groups)) {
    # if phy = NULL, use no grouping at all
    macro_phys <- phy
  } else {
    # Do grouping
    macro_phys <-
      lapply(macro_groups, function(m) {
        tree_set <- which(main_macro %in% m)
        macro_tree <- 
          suppressWarnings(assemble_rake(phy[tree_set]))
        macro_label <- str_c(.name_to_label(m), collapse = "-")
        macro_tree$node.label[1] <- macro_label
        macro_tree
      }) %>%
      do.call(c, .)
  }
  
  # Group globally
  super_phy <- suppressWarnings(assemble_rake(macro_phys))
  super_phy$node.label[1] <- "World"
  super_phy$edge.length <- rep(1, Nedge(super_phy))
  
  super_phy
}


#' Keep tips
#'
#' Nominate which tips are to be kept. Others are removed.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip labels.
#' @return A phylo object containing the modified tree.
#' @examples 
#' 
#' library(ape)
#' 
#' tree <- abridge_labels(get_glottolog_trees("Tangkic"))
#' plot(tree)
#' nodelabels(tree$node.label)
#' 
#' tree2 <- keep_tip(tree, c("lard1243", "kang1283", "kaya1319"))
#' plot(tree2)
#' nodelabels(tree2$node.label)
keep_tip = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  # Keep only the named tips
  drop_tips <- setdiff(phy$tip.label, label)
  phy <- drop.tip(phy, drop_tips, collapse.singles = FALSE)
  
  phy
}


#' Keep tips and convert nodes to tips
#'
#' Nominate which tips and nodes are to be kept as tips. Others tips are
#' removed. If any element of \code{label} is both a node label and tip label,
#' it will be treated as referring to the tip, not the node.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip and node labels.
#' @return A phylo object containing the modified tree.
#' @examples
#'
#' library(ape)
#'
#' tree <- abridge_labels(get_glottolog_trees("Tangkic"))
#' plot_glotto(tree)
#'
#' tree2 <- keep_as_tip(tree, c("lard1243", "kaya1319", "nyan1300", "gang1267"))
#' plot_glotto(tree2)
keep_as_tip = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "both")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  # Call convert_to_tip() for node labels which are not also tip labels
  node_label <- label %>% setdiff(phy$tip.label)
  if (length(node_label) > 0) {
    phy <- phy %>% convert_to_tip(node_label)
  }

  phy %>% keep_tip(label)
}


#' Remove tips
#'
#' From a tree, remove tips identified by their node labels.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip labels.
#' @return A phylo object containing the modified tree.
#' @examples 
#' 
#' library(ape)
#' 
#' tree <- abridge_labels(get_glottolog_trees("Tangkic"))
#' plot_glotto(tree)
#' 
#' tree2 <- remove_tip(tree, c("kang1283", "kaya1319"))
#' plot_glotto(tree2)
remove_tip = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy

  # Check label
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }

  # Drop tips
  drop.tip(phy, label, collapse.singles = FALSE)
}


#' Add tips to a tree
#'
#' Add one or more tips below a parent node specified by its label.
#'
#' The length of the branches, between the added tips and their parent node
#' is set equal to the longest of the original branches directly below node n.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip labels.
#' @param parent_label A character string containing the label of the parent
#'   node.
#' @return A phylo object containing the modified tree.
#' @examples 
#' 
#' library(ape)
#' 
#' tree <- abridge_labels(get_glottolog_trees("LeftMay"))
#' tree <- ultrametricize(rescale_branches_exp(tree))
#' plot_glotto(tree)
#' 
#' # Attach one or more new tips to a tree:
#' tree2 <- add_tip(tree, label = "rockypeak", parent_label = "iter1240")
#' plot_glotto(tree2)
#' tree3 <- add_tip(tree, label = c("bo", "kaumifi"),  parent_label = "bopa1235")
#' plot(tree3)
#' nodelabels(tree3$node.label)
#' 
#' # Move tips by using remove_tip() and add_tip():
#' tree4 <- remove_tip(tree, "amap1240")
#' tree4a <- add_tip(tree4, "amap1240", parent_label = "left1242")
#' plot(tree4a)
#' nodelabels(tree4a$node.label)
add_tip = function(phy, label, parent_label) {
  
  # Note: TreeTools::AddTip doesn't suffice, because it scrambles node labels
  # and introduces new tips in binary branches rather than multichotomies.
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  if (!is.character(label)) {
    cl <- class(label)
    return(list(
      error_msg = str_c("`label` must be a character vector.\n",
                        "You supplied an object of class ", cl, "."),
      warning_msg = NA
    ))
  }
  n_label <- length(label)
  if (n_label == 0) {
    stop(str_c("`label` should be length 1 or more.\n",
               "You supplied a vector length 0."))
  }
  
  # Check parent label
  check_result <- .check_labels(phy, parent_label, type = "parent")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  if (length(parent_label) != 1) {
    stop(str_c("`parent_label` should be length 1.\n",
               "You supplied a vector length ", length(parent_label), "."))
  }
  
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)
  n_vertex <- n_tip + n_node
  n_edge <- Nedge(phy)
  e <- phy$edge
  parent_node <- n_tip + which(phy$node.label == parent_label)
  if (length(parent_node) > 1) {
    stop(str_c("`parent_label` should uniquely identify one node.\n",
               "You provided the value `", parent_label, "`, which matches ",
               length(parent_node), " nodes in `phy`.\n",
               "Consider using apply_duplicate_suffixes() to make node ",
               "labels unique before using add_tip()."))
  }
  
  # Get new edge length
  sibling_edges <- which(phy$edge[, 1] == parent_node)
  new_edgelength <- max(phy$edge.length[sibling_edges])
  
  # # Increase node indices by 1
  # e <- e + (e > n_tip)
  # 
  # # Add new pendant, placing the new edge & tip last among the edges & tips
  # e <- rbind(e, c(parent_node + 1, n_tip + 1))
  
  # Increase node indices by n_label
  e <- e + (e > n_tip) * n_label
  parent_node <- parent_node + n_label
  
  # Add new pendants, placing the new edges & tips last among the edges & tips
  new_pend <- matrix(c(rep(parent_node, n_label), n_tip + 1:n_label), ncol = 2)
  e <- rbind(e, new_pend)
  
  # Rebuild the tree
  new_tree <- list(
    edge = e,
    Nnode = phy$Nnode,
    node.label = phy$node.label,
    tip.label = c(phy$tip.label, label),
    edge.length = c(phy$edge.length, rep(new_edgelength, n_label))
  )
  class(new_tree) <- "phylo"
  
  # Reorder edges and reassign indices
  .reindex(new_tree)
}


#' Move a tip
#'
#' Move one tip to a new parent node.
#'
#' In the tip's new position, the length of the branch above it is the same
#' as the longest branch above any new sister of the tip.
#'
#' If moving the tip would result in any node(s) having no descendant tips, then
#' those nodes are removed.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character string containing the tip label.
#' @param parent_label A character string containing the label of the parent
#'   node.
#' @return A phylo object containing the modified tree.
#' @examples
#'
#' library(ape)
#'
#' tree <- abridge_labels(get_glottolog_trees("LeftMay"))
#' tree <- ultrametricize(rescale_branches_exp(tree))
#' plot(tree)
#' nodelabels(tree$node.label)

#' tree2 <- move_tip(tree, "amap1240", parent_label = "left1242")
#' plot(tree2)
#' nodelabels(tree2$node.label)
move_tip = function(phy, label, parent_label) {
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  if (length(label) != 1) {
    stop(str_c("`label` should be length 1.\n",
               "You supplied a vector length ", length(label), "."))
  }
  
  # Check parent label
  check_result <- .check_labels(phy, parent_label, type = "parent")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  if (length(parent_label) != 1) {
    stop(str_c("`parent_label` should be length 1.\n",
               "You supplied a vector length ", length(parent_label), "."))
  }
  
  phy %>%
    relabel_tip(label, "##REMOVE_ME##") %>%
    add_tip(label, parent_label) %>%
    remove_tip("##REMOVE_ME##")
    
}



#' Clone tips
#'
#' Clones tips as sisters of the original. Optionally, places the new clones and
#' the original in their own subgroup, in which case the node for the new
#' subgroup is assigned the same label as the original tip.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip labels.
#' @param n A numeric vector. The number of clones to make.
#' @param subgroup A logical. Whether to create a subgroup containing the new
#'   clones and their original.
#' @return A phylo object containing the modified tree.
#' @examples 
#' 
#' library(ape)
#' 
#' tree <- 
#'   rescale_branches_exp(abridge_labels(get_glottolog_trees("Tangkic")))
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- clone_tip(tree, "nyan1300")
#' plot(tree2)
#' nodelabels(tree2$node.label)
#' 
#' tree3 <- clone_tip(tree, "nyan1300", subgroup = TRUE)
#' plot(tree3)
#' nodelabels(tree3$node.label)
#' # Add suffixes to labels, to keep all labels distinct
#' tree3a <- apply_duplicate_suffixes(tree3)
#' plot(tree3a)
#' nodelabels(tree3a$node.label)
#' 
#' tree4 <- clone_tip(tree, "lard1243", n = 3)
#' plot(tree4)
#' nodelabels(tree4$node.label)
#' 
#' tree5 <- clone_tip(tree, "lard1243", n = 3, subgroup = TRUE)
#' plot(tree5)
#' nodelabels(tree5$node.label)
#' 
#' tree6 <- clone_tip(tree, c("lard1243", "nyan1300"), n = 2, subgroup = TRUE)
#' plot(tree6)
#' nodelabels(tree6$node.label)
#' tree6a <- apply_duplicate_suffixes(tree6)
#' plot(tree6a)
#' nodelabels(tree6a$node.label)
#' 
#' \dontrun{
#' # Returns error if any element of `label` is not in `phy`
#' tree7 <-  clone_tip(tree, c("lard1243", "xxxx1234"))
#' }
clone_tip = function(
  phy, 
  label, 
  n = 1, 
  subgroup = FALSE
) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  phy$node.label[phy$node.label == ""] <- "##NOLABEL##"
  phy$tip.label[phy$tip.label == ""] <- "##NOLABEL##"
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  # Check for duplicate tips that match label
  l_tips <- which(phy$tip.label %in% label)
  dup_l_tips <- duplicated(phy$tip.label[l_tips])
  dup_tip_labels <- unique(phy$tip.label[l_tips[dup_l_tips]])
  if (any(dup_l_tips)) {
    stop(str_c("Cannot clone already-duplicated tips: ",
               "`label` must not contain names of tip labels ",
               "that occur more than once in `phy`.\n",
               "You supplied one or more labels that match ",
               "more than one tip: ",
               str_c(head(dup_tip_labels, 4), collapse = ","),
               ifelse(length(dup_tip_labels) > 4, "..", ""), ".\n",
               "Use `apply_duplicate_suffixes()` to add suffixes ",
               "to labels, to ensure that none are duplicates."
               ))
  }
  
  # Check n
  if (!is.numeric(n)) {
    cn <- class(n)
    stop(str_c("`n` must be numeric.\n",
               "You supplied an object of class ", cn, "."))
  }
  n_n <- length(n)
  n_label <- length(label)
  if (n_n != 1 & n_n != n_label) {
    stop(str_c("`n` must be length 1 or the same length as `label`.\n",
               "`label` is length ", n_label, " but `n` is length ", n_n, "."))
  } else if (n_n == 1) {
    n <- rep(n, n_label)
  }
  
  n_label <- length(label)
  is_dupl <- duplicated(label)
  
  for (i in (1:n_label)[!is_dupl]) {
    l <- label[i]
    l_tip <- which(phy$tip.label == l)[1]
    orig_edge_length <- phy$edge.length[which(phy$edge[,2] == l_tip)]
    
    for(j in 1:n[i]) {
      
      n_tips <- length(phy$tip.label)
      l_tip <- which(phy$tip.label == l)[1]
      l_parent <- 
        ifelse(subgroup, 
               which(phy$node.label == "NA") + n_tips,
               phy$edge[phy$edge[,2] == l_tip, 1])
      
      # Add a tip at the destination
      phy <-
        phy %>%
        bind.tip(
          edge.length = 1,
          # Temporary label avoids unexpected behaviour from bind.tip():
          tip.label = str_c("#TEMP#"),
          where = ifelse(subgroup & (j == 1), l_tip, l_parent)
          )
      phy$tip.label[phy$tip.label == "#TEMP#"] <- l
    }
    
    # Set total edge lengths the same as the original:
    n_tips <- length(phy$tip.label)
    is_target_daughter <- c(phy$tip.label == l, phy$node.label == "NA")
    target_edges <- which(phy$edge[,2] %in% which(is_target_daughter))
    if (subgroup) {
      phy$edge.length[target_edges] <- orig_edge_length / 2
    } else {
      phy$edge.length[target_edges] <- orig_edge_length
    }
    
    # label the newly created node
    phy$node.label[phy$node.label == "NA"] <- l
  }
  
  phy$node.label[phy$node.label == "##NOLABEL##"] <- ""
  phy$tip.label[phy$tip.label == "##NOLABEL##"] <- ""
  
  # Place clones adjacent to one another
  .group_cloned_sisters(phy)
}


#' Convert nodes to tips
#'
#' For nodes identified by their label, remove the clade that they dominate and
#' replace it with a tip of the same name.
#'
#' Any labels included in \code{label} which are both node labels and tip labels
#' are regarded as node label, and that node will be removed and replaced by a
#' tip.
#'
#' Any labels included in \code{label} which are tip labels only are ignored.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing node or tip labels.
#' @param warn A logical, whether to issue warning when \code{label} cantains
#'   tip labels.
#' @return A phylo object containing the modified tree.
#' @examples 
#' 
#' tree1 <- abridge_labels(get_glottolog_trees("GreatAndamanese"))
#' plot_glotto(tree1)
#' tree2 <- convert_to_tip(tree1, label = c("okol1242", "sout2683"))
#' plot_glotto(tree2)
convert_to_tip = function(phy, label, warn = TRUE) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check label
  check_result <- .check_labels(phy, label, type = "both")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  
  # Ascertain which elements of label are node labels
  n_tip <- Ntip(phy)
  e <- phy$edge
  node_labels <- intersect(label, phy$node.label)
  ignored_labels <- setdiff(label, node_labels)
  nodes <- match(node_labels, phy$node.label) + n_tip
  ignoreds <- match(ignored_labels, phy$tip.label)
  
  # Return error if any nodes targeted for replacement dominate any other
  # node or tip in label
  for (nd in nodes) {
    descendants <- getDescendants(phy, nd)
    dom_nd <- nodes[nodes %in% descendants]
    dom_ig <- ignoreds[ignoreds %in% descendants]
    n_dom_nd <- length(dom_nd)
    n_dom_ig <- length(dom_ig)
    if (n_dom_nd > 0) {
      stop(str_c("It is not possible to convert two nodes to tips if ",
                 "one of them dominates the other.\n",
                 "You provided node ", phy$node.label[nd - n_tip],
                 ", as well as ", phy$node.label[dom_nd[1] - n_tip],
                 ", which it dominates."
                 ))
    }
    if (n_dom_ig > 0) {
      stop(str_c("It is not possible to convert a node to a tip and ",
                 "also preserve a tip that is dominates.\n",
                 "You provided node ", phy$node.label[nd - n_tip],
                 ", as well as ", phy$tip.label[dom_ig[1]],
                 ", which is a tip that it dominates."
      ))
    }
  }
  
  # Warn if any labels are not node labels
  n_ignored <- length(ignored_labels)
  if (warn & n_ignored > 0) {
    warning(str_c("Any elements in `label` which are tip labels and ",
                  "not node labels are ignored.\n", 
                  "You provided ", n_ignored, " of this kind: ",
                  str_c(head(ignored_labels, 4), collapse = ", "),
                  ifelse(n_ignored > 4, "..", ""), "."))
  }
  
  # Warn if any targeted node label is also a tip label in the tree
  double_labels <- intersect(node_labels, phy$tip.label)
  n_double <- length(double_labels)
  if (n_double > 0) {
    warning(str_c("Requested nodes are converted tips, even if a tip with ",
                  "the same label is present elsewhere in the tree.\n",
                  "You provided ", n_double, " label(s) which match both ",
                  "a node and a tip, in which case the node has been ",
                  "converted and the tip left unchanged: ",
                  str_c(head(double_labels, 4), collapse = ", "),
                  ifelse(n_double > 4, "..", ""), "."))
  }
  
  # Get the parents of the nodes to be converted
  parents <- e[match(nodes, e[,2]), 1]
  parent_labels <- phy$node.label[parents - n_tip]

  # "Convert" by adding a tip and removing the node
  for (i in 1:length(node_labels)) {
    phy <-
      phy %>%
      add_tip(node_labels[i], parent_labels[i]) %>%
      remove_clade(node_labels[i])
  }
  phy
}


#' Remove clades
#'
#' From a tree, remove clacdes, identified by the labels of their deepest node.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing node labels.
#' @return A phylo object containing the modified tree.
#' @examples 
remove_clade = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check label
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  
  # Get all descendant tips of target nodes
  n_tip <- Ntip(phy)
  target_nodes <- match(label, phy$node.label) + n_tip
  descendants <- getDescendants(phy, target_nodes)
  target_tips <- descendants[descendants <= n_tip]
  
  # Drop tips
  drop.tip(phy, target_tips, collapse.singles = FALSE)
}


#' Collapse one or more nodes rootwards
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing node labels.
#' @return A phylo object containing the modified tree.
#' @examples 
#' 
#' library(ape)
#' tree <- 
#'   rescale_branches_exp(abridge_labels(get_glottolog_trees("Tangkic")))
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- collapse_node(tree, "gang1267")
#' plot(tree2)
#' nodelabels(tree2$node.label)
collapse_node = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)
  vertices <- 1:(n_tip + n_node)
  e <- phy$edge
  root <- setdiff(e[,1], e[,2])
  
  # Identify the target nodes and their parents
  target_nodes <- match(label, phy$node.label) + n_tip
  if (root %in% target_nodes) {
    stop(str_c("The root node cannot be collapsed rootwards.\n",
               "In `label` you supplied the root node: ", 
               phy$node.label[root - n_tip], "."))
  }
  n_target <- length(target_nodes)
  target_edges <- match(target_nodes, e[,2])
  target_parents <- e[target_edges, 1]
  
  # So long as any target's parent is itself a target node, then get the next
  # deeper ancestor
  while (any(target_parents %in% target_nodes)) {
    is_bad_p <- target_parents %in% target_nodes
    bad_p_edges <- match(target_parents[is_bad_p], e[,2])
    target_parents[is_bad_p] <- e[bad_p_edges, 1]
  }
  
  ## Remappings of vertex indices
  # Map targets to parents
  v_remap <- vertices
  v_remap[target_nodes] <- target_parents
  # Adjust for removal of targets
  adjustment <- cumsum(vertices %in% target_nodes)
  v_remap <- (v_remap - adjustment)[v_remap]
    
  # Rebuild the tree
  new_tree <- list(
    edge = matrix(v_remap[e], ncol = 2)[-target_edges, ],
    Nnode = n_node - n_target,
    node.label = phy$node.label[-(target_nodes - n_tip)],
    tip.label = phy$tip.label,
    edge.length = phy$edge.length[-target_edges]
  )
  class(new_tree) <- "phylo"
  .reindex(new_tree)
}


#' Find nodes with a single child
#' 
#' @param phy A phylo object
#' @return A vector of node labels
nonbranching_nodes = function(phy) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  s <- as.data.frame(phy$edge) %>%
    group_by(V1) %>%
    summarise(n_child = n()) %>%
    filter(n_child == 1)
  
  phy$node.label[s$V1 - Ntip(phy)]
}

  
#' Move a node
#'
#' Move one node to a position dominated by a new parent node.
#'
#' The branch length above the moved node remains unchanged.
#'
#' If moving the node would result in any other node(s) having no descendant 
#' tips, then those other nodes are removed.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character string containing the label of the node to move.
#' @param parent_label A character string containing the label of the new parent
#'   node.
#' @return A phylo object containing the modified tree.
#' @examples
#'
#' library(ape)
#'
#' tree <- abridge_labels(get_glottolog_trees("LeftMay"))
#' tree <- ultrametricize(rescale_branches_exp(tree))
#' plot(tree)
#' nodelabels(tree$node.label)

#' tree2 <- move_node(tree, "iter1240", parent_label = "left1242")
#' plot(tree2)
#' nodelabels(tree2$node.label)
move_node = function(phy, label, parent_label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  if (length(label) != 1) {
    stop(str_c("`label` should be length 1.\n",
               "You supplied a vector length ", length(label), "."))
  }
  
  # Check parent label
  check_result <- .check_labels(phy, parent_label, type = "parent")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  if (length(parent_label) != 1) {
    stop(str_c("`parent_label` should be length 1.\n",
               "You supplied a vector length ", length(parent_label), "."))
  }
  
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)
  n_edge <- Nedge(phy)
  e <- phy$edge
  
  # Identify the target node and related nodes and tips
  target_node <- n_tip + which(phy$node.label == label)
  target_edge <- which(e[,2] == target_node)
  old_parent_node <- e[target_edge, 1]
  new_parent_node <- n_tip + which(phy$node.label == parent_label)
  new_sibling_edges <- which(e[, 1] == new_parent_node)
  descendants <- phytools::getDescendants(phy, target_node)
  descendant_edges <- which(e[, 2] %in% descendants)
  
  # Check position of new parent
  if (new_parent_node %in% descendants) {
    stop(str_c("A node cannot be moved to a position which it dominates.\n",
               "`label` ", label, " dominates `parent_label` ", 
               parent_label, "."))
  }
  if (new_parent_node == old_parent_node) {
    warning(str_c("Node ", label, " was not moved, because ",
                  "`parent_label` ", parent_label, " is already its ",
                  "parent node."))
    return(phy)
  }
  
  # Change the parent of the target node
  e[target_edge, 1] <- new_parent_node
  
  # In $edge, place the moved clade next to its sisters
  # This is necessary to prevent some ape functions crashing.
  layer1 <- (1:max(new_sibling_edges)) %>% setdiff(descendant_edges)
  layer2 <- descendant_edges
  layer3 <- (1:n_edge)[-c(layer1, layer2, target_edge)]
  new_order <- c(layer1, target_edge, layer2, layer3)
  e <- e[new_order, ]
  phy$edge.length <- phy$edge.length[new_order]
  
  # Tidy up any childless nodes
  while (!(old_parent_node %in% e[, 1])) {
    remove_edge <- which(e[, 2] == old_parent_node)
    old_grandparent_node <- e[remove_edge, 1]
    # Remove edge and edge.length
    e <- e[-remove_edge, ]
    phy$edge.length <- phy$edge.length[-remove_edge]
    # Reindex nodes
    e <- e - (e > old_parent_node)
    old_grandparent_node <- 
      old_grandparent_node - (old_grandparent_node > old_parent_node)
    # Remove node
    phy$Nnode <- phy$Nnode - 1
    phy$node.label <- phy$node.label[-(old_parent_node - n_tip)]
    # Shift focus onto the grandparent
    old_parent_node <- old_grandparent_node
  }
  phy$edge <- e
  
  # Reorder edges and reassign indices
  .reindex(phy)
}


#' #' Safely bind tips
#' #'
#' #' \code{ape::bind.tip()} deletes tip and node label substrings enclosed in
#' #' square brackets. This is a safe version which doesn't.
#' #'
#' #' @noRd
#' .bind_tip = function(phy, tip.label, ...) {
#'   
#'   tip.label <- 
#'     tip.label %>%
#'     str_replace_all("\\[", "〔") %>%
#'     str_replace_all("\\]", "〕")
#'   
#'   phy$tip.label <- 
#'     phy$tip.label %>%
#'     str_replace_all("\\[", "〔") %>%
#'     str_replace_all("\\]", "〕")
#'   
#'   phy$node.label <- 
#'     phy$node.label %>%
#'     str_replace_all("\\[", "〔") %>%
#'     str_replace_all("\\]", "〕")
#'   
#'   phy <- bind.tip(phy, tip.label, ...)
#'   
#'   phy$tip.label <- 
#'     phy$tip.label %>%
#'     str_replace_all("〔", "\\[") %>%
#'     str_replace_all("〕", "\\]")
#'   
#'   phy$node.label <- 
#'     phy$node.label %>%
#'     str_replace_all("〔", "\\[") %>%
#'     str_replace_all("〕", "\\]")
#'   
#'   phy
#' }


#' Group cloned sisters in the tree ordering
#'
#' In Glottolog trees, nodes and tips are ordered by their glottocode.
#' clone_tip() places tips either in the right order or too early, as the first
#' sibling. This function groups sisters by moving early ones to a position
#' adjacent to their last-ordered sibling(s).
#'
#' @param phy A phylo object.
#' @return A phylo object, in which the edge orders potentially are changed.
#' @noRd
.group_cloned_sisters = function(phy) {
  
  # helper function
  move_row_fwd <- function(df, from, to) {
    slice(df, c((1:to)[-from], from, (to:nrow(df))[-1]))
  }
  
  labels <- str_c(phy$tip.label, phy$edge.label)
  edges <- phy$edge
  df <- data.frame(orig_row = 1:nrow(edges),
                   parent = edges[, 1],
                   child_label = labels[edges[, 2]])
  
  # Cycle through parents and the labels of their children
  for (p in unique(df$parent)) {
    clabels <- filter(df, parent == p, !is.na(child_label))$child_label
    for (l in unique(clabels)) {
      lp_row <- which(df$parent == p & df$child_label == l)
      n_lp_row = length(lp_row)
      # Skip ahead if there's only one match to this parent & label
      if (n_lp_row == 1) next
      first_lp_row <- lp_row[1]
      last_lp_row <- lp_row[n_lp_row]
      # Skip ahead if all matches are in adjacent rows
      if (last_lp_row == (first_lp_row + n_lp_row - 1)) next
      # At least one not adjacent, so cycle through, moving each to last
      # This is brute-force but it works
      for (i in (n_lp_row - 1):1) {
        df <- df %>% move_row_fwd(from = lp_row[i], to = last_lp_row)
      }
    }
  }
  
  # Change order of edges, leaving the numbering of nodes and tips unchanged
  phy$edge <- phy$edge[df$orig_row, ]
  phy$edge.length <- phy$edge.length[df$orig_row]
  phy
}


.reindex = function(phy) {
  
  # Reorder one way then another, because just applying cladewise reordering
  # can -- for whatever reason -- sometimes produce no change.
  # This ensures edges are in the desired, cladewise order
  phy <- phy %>% reorder(order = "pruning") %>% reorder(order = "cladewise")
  
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)
  e <- phy$edge
  pendants <- which(e[,2] <= n_tip)
  
  # Observe the current numbering sequence of nodes and tips.
  nd_seq <- c(e[1,1], e[-pendants, 2])
  t_seq <- e[pendants, 2]
  
  # Remapping of current order to the desired order, i.e., according to first
  # appearance in phy$edge
  nd_remap <- match((1:n_node) + n_tip, nd_seq)
  t_remap <- match(1:n_tip, t_seq)
  v_remap <- c(t_remap, nd_remap + n_tip)
  
  # Rebuild the tree
  new_tree <- list(
    edge = matrix(v_remap[e], ncol = 2),
    Nnode = n_node,
    node.label = phy$node.label[nd_seq - n_tip],
    tip.label = phy$tip.label[t_seq],
    edge.length = phy$edge.length
  )
  class(new_tree) <- "phylo"
  new_tree
}
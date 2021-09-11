# Setting branch lengths

#' Set length of deepest branches
#'
#' Sets lengths of branches immediately below the root to the same,
#' user-specified length.
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @param length A numeric stating the branch length.
#' @return A phlyo object, the manipulated tree.
#' @examples 
#' 
#' library(ape)
#' arnhem_hypothesis <- 
#'   c("Gunwinyguan", "Mangarrayi-Maran", "Maningrida",
#'     "Kungarakany", "Gaagudju")
#' tree <- assemble_rake(abridge_labels(get_glottolog_trees(arnhem_hypothesis)))
#' plot(tree)
#' nodelabels(tree$node.label)
#' # tree now contains five language families. All branch lengths are 1.
#' # Set the deepest branch lengths to 5, implying a great genealogical
#' # distance between the families within the tree.
#' tree2 <- set_deepest_branch_lengths(tree, length = 5)
#' plot(tree2)
#' nodelabels(tree2$node.label)
rescale_deepest_branches = function(
  phy, 
  length = 1
) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check length
  if (!is.numeric(length)) {
    cl <- class(length)
    stop(str_c("`length` must be numeric.\n",
               "You supplied an object of class ", cl, "."))
  }
  if (length(length) != 1) {
    stop(str_c("`length` must be a vector length 1.\n",
               "You supplied a vector length ", length(length), "."))
  }
  if (length < 0) {
    stop(str_c("`length` must be non-negative.\n",
               "You supplied the negative value ", length, "."))
  }
  
  root <- Ntip(phy) + 1
  first_edges <- which(phy$edge[,1] == root)
  phy$edge.length[first_edges] <- length
  phy
}


#' Set all branch length to 1
#' 
#' Sets all branch lengths in a tree to the same length.
#' 
#' By default, sets all branch lengths to 1.
#' 
#' @param phy A phlyo object, the tree to manipulate.
#' @param length A numeric stating the branch length.
#' @return A phlyo object, the manipulated tree.
#' @examples 
#' 
#' library(ape)
#' tree <- abridge_labels(get_glottolog_trees("Tangkic"))
#' tree2 <- clone_tip(tree, "nyan1300", n = 2, subgroup = TRUE)
#' plot(tree2)
#' nodelabels(tree2$node.label)
#' tree3 <- rescale_branches(tree2)
#' plot(tree3)
#' nodelabels(tree3$node.label)
rescale_branches = function(
  phy,
  length = 1
) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check length
  if (!is.numeric(length)) {
    cl <- class(length)
    stop(str_c("`length` must be numeric.\n",
               "You supplied an object of class ", cl, "."))
  }
  if (length(length) != 1) {
    stop(str_c("`length` must be a vector length 1.\n",
               "You supplied a vector length ", length(length), "."))
  }
  if (length < 0) {
    stop(str_c("`length` must be non-negative.\n",
               "You supplied the negative value ", length, "."))
  }
  
  phy$edge.length <- rep(length, Nedge(phy))
  phy
}


#' Exponentialize branch lengths
#'
#' Sets the deepest branches to length 1/2, the next deepest to 1/4, the next to
#' 1/8, etc., all multiplied by the parameter \code{length}.
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @param length A positive numeric, a multiplier for the exponential branch
#'   lengths 1/2, 1/4, 1/8...
#' @return A phlyo object, the manipulated tree.
#' @examples
#'
#' library(ape)
#' tree <- abridge_labels(get_glottolog_trees("Siouan"))
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- rescale_branches_exp(tree)
#' plot(tree2)
#' nodelabels(tree2$node.label)
rescale_branches_exp = function(
  phy,
  length = 1
) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check length
  if (!is.numeric(length)) {
    cl <- class(length)
    stop(str_c("`length` must be numeric.\n",
               "You supplied an object of class ", cl, "."))
  }
  if (length(length) != 1) {
    stop(str_c("`length` must be a vector length 1.\n",
               "You supplied a vector length ", length(length), "."))
  }
  if (length < 0) {
    stop(str_c("`length` must be non-negative.\n",
               "You supplied the negative value ", length, "."))
  }
  
  nonroot <- phy$edge[,2]
  
  # Get node depths for non-roots
  phy$edge.length <- rep(1, Nedge(phy))
  depth <- node.depth.edgelength(phy)[nonroot]
  
  # Assign new branch length according to depth
  phy$edge.length <- length * (0.5 ^ depth)
  phy
}


#' Ultrametricize tree by stretching final edges
#'
#' Alters branches ending in a tip in such a way that all tips are equidistant
#' from the root. Does this by lengthening branches above all but the existing,
#' most-distance tip(s).
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @return A phlyo object, the manipulated tree.
#' @examples 
#' 
#' library(ape)
#' tree <- rescale_branches_exp(abridge_labels(get_glottolog_trees("Siouan")))
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- ultrametricize(tree)
#' plot(tree2)
#' nodelabels(tree2$node.label)
ultrametricize = function(phy) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  if (any(is.na(phy$edge.length))) {
    warning("Converting NA branch lengths to 1")
    phy$edge.length[is.na(phy$edge.length)] <- 1
  }
  
  n_tips <- Ntip(phy)
  tip_depths <- node.depth.edgelength(phy)[1:n_tips]
  added_depth <- max(tip_depths) - tip_depths
  tip_edges <- match(1:n_tips, phy$edge[,2])
  phy$edge.length[tip_edges] <- 
    phy$edge.length[tip_edges] + added_depth
  
  phy
}
# Setting branch lengths

#' Set length of deepest branches
#'
#' Sets lengths of branches immediately below the root to the same,
#' user-specified length.
#'
#' @param phy A phylo object, the tree to manipulate.
#' @param length A numeric stating the branch length.
#' @return A phylo object, the manipulated tree.
#' @examples 
#' 
#' library(ape)
#' arnhem_hypothesis <- 
#'   c("Gunwinyguan", "Mangarrayi-Maran", "Maningrida",
#'     "Kungarakany", "Gaagudju")
#' tree <- assemble_rake(abridge_labels(get_glottolog_trees(arnhem_hypothesis)))
#' plot_glotto(tree)
#' # tree now contains five language families. All branch lengths are 1.
#' # Set the deepest branch lengths to 5, implying a great genealogical
#' # distance between the families within the tree.
#' tree2 <- rescale_deepest_branches(tree, length = 5)
#' plot_glotto(tree2)
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
#' @param phy A phylo object, the tree to manipulate.
#' @param length A numeric stating the branch length.
#' @return A phylo object, the manipulated tree.
#' @examples 
#' 
#' library(ape)
#' tree <- abridge_labels(get_glottolog_trees("Tangkic"))
#' tree2 <- clone_tip(tree, "nyan1300", n = 2, subgroup = TRUE)
#' plot_glotto(tree2)
#' tree3 <- rescale_branches(tree2)
#' plot_glotto(tree3)
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
#' @param phy A phylo object, the tree to manipulate.
#' @param length A positive numeric, a multiplier for the exponential branch
#'   lengths 1/2, 1/4, 1/8...
#' @return A phylo object, the manipulated tree.
#' @examples
#'
#' library(ape)
#' tree <- abridge_labels(get_glottolog_trees("Siouan"))
#' plot_glotto(tree)
#' tree2 <- rescale_branches_exp(tree)
#' plot_glotto(tree2)
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


#' Exponentialize branch lengths using constraints
#'
#' Sets the branches furthest from the root to \code{shortest_length} and then
#' scales other branches exponentially in such a way that the total tree has a
#' height of \code{tree_height}.
#'
#' @param phy A phylo object, the tree to manipulate.
#' @param tree_height A positive numeric, the desired height of the whole tree.
#' @param shortest_length A positive numeric, the desired length of the shortest
#'   branches, i.e., the branches farthest from the root.
#' @return A phylo object, the manipulated tree.
#' @examples
#'
#' library(ape)
#' tree <- abridge_labels(get_glottolog_trees("Siouan"))
#' plot_glotto(tree)
#' tree2 <- rescale_branches_constrained(tree, tree_height = 2000, shortest_length = 200)
#' plot_glotto(tree2)
#' tree3 <- rescale_branches_constrained(tree, tree_height = 2000, shortest_length = 100)
#' plot_glotto(tree3)
rescale_branches_constrained = function(
    phy,
    tree_height,
    shortest_length
) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check tree_height
  if (!is.numeric(tree_height)) {
    cl <- class(tree_height)
    stop(str_c("`tree_height` must be numeric.\n",
               "You supplied an object of class ", cl, "."))
  }
  if (length(tree_height) != 1) {
    stop(str_c("`shortest_length` must be a vector length 1.\n",
               "You supplied a vector length ", length(length), "."))
  }
  if (tree_height < 0) {
    stop(str_c("`shortest_length` must be non-negative.\n",
               "You supplied the negative value ", length, "."))
  }
  
  # Check shortest_length
  if (!is.numeric(shortest_length)) {
    cl <- class(shortest_length)
    stop(str_c("`shortest_length` must be numeric.\n",
               "You supplied an object of class ", cl, "."))
  }
  if (length(shortest_length) != 1) {
    stop(str_c("`shortest_length` must be a vector length 1.\n",
               "You supplied a vector length ", length(length), "."))
  }
  if (shortest_length < 0) {
    stop(str_c("`shortest_length` must be non-negative.\n",
               "You supplied the negative value ", length, "."))
  }
  
  nonroot <- phy$edge[,2]
  
  # Get node depths for non-roots
  phy$edge.length <- rep(1, Nedge(phy))
  depth <- node.depth.edgelength(phy)[nonroot]
  tree_depth <- max(depth)
  
  # Check if tree is has depth > 1
  if (tree_depth == 1) {
    stop(str_c("`phy` must be a tree with hierarchical branching. You ",
               "supplied a tree without hierarchical branching."))
  }
  
  # Check if the constraints are inconsistent
  if (shortest_length * tree_depth > tree_height) {
    stop(str_c("For a tree with n hierarchical levels of branches, the ",
               "shortest branch length must be no greater than 1/n times the ",
               "tree height. However, `phy` has ", tree_depth, " hierarchical ",
               "levels and you have specified `shortest_length` as ", 
               shortest_length, ", which is greater than 1/n times ",
               "`tree_height`, which you specified as ", tree_height, "."))
  }
  
  # Calculate the branch length ratio, by solving: 
  # (1-r^n) / (1-r) = height / shortest
  if (tree_depth == 2) {
    r <- (tree_height - shortest_length) / shortest_length
  } else {
    p <- tree_height / shortest_length
    coeffs <- c((p - 1), -p, rep(0, tree_depth - 2), 1)
    rts <- polyroot(coeffs)
    # Take the real solution which is not 1
    rt <- rts[which(Im(rts) < 1e-10 & Re(rts) > 1)]
    r <- Re(rt)
  }
  
  # Assign new branch length according to depth
  phy$edge.length <- shortest_length * (r ^ (tree_depth - depth))
  phy
}


#' Ultrametricize tree by stretching final edges
#'
#' Alters branches ending in a tip in such a way that all tips are equidistant
#' from the root. Does this by lengthening branches above all but the existing,
#' most-distance tip(s).
#'
#' @param phy A phylo object, the tree to manipulate.
#' @return A phylo object, the manipulated tree.
#' @examples 
#' 
#' library(ape)
#' tree <- rescale_branches_exp(abridge_labels(get_glottolog_trees("Siouan")))
#' plot_glotto(tree)
#' tree2 <- ultrametricize(tree)
#' plot_glotto(tree2)
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
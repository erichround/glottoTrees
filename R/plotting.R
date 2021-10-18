# Plotting

#' Plot a tree in downward, linguistic style
#'
#' An attempt is made to choose a reasonable plot width, height and offset of
#' the labels from the tree's tips. If the choices are not satisfactory,
#' manipulate them using positive or negative values of \code{extra_width},
#' \code{extra_height} and \code{extra_offset}.
#' 
#' By default, tip labels will run vertically, except in the case that all tip
#' labels are a maximum of one character long, in which case they are oriented
#' horizontally.
#'
#' @param phy A phylo object.
#' @param nodelabels A logical, whether to plot node labels.
#' @param extra_width A numeric, extra width to add beyond the default.
#' @param extra_height A numeric, extra height to add beyond the default.
#' @param extra_offset A numeric, extra offset to add beyond the default.
#' @param srt A numeric giving how much the labels are rotated in degrees.
#' @param adj A numeric specifying the justification of the text strings of the
#'   labels: 0 (left-justification), 0.5 (centering), or 1
#'   (right-justification).
#' @param ... Additional arguments passed to plot().
#' @examples
#'
#' library(ape)
#'
#' tree <- abridge_labels(get_glottolog_trees("Tangkic"))
#' plot_glotto(tree)
#' tree2 <- rescale_deepest_branches(tree, 7)
#' plot_glotto(tree2)
#' plot_glotto(tree2, srt = 90, adj = 0.5)
plot_glotto = function(phy, nodelabels = TRUE,
                       extra_width = 0, extra_height = 0,
                       extra_offset = 0, srt = 0, adj = NULL,
                       ...) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  has_nodelabels <- "node.label" %in% names(phy)
  has_lengths <- "edge.length" %in% names(phy)
  max_tip_len <- max(str_length(phy$tip.label))
  height_mult <- 1.01 + 0.025 * max_tip_len
  
  default_width <- Ntip(phy) + 1

  if (has_lengths) {
    default_height <- max(node.depth.edgelength(phy)) * height_mult
  } else {
    default_height <- (max(node.depth(phy, method = 2)) - 1) * height_mult
  }
  
  if (max_tip_len == 1 & srt == 0 & is.null(adj)) {
    # All tip labels are single-character. Rotate them to run horizontally.
    srt <- 90
    adj <- 0.5
  }
  
  default_offset <- default_height * 0.01
  
  plot(phy, direction = "downwards", no.margin = TRUE,
       label.offset = default_offset + extra_offset, 
       x.lim = c(0, default_width + extra_width),
       y.lim = c(0, default_height + extra_height),
       srt = srt, adj= adj, ...)
  if (nodelabels & has_nodelabels) { nodelabels(text = phy$node.label) }
}
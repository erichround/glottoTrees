% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/branch_lengths.R
\name{rescale_deepest_branches}
\alias{rescale_deepest_branches}
\title{Set length of deepest branches}
\usage{
rescale_deepest_branches(phy, length = 1)
}
\arguments{
\item{phy}{A phylo object, the tree to manipulate.}

\item{length}{A numeric stating the branch length.}
}
\value{
A phylo object, the manipulated tree.
}
\description{
Sets lengths of branches immediately below the root to the same,
user-specified length.
}
\examples{

library(ape)
arnhem_hypothesis <- 
  c("Gunwinyguan", "Mangarrayi-Maran", "Maningrida",
    "Kungarakany", "Gaagudju")
tree <- assemble_rake(abridge_labels(get_glottolog_trees(arnhem_hypothesis)))
plot_glotto(tree)
# tree now contains five language families. All branch lengths are 1.
# Set the deepest branch lengths to 5, implying a great genealogical
# distance between the families within the tree.
tree2 <- rescale_deepest_branches(tree, length = 5)
plot_glotto(tree2)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topology.R
\name{keep_tip}
\alias{keep_tip}
\title{Keep tips}
\usage{
keep_tip(phy, label)
}
\arguments{
\item{phy}{A phylo object. The tree to manipulate.}

\item{label}{A character vector containing tip labels.}
}
\value{
A phylo object containing the modified tree.
}
\description{
Nominate which tips are to be kept. Others are removed.
}
\examples{

library(ape)

tree <- abridge_labels(get_glottolog_trees("Tangkic"))
plot_glotto(tree)

tree2 <- keep_tip(tree, c("lard1243", "kang1283", "kaya1319"))
plot_glotto(tree2)
}

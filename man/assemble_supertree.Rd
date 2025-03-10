% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topology.R
\name{assemble_supertree}
\alias{assemble_supertree}
\title{Create a glottolog super-tree}
\usage{
assemble_supertree(macro_groups, glottolog_version)
}
\arguments{
\item{macro_groups}{A list of character vectors, in which each vector
contains the names of one or more macroareas which define a group.
Alternatively, setting \code{macro_groups} to \code{NULL} causes the tree
to be assembled without groups.}

\item{glottolog_version}{A character string, specifying which glottolog
version to use. Currently available options are \code{'4.0'} through to
\code{'4.8'} and \code{'5.0'}. If no value is specified then the newest
available version is used.}
}
\description{
Combining glottolog family trees into one large tree. Families can be
assembled directly below a rake structure at the root, or can be grouped, so
that the root first branches into groups, and the families then branch out
below the group nodes.
}
\details{
Grouping is controlled by the \code{macro_groups} parameter. Groups can
comprise a single glottolog macroarea, or multiple macroareas. Current
macroareas are \code{Africa}, \code{Australia}, \code{Eurasia}, \code{North
America}, \code{Papunesia} and \code{South America}. Setting
\code{macro_groups} to \code{NULL} causes the tree to be assembled without
groups.
}
\examples{

# Supertree whose first order branches are the glottolog macroareas
supertree <- assemble_supertree()
supertree_v.4.3 <- assemble_supertree(glottolog_version = "4.3")

# Supertree whose first order branches are glottolog families
supertree <- assemble_supertree(macro_groups = NULL)

# Supertree whose first order branches are the African & Eurasian macroareas
supertree <- assemble_supertree(macro_groups = list("Africa", "Eurasia"))

# Supertree whose first order branches are the glottolog macroareas, but
# with the Americas combined:
supertree <- assemble_supertree(
  macro_groups = list("Africa", "Australia", "Eurasia", "Papunesia", 
                      c("South America", "North America"))
  )
}

#' Trees from glottolog, versions 4.0 - 4.4
#'
#' multiPhylo objects, which provide a representation of the phylogenetic
#' relationships of the languages in glottolog.
#'
#' @source \url{https://glottolog.org/meta/downloads}
#' @name glottolog_trees
NULL

#' @rdname glottolog_trees
#' @format v4.0: A multiPhylo object of 421 trees
"glottolog_trees_v4.0"

#' @rdname glottolog_trees
#' @format v4.1: A multiPhylo object of 421 trees
"glottolog_trees_v4.1"

#' @rdname glottolog_trees
#' @format v4.2: A multiPhylo object of 422 trees
"glottolog_trees_v4.2"

#' @rdname glottolog_trees
#' @format v4.3: A multiPhylo object of 418 trees
"glottolog_trees_v4.3"

#' @rdname glottolog_trees
#' @format v4.4: A multiPhylo object of 420 trees
"glottolog_trees_v4.4"


#' Geographical data from glottolog, versions 4.0 - 4.4
#'
#' Datasets of geographical information about the languages in glottolog.
#'
#' @format A dataframe:
#' \describe{
#'   \item{glottocode}{glottocode of the lect}
#'   \item{name}{name of the lect}
#'   \item{isocodes}{the ISO-639-3 code of the lect}
#'   \item{level}{"language" or "dialect"}
#'   \item{macroarea}{glottolog's geographical macroarea}
#'   \item{latitude}{}
#'   \item{longitude}{}
#' }
#' @source \url{https://glottolog.org/meta/downloads}
#' @name glottolog_geography
NULL

#' @rdname glottolog_geography
#' @format v4.0: A dataframe of 20,049 rows
"glottolog_geography_v4.0"

#' @rdname glottolog_geography
#' @format v4.1: A dataframe of 20,290 rows
"glottolog_geography_v4.1"

#' @rdname glottolog_geography
#' @format v4.2: A dataframe of 20,752 rows
"glottolog_geography_v4.2"

#' @rdname glottolog_geography
#' @format v4.3: A dataframe of 20,930 rows
"glottolog_geography_v4.3"

#' @rdname glottolog_geography
#' @format v4.4: A dataframe of 21,329 rows
"glottolog_geography_v4.4"
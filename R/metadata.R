# Glottolog metadata


#' Simple language metadata
#'
#' Returns a dataframe of metadata on glottolog's languages.
#'
#' Returned columns are: \code{glottocode}, \code{isocodes}, \code{name},
#' \code{name_in_tree}, \code{position}, \code{tree} and \code{tree_name}.
#'
#' @param glottolog_version A character string, specifying which glottolog
#'   version to use. Currently available options are \code{'4.3'} and
#'   \code{'4.4'}. If no value is specified then the newest available version is
#'   used.
#' @examples
#' head(get_glottolog_languages())
#' head(get_glottolog_languages(glottolog_version = "4.3"))
get_glottolog_languages = function(
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
  
  get_glottolog_phylo_geo(glottolog_version) %>%
    select(
      glottocode,
      isocodes,
      name,
      name_in_tree = vertex_name,
      position = vertex_type,
      tree = tree,
      tree_name = family_name
    )
}


#' Simple family metadata
#'
#' Returns a dataframe of metadata on glottolog's language families.
#'
#' Returned columns are: \code{tree}, \code{tree_name}, \code{n_tips},
#' \code{n_nodes} and \code{main_macroarea}.
#'
#' @inheritParams get_glottolog_languages
#' @examples
#' head(get_glottolog_families())
#' head(get_glottolog_families(glottolog_version = "4.3"))
get_glottolog_families = function(
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
  
  main_macro <-
    get_glottolog_phylo_geo(glottolog_version) %>%
    filter(!is.na(tree)) %>%
    group_by(tree, macroarea) %>%
    summarise(n = n()) %>%
    group_by(tree) %>%
    arrange(-n) %>%
    slice(1) %>%
    select(tree, main_macroarea = macroarea)
  
  get_glottolog_phylo_geo(glottolog_version) %>%
    select(
      tree = tree,
      tree_name = family_name,
      vertex_type
    ) %>%
    filter(!is.na(tree)) %>%
    group_by(tree, tree_name) %>%
    summarise(
      n_tips = sum(vertex_type == "tip"),
      n_nodes = sum(vertex_type == "node")
      ) %>%
    arrange(tree) %>%
    left_join(main_macro, by = "tree") %>%
    as.data.frame()
}


#' Extended glottolog metadata
#'
#' Returns a dataframe of glottolog geographical and phylogenetic metadata.
#' 
#' Returned columns are: \code{glottocode}, \code{isocodes}, \code{name}, \code{level},
#' \code{vertex_type}, \code{vertex_label}, \code{vertex_name},
#' \code{macroarea}, \code{latitude}, \code{longitude},
#' \code{family_glottocode}, \code{family_name} and \code{tree}.
#'
#' @inheritParams get_glottolog_languages
#' @examples
#' head(get_glottolog_phylo_geo())
#' head(get_glottolog_phylo_geo(glottolog_version = "4.3"))
get_glottolog_phylo_geo = function(
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
  if (glottolog_version == "4.3") {
    phylo_geo <- glottolog_phylo_geo_v4.3
  } else if (glottolog_version == "4.4") {
    phylo_geo <- glottolog_phylo_geo_v4.4
  }
  
  phylo_geo
}


#' Tree numbers of glottolog families
#' 
#' Returns the tree number of one or more glottolog families.
#' 
#' @inheritParams get_glottolog_trees
#' @return A named vector of integers, giving the tree numbers and the
#'   family names as the vector names.
#' @examples
#' which_tree("Caddoan")
#' which_tree(c("Caddoan", "Tangkic"), glottolog_version = "4.3")
#' # If some family names are unrecognized, a warning is issued
#' which_tree(c("Caddoan", "Zzz"), glottolog_version = "4.4")
#' \dontrun{
#' # If no family names are recognized, an error results
#' which_tree()
#' which_tree("Zzz")
#' }
which_tree = function(
  family,
  glottolog_version
) {
  
  if (missing(family)) {
    stop(str_c("A value for `family` needs to be supplied.\n",
               "You didn't supply one."))
  }
  
  # Check glottolog_version
  if (missing(glottolog_version)) {
    glottolog_version <- .get_newest_version()
  } else {
    error_msg <- .check_glottolog_version(glottolog_version)
    if (!is.na(error_msg)) { stop(error_msg) }
    glottolog_version <- as.character(glottolog_version)
  }
  
  # Check family
  check_result <- 
    .check_family(family, glottolog_version)
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  f <- get_glottolog_families(glottolog_version)
  tree_nums <- f$tree[match(family, f$tree_name)]
  names(tree_nums) <- family
  tree_nums
}


#' The current glottolog version
#'
#' Returns the newest version of glottolog for which data is included in this
#' package, which is \code{'4.4'}.
#'
#' @return A character string.
#' @noRd
.get_newest_version = function() {
  "4.4"
}


#' Vector of glottolog families
#' 
#' @inheritParams get_glottolog_languages
#' @return A vector.
#' @noRd
.get_family_vector = function(
  glottolog_version
) {
  get_glottolog_phylo_geo(glottolog_version) %>%
    .$family_name %>%
    unique()
}

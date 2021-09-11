# Input checking


#' Check glottolog versions
#' 
#' Check the glottolog_version parameter 
#' and return informative error messages
#' 
#' @param glottolog_version An object to be 
#'   checked. Should be a character vector,
#'   (or numeric or factor) of length 1,
#'   equal to '4.3' or '4.4' (or 4.3 or 4.4).
#' @return An error message, or else NA.
#' @noRd
.check_glottolog_version = function(glottolog_version) {
  
  available_versions <- c("4.3", "4.4")
  
  if (!is.character(glottolog_version)) {
    if (is.numeric(glottolog_version) | is.factor(glottolog_version)) {
      glottolog_version <- as.character(glottolog_version)
    } else {
      gc <- class(glottolog_version)
      return(str_c("`glottolog_version` must be a character string.\n",
                   "You have supplied an object of class ", gc, "."))
    }
  }
  
  if (length(glottolog_version) != 1) {
    return(str_c("`glottolog_version` must have length 1.\n",
                 "You have supplied a vector of length ", 
                 length(glottolog_version), "."))
  }
  
  if (!(glottolog_version %in% available_versions)) {
    return(str_c("`glottolog_version` must be one the following values: '",
                 str_c(available_versions, collapse = "', '"), "'.\n",
                 "You have supplied ", glottolog_version, "."))
  }
  
  NA
}


#' Check families
#'
#' Check non-NULL values of the family parameter and return informative error
#' and warning messages.
#'
#' @param family Should be a character vector whose elements are glottolog
#'   family names.
#' @param glottolog_version Should be a character string.
#' @return A list with elements error_msg and warning_msg. These are NA if
#'   there's no message.
#' @noRd
.check_family = function(
  family,
  glottolog_version,
  check_misplaced_gversion = TRUE
) {
  
  # Assumes glottolog_version is okay!
  fvec <- .get_family_vector(glottolog_version)
  is_good <- (family %in% fvec)
  n_good <- sum(is_good)
  n_bad <- sum(!is_good)
  
  # Check family
  if (is.factor(family)) {
    family <- as.character(family)
  }
  if (!is.character(family)) {
    fc <- class(family)
    return(list(
      error_msg = str_c("`family` must be a character vector.\n",
                        "You have supplied on object of class ", fc, "."),
      warning_msg = NA
    ))
  }
  if (length(family) == 0) {
    return(list(
      error_msg = str_c("`family` be length 1 or more.\n",
                        "You have supplied a value of length 0."),
      warning_msg = NA
    ))
  }
  if (check_misplaced_gversion &
      length(family) == 1 &&
      str_detect(family, "^[0-9]+\\.[0-9]+$")) {
    return(list(
      error_msg = str_c("`family` must contain at least one glottolog family name.\n",
                        "You supplied the value '", family, "'.\n",
                        "Did you mean `glottolog_version = '", family, "'`?"),
      warning_msg = NA
    ))
  }
  if (n_bad > 0) {
    bad_family <- family[!is_good]
    if (n_bad > 4) { bad_family <- c(bad_family[1:4], "..") }
    if (n_good == 0) {
      return(list(
        error_msg =str_c("`family` must contain at least one glottolog family name.\n",
                         "The value(s) you supplied are not glottolog family names.\n",
                         "Did you use the right capitalization and punctuation?\n"),
        warning_msg = NA
        ))
    } else {
      return(list(
        error_msg = NA,
        warning_msg = str_c("Elements of `family` should be glottolog family names.\n",
                            "You supplied one or more values which are not: ",
                            str_c(bad_family, collapse = ", "), ".\n",
                            "Did you use the right capitalization and punctuation?\n")
        ))   
    }
  }
  list(error_msg = NA, warning_msg = NA)
}


#' Check labels
#'
#' Check the label parameter and return informative error and warning messages
#'
#' @param phy A phylo object.
#' @param label An object to be checked. Should be a character vector of length
#'   > 0, comprised of elements also found in phy$tip.label and/or
#'   phy$node.label, with no duplicates
#' @param type A string. Either 'tip', 'node' or 'both'.
#' @return A list with elements `error_msg` and `warning_msg`. These are `NA`` if
#'   there's no message.
#' @noRd
.check_labels = function(phy, label, type) {
  
  if (!is.character(label)) {
    cl <- class(label)
    return(list(
      error_msg = str_c("`label` must be a character vector.\n",
                        "You supplied an object of class ", cl, "."),
      warning_msg = NA
    ))
  }
  
  if (length(label) == 0) {
    return(list(
      error_msg = str_c("`label` must be length 1 or more.\n",
                        "You supplied a vector length 0."),
      warning_msg = NA
    ))
  }
  
  if (any(label == "")) {
    return(list(
      error_msg = str_c("Elements of `label` cannot be an empty string.\n",
                        "You provided an empty string as element ", 
                        which(label == "")[1], " of `label`."),
      warning_msg = NA
    ))
  }
  
  # This assumes phy is okay!
  param <- "`label`"
  if (type == "tip") {
    extra_label <- setdiff(label, phy$tip.label)
    type_str <- "tip"
  } else if (type == "node") {
    extra_label <- setdiff(label, phy$node.label)
    type_str <- "node"
  } else if (type == "parent") {
    extra_label <- setdiff(label, phy$node.label)
    type_str <- "node"
    param <- "`parent_label`"
  } else if (type == "both") {
    extra_label <- setdiff(label, c(phy$tip.label, phy$node.label))
    type_str <- "tip and/or node"
  }
  n_extra <- length(extra_label)
  if (n_extra != 0) {
    return(list(
      error_msg = 
        str_c("Elements of ", param, " should match ", type_str, 
              " labels in `phy`.\n",
              "In the values you supplied, there are no matches for: ",
              str_c(head(extra_label, 4), collapse = ", "),
              ifelse(n_extra > 4, "..", ""), "."),
      warning_msg = NA
    ))
  }
  is_dupl <- duplicated(label)
  dupl <- unique(label[is_dupl])
  if (any(is_dupl)) {
    return(list(
      error_msg = NA,
      warning_msg =
        str_c(param, " contained duplicate entries for: ",
              str_c(head(dupl, 4), collapse = ", "),
              ifelse(length(dupl) > 4, "..", ""), ".\n",
              "These were treated as if just one copy had been provided.")
    ))
  }
  
  list(error_msg = NA, warning_msg = NA)
}


#' Check phylo
#'
#' Checks if object phy has is class phylo. If phy is missing tip.label
#' and/or node.label, then these are added, and the labels are \code{''}.
#'
#' @param phy A phlyo object
#' @return A list with elements error_msg and warning_msg (which are NA is
#'   there's no message) and phy, a possibly modified version of phy
#' @noRd
.check_phy = function(phy) {
  
  # Check class
  if (class(phy) != "phylo") {
    cp <- class(phy)
    return(list(
      error_msg = str_c("`phy` must be of class phylo.\n",
                        "You supplied an object of class ", cp, "."),
      warning_msg = NA,
      phy = phy
      ))
  }
  
  # Check for node labels. Add "" if missing.
  if (!("node.label") %in% names(phy)) {
    phy$node.label <- rep("", phy$Nnode)
  }
  
  # Check for tip labels. Add "" if missing.
  n_tip <- length(setdiff(phy$edge[,2], phy$edge[,1]))
  if (!("tip.label") %in% names(phy)) {
    phy$tip.label <- rep("", n_tip)
  }
  
  # Check for edge lengths. Add 1 if missing.
  n_edge <- nrow(phy$edge)
  if (!("edge.length") %in% names(phy)) {
    phy$edge.length <- rep(1, n_edge)
  }
  
  return(list(error_msg = NA, warning_msg = NA, phy = phy))
}
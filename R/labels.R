# Labeling

#' Change tip labels
#' 
#' Replaces existing tip labels with new labels.
#' 
#' @param phy A phylo object.
#' @param label A vector of strings, the tip labels to be replaced.
#' @param new_label A vector of strings, the same length as \code{label}. The
#'   corresponding replacement labels.
#' @examples
#' 
#' library(ape)
#' 
#' tree <- relabel_with_names(get_glottolog_trees("Tangkic"))
#' plot(tree)
#' tree2 <- relabel_tip(tree, "Kangkalita", "Gangalidda")
#' plot(tree2)
#' 
#' tree3 <- clone_tip(tree, "Yangkaal", subgroup = TRUE)
#' plot(tree3)
#' tree4 <- apply_duplicate_suffixes(tree3)
#' plot(tree4)
#' tree5 <- relabel_tip(tree4, c("Kangkalita", "Yangkaal-1", "Yangkaal-2"), 
#'                      c("Gangalidda", "Yangkaal", "Yangarella"))
#' plot(tree5)                
relabel_tip = function(phy, label, new_label) {
  
  # Check phy
  if (class(phy) != "phylo" & 
      class(phy) != "multiPhylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo or multiPhylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check label
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  
  # Check new_label
  if (length(label) != length(new_label)) {
    stop(str_c(
      "`label` and `new_label` must contain the same number of elements.\n",
      "You provided `label` with ", length(label),
      " and `new_label` with ", length(new_label), "."
      ))
  }
  
  label_indices <- match(label, phy$tip.label)
  phy$tip.label[label_indices] <- new_label
  
  phy
}


#' Change node labels
#' 
#' Replaces existing tip labels with new labels.
#' 
#' @param phy A phylo object.
#' @param label A vector of strings, the node labels to be replaced.
#' @param new_label A vector of strings, the same length as \code{label}. The
#'   corresponding replacement labels.
#' @examples
#' 
#' library(ape)
#' 
#' tree <- relabel_with_names(get_glottolog_trees("Tangkic"))
#' tree2 <- apply_duplicate_suffixes(clone_tip(tree, "Yangkaal", subgroup = TRUE))
#' plot(tree2)
#' tree3 <- relabel_tip(tree2, c("Kangkalita", "Yangkaal-1", "Yangkaal-2"), 
#'                      c("Gangalidda", "Yangkaal", "Yangarella"))
#' plot(tree3)
#' nodelabels(tree3$node.label)
#' tree4 <- relabel_node(tree3, c("Yangkaal", "Kayardild-Yangkaal"),
#'                       c("Yangkaralda", "Kayardild-Yangkaralda"))
#' plot(tree4)
#' nodelabels(tree4$node.label)
relabel_node = function(phy, label, new_label) {
  
  # Check phy
  if (class(phy) != "phylo" & 
      class(phy) != "multiPhylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo or multiPhylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check label
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg)
  }
  
  # Check new_label
  if (length(label) != length(new_label)) {
    stop(str_c(
      "`label` and `new_label` must contain the same number of elements.\n",
      "You provided `label` with ", length(label),
      " and `new_label` with ", length(new_label), "."
    ))
  }
  
  label_indices <- match(label, phy$node.label)
  phy$node.label[label_indices] <- new_label
  
  phy
}


#' Shorten labels to a glottocode substring
#'
#' Shortens tip and node labels to a glottocode substring within them. If any
#' labels lack a glottocode substring, a warning is given.
#'
#' Glottocodes comprise four lowercase letters (or b10b or 3adt) followed by
#' four numbers, and are only identified if they are initial in the string or
#' are preceded by [.
#' 
#' Also recognizes and retains duplicate suffixes, i.e., a hyphen followed by
#' one or more numerals at the end of the string (see
#' \link{apply_duplicate_suffixes}).
#'
#' @param phy A phylo or multiPhylo object, containing one or more trees to
#'   manipulate.
#' @return A phylo or multiPhylo object, the manipulated tree(s).
#' @examples 
#' 
#' library(ape)
#' 
#' tree <- get_glottolog_trees("Koreanic")
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- abridge_labels(tree)
#' plot(tree2)
#' nodelabels(tree2$node.label)
#' 
#' # Retain duplicate suffixes:
#' tree3 <- clone_tip(tree, "Jollado<chol1278>", n = 2, subgroup = TRUE)
#' tree3a <- apply_duplicate_suffixes(tree3)
#' plot(tree3a)
#' nodelabels(tree3a$node.label)
#' tree4 <- abridge_labels(tree3a)
#' plot(tree4)
#' nodelabels(tree4$node.label)
#' 
#' # A warning is issued if any label does not contain a glottocode
#' supertree <- assemble_supertree() # contains nodes without glottocodes
#' supertree2 <- abridge_labels(assemble_supertree())
#' 
#' # Applied to a `multiPhylo` object:
#' trees <- get_glottolog_trees(c("Kartvelian", "Basque"))
#' trees2 <- abridge_labels(trees)
#' plot(trees[[1]])
#' nodelabels(trees[[1]]$node.label)
#' plot(trees2[[1]])
#' nodelabels(trees2[[1]]$node.label)
#' plot(trees[[2]])
#' nodelabels(trees[[2]]$node.label)
#' plot(trees2[[2]])
#' nodelabels(trees2[[2]]$node.label) 
abridge_labels = function(phy) {
  
  # Check phy
  if (class(phy) != "phylo" & 
      class(phy) != "multiPhylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo or multiPhylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  tip_orig <- .get_tip_labels(phy)
  node_orig <- .get_node_labels(phy)
  
  is_multi <- (class(phy) == "multiPhylo")
  if (!is_multi) { phy <- c(phy) }
  
  tip_fail <- node_fail <- logical(0)
  
  for (i in 1:length(phy)) {
    p <- phy[[i]]
    tip_glotto  <- extract_glottocode(p$tip.label)
    node_glotto <- extract_glottocode(p$node.label)
    tip_matched  <- !is.na(tip_glotto)
    node_matched <- !is.na(node_glotto)
    phy[[i]]$tip.label[tip_matched]  <- tip_glotto[tip_matched]
    phy[[i]]$node.label[node_matched] <- node_glotto[node_matched]
    tip_fail <- c(tip_fail, !tip_matched)
    node_fail <- c(node_fail, !node_matched)
  }
  
  if (!is_multi) { phy <- phy[[1]] }
  
  if (any(c(tip_fail, node_fail))) {
    warning(
      str_c("Labels without glottocodes were detected and left unchanged for: ",
            sum(tip_fail), " tip(s)",
            ifelse(!any(tip_fail), "", 
                   str_c(": ",
                         str_c(str_replace_na(head(tip_orig[tip_fail])),
                               collapse = ", "))),
            ifelse(sum(tip_fail) > 6, "...", ""),
            "; ", sum(node_fail), " node(s)",
            ifelse(!any(node_fail), "", 
                   str_c(": ",
                         str_c(str_replace_na(head(node_orig[node_fail])),
                               collapse = ", "))),
            ifelse(sum(node_fail) > 6, "...", "")
      ))    
  }
  
  phy
}


#' Change labels from glottocodes to names
#'
#' Looks up glottolog language names corresponding to glottocodes and replaces
#' tip and node labels, which contain glottocodes, with the appropriate names.
#'
#' Also recognizes and retains duplicate suffixes, i.e., a hyphen followed by
#' one or more numerals at the end of the string (see
#' \link{apply_duplicate_suffixes}).
#'
#' Labels without glottocodes are left unchanged and a warning is given. The
#' version of glottolog to use for look-up can be controlled with
#' \code{glottolog_version}.
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @inheritParams get_glottolog_languages
#' @return A phlyo object, the manipulated tree.
#'
#' @examples
#'
#' library(ape)
#'
#' # Replace full glottolog labels with names
#' tree <- get_glottolog_trees("Kresh-Aja")
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- relabel_with_names(tree)
#' plot(tree2)
#' nodelabels(tree2$node.label)
#'
#' # Replace abridged labels with names
#' tree3 <- abridge_labels(tree)
#' plot(tree3)
#' nodelabels(tree3$node.label)
#' tree4 <- relabel_with_names(tree3)
#' plot(tree4)
#' nodelabels(tree4$node.label)
#'
#' # Use names from earlier glottolog version:
#' tree5 <- relabel_with_names(tree, glottolog_version = "4.3")
#' plot(tree5)
#' nodelabels(tree5$node.label)
relabel_with_names = function(
  phy,
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
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check glottolog_version
  error_msg <- .check_glottolog_version(glottolog_version)
  if (!is.na(error_msg)) { stop(error_msg) }
  
  tip_labels <- .get_tip_labels(phy)
  node_labels <- .get_node_labels(phy)
  labels <- c(tip_labels, node_labels)
  labels_df <-
    data.frame(
      label = labels,
      vertex = c(rep("tip", length(tip_labels)),
                 rep("node", length(node_labels))),
      stringsAsFactors = FALSE
    ) %>%
    mutate(
      glottocode = extract_glottocode(labels),
      copy_suffix = .extract_copy_suffix(glottocode),
      glotto_base = .remove_copy_suffix(glottocode)
    ) %>%
    left_join(
      get_glottolog_phylo_geo(glottolog_version) %>%
        select(glotto_base = glottocode, 
               name_base = vertex_name),
      by = "glotto_base"
    ) %>%
    mutate(name = str_c(name_base, .hyphenate(copy_suffix)))
  
  missing_tips <- 
    filter(labels_df, is.na(name) & vertex == "tip")$label
  missing_nodes <- 
    filter(labels_df, is.na(name) & vertex == "node")$label
  n_missing_tips <- length(missing_tips)
  n_missing_nodes <- length(missing_nodes)
  
  # If a name is missing, give warning and retain the glottocode
  if (n_missing_tips + n_missing_nodes > 0) {
    warning(
      str_c("Labels without glottocodes were detected and left unchanged for: ",
            n_missing_tips, " tip(s)",
            ifelse(n_missing_tips == 0, "", 
                   str_c(": ",
                         str_c(str_replace_na(head(missing_tips)),
                               collapse = ", "))),
            ifelse(n_missing_tips > 6, "...", ""),
            "; ", n_missing_nodes, " node(s)",
            ifelse(n_missing_nodes == 0, "", 
                   str_c(": ",
                         str_c(str_replace_na(head(missing_nodes)),
                               collapse = ", "))),
            ifelse(n_missing_nodes > 6, "...", "")
      ))   
    labels_df <- labels_df %>%
      mutate(name = ifelse(is.na(name), label, name))
  }
  
  is_multi <- (class(phy) == "multiPhylo")
  if (!is_multi) { phy <- c(phy) }
  for (i in 1:length(phy)) {
    tip_matches <- match(phy[[i]]$tip.label, labels_df$label)
    node_matches <- match(phy[[i]]$node.label, labels_df$label)
    phy[[i]]$tip.label <- labels_df$name[tip_matches]
    phy[[i]]$node.label <- labels_df$name[node_matches]
  }
  if (!is_multi) { phy <- phy[[1]] }
  
  phy
}


#' Extract glottocode substrings
#'
#' From a character vector, extracts the first glottocode from each element.
#'
#' Glottocodes comprise four lowercase letters (or b10b or 3adt) followed by
#' four numbers, and are only identified if they are initial in the string or
#' are preceded by <.
#'
#' Also recognizes and retains duplicate suffixes, i.e., a hyphen followed by
#' one or more numerals at the end of the string (see
#' \link{apply_duplicate_suffixes}).
#'
#' @param label A string
#' @return A string
#'
#' @examples
#'
#' extract_glottocode("DongoKresh<dong1296>-l-")
#' extract_glottocode(c("DongoKresh<dong1296>-l-", "Goro-Golo<orlo1238>"))
#'
#' # Duplicate suffixes are recognised and retained
#' extract_glottocode(c("Goro-Golo<orlo1238>-1", "Goro-Golo<orlo1238>-2"))
extract_glottocode = function(label) {
  # note: two glottocodes contain numbers in the initial
  #       four characters: b10b and 3adt
  regex <- "(?<=(^|<))([a-z]{4}|b10b|3adt)[0-9]{4}(?=(>|-[0-9]{1,3}$|$))"
  g <- str_extract(label, regex)
  copy_suffix <- .extract_copy_suffix(label)
  str_c(g, .hyphenate(copy_suffix))
}


#' Apply duplicate suffixes to tips and nodes
#'
#' Suffixes are applied to ensure tip labels and node labels are not duplicates.
#' Suffixes have the form -1, -2, -3, ...
#'
#' The function recognizes existing duplicate suffixes and deals with them in
#' one of two ways. If a label has n duplicates that are already suffixed -1,
#' -2, ... -n, then the suffixes are not changed. Under any other conditions,
#' old suffixes are removed and new ones applied.
#'
#' Suffixation of tips and of nodes are handled independently of one another.
#'
#' @param phy A phylo object, the tree whose labels are to have copy suffixes
#'   applied.
#' @return A phylo object, the same tree but with suffixes applied to the
#'   labels.
#' @examples 
#' 
#' library(ape)
#' 
#' tree <- abridge_labels(get_glottolog_trees("Koreanic"))
#' plot(tree)
#' nodelabels(tree$node.label)
#' tree2 <- clone_tip(tree, "chol1278", n = 2, subgroup = TRUE)
#' plot(tree2)
#' nodelabels(tree2$node.label)
#' # Technically, tree2 is ill-formed because it has duplicate tip labels.
#' # Note how this causes problems if we try to clone one of them, since it
#' # is unclear which should be cloned:
#' \dontrun{
#' tree2a <- clone_tip(tree2, "chol1278")
#' }
#' 
#' # Suffixation of duplicate tips
#' tree3 <- apply_duplicate_suffixes(tree2)
#' plot(tree3)
#' nodelabels(tree3$node.label)
#' 
#' # Once they are suffixed, these tips can be cloned successfully:
#' tree4 <- clone_tip(tree3, c("chol1278-2", "chol1278-3"), subgroup = TRUE)
#' plot(tree4)
#' nodelabels(tree4$node.label)
#' # Suffixing is applied across all tips that share a glottocode, and
#' # separately, across all nodes that share a glottocode:
#' tree5 <- apply_duplicate_suffixes(tree4)
#' plot(tree5)
#' nodelabels(tree5$node.label)
apply_duplicate_suffixes = function(
  phy
) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing phy$node.label, phy$tip.label
  phy <- check$phy
  
  # Helper function
  reassign_sufnum = function(x) {
    x <- as.numeric(x)
    nums <- seq(length(x))
    # If all numbers are present, return them as is
    if (all(nums %in% x)) { return(x) }
    # Else, reassign
    nums
  }
  
  # # Change empty labels
  # phy$node.label[phy$node.label == ""] <- label_missing_as[1]
  # phy$tip.label[phy$tip.label == ""] <- label_missing_as[1]
  
  n_node <- phy$Nnode
  n_tip <- length(phy$tip.label)
  label_df <- 
    data.frame(is_node = c(rep(FALSE, n_tip), rep(TRUE, n_node)),
               orig = c(phy$tip.label, phy$node.label),
               stringsAsFactors = FALSE) %>%
    mutate(
      orig_suf = .extract_copy_suffix(orig),
      orig_base = .remove_copy_suffix(orig)
      ) %>%
    group_by(orig_base, is_node) %>%
    mutate(
      n_type = n(),
      new_suf = .hyphenate(reassign_sufnum(orig_suf)),
      is_suffixable = (n_type > 1) & (str_length(orig_base) != 0),
      new = str_c(orig_base, ifelse(is_suffixable, new_suf, ""))
    )
  phy$tip.label <- filter(label_df, !is_node)$new
  phy$node.label <- filter(label_df, is_node)$new
  
  phy
}


#' Add copy suffixes
#'
#' Adds copy suffixes to duplicated items in a vector of character strings. Copy
#' suffixes are \code{-1}, \code{-2}, \code{-3}, ...
#'
#' @param x A character vector.
#' @param strip_first A logical, whether to remove existing copy suffixes first.
#' @return A character vector, with suffixes added.
#' @noRd
.add_copy_suffix = function(
  x,
  strip_first = TRUE
) {
  if (strip_first) {
    x <- .remove_copy_suffix(x)
  }
  data.frame(x = x) %>%
    group_by(x) %>%
    mutate(
      n = n(), 
      suffix = ifelse(n > 1, .hyphenate(row_number()), ""),
      x = str_c(x, suffix)
      ) %>%
    .$x
}


#' Remove copy suffixes
#'
#' Removes copy suffixes from a vector of character strings. Copy suffixes are
#' string-final and are comprised of a hyphen followed by numbers.
#'
#' @param x A character vector.
#' @return A character vector, with suffixes removed.
#' @noRd
.remove_copy_suffix = function(x) {
  x <- str_remove(x, "-[0-9]+$")
}


#' Extract copy suffixes
#'
#' Extracts copy suffixes from a vector of character strings. Copy suffixes are
#' string-final and are comprised of a hyphen followed by numbers.
#' 
#' If there is no copy suffix, the return is a zero-length string, \code{""}.
#'
#' @param label A string
#' @return A string
#' @noRd
.extract_copy_suffix = function(label) {
  regex <- "(?<=-)[0-9]+$"
  e <- str_extract(label, regex)
  e[is.na(e)] <- ""
  e
}


#' Get all tip labels from phylo or multiPhylo
#' 
#' Returns all tip labels from one or more trees.
#' 
#' Does this by concatenating the \code{tip.label} vector from each tree.
#' 
#' @param phy A phlyo or multiPhlyo, containing the tree(s).
#' @return A character vector, the tip labels.
#' @noRd
.get_tip_labels = function(phy) {
  if (class(phy) == "phylo") {
    phy$tip.label
  } else if (class(phy) == "multiPhylo") {
    lapply(phy, function(x) x$tip.label) %>%
      unlist()
  }
}


#' Get all node labels from phylo or multiPhylo
#' 
#' Returns all node labels from one or more trees.
#' 
#' Does this by concatenating the \code{node.label} vector from each tree.
#' 
#' @param phy A phlyo or multiPhlyo, containing the tree(s).
#' @return A character vector, the tip labels.
#' @noRd
.get_node_labels = function(phy) {
  if (class(phy) == "phylo") {
    phy$node.label
  } else if (class(phy) == "multiPhylo") {
    lapply(phy, function(x) x$node.label) %>%
      unlist()
  }
}


#' Convert glottolog names to tree labels
#'
#' Converts a vector of glottolog lect names to glottolog tree labels.
#'
#' Comma is converted to forward slash, parentheses are converted to braces, and
#' spaces and apostrophes are removed. The resulting strings are then acceptable
#' labels for nodes and tip in the Newick format.
#'
#' @param name A character vector, the names.
#' @return A character vector, the corresponding tree labels.
#' @noRd
.name_to_label = function(name) {
  name %>% 
    str_replace_all(", ", "/") %>%
    str_replace_all("\\(", "{") %>%
    str_replace_all("\\)", "}") %>%
    str_remove_all("[' ]")
}

#' Add a hyphen before a nonzero-length string
#' 
#' @param x A character vector
#' @return A character vector
#' @noRd
.hyphenate = function(x) {
  str_c(ifelse(str_length(x) == 0, "", "-"), x)
}
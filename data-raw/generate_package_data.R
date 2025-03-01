# This script is for updating with version 4.8

# Note: it uses the package functions extract_glottocode() and .name_to_label()

library(tidyverse)


#### FYI -- what to update in the package for a new glottolog version:

# data.R
#   [ ] title line "Trees from glottolog, versions 4.0 - 5.X"
#   [ ] list of abjects in @rdname glottolog_trees (in @format and ",,,")
#   [ ] the number of elements in these
#   [ ] title line "Geographical data from glottolog, versions 4.0 - 5.X"
#   [ ] list of abjects in @rdname glottolog_geography (in @format and ",,,")
#   [ ] the number of elements in these
# input_checks.R
#   [ ] .check_glottolog_version()
# metadata.R
#   [ ] get_glottolog_languages()  -- in the Roxygen metadata for @param
#   [ ] get_glottolog_phylo_geo() -- in the if/else list
#   [ ] .get_newest_version() -- in Roxygen metadata and the function code
# topology.R
#   [ ] get_glottolog_trees() -- in the if/else list
# zzz.R
#   [ ] .onAttach()
# DESCRIPTION
#   [ ] version number
# inst/CITATION
#   [ ] reference to glottolog


#### EXTERNAL DATA

# Read the new glottolog data, which has been put in the data-raw directory with
# the appropriate name (you'll need to add a final "_X.X") Note that glottolog
# uses square brackets, in violation of the Newick standard. A recent update to
# ape::read.tree() now removes substrings between square brackets, making the
# reading of glottolog's files more tricky. Here, we read the text and change
# square brackets to < and >.

read_newick = function(version) {
  scan(file = str_c("data-raw/tree_glottolog_newick_", version, ".txt"),
       what = "character", sep = "\n", quiet = TRUE) %>%
    str_replace_all("\\[", "<") %>%
    str_replace_all("\\]", ">") %>%
    # New step needed due to the change in behavior of ape::read.tree() as
    # of ape v4.1
    str_remove_all("[' ]") 
}

read_geo = function(version) {
  read.csv(str_c("data-raw/languages_and_dialects_geo_", version, ".csv"),
           stringsAsFactors = FALSE)
}


newick_v4.0 <- read_newick(version = "4-0")
newick_v4.1 <- read_newick(version = "4-1")
newick_v4.2 <- read_newick(version = "4-2")
newick_v4.3 <- read_newick(version = "4-3")
newick_v4.4 <- read_newick(version = "4-4")
newick_v4.5 <- read_newick(version = "4-5")
newick_v4.6 <- read_newick(version = "4-6")
newick_v4.7 <- read_newick(version = "4-7")
newick_v4.8 <- read_newick(version = "4-8")
newick_v5.0 <- read_newick(version = "5-0")

glottolog_trees_v4.0 <- ape::read.tree(text = newick_v4.0)
glottolog_trees_v4.1 <- ape::read.tree(text = newick_v4.1)
glottolog_trees_v4.2 <- ape::read.tree(text = newick_v4.2)
glottolog_trees_v4.3 <- ape::read.tree(text = newick_v4.3)
glottolog_trees_v4.4 <- ape::read.tree(text = newick_v4.4)
glottolog_trees_v4.5 <- ape::read.tree(text = newick_v4.5)
glottolog_trees_v4.6 <- ape::read.tree(text = newick_v4.6)
glottolog_trees_v4.7 <- ape::read.tree(text = newick_v4.7)
glottolog_trees_v4.8 <- ape::read.tree(text = newick_v4.8)
glottolog_trees_v5.0 <- ape::read.tree(text = newick_v5.0)

glottolog_geography_v4.0 <- read_geo(version = "4-0")
glottolog_geography_v4.1 <- read_geo(version = "4-1")
glottolog_geography_v4.2 <- read_geo(version = "4-2")
glottolog_geography_v4.3 <- read_geo(version = "4-3")
glottolog_geography_v4.4 <- read_geo(version = "4-4")
glottolog_geography_v4.5 <- read_geo(version = "4-5")
glottolog_geography_v4.6 <- read_geo(version = "4-6")
glottolog_geography_v4.7 <- read_geo(version = "4-7")
glottolog_geography_v4.8 <- read_geo(version = "4-8")
glottolog_geography_v5.0 <- read_geo(version = "5-0")


# Add the new external datasets 

usethis::use_data(
  internal = FALSE, 
  overwrite = TRUE,
  glottolog_trees_v4.0, 
  glottolog_trees_v4.1, 
  glottolog_trees_v4.2, 
  glottolog_trees_v4.3, 
  glottolog_trees_v4.4, 
  glottolog_trees_v4.5, 
  glottolog_trees_v4.6, 
  glottolog_trees_v4.7, 
  glottolog_trees_v4.8, 
  glottolog_trees_v5.0, 
  glottolog_geography_v4.0,
  glottolog_geography_v4.1,
  glottolog_geography_v4.2,
  glottolog_geography_v4.3,
  glottolog_geography_v4.4,
  glottolog_geography_v4.5,
  glottolog_geography_v4.6,
  glottolog_geography_v4.7,
  glottolog_geography_v4.8,
  glottolog_geography_v5.0
)


#### INTERNAL DATA

extract_name = function(labels) {
  regex <- "^[^<]+"
  str_extract(labels, regex)
}

phy4.0 <- glottolog_trees_v4.0
phy4.1 <- glottolog_trees_v4.1
phy4.2 <- glottolog_trees_v4.2
phy4.3 <- glottolog_trees_v4.3
phy4.4 <- glottolog_trees_v4.4
phy4.5 <- glottolog_trees_v4.5
phy4.6 <- glottolog_trees_v4.6
phy4.7 <- glottolog_trees_v4.7
phy4.8 <- glottolog_trees_v4.8
phy5.0 <- glottolog_trees_v5.0

geo4.0 <- glottolog_geography_v4.0
geo4.1 <- glottolog_geography_v4.1
geo4.2 <- glottolog_geography_v4.2
geo4.3 <- glottolog_geography_v4.3
geo4.4 <- glottolog_geography_v4.4
geo4.5 <- glottolog_geography_v4.5
geo4.6 <- glottolog_geography_v4.6
geo4.7 <- glottolog_geography_v4.7
geo4.8 <- glottolog_geography_v4.8
geo5.0 <- glottolog_geography_v5.0

root_labels4.0 <- lapply(phy4.0, function(p) p$node.label[1]) %>% unlist()
root_labels4.1 <- lapply(phy4.1, function(p) p$node.label[1]) %>% unlist()
root_labels4.2 <- lapply(phy4.2, function(p) p$node.label[1]) %>% unlist()
root_labels4.3 <- lapply(phy4.3, function(p) p$node.label[1]) %>% unlist()
root_labels4.4 <- lapply(phy4.4, function(p) p$node.label[1]) %>% unlist()
root_labels4.5 <- lapply(phy4.5, function(p) p$node.label[1]) %>% unlist()
root_labels4.6 <- lapply(phy4.6, function(p) p$node.label[1]) %>% unlist()
root_labels4.7 <- lapply(phy4.7, function(p) p$node.label[1]) %>% unlist()
root_labels4.8 <- lapply(phy4.8, function(p) p$node.label[1]) %>% unlist()
root_labels5.0 <- lapply(phy5.0, function(p) p$node.label[1]) %>% unlist()


##### Tabulate the labels of families' trees

tabulate_fam_labs = function(phy, root_labels) {
  data.frame(
    tree = 1:length(phy),
    family_name = extract_name(root_labels),
    family_glottocode = extract_glottocode(root_labels),
    stringsAsFactors = FALSE
  )
}

glottolog_family_labels_v4.0 <- tabulate_fam_labs(phy4.0, root_labels4.0)
glottolog_family_labels_v4.1 <- tabulate_fam_labs(phy4.1, root_labels4.1)
glottolog_family_labels_v4.2 <- tabulate_fam_labs(phy4.2, root_labels4.2)
glottolog_family_labels_v4.3 <- tabulate_fam_labs(phy4.3, root_labels4.3)
glottolog_family_labels_v4.4 <- tabulate_fam_labs(phy4.4, root_labels4.4)
glottolog_family_labels_v4.5 <- tabulate_fam_labs(phy4.5, root_labels4.5)
glottolog_family_labels_v4.6 <- tabulate_fam_labs(phy4.6, root_labels4.6)
glottolog_family_labels_v4.7 <- tabulate_fam_labs(phy4.7, root_labels4.7)
glottolog_family_labels_v4.8 <- tabulate_fam_labs(phy4.8, root_labels4.8)
glottolog_family_labels_v5.0 <- tabulate_fam_labs(phy5.0, root_labels5.0)


#### Tabulate gottolog tree vertices and geo data

tabulate_phylo_geo = function(phy, geo, family_labels) {
  # Compile a table of vertices in trees
  lapply(
    1:length(phy),
    function(i) {
      bind_rows(
        data.frame(vertex_type = "node",
                   vertex_label = phy[[i]]$node.label,
                   stringsAsFactors = FALSE),
        data.frame(vertex_type = "tip",
                   vertex_label = phy[[i]]$tip.label,
                   stringsAsFactors = FALSE)
      ) %>% 
        mutate(tree = i) }
  ) %>% 
    bind_rows() %>%
    mutate(
      glottocode = extract_glottocode(vertex_label),
      vertex_name = extract_name(vertex_label)
    ) %>%
    
    # Combine the geo and vertex info
    full_join(
      
      # Add a column vertex_name, the equivalent of name
      # that we'd expect to see in a vertex label
      geo %>% mutate(vertex_name = .name_to_label(name)), 
      ., 
      by = c("vertex_name", "glottocode")
    ) %>% 
    
    # Add family names
    left_join(family_labels, by = "tree") %>%
    select(glottocode, isocodes, name, level,
           vertex_type, vertex_label, vertex_name, 
           macroarea, latitude, longitude,
           family_glottocode, family_name, tree) %>%
    arrange(glottocode) %>%
    as.data.frame()
}

glottolog_phylo_geo_v4.0 <-
  tabulate_phylo_geo(phy4.0, geo4.0, glottolog_family_labels_v4.0)
glottolog_phylo_geo_v4.1 <-
  tabulate_phylo_geo(phy4.1, geo4.1, glottolog_family_labels_v4.1)
glottolog_phylo_geo_v4.2 <-
  tabulate_phylo_geo(phy4.2, geo4.2, glottolog_family_labels_v4.2)
glottolog_phylo_geo_v4.3 <-
  tabulate_phylo_geo(phy4.3, geo4.3, glottolog_family_labels_v4.3)
glottolog_phylo_geo_v4.4 <-
  tabulate_phylo_geo(phy4.4, geo4.4, glottolog_family_labels_v4.4)
glottolog_phylo_geo_v4.5 <-
  tabulate_phylo_geo(phy4.5, geo4.5, glottolog_family_labels_v4.5)
glottolog_phylo_geo_v4.6 <-
  tabulate_phylo_geo(phy4.6, geo4.6, glottolog_family_labels_v4.6)
glottolog_phylo_geo_v4.7 <-
  tabulate_phylo_geo(phy4.7, geo4.7, glottolog_family_labels_v4.7)
glottolog_phylo_geo_v4.8 <-
  tabulate_phylo_geo(phy4.8, geo4.8, glottolog_family_labels_v4.8)
glottolog_phylo_geo_v5.0 <-
  tabulate_phylo_geo(phy5.0, geo5.0, glottolog_family_labels_v5.0)


# Tabulate glottolog families and macroareas

tabulate_family_geo = function(phylo_geo) {
  phylo_geo %>%
    filter(!is.na(macroarea), !is.na(tree)) %>% 
    group_by(tree, family_name, family_glottocode, macroarea) %>% 
    summarise(n = n()) %>% 
    group_by(tree, family_name, family_glottocode) %>%
    arrange(-n) %>%
    mutate(macroarea_n = str_c(macroarea, ":", n)) %>% 
    summarise(
      main_macroarea = macroarea[1],
      all_macroareas = str_c(macroarea_n, collapse = ", ")
    ) %>% 
    arrange(tree) %>%
    as.data.frame()
}

glottolog_family_geo_v4.0 <- tabulate_family_geo(glottolog_phylo_geo_v4.0)
glottolog_family_geo_v4.1 <- tabulate_family_geo(glottolog_phylo_geo_v4.1)
glottolog_family_geo_v4.2 <- tabulate_family_geo(glottolog_phylo_geo_v4.2)
glottolog_family_geo_v4.3 <- tabulate_family_geo(glottolog_phylo_geo_v4.3)
glottolog_family_geo_v4.4 <- tabulate_family_geo(glottolog_phylo_geo_v4.4)
glottolog_family_geo_v4.5 <- tabulate_family_geo(glottolog_phylo_geo_v4.5)
glottolog_family_geo_v4.6 <- tabulate_family_geo(glottolog_phylo_geo_v4.6)
glottolog_family_geo_v4.7 <- tabulate_family_geo(glottolog_phylo_geo_v4.7)
glottolog_family_geo_v4.8 <- tabulate_family_geo(glottolog_phylo_geo_v4.8)
glottolog_family_geo_v5.0 <- tabulate_family_geo(glottolog_phylo_geo_v5.0)


# Add the new internal datasets

usethis::use_data(
  internal = TRUE, 
  overwrite = TRUE,
  
  glottolog_family_labels_v4.0,
  glottolog_family_labels_v4.1,
  glottolog_family_labels_v4.2,
  glottolog_family_labels_v4.3,
  glottolog_family_labels_v4.4,
  glottolog_family_labels_v4.5,
  glottolog_family_labels_v4.6,
  glottolog_family_labels_v4.7,
  glottolog_family_labels_v4.8,
  glottolog_family_labels_v5.0,
  
  glottolog_phylo_geo_v4.0,
  glottolog_phylo_geo_v4.1,
  glottolog_phylo_geo_v4.2,
  glottolog_phylo_geo_v4.3,
  glottolog_phylo_geo_v4.4,
  glottolog_phylo_geo_v4.5,
  glottolog_phylo_geo_v4.6,
  glottolog_phylo_geo_v4.7,
  glottolog_phylo_geo_v4.8,
  glottolog_phylo_geo_v5.0,
  
  glottolog_family_geo_v4.0,
  glottolog_family_geo_v4.1,
  glottolog_family_geo_v4.2,
  glottolog_family_geo_v4.3,
  glottolog_family_geo_v4.4,
  glottolog_family_geo_v4.5,
  glottolog_family_geo_v4.6,
  glottolog_family_geo_v4.7,
  glottolog_family_geo_v4.8,
  glottolog_family_geo_v5.0
)

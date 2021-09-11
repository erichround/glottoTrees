# This script is for updating with version 4.4

# Note: it uses the package functions extract_glottocode() and .name_to_label()

#### EXTERNAL DATA

# Read the new glottolog data, which has been put in the data-raw directory with
# the appropriate name (you'll need to add a final "_vX.X")

glottolog_trees_v4.4 <- ape::read.tree("data-raw/tree_glottolog_newick_4-4.txt")
glottolog_geography_v4.4 <- 
  read.csv("data-raw/languages_and_dialects_geo_4-4.csv",
           stringsAsFactors = FALSE)

# Add the new external datasets 

usethis::use_data(
  internal = FALSE, 
  overwrite = FALSE,
  glottolog_trees_v4.4, 
  glottolog_geography_v4.4
)


#### INTERNAL DATA

extract_name = function(labels) {
  regex <- "^[^\\[]+"
  str_extract(labels, regex)
}

phy <- glottolog_trees_v4.4
geo <- glottolog_geography_v4.4
root_labels <- lapply(phy, function(p) p$node.label[1]) %>% unlist()


##### Tabulate the labels of families' trees

glottolog_family_labels_v4.4 <-
  data.frame(
    tree = 1:length(phy),
    family_name = extract_name(root_labels),
    family_glottocode = extract_glottocode(root_labels),
    stringsAsFactors = FALSE
  )


#### Tabulate gottolog tree vertices and geo data

glottolog_phylo_geo_v4.4 <-
  
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
  left_join(glottolog_family_labels_v4.4, by = "tree") %>%
  select(glottocode, isocodes, name, level,
         vertex_type, vertex_label, vertex_name, 
         macroarea, latitude, longitude,
         family_glottocode, family_name, tree) %>%
  arrange(glottocode) %>%
  as.data.frame()



# Tabulate glottolog families and macroareas
  
glottolog_family_geo_v4.4 <-
  glottolog_phylo_geo_v4.4 %>%
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


# Add the new internal datasets 

usethis::use_data(
  internal = TRUE, 
  overwrite = TRUE,
  
  # Older versions
  glottolog_family_labels_v4.3,
  glottolog_phylo_geo_v4.3,
  glottolog_family_geo_v4.3,
  
  # New version beign added now
  glottolog_family_labels_v4.4,
  glottolog_phylo_geo_v4.4,
  glottolog_family_geo_v4.4
)
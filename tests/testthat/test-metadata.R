context("Glottolog metadata")
library(PhyloWeights)

test_that("languages metadata v4.3 has the expected contents", {
  m <- get_glottolog_languages("4.3")
  expect_equal(class(m), "data.frame")
  expect_equal(colnames(m), 
               c("glottocode", "isocodes", "name", "name_in_tree", 
                 "position", "tree", "tree_name"))
  expect_equal(class(m$glottocode), "character")
  expect_equal(class(m$isocodes), "character")
  expect_equal(class(m$name), "character")
  expect_equal(class(m$name_in_tree), "character")
  expect_equal(class(m$position), "character")
  expect_equal(class(m$tree), "integer")
  expect_equal(class(m$tree_name), "character")
})

test_that("languages metadata v4.4 has the expected contents", {
  m <- get_glottolog_languages("4.4")
  expect_equal(class(m), "data.frame")
  expect_equal(colnames(m), 
               c("glottocode", "isocodes", "name", "name_in_tree", 
                 "position", "tree", "tree_name"))
  expect_equal(class(m$glottocode), "character")
  expect_equal(class(m$isocodes), "character")
  expect_equal(class(m$name), "character")
  expect_equal(class(m$name_in_tree), "character")
  expect_equal(class(m$position), "character")
  expect_equal(sort(unique(m$position)), c("node", "tip"))
  expect_equal(class(m$tree), "integer")
  expect_equal(class(m$tree_name), "character")
})

test_that("families metadata v4.3 has the expected contents", {
  m <- get_glottolog_families("4.3")
  expect_equal(class(m), "data.frame")
  expect_equal(colnames(m), 
               c("tree", "tree_name", "n_tips", "n_nodes", 
                 "main_macroarea"))
  expect_equal(class(m$tree), "integer")
  expect_equal(class(m$tree_name), "character")
  expect_equal(class(m$n_tips), "integer")
  expect_equal(class(m$n_nodes), "integer")
  expect_equal(class(m$main_macroarea), "character")
  expect_equal(sort(unique(m$main_macroarea)),
               c("Africa", "Australia", "Eurasia", "North America", 
                 "Papunesia", "South America"))
})

test_that("families metadata v4.4 has the expected contents", {
  m <- get_glottolog_families("4.3")
  expect_equal(class(m), "data.frame")
  expect_equal(colnames(m), 
               c("tree", "tree_name", "n_tips", "n_nodes", 
                 "main_macroarea"))
  expect_equal(class(m$tree), "integer")
  expect_equal(class(m$tree_name), "character")
  expect_equal(class(m$n_tips), "integer")
  expect_equal(class(m$n_nodes), "integer")
  expect_equal(class(m$main_macroarea), "character")
  expect_equal(sort(unique(m$main_macroarea)),
               c("Africa", "Australia", "Eurasia", "North America", 
                 "Papunesia", "South America"))
})

test_that("phylo_geo metadata v4.3 has the expected contents", {
  m <- get_glottolog_phylo_geo("4.3")
  expect_equal(class(m), "data.frame")
  expect_equal(colnames(m),
               c("glottocode", "isocodes", "name", "level", 
                 "vertex_type", "vertex_label", "vertex_name", 
                 "macroarea", "latitude", "longitude", 
                 "family_glottocode", "family_name", "tree"))
  
  expect_equal(class(m$glottocode), "character")
  expect_equal(class(m$isocodes), "character")
  expect_equal(class(m$name), "character")
  expect_equal(class(m$level), "character")
  expect_equal(sort(unique(m$level)), c("dialect", "language"))
  expect_equal(class(m$vertex_type), "character")
  expect_equal(sort(unique(m$vertex_type)), c("node", "tip"))
  expect_equal(class(m$vertex_label), "character")
  expect_equal(class(m$vertex_name), "character")
  expect_equal(class(m$macroarea), "character")
  expect_equal(sort(unique(m$macroarea)),
               c("", "Africa", "Australia", "Eurasia", "North America", 
                 "Papunesia", "South America"))
  expect_equal(class(m$latitude), "numeric")
  expect_equal(class(m$longitude), "numeric")
  expect_equal(class(m$family_glottocode), "character")
  expect_equal(class(m$family_name), "character")
  expect_equal(class(m$tree), "integer")
})

test_that("phylo_geo metadata v4.4 has the expected contents", {
  m <- get_glottolog_phylo_geo("4.4")
  expect_equal(class(m), "data.frame")
  expect_equal(colnames(m),
               c("glottocode", "isocodes", "name", "level", 
                 "vertex_type", "vertex_label", "vertex_name", 
                 "macroarea", "latitude", "longitude", 
                 "family_glottocode", "family_name", "tree"))
  
  expect_equal(class(m$glottocode), "character")
  expect_equal(class(m$isocodes), "character")
  expect_equal(class(m$name), "character")
  expect_equal(class(m$level), "character")
  expect_equal(sort(unique(m$level)), c("dialect", "language"))
  expect_equal(class(m$vertex_type), "character")
  expect_equal(sort(unique(m$vertex_type)), c("node", "tip"))
  expect_equal(class(m$vertex_label), "character")
  expect_equal(class(m$vertex_name), "character")
  expect_equal(class(m$macroarea), "character")
  expect_equal(sort(unique(m$macroarea)),
               c("", "Africa", "Australia", "Eurasia", "North America", 
                 "Papunesia", "South America"))
  expect_equal(class(m$latitude), "numeric")
  expect_equal(class(m$longitude), "numeric")
  expect_equal(class(m$family_glottocode), "character")
  expect_equal(class(m$family_name), "character")
  expect_equal(class(m$tree), "integer")
})

test_that("the glottolog version parameter is handled correctly", {
  expect_equal(get_glottolog_languages("4.3"),
               get_glottolog_languages(4.3))
  expect_equal(get_glottolog_families("4.3"),
               get_glottolog_families(4.3))
  expect_equal(get_glottolog_phylo_geo("4.3"),
               get_glottolog_phylo_geo(4.3))
  expect_equal(which_tree("Tangkic", "4.3"),
               which_tree("Tangkic", 4.3))
  expect_error(get_glottolog_languages("4.2"))
  expect_error(get_glottolog_languages(c("4.3", "4.4")))
})

test_that("which_tree returns named integers", {
  w <- which_tree("Tangkic", "4.4")
  expect_equal(class(w), "integer") 
  expect_equal(class(w), "integer") 
  expect_equal(is.null(names(w)), FALSE)
  expect_equal(w, c("Tangkic" = 300))
})

test_that("which_tree with no recognised family returns an error", {
  expect_error(which_tree())
  expect_error(which_tree(glottolog_version = "4.4"))
  expect_error(which_tree(c("caddoan", "tangkic")))
})
  
test_that("which_tree with some unrecognised families returns warning and NA", {
  expect_warning(which_tree(c("caddoan", "Tangkic"), glottolog_version = "4.4"))
  expect_equal(
    sum(is.na(suppressWarnings(
      which_tree(c("caddoan", "Tangkic"), glottolog_version = "4.4")))), 
    1)
})        
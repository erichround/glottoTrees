.onAttach = function(libname, pkgname) {
  packageStartupMessage(
    str_c(
    "glottoTrees contains data from glottolog.com under glottolog's ",
    "CC-BY 4.0 licence.\n",
    "Available versions are glottolog 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, ",
    "4.6, 4.7, 4.8 and 5.0.\n",
    "The newest version is used by default."
  ))
}
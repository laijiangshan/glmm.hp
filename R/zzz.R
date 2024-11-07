
.onAttach <- function(libname, pkgname) {

  cite_info <- citation(pkgname)
  

  packageStartupMessage("Thank you for using this package! If you use this package in your research, please cite the following references：\n")

  packageStartupMessage(cite_info)
}


.onAttach <- function(libname, pkgname) {
  cite_info <- utils::citation(pkgname)
  cite_info <- "Jiangshan Lai, Weijie Zhu, Dongfang Cui, Lingfeng Mao(2023). Extension of the glmm.hp package to Zero-Inflated generalized linear mixed models and multiple regression. Journal of Plant Ecology,16(6):rtad038"
  packageStartupMessage("Thank you for using this package! If you use this package in your research, please cite the following references \n")
  packageStartupMessage(cite_info)
}

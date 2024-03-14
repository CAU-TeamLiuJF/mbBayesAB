#!/public/home/liujf/software/program/R-4.3.1/bin/Rscript
#!

## 检查需要的R包是否都已安装  ##
## liwn 2023-09-07 ##

# 要检查的 R 包列表
packages <- c("data.table", "getopt", "MASS", "RcppEigen", "Matrix", "ggplot2", "reshape2",
  "pedigree", "RcppArmadillo", "dplyr", "coda", "lattice", "stringr", "Cairo", "Rcpp", "reshape")

# 检查并安装缺失的包
missing_packages <- packages[!packages %in% installed.packages()]
if (length(missing_packages) > 0) {
  message("以下R包未安装: ")
  cat(missing_packages, sep = "\n")
  message("请在R语言中运行以下命令安装: ")
  command <- paste0("install.packages(c(\"", paste(missing_packages, collapse = "\", \""), "\"))")
  cat(command, "\n")
}

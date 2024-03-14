#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 将两个群体的dmu表型文件进行合并，以便进行联合评估  ###
## debug
# opt = list(fold=5, pheno='LineA_LineB_pheno.txt', iyse='2', phe='10')

## command parameters
spec <- matrix(
  c("phef",    "P", 1, "character", "[Required] phenotype file name",
    "pops",    "N", 1, "character", "[Required] population names",
    "nInt",    "n", 1, "integer",   "[Optional] number of integer columns",
    "rep",     "r", 1, "integer",   "[Optional] repeat times [1]",
    "fold",    "f", 1, "integer",   "[Optional] cross validation fold [1]",
    "type",    "m", 1, "character", "[Optional] blend / union/multi(two trait model) [blend]",
    "pheCol",  "p", 1, "integer",   "[Optional] Number of columns of phenotypic in real columns [1]",
    "miss",    "M", 1, "double",    "[Optional] Identifier of the missing phenotype [-99]",
    "fixPop",  "F", 2, "character", "[Optional] Add fixed breed effects to the model [NULL]",
    "addpop",  "a", 1, "character", "[Optional] Add breed labels on some effects, like:'2:3', e.g. 2-breed, 3-breed",
    "header",  "H", 1, "logical",   "[Optional] whether include header in output phenotype files/NULL",
    "out",     "o", 1, "character", "[Optional] Output file prefix\n",
    "overwri", "O", 2, "character", "[Optional] whether overwrite file if pheno.txt already exist\n",
    "help",    "h", 0, "logical",   "This is Help!"),
    byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)

## 参数检查
if (!is.null(opt$help) || is.null(opt$phef) || is.null(opt$pops)) {
  ## print help message
  cat(paste(getopt::getopt(spec = spec, usage = TRUERUE), "\n"))
  quit(status = -1)
}

# Load packages
# cat('Loading packages needed...\n')
suppressPackageStartupMessages(library("data.table"))

## 默认参数
if (is.null(opt$rep)) opt$rep <- 1
if (is.null(opt$fold)) opt$fold <- 1
if (is.null(opt$out)) opt$out <- "pheno.txt"
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$pheCol)) opt$pheCol <- 1
if (is.null(opt$type)) opt$type <- "blend"
if (is.null(opt$header)) header <- FALSE else header <- TRUE

## 解析参数
pops <- strsplit(opt$pops, " ")[[1]]
np <- length(pops)

## 合并表型
for (r in 1:opt$rep) { # nolint
  for (f in 1:opt$fold) {
    ## 输出表型文件
    out <- gsub("#rep#", r, opt$out)
    out <- gsub("#val#", f, out)

    ## 确保目标文件夹存在
    outdir <- dirname(out)
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    ## 查看输出文件是否存在
    if (file.exists(out) && opt$overwri != "true") next

    ## 合并各个品种的表型文件
    phe <- data.table()
    for (i in seq_len(np)) {
      nint <- opt$nInt
      ## 文件名准备
      phefi <- gsub("#rep#", r, opt$phef)
      phefi <- gsub("#val#", f, phefi)
      phefi <- gsub("#breed#", pops[i], phefi)

      ## 检查品种内表型文件是否存在
      if (!file.exists(phefi)) {
        cat(paste("Error: ", phefi, " not exist!\n"))
        quit(status = -1)
      }

      ## 读取表型文件
      phei <- fread(phefi)

      ## 整型变量列数
      if (is.null(nint)) nint <- sum(sapply(phei, class) == "integer")

      ## 提取需要的列
      phei2 <- subset(phei, select = c(1:nint, nint + opt$pheCol))
      names(phei2) <- c(paste0("I", 1:nint), "P")

      ## 添加品种效应
      if (!is.null(opt$fixPop)) {
        phei2$pop <- i
        phei2 <- subset(phei2, select = c(1:nint, ncol(phei2), ncol(phei2) - 1))
        nint <- nint + 1
      }

      ## 合并
      if (opt$type == "union" || opt$type == "multi") {
        ## 多性状模型
        phei2[, paste0("P", seq_len(np))] <- opt$miss
        phei2[[paste0("P", i)]] <- phei2[[nint + 1]]
        phei2 <- subset(phei2, select = names(phei2) != "P") # drop the original phenotype
      } else if (opt$type != "blend") {
        cat("The type parameter is incorrect. It can only be union/blend/multi\n")
        quit(status = -1)
      }

      ## 合并
      phe <- rbind(phe, phei2)
    }

    if (!is.null(opt$addpop)) {
      ## 在指定原固定效应分组添加品种标识符以区别
      cols <- as.numeric(unlist(strsplit(opt$addpop, ":")))
      for (c in cols) {
        if (c > opt$nInt) {
          cat("Non-integer columns cannot add breed identifiers.\n")
          quit(status = -1)
        }
        phe[[paste0("V", c)]] <- paste0(phe$pop, phe[[paste0("V", c)]])
      }
    }

    ## 写出表型文件
    fwrite(phe, out, col.names = header, sep = " ")

    ## 报告
    if (r == 1 && opt$fold > 1) {
      ## 参考群大小
      if (opt$type == "blend") {
        nind_ref <- sum(subset(phe, select = ncol(phe)) != opt$miss)
      } else {
        nind_val <- sum(apply(phe[, (ncol(phe) - np + 1):ncol(phe)], 1, function(x) all(x == opt$miss)))
        nind_ref <- nrow(phe) - nind_val
      }

      cat("records in fold ", f, " reference populations: ", nind_ref, "\n", sep = "")
    }
  }
}

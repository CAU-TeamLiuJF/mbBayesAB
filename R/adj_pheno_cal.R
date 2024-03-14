#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 根据dmu结果文件和表型文件计算校正表型(校正固定效应和非加性效应)

# 加载需要的程序包
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## 命令行参数
spec <- matrix(
  c("DIR",     "d", 1, "character", "[Required] Full dataset parameter card prefix",
    "phe",     "p", 1, "character", "[Required] Full dataset phenotype file name",
    "idc",     "c", 1, "integer",   "[Optional] id columns in phenotype file [1]",
    "traiti",  "i", 1, "integer",   "[Optional] trait rank [1]",
    "nTrait",  "n", 1, "integer",   "[Optional] number of traits [1]",
    "add_sol", "a", 1, "integer",   "[Optional] number of traits [1]",
    "append",  "A", 1, "character", "[Optional] append to the result file\n",
    "out",     "o", 1, "character", "[Optional] output file name/phe_adj.txt\n",
    "help",    "h", 0, "logical",   "This is Help!"),
    byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$phe) || is.null(opt$DIR)) {
    cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
}

## 检查需要的文件是否存在
ex <- TRUE
if (!file.exists(opt$phe)) {
    file_name <- opt$phe
} else if (!file.exists(paste0(opt$DIR, ".SOL"))) {
    file_name <- paste0(opt$DIR, ".SOL")
} else if (!file.exists(paste0(opt$DIR, ".RESIDUAL"))) {
    file_name <- paste0(opt$DIR, ".RESIDUAL")
} else {
    ex <- FALSE
}
if (ex) {
    cat(file_name, "not found!\n")
    quit(status = 1)
}

## 默认参数
if (is.null(opt$traiti)) opt$traiti <- 1
if (is.null(opt$nTrait)) opt$nTrait <- 1
if (is.null(opt$idc)) opt$idc <- 1
if (is.null(opt$out)) opt$out <- "phe_adj.txt"

## 表型文件，用于匹配残差id
phe <- fread(opt$phe)
names(phe)[opt$idc] <- "id"

## ebv
sol <- fread(paste0(opt$DIR, ".SOL"))
if (is.null(opt$add_sol)) opt$add_sol <- max(sol$V1)
if (!opt$add_sol %in% unique(sol$V1)) {
    cat("add_sol cant be ", opt$add_sol, "!\n")
    quit(status = 1)
}
ebv_all <- subset(sol, V1 == opt$add_sol & V2 == opt$traiti, c(5, 8))
names(ebv_all) <- c("id", "ebv")

## 提取残差（暂时只能用于单性状模型）
res <- fread(paste0(opt$DIR, ".RESIDUAL"))
re <- subset(res, select = c(1, 4))
names(re) <- c("rows", "re")
if (max(re$rows) > nrow(phe)) {
    cat("sol file does not match RESIDUAL file!\n")
    quit(status = 1)
}
re$id <- phe$id[re$rows]

## 计算校正表型
y_adj <- left_join(ebv_all, re, by = "id")
y_adj$Yhat <- y_adj$ebv + y_adj$re

## 去除缺失值
adj <- subset(y_adj, !is.na(Yhat), select = c("id", "Yhat"))
nadj <- nrow(adj)

## 是否追加
if (opt$append == "true") {
    append <- TRUE
} else {
    append <- FALSE
}

if (nadj > 0) {
    ## 输出校正表型
    fwrite(adj, opt$out, col.names = FALSE, append = append, sep = " ")
    cat("phenotypes (", nrow(adj), ") corrected for all other effects output to:", opt$out, "\n")
} else {
    cat("error! no records in results file.\n")
}

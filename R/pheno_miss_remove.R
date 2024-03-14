#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 将输入的表型文件按照基因型文件个体顺序进行排序
## opt <- list(fam='DD_num.fam', pheno='gwas_pca_DD.fam')

# 加载需要的程序包
suppressPackageStartupMessages(library("getopt"))

## 命令行参数
spec <- matrix(
  c("file",     "I",  1, "character", "[Required] input file name",
    "col",      "C",  1, "integer",   "[Required] column to find missing value",
    "idC",      "i",  1, "integer",   "[Required] id column/1",
    "miss",     "M",  1, "double",    "[Optional] value to be treat as missing value/-99",
    "backName", "B",  1, "character", "[Optional] new file names to save origin file/*.bc",
    "map",      "m",  1, "character", "[Optional] plink map/fam file",
    "missid",   "s",  1, "character", "[Optional] file contain missing phenotype id/No_phe.id",
    "out",      "o",  1, "character", "[Optional] file contain missing phenotype id/No_phe.id",
    "NotRm",    "n",  0, "logical",   "[Optional] do not remove id with missing phenotype",
    "help",     "h",  0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$file) || is.null(opt$col)) {
    cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
    quit()
}

## 默认参数
if (is.null(opt$idC)) opt$idC <- 1
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$backName)) opt$backName <- paste0(opt$file, ".bc")
if (is.null(opt$missid)) opt$missid <- "No_phe.id"

## 读取文件
phe <- read.table(opt$file)

## 命名
if (opt$col > ncol(phe) || opt$idC > ncol(phe)) {
    cat("parameter --col (", opt$col, ") or --idC (", opt$idC, ") bigger than phenotype file columns!\n")
    quit()
} else {
    names(phe)[opt$col] <- "value"
    names(phe)[opt$idC] <- "id"
}


## 查找缺失表型个体
missid <- phe$id[phe$value <= opt$miss]
if (length(missid) > 0) {
    if (!is.null(opt$map)) {
        ## 输出缺失表型的基因型个体
        map <- read.table(opt$map)
        fid_iid <- map[map$V2 %in% missid, 1:2]
        write.table(fid_iid, opt$missid, row.names = FALSE, col.names = FALSE, quote = FALSE)
    } else {
        write.table(missid, opt$missid, row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    cat("individuals ID with the missing phenotype output to:", opt$missid, "\n")

    ## 在表型文件中删除缺失表型的个体
    if (is.null(opt$NotRm)) {
        phe_new <- phe[!phe$id %in% missid, ]
        write.table(phe_new, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
        # write.table(phe, opt$backName, row.names = FALSE, col.names = FALSE, quote = FALSE)
        cat("phenotype (", nrow(phe_new), ") that no missing values in column", opt$col,
            "has been output to", opt$out, "\n")
        cat("remove", length(missid), "individuals with the missing value in input files\n")
    }
} else {
    write.table(phe, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)
    cat("No individuals with the missing phenotype found, output file copy.\n")
}

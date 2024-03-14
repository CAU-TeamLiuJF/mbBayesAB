#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 用于跑完完整数据集和剔除验证群的缩减数据集dmu4后进行准确性计算

# 加载需要的程序包
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## 命令行参数
spec <- matrix(
  c("famf",        "G", 1, "character", "[Required] plink fam file",
    "phef",        "P", 1, "character", "[Required] phenotype file",
    "phec",        "p", 1, "integer",   "[Optional] phenotype column in phenotype file",
    "phe_idc",     "I", 1, "integer",   "[Optional] id column in phenotype file",
    "miss",        "m", 1, "double",    "[Optional] miss value for real value",
    "num_int",     "i", 1, "integer",   "[Optional] number of integer variables",
    "mean",        "M", 1, "integer",   "[Optional] population mean columns",
    "breed",       "b", 1, "integer",   "[Optional] breed effect column",
    "out",         "O", 1, "character", "[Optional] output file name prefix/phef\n",
    "rm",          "r", 1, "character", "[Optional] remove all/single trait missing recodes/NULL\n",
    "rmOut",       "R", 1, "character", "[Optional] remove ids\n",
    "header",      "H", 0, "logical",   "[Optional] add header in output phenotype file",
    "sort",        "d", 0, "logical",   "[Optional] sort individuals in phenotype file (consistent with famf)",
    "help",        "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$phef) || is.null(opt$famf)) {
    cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
    quit()
}

## 默认参数
if (is.null(opt$header)) header <- FALSE else header <- TRUE
if (is.null(opt$miss)) opt$miss <- -99
if (is.null(opt$phe_idc)) opt$phe_idc <- 1
if (is.null(opt$rmOut)) opt$rmOut <- "miss_phe.id"
if (is.null(opt$out)) opt$out <- opt$phef  ## 注意：会覆盖输入文件

## 表型文件
if (file.exists(opt$phef)) {
    phe <- fread(opt$phef)
} else {
    cat(opt$phef, "not found!\n")
    quit()
}


## 整型变量列数
if (is.null(opt$num_int)) {
    index <- apply(phe[1, ], 1, function(x) as.integer(x) == as.numeric(x))
    opt$num_int <- which(!index)[1] - 1
}

## 增加一列辅助的均值列
if (ncol(phe) < 2) phe$mean <- 1

## 命名
names(phe)[opt$phe_idc] <- "id"
if (!is.null(opt$mean)) names(phe)[opt$mean] <- "mean"
if (!is.null(opt$breed)) names(phe)[opt$breed] <- "breed"
names(phe)[(opt$num_int + 1):ncol(phe)] <- paste0("R", 1:(ncol(phe) - opt$num_int))
phe_names <- names(phe)

## 基因型个体id
fam <- fread(opt$famf)
names(fam) <- c("fid", "id", "sire", "dam", "sex", "phe")
cat("genotyped individuals in plink:", nrow(fam), "\n")


if (all(fam$id %in% phe$id) && is.null(opt$rm)) {
    ## 所有基因型个体均出现在表型文件中
    phe <- phe[phe$id %in% fam$id, ]

    if (!is.null(opt$out)) {
        cat("phenotypes of all genotyped individuals (", nrow(fam), ") output to:", opt$out, "\n")
        write.table(phe, opt$out, row.names = FALSE, col.names = header, quote = FALSE)
    }
    quit(status = 1)
} else if (all(!fam$id %in% phe$id)) {
    ## 所有基因型个体都不出现在表型文件中
    cat("ID in genotype file and phenotype file are all inconsistent.\n")
    quit(status = 1)
}

## 合并文件
phe_g <- left_join(fam, phe, by = "id")

## 填充群体均值
if (!is.null(opt$mean)) {
    phe_g$mean <- 1
}

## 品种固定效应变量
if (!is.null(opt$breed)) {
    phe_breed <- unique(phe_g$breed)
    if (any(is.na(phe_breed))) {
        if (length(unique(phe_g$fid)) > 10) {
            cat("Please provide the breed label in the first column of the fam file!\n")
            print(head(fam, n = 3))
            quit(status = 1)
        } else {
            phe_g$breed <- as.numeric(as.factor(phe_g$fid))
        }
    }
}

## 选取原来表型的列
phe_out <- subset(phe_g, select = phe_names)

## 填充其他列缺失值
for (i in seq_len(ncol(phe_out))) {
    ## 整型缺失为0，实型为设定值
    if (i <= opt$num_int) {
        miss <- 0
    } else {
        miss <- opt$miss
    }

    ## 填充
    na_index <- is.na(phe_out[[i]])
    phe_out[[i]][na_index] <- miss
}

## 删除所有表型都缺失的个体记录
if (!is.null(opt$rm)) {
    if (opt$rm == "all") {
        all_miss_id <- apply(phe_out[, (opt$num_int + 1):ncol(phe_out)], 1, function(x) all(x == opt$miss))
        all_miss_id <- phe_out$id[all_miss_id]
        phe_out <- subset(phe_out, !id %in% all_miss_id)
        mv_id <- fam[fam$id %in% all_miss_id, 1:2]
        if (nrow(mv_id) > 0) {
            cat(nrow(mv_id), "individual IDs missing all phenotypes output to:", opt$rmOut, "\n")
            write.table(mv_id, opt$rmOut, row.names = FALSE, col.names = header, quote = FALSE)
        }
    } else if (opt$rm == "single") {
        if (is.null(opt$phec)) opt$phec <- opt$num_int + 1
        names(phe_out)[opt$phec] <- "phe"
        miss_index <- phe_out$phe == opt$miss

        if (sum(miss_index) > 0) {
            ## 表型缺失个体
            mv_id <- fam[fam$id %in% phe_out$id[miss_index], c("fid", "id")]

            ## 筛选出有表型的个体
            phe_out <- phe_out[!miss_index, ]

            cat(nrow(mv_id), "individuals missing specified phenotypes output to:", opt$rmOut, "\n")
            write.table(mv_id, opt$rmOut, row.names = FALSE, col.names = header, quote = FALSE)
        }
    }
}

## 保证表型文件中的个体排序与基因型文件famf中个体顺序一致
if (!is.null(opt$sort)) {
    phe_out <- phe_out[phe_out$id %in% fam$id, ]
}

## 输出文件
if (!is.null(opt$out)) {
    write.table(phe_out, opt$out, row.names = FALSE, col.names = header, quote = FALSE)
    cat("phenotypes of all genotyped individuals (", nrow(phe_out), ") output to: ", opt$out, "\n", sep = "")
}

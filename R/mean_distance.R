#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 计算群体间的distance均值

## 加载需要的程序包
cat("Loading required packages... \n\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))

## 参数列表
spec <- matrix(
  c(
    "prefix", "f", 1, "character", "[Required] mdist file path\n",
    "indexf", "i", 1, "character", "[Optional] family id\n",
    "out",    "o", 1, "character", "[Optional] output file name prefix [mdist_summary]\n",
    "help",   "h", 0, "logical", "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$prefix)) {
  cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## 默认参数
if (is.null(opt$out)) opt$out <- "mdist_summary"

## 文件名
dist_file <- paste0(opt$prefix, ".mdist")
if (is.null(opt$indexf)) {
  id_file <- paste0(opt$prefix, ".mdist.id")
} else {
  id_file <- opt$indexf
}

## distance距离矩阵
if (file.exists(dist_file)) {
  mdist <- fread(dist_file)
  mdist <- as.matrix(mdist)
} else {
  cat(dist_file, "not found!\n")
  quit()
}

## family ID
if (file.exists(id_file)) {
  id <- fread(id_file)
  names(id)[1] <- "fid"
} else {
  cat(id_file, "not found!\n")
  quit()
}

## 检查品种id合理性
fids <- unique(id$fid)
if (nrow(id) != nrow(mdist)) {
  cat("The number of id is not equal to the number of rows in mdist file!\n")
  quit()
} else if (length(fids) == nrow(mdist)) {
  cat("The number of family IDS is equal to the individual IDS, please check\n")
  quit()
} else if (length(fids) == 1) {
  cat("Only one family id:", fids, "\n")
  quit()
}

## 计算均值
mean <- matrix(NA, nrow = length(fids), ncol = length(fids))
for (i in seq_len(length(fids))) {
  for (j in 1:i) {
    if (i == j) {
      mean[i, j] <- NA
    } else {
      fid_i <- id$fid == fids[i]
      fid_j <- id$fid == fids[j]
      mean[i, j] <- mean(mdist[fid_i, fid_j])
    }
  }
}

## 命名
colnames(mean) <- fids
rownames(mean) <- fids

## 输出
output_name <- paste0(opt$out, ".csv")
write.csv(mean, output_name, row.names = TRUE, na = "")

cat("average Genetic distances output to:", output_name, "\n")

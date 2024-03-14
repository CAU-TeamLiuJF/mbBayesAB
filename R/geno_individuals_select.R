#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 根据QMSim输出的表型文件_data_***.txt选择基因型个体

# Load packages
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## command parameters
spec <- matrix(
  c("dataf",    "I", 1, "character", "[Required] phenotype file name/path",
    "nlitter",  "l", 1, "double",    "[Optional] numbers individuals select/litter [2]",
    "nsel",     "N", 1, "double",    "[Optional] total selected individuals number [2]",
    "gen_all",  "g", 1, "character", "[Optional] generations in genotype file. 8-10 [last one]",
    "gen_sel",  "G", 1, "character", "[Optional] generation to selected eg. 10 [last one]",
    "fid",      "F", 1, "character", "[Optional] Output ID file with FID column",
    "out",      "O", 1, "character", "[Optional] Output file name prefix [keep_geno_ids.txt]\n",
    "seed",     "s", 2, "integer",   "[Optional] [NULL]",
    "outIndex", "i", 0, "logical",   "[Optional] selected Ind index in dataf \n",    
    "help",     "h", 0, "logical",   "This is Help!"),
  byrow = TRUE, ncol = 5)
opt <- getopt::getopt(spec = spec)


## check parameters
if (!is.null(opt$help) || is.null(opt$dataf)) {
  ## print help message
  cat(paste(getopt::getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## default parameters
if (is.null(opt$out)) opt$out <- "keep_geno_ids.txt"
if (is.null(opt$nlitter)) opt$nlitter <- 2

## 随机数种子
if (is.null(opt$seed)) {
  opt$seed <- Sys.time()
  write.table(opt$seed, paste0(out, ".seed"))
}
set.seed(opt$seed)

## 读取表型文件
if (file.exists(opt$dataf)) {
  data <- fread(opt$dataf)
} else {
  cat(opt$dataf, "not found.\n")
  quit()
}

## QMSim中输出基因型个体的世代
if (!is.null(opt$gen_all)) {
  gen <- as.numeric(unlist(strsplit(opt$gen_all, "-")))
  if (length(gen) > 1) gen <- gen[1]:gen[2]
  data <- subset(data, G %in% gen)
}

## 选择的世代
if (!is.null(opt$gen_sel)) {
  gen <- as.numeric(unlist(strsplit(opt$gen_sel, "-")))
  if (length(gen) > 1) gen <- gen[1]:gen[2]
  data_gen <- subset(data, G %in% gen)
  cat("selected generation:", gen, "\n")
} else {
  gen <- unique(data$G)
  gen <- gen[length(gen)]
  data_gen <- subset(data, G %in% gen)
}
cat("number of candidates:", nrow(data_gen), "\n")
cat("number of litter:", length(unique(data_gen$Dam)), "\n")

## 每窝中选出nlitter个个体
sample_litter <- data_gen %>% group_by(G, Dam) %>% slice_sample(n = opt$nlitter)

## 只选择部分窝
if (!is.null(opt$nsel)) {
  ## 选择合适的窝数
  nlitter_per_gen <- opt$nsel %/% opt$nlitter %/% length(gen)
  dam_sel <- sample_litter %>% select(G, Dam) %>% unique %>% group_by(G) %>% slice_sample(n = nlitter_per_gen)

  ## 筛选指定窝的个体
  sample_litter <- subset(sample_litter, Dam %in% dam_sel$Dam)
}

if (is.null(opt$outIndex)) {
  ## 输出文件中是否包含fid(供plink提取个体基因型用)
  if (!is.null(opt$fid)) {
    keep <- data.frame(iid = opt$fid, iid = sample_litter$Progeny)
  } else {
    keep <- sample_litter$Progeny
  }
  cat("output selected individuals ID\n")
  nsel_final <- nrow(keep)
} else {
  ## 输出选择个体在原有文件中的位置(用一列0/1变量指示)
  keep <- as.integer(data$Progeny %in% sample_litter$Progeny)
  nsel_final <- sum(keep)
  cat("output index of selected individuals in input file\n")
}

## 输出文件
write.table(keep, opt$out, row.names = FALSE, col.names = FALSE, quote = FALSE)

cat("Total selected:", nsel_final, "\n")
cat("ids selected output to:", opt$out, "\n")

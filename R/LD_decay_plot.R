#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 根据plink的LD结果进行 LD decay 和 LD phase correlation 作图,LD结果文件过大可能会运行慢
# debug
# opt <- list(files="breedAsq breedBsq", mapf = "breedAsq.map")

## 加载需要的程序包
cat("Loading required packages... \n\n")
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("Cairo"))

## 参数列表
spec <- matrix(
  c(
    "files", "I", 1, "character", "[Required] LD results files prefix of populations\n",
    "mapf", "m", 1, "character", "[Optional] plink map file [files[1].bim]\n",
    "popN", "N", 1, "character", "[Optional] populations names [files]\n",
    "bin1", "1", 1, "double", "[Optional] The size bin for mean r^2 of Short Dist/kb [10]\n",
    "bin2", "2", 1, "double", "[Optional] The size bin for mean r^2 of Long  Dist/kb [100]\n",
    "type", "t", 1, "character", "[Optional] ld/frq [ld]\n",
    "breaks", "B", 1, "double", "[Optional] Break point to distinguish Short or Long Dist/kb [100]\n",
    "max", "M", 1, "double", "[Optional] Max distance/kb [10000]\n",
    "decayU", "U", 1, "double", "[Optional] Y-axis upper limit of LD decay/[0.8]\n",
    "decayL", "L", 1, "double", "[Optional] Y-axis lower limit of LD decay/[0]\n",
    "corU", "u", 1, "double", "[Optional] Y-axis upper limit of LD correlation/[0.6]\n",
    "corL", "l", 1, "double", "[Optional] Y-axis lower limit of LD correlation/[0]\n",
    "out", "O", 1, "character", "[Optional] output file names prefix/[files]\n",
    "dir", "D", 1, "character", "[Optional] Working directory/NULL\n",
    "plot", "P", 1, "character", "[Optional] Output plot [T/F]/T\n",
    "help", "h", 0, "logical", "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || is.null(opt$files)) {
  cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## 文件名前缀
prefixs <- unlist(strsplit(opt$files, " "))

## 默认参数
if (is.null(opt$popN)) opt$popN <- opt$files
if (is.null(opt$bin1)) opt$bin1 <- 10
if (is.null(opt$bin2)) opt$bin2 <- 100
if (is.null(opt$breaks)) opt$breaks <- 100
if (is.null(opt$max)) opt$max <- 10000
if (is.null(opt$out)) opt$out <- paste(prefixs, collapse = "_")
if (is.null(opt$plot)) opt$plot <- "T"
if (is.null(opt$type)) opt$type <- "ld"

## 工作文件夹
if (!is.null(opt$dir)) setwd(opt$dir)

## 文件名
files <- paste0(prefixs, ".", opt$type)

## 群体个数
np <- length(files)

## 群体标识符
pop_name <- unlist(strsplit(opt$popN, " "))

## 确定读取plink的ld还是frq结果
if (opt$type == "ld") {
  cols <- c(1, 2, 5, 7)
  value_cols <- 4
  by_cols <- c("BP_A", "BP_B")
} else {
  cols <- 1:5
  value_cols <- 5
  by_cols <- c("CHR", "SNP")
}

cat("Reading and merging results...\n")
for (i in 1:np) {
  ## 检查文件是否存在
  if (!file.exists(files[i])) {
    cat("file ", i, ": ", files[i], "not found.\n")
    quit()
  }

  ## 读取处理文件
  datai <- fread(files[i], select = cols)
  names(datai)[value_cols] <- pop_name[i]
  
  ## 去除pos为0的标记位点
  datai <- subset(datai, !(BP_A == 0 | BP_B == 0))

  datai$dist <- abs(datai$BP_A - datai$BP_B) / 1000
  datai <- datai[!datai$dist > opt$max, ]

  if (i == 1) {
    data <- datai

    ## 剔除大于设定值的SNP对
    if (opt$type == "ld") {
      data$dist <- abs(data$BP_A - data$BP_B) / 1000
      data <- data[!data$dist > opt$max, ]
    }
  } else {
    datai <- datai[, -c(1, 5)]
    data <- inner_join(data, datai, by = by_cols, suffix = c("_A", "_B"))
  }
}

## 去掉染色体列
data <- data[, -1]

## 生成LD作图区间
bins <- c(
  seq(0, opt$breaks, opt$bin1),
  seq(opt$breaks + opt$bin2, opt$max, opt$bin2)
)

# ## 删除R2为0的标记对
# zero <- apply(data, 1, function(x) any(x <= 0))
# R2_all <- subset(data, !zero)

## 列重排序
cols <- c("BP_A", "BP_B", "dist", pop_name)
data <- subset(data, select = cols)

## 统计两两群体之间的相关
cat("Performing descriptive statistics on LD results...\n")
mean_cor <- matrix(NA, nrow = np, ncol = np)
rownames(mean_cor) <- colnames(mean_cor) <- pop_name
for (i in 4:(ncol(data) - 1)) {
  for (j in (i + 1):ncol(data)) {
    ## 提取i、j群体数据
    colsi <- c("dist", pop_name[i - 3], pop_name[j - 3])
    cor_data <- subset(data, select = colsi)
    names(cor_data) <- c("dist", "R2_1", "R2_2")

    ## 报告总体相关
    corij <- cor(cor_data$R2_1, cor_data$R2_2)
    mean_cor[i - 3, j - 3] <- corij
    # cat("LD correlation of", colsi[2], "and", colsi[3], "is:", corij, "\n")

    ## 统计R2均值，不同群体间相关系数
    cori <- mutate(cor_data, dist = cut(dist, breaks = bins)) %>%
      group_by(dist) %>%
      summarise(R1 = mean(R2_1), R2 = mean(R2_2), cor = cor(R2_1, R2_2))

    ## 命名
    cori_names <- paste(colsi[2], colsi[3], sep = "_")
    names(cori) <- c(colsi, cori_names)

    ## 合并结果
    if (i == 4 && j == i + 1) {
      ## 第一个组，不需要合并
      cor_all <- cori

      ## 将区间表示列更换
      cor_all$dist <- bins[2:(nrow(cor_all) + 1)]

      ## 组的名字
      cor_names <- cori_names
    } else {
      ## 将新的列合并到已有结果
      cols_new <- !names(cori) %in% names(cor_all)
      cori <- subset(cori, select = cols_new)
      cor_all <- cbind(cor_all, cori)

      ## 合并组名
      cor_names <- c(cor_names, cori_names)
    }
  }
}

## 删除na
cor_all <- subset(cor_all, !is.na(dist))

## 输出统计结果
print(mean_cor)
narm <- !apply(cor_all, 1, function(x) any(is.na(x)))
write.table(cor_all[narm, ], paste0(opt$out, "_ld.txt"), row.names = FALSE, quote = FALSE)
write.table(mean_cor, paste0(opt$out, "_cor_ld.txt"), quote = FALSE)

cat("Ploting...\n")

## LD衰减作图
## 提取数据
ld_decay <- subset(cor_all, select = c("dist", pop_name))
ld_decay_p <- melt(ld_decay,
  id.vars = "dist",
  variable.name = "populations",
  value.name = "R2"
)
## 作图
plot_decay <- ggplot(ld_decay_p, aes(x = dist, y = R2, color = populations)) +
  geom_line(size = 1.2) +
  theme(
    axis.text.x = element_text(colour = "black", size = 14),
    axis.text.y = element_text(size = 14, face = "plain"),
    axis.title.x = element_text(size = 14, face = "plain"),
    axis.title.y = element_text(size = 14, face = "plain"),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 0.9),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Distance/kb") +
  ylab("LD/R2")

## y轴约束
if (!is.null(opt$decayL) && !is.null(opt$decayU)) {
  plot_decay <- plot_decay +
    scale_y_continuous(breaks = seq(0, 0.9, by = 0.1), limit = c(opt$decayL, opt$decayU))
}

## LD一致性图
## 提取数据
ld_cor <- subset(cor_all, select = c("dist", cor_names))
ld_cor_p <- melt(ld_cor,
  id.vars = "dist",
  variable.name = "populations",
  value.name = "cor"
)
ld_cor_p <- ld_cor_p[!is.na(ld_cor_p$cor), ]

## 作图
if (np > 2) {
  plot_cor_pre <- ggplot(ld_cor_p, aes(x = dist, y = cor, color = populations))
} else {
  plot_cor_pre <- ggplot(ld_cor_p, aes(x = dist, y = cor))
}

plot_cor <- plot_cor_pre +
  geom_line(size = 1.2) +
  theme(
    axis.text.x = element_text(colour = "black", size = 20),
    axis.text.y = element_text(size = 20, face = "plain"),
    axis.title.x = element_text(size = 20, face = "plain"),
    axis.title.y = element_text(size = 20, face = "plain"),
    legend.title = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", size = 0.9),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.1), limit = c(opt$corL, opt$corU)) +
  xlab("Distance/kb") +
  ylab("Pearson correlation of R2")


## y轴约束
if (!is.null(opt$corL) && !is.null(opt$corU)) {
  plot_cor <- plot_cor +
    scale_y_continuous(breaks = seq(0, 0.9, by = 0.1), limit = c(opt$corL, opt$corU))
}

## 输出
# widths = 600 * (opt$max / 1000 / 2.5)
if (opt$plot == "F") {
  theme_p
} else {
  widths <- 600 * 1.5

  prefix <- paste0(opt$out, "_decay")
  CairoPNG(paste0(prefix, ".png"), width = widths, height = 600, bg = "white")
  # png(paste0(prefix, ".png"), width = widths, height = 600, bg = "white")
  print(plot_decay)
  hide_message <- dev.off()

  prefix <- paste0(opt$out, "_correlation")
  CairoPNG(paste0(prefix, ".png"), width = widths, height = 600, bg = "white")
  # png(paste0(prefix, ".png"), width = widths, height = 600, bg = "white")
  print(plot_cor)
  hide_message <- dev.off()
}

cat("LD decay calculation and plot successfully.\n")

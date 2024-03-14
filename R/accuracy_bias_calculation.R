#!/public/home/liujf/software/program/R-4.3.1-no-dev/bin/Rscript
#!
## 用于跑完完整数据集和剔除验证群的缩减数据集dmu4后进行准确性计算

# 加载需要的程序包
suppressPackageStartupMessages(library("getopt"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))

## 命令行参数
spec <- matrix(
  c(
    "dir_all", "A", 1, "character", "[Optional] Full dataset parameter card prefix",
    "phe_all", "P", 1, "character", "[Optional] Full dataset phenotype file name",
    "tbvf", "B", 1, "character", "[Optional] file contain true breeding value",
    "ebvf", "e", 1, "character", "[Optional] file contain Estimated breeding value",
    "label", "l", 1, "character", "[Optional] breed label",
    "famf", "f", 1, "character", "[Optional] fam file whose ID order is consistent with that in EBV file",
    "ebvidc", "i", 1, "integer", "[Optional] id columns in Estimated breeding value file",
    "id_col", "C", 1, "character", "[Optional] id columns in phenotype file [1]",
    "ebv_id", "c", 1, "integer", "[Optional] id columns in true breeding value file [1]",
    "dir_val", "R", 1, "character", "[Optional] reduce(val) dataset parameter card prefix/valX",
    "fold", "F", 1, "integer", "[Optional] fold/1",
    "val_idf", "W", 1, "character", "[Optional] Validation ids file/valX_id_2col.txt",
    "ebv_col", "t", 1, "integer", "[Optional] trait rank [1]",
    "add_rnd", "n", 1, "integer", "[Optional] Random effects group number for additive genetic effects [1]",
    "tbv_col", "T", 1, "integer", "[Optional] Position of TBV in real variable columns of phe_all",
    "add_sol", "a", 1, "integer", "[Optional] code for type of additive effect in sol first column",
    "rep", "r", 1, "integer", "[Optional] number of traits [4]",
    "digit", "d", 1, "integer", "[Optional] digit of output accuracy [4]",
    "fid", "I", 1, "character", "[Optional] use to index valid in column 1 of vidf  [4]",
    "out", "O", 1, "character", "[Optional] output file name prefix/NULL\n",
    "rmNeg", "g", 0, "logical", "[Optional] rm negative value in accuracy",
    "mean", "m", 0, "logical", "[Optional] output mean of cor",
    "help", "h", 0, "logical", "This is Help!"
  ),
  byrow = TRUE, ncol = 5
)
opt <- getopt(spec = spec)

## 检查参数
if (!is.null(opt$help) || (is.null(opt$phe_all) && is.null(opt$tbvf)) ||
  (is.null(opt$dir_val) && is.null(opt$ebvf))) {
  cat(paste(getopt(spec = spec, usage = TRUERUE), "\n"))
  quit()
}

## 默认参数
if (is.null(opt$add_sol)) opt$add_sol <- 4
if (is.null(opt$fold)) opt$fold <- 1
if (is.null(opt$rep)) opt$rep <- 1
if (is.null(opt$digit)) opt$digit <- 5
if (is.null(opt$dir_val)) opt$dir_val <- "./val#val#/rep#rep#/dmu"
if (is.null(opt$val_idf)) opt$val_idf <- "./val#val#/rep#rep#/val.id"
if (is.null(opt$tbvf) && !is.null(opt$tbv_col)) opt$tbvf <- opt$phe_all
if (is.null(opt$id_col)) {
  opt$id_col <- 1
} else {
  opt$id_col <- as.numeric(opt$id_col)
}


#########################################################
##################  真实育种值/校正表型  #################
#########################################################
## 读取真实育种值
if (!is.null(opt$tbvf)) {
  if (file.exists(opt$tbvf)) {
    if (is.null(opt$tbv_col)) opt$tbv_col <- 2
    y_adj <- fread(opt$tbvf)
    names(y_adj)[opt$id_col] <- "id"
    names(y_adj)[opt$tbv_col] <- "yhat"
  } else {
    cat(opt$tbvf, "not found!\n")
    quit()
  }
} else if (!is.null(opt$dir_all)) {
  ## 从完整数据集估计ebv和residual中获得校正表型
  phe_all <- fread(opt$phe_all)
  names(phe_all)[opt$id_col] <- "id"
  sol_all <- fread(paste0(opt$dir_all, ".SOL"))
  res_all <- fread(paste0(opt$dir_all, ".RESIDUAL"))
  ## 提取ebv(群体内评估)
  ebv_all <- subset(sol_all, V1 == opt$add_sol & V2 == 1, c(5, 8))
  names(ebv_all) <- c("id", "ebv")

  ## 提取残差(多性状模型用的校正表型是群体内评估计算的，所以这里为单性状模型残差文件)
  re <- subset(res_all, select = c(1, 4))
  names(re) <- c("rows", "re")
  re$id <- phe_all$id[re$rows]

  ## 计算校正表型
  y_adj <- left_join(re, ebv_all, by = "id")
  y_adj$yhat <- y_adj$ebv + y_adj$re

  ## 输出校正表型
  if (!is.null(opt$tbvf)) {
    fwrite(y_adj[, c("id", "yhat")], opt$tbvf, col.names = FALSE, sep = "\t")
  }
}


## 结果存放
accs <- data.frame(
  rep = rep(1:opt$rep, each = opt$fold),
  fold = rep(1:opt$fold, times = opt$rep),
  cor = NA, bias = NA, rank_cor = NA, num = NA
)

## 计算准确性
for (r in 1:opt$rep) { # nolint r=5;f=4
  for (f in 1:opt$fold) {
    if (!is.null(opt$ebvf)) {
      ## 文件名
      ebvf <- gsub("#val#", f, opt$ebvf)
      ebvf <- gsub("#rep#", r, ebvf)

      ## 检查文件是否存在
      if (!file.exists(ebvf)) {
        ## 没有这种情形的结果文件，跳过
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
      } else if (file.info(ebvf)$size == 0) {
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
      }

      ## 读取估计育种值
      ebv_val <- fread(ebvf)
      if (!is.null(opt$famf)) {
        ## ebv所在列
        if (is.null(opt$ebv_col)) opt$ebv_col <- 1
        names(ebv_val)[opt$ebv_col] <- "ebv_val"

        ## 文件名
        famf <- gsub("#val#", f, opt$famf)
        famf <- gsub("#rep#", r, famf)

        ## bayes估计育种值id与fam文件中id一致
        ebv_id <- fread(famf)
        ebv_id <- ebv_id[, 1:2]
        names(ebv_id) <- c("fid", "id")
        ebv_val <- cbind(ebv_id, ebv_val)
      } else {
        if (is.null(opt$ebv_id)) opt$ebv_id <- 1
        if (is.null(opt$ebv_col)) opt$ebv_col <- 2
        names(ebv_val)[opt$ebv_id] <- "id"
        names(ebv_val)[opt$ebv_col] <- "ebv_val"
      }
    } else {
      ## dmu结果文件名
      dir_val <- gsub("#val#", f, opt$dir_val)
      dir_val <- gsub("#rep#", r, dir_val)
      solf <- paste0(dir_val, ".SOL")

      ## 检查dmu结果文件是否正常
      if (!file.exists(solf)) {
        cat(solf, "not found!\n")
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
      } else if (file.info(solf)$size == 0) {
        cat(solf, "zero size!\n")
        accs <- accs[!(accs$rep == r & accs$fold == f), ]
        next
      }

      ## SOL中性状代码
      if (is.null(opt$ebv_col)) opt$ebv_col <- 1

      ## SOL中遗传效应在随机效应第几组
      if (is.null(opt$add_rnd)) opt$add_rnd <- 1

      ## 提取验证群的ebv
      sol_val <- fread(solf)
      ebv_val <- subset(sol_val, V1 == opt$add_sol & V2 == opt$ebv_col & V4 == opt$add_rnd, c(5, 8))
      names(ebv_val) <- c("id", "ebv_val")
    }

    ## 匹配验证群ebv
    y_adj_ebv <- inner_join(y_adj, ebv_val, by = "id")

    ## 验证群id
    val_idf <- gsub("#val#", f, opt$val_idf)
    val_idf <- gsub("#rep#", r, val_idf)
    if (file.exists(val_idf)) {
      val_id <- fread(val_idf)
    } else {
      cat(val_idf, " not found!\n")
      next
    }
    ## 命名
    if (ncol(val_id) > 1) {
      names(val_id)[1:2] <- c("fid", "id")
      if (!is.null(opt$fid)) {
        val_id <- subset(val_id, fid == opt$fid)
      }
    } else {
      names(val_id) <- "id"
    }

    ## 筛选出验证群个体
    ids_index <- y_adj_ebv$id %in% val_id$id
    if (sum(ids_index) < 1) {
      cat("The number of validation individuals in val", f, "rep", r, "is 0, please check\n")
      next
    }
    y_adj_ebv_val <- subset(y_adj_ebv, ids_index)


    ## 准确性(pearson相关系数)
    accur_rf <- cor(y_adj_ebv_val$yhat, y_adj_ebv_val$ebv_val)

    ## 无偏性(回归系数系数)
    bias_rf <- lm(yhat ~ ebv_val, data = y_adj_ebv_val)$coefficients[2]

    ## 秩相关(spearman相关系数)
    rank_cor <- cor(y_adj_ebv_val$yhat, y_adj_ebv_val$ebv_val, method = "spearman")

    ## 保存结果
    accs$cor[accs$rep == r & accs$fold == f] <- round(rank_cor, opt$digit)
    accs$bias[accs$rep == r & accs$fold == f] <- round(bias_rf, opt$digit)
    accs$rank_cor[accs$rep == r & accs$fold == f] <- round(accur_rf, opt$digit)
    accs$num[accs$rep == r & accs$fold == f] <- sum(ids_index)
  }
}

if (any(!is.na(accs$cor))) {
  ## 检查结果中是否有负数
  if (any(accs$cor < 0)) {
    if (is.null(opt$rmNeg)) {
      cat("Negative values founded in result. Provide the --rmNeg parameter can delete these values\n")
    } else {
      cat("delete", sum(accs$cor < 0), "values in results\n")
      accs <- subset(accs, cor >= 0)
    }
  }

  ## 均值
  mean_cor <- paste(round(mean(accs$cor, na.rm = TRUE), opt$digit),
    round(sd(accs$cor, na.rm = TRUE), opt$digit),
    sep = "±"
  )
  mean_bias <- paste(round(mean(accs$bias, na.rm = TRUE), opt$digit),
    round(sd(accs$bias, na.rm = TRUE), opt$digit),
    sep = "±"
  )
  mean_rank_cor <- paste(round(mean(accs$rank_cor, na.rm = TRUE), opt$digit),
    round(sd(accs$rank_cor, na.rm = TRUE), opt$digit),
    sep = "±"
  )

  mean_cor <- data.frame(
    rep = "mean",
    fold = NA,
    cor = mean_cor,
    bias = mean_bias,
    rank_cor = mean_rank_cor,
    num = mean(accs$num, na.rm = TRUE)
  )
  accs$cor <- as.character(accs$cor)
  accs$bias <- as.character(accs$bias)
  accs$rep <- as.character(accs$rep)
  accs$rank_cor <- as.character(accs$rank_cor)

  accs <- bind_rows(accs, mean_cor)

  ## 打印在屏幕
  print(accs)

  ## 输出结果文件
  if (!is.null(opt$out)) fwrite(accs, opt$out, sep = "\t")
} else {
  cat("No outcome document in all cases!\n")
}

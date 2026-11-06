source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 02_build_herb_library.R
# 从平台内置 TCMSP 全库构建 药材→靶点 表（L3 伪方池的数据基础）
#
# 输入（平台自带，只读）:
#   LKFlow2/0.data/TCMSP_merged_targets.rds   31 万条 herb-compound-target
#   LKFlow2/0.data/TCMSP_merged_ingredients.rds  成分理化性质（ob/dl）
#   LKFlow2/0.data/TCMSP-uniport.rds           target → 基因 Symbol
#
# 输出:
#   data/processed/herb_targets.tsv     herb \t gene
#   data/processed/herb_stats.tsv       herb \t n_compound \t n_target
# =============================================================

suppressWarnings({
  dat_dir <- file.path(ROOT, "data", "raw")
  out_dir <- file.path(ROOT, "data", "raw")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  OB_CUT <- 30
  DL_CUT <- 0.18

  message("读取平台内置 TCMSP 数据 ...")
  tg  <- readRDS(file.path(dat_dir, "TCMSP_merged_targets.rds"))
  ing <- readRDS(file.path(dat_dir, "TCMSP_merged_ingredients.rds"))
  uni <- readRDS(file.path(dat_dir, "TCMSP-uniport.rds"))

  message(sprintf("  targets: %d 行 / %d 味药 / %d 成分",
                  nrow(tg), length(unique(tg[["中药名"]])), length(unique(tg$MOL_ID))))
  message(sprintf("  ingredients: %d 行 | uniprot 映射: %d 行", nrow(ing), nrow(uni)))

  # ---- 1. target_name → gene symbol ----
  map <- data.frame(
    target_name = trimws(as.character(uni$target_name)),
    gene        = trimws(as.character(uni[["GENE name"]])),
    stringsAsFactors = FALSE)
  map <- map[!is.na(map$gene) & map$gene != "" & !is.na(map$target_name), ]
  # 多别名取第一个
  map$gene <- toupper(sub("[,; ].*$", "", map$gene))
  map <- map[!duplicated(map$target_name), ]
  message(sprintf("  基因映射: %d 个 target_name", nrow(map)))

  # ---- 2. 成分 OB / DL ----
  phys <- data.frame(
    MOL_ID = as.character(ing$MOL_ID),
    ob     = suppressWarnings(as.numeric(ing$ob)),
    dl     = suppressWarnings(as.numeric(ing$dl)),
    stringsAsFactors = FALSE)
  phys <- phys[!duplicated(phys$MOL_ID), ]

  # ---- 3. 合并 + 过滤 ----
  x <- data.frame(
    herb = trimws(as.character(tg[["中药名"]])),
    MOL_ID = as.character(tg$MOL_ID),
    target_name = trimws(as.character(tg$target_name)),
    stringsAsFactors = FALSE)

  x <- merge(x, phys, by = "MOL_ID", all.x = TRUE)
  message(sprintf("  合并后 %d 行；OB/DL 缺失 %d 行",
                  nrow(x), sum(is.na(x$ob) | is.na(x$dl))))

  # 无 OB/DL 的成分：保守起见剔除（与平台 3.1 默认阈值一致）
  x <- x[!is.na(x$ob) & !is.na(x$dl) & x$ob >= OB_CUT & x$dl >= DL_CUT, ]
  message(sprintf("  OB>=%g & DL>=%g 过滤后: %d 行 / %d 成分 / %d 味药",
                  OB_CUT, DL_CUT, nrow(x), length(unique(x$MOL_ID)), length(unique(x$herb))))

  x <- merge(x, map, by = "target_name", all.x = FALSE)
  message(sprintf("  映射到基因后: %d 行 / %d 基因", nrow(x), length(unique(x$gene))))

  # ---- 4. 药材 → 基因集 ----
  ht <- unique(x[, c("herb", "gene")])
  ht <- ht[!is.na(ht$gene) & ht$gene != "", ]

  # 过滤太小或太大的药材（靶点 <5 无统计意义；>800 可能是录入异常）
  n_per <- table(ht$herb)
  keep <- names(n_per)[n_per >= 5 & n_per <= 800]
  dropped <- setdiff(names(n_per), keep)
  ht <- ht[ht$herb %in% keep, ]

  write.table(ht, file.path(out_dir, "herb_targets.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  st <- data.frame(herb = names(n_per), n_target = as.integer(n_per))
  st <- st[st$herb %in% keep, ]
  st <- st[order(-st$n_target), ]
  write.table(st, file.path(out_dir, "herb_stats.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message(sprintf("\n=== 完成 ==="))
  message(sprintf("药材 %d 味（剔除 %d 味：靶点<5 或 >800）", length(keep), length(dropped)))
  message(sprintf("药材-基因对 %d 条", nrow(ht)))
  message(sprintf("每味药靶点数：中位 %.0f，四分位 %.0f–%.0f，范围 %d–%d",
                  median(st$n_target), quantile(st$n_target, .25), quantile(st$n_target, .75),
                  min(st$n_target), max(st$n_target)))

  # 抽 8 味的伪方预期靶点数
  set.seed(1)
  sim <- replicate(2000, length(unique(ht$gene[ht$herb %in% sample(st$herb, 8)])))
  message(sprintf("随机 8 味伪方的靶点数：中位 %d（四分位 %d–%d）",
                  median(sim), quantile(sim, .25), quantile(sim, .75)))

  cat("\n前 10 味（按靶点数）:\n")
  print(head(st, 10), row.names = FALSE)
})

# 校验：六首方的药材是否都在库里
cat("\n=== 六方药材覆盖检查 ===\n")
f <- file.path(ROOT, "data", "raw")
if (file.exists(f)) {
  fh <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  have <- unique(ht$herb)
  miss <- setdiff(unique(fh$V2), have)
  cat("缺失药材:", if (length(miss)) paste(miss, collapse = ", ") else "无", "\n")
} else {
  cat("（formula_herbs.tsv 尚未创建，跳过）\n")
}

source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 04_build_formula_targets.R
# 从本地合并库计算六方 + 单味药的靶点集（替代平台 3.1 作为主流程）
#
# 设计：
#   1. 药材名对照表（含别名），报告未匹配药材
#   2. 每方靶点 = 其药材靶点并集
#   3. 统计：靶点数、GC 轴命中、HT 轴命中
#   4. 输出 formulas.tsv / formula_herbs.tsv / formula_summary.tsv
# =============================================================

root <- ROOT
lib     <- file.path(root, "data/processed/herb_targets_v2.tsv")
out     <- file.path(root, "data/processed")

GC <- c("NR3C1","FKBP4","HSP90AA1","HSP90AB1","NR3C2","HSD11B1","HSD11B2","CRH",
        "CRHR1","CRHR2","POMC","MC2R","SERPINA6","TSC22D3","SGK1","PER1","TDO2","TPH1","TPH2")
HT <- c("SLC6A4","TPH1","TPH2","HTR1A","HTR1B","HTR2A","HTR2C","HTR4","SLC6A3","MAOA","MAOB","DDC")

ht <- read.delim(lib, stringsAsFactors = FALSE)
message(sprintf("库: %d 味药 / %d 基因", length(unique(ht$herb)), length(unique(ht$gene))))

# ---- 六方 ----
FORMULAS <- list(
  "逍遥散"         = c("柴胡","当归","白芍","白术","茯苓","甘草","薄荷","生姜"),
  "柴胡疏肝散"     = c("柴胡","陈皮","川芎","香附","枳壳","白芍","甘草"),
  "越鞠丸"         = c("香附","川芎","苍术","栀子","神曲"),
  "酸枣仁汤"       = c("酸枣仁","甘草","知母","茯苓","川芎"),
  "天王补心丹"     = c("生地黄","人参","玄参","丹参","茯苓","远志","桔梗","五味子",
                       "当归","天冬","麦冬","柏子仁","酸枣仁"),
  "温胆汤"         = c("半夏","竹茹","枳实","陈皮","甘草","茯苓"),
  "半夏白术天麻汤" = c("半夏","白术","天麻","茯苓","陈皮","甘草","生姜","大枣")
)

# 别名 → 库内名（按需扩充）
ALIAS <- c("生地"="生地黄","生地黄"="生地黄","天门冬"="天冬","麦门冬"="麦冬",
           "酸枣仁"="酸枣仁","枣仁"="酸枣仁","制半夏"="半夏","法半夏"="半夏",
           "炙甘草"="甘草","生甘草"="甘草","云苓"="茯苓","白茯苓"="茯苓",
           "广陈皮"="陈皮","橘皮"="陈皮","川穹"="川芎","元参"="玄参")

lib_herbs <- unique(ht$herb)
norm <- function(x) { y <- ifelse(x %in% names(ALIAS), ALIAS[x], x); y }

miss_all <- character(0)
rows <- list(); fh <- list()
for (f in names(FORMULAS)) {
  hs <- norm(FORMULAS[[f]])
  ok <- hs %in% lib_herbs
  miss <- hs[!ok]
  if (length(miss)) miss_all <- c(miss_all, miss)
  hs_ok <- hs[ok]
  fh[[f]] <- data.frame(formula = f, herb = FORMULAS[[f]], stringsAsFactors = FALSE)
  tg <- unique(ht$gene[ht$herb %in% hs_ok])
  rows[[f]] <- data.frame(
    formula = f,
    n_herb = length(FORMULAS[[f]]),
    n_herb_matched = length(hs_ok),
    n_target = length(tg),
    GC_hit = sum(GC %in% tg),
    HT_hit = sum(HT %in% tg),
    GC_genes = paste(intersect(GC, tg), collapse = ","),
    HT_genes = paste(intersect(HT, tg), collapse = ","),
    stringsAsFactors = FALSE)
  attr(rows[[f]], "targets") <- tg
}

summ <- do.call(rbind, rows)
# 组装 formulas.tsv（长表）
flong <- do.call(rbind, lapply(names(rows), function(f)
  data.frame(formula = f, gene = attr(rows[[f]], "targets"), stringsAsFactors = FALSE)))
fhtab <- do.call(rbind, fh)

write.table(flong, file.path(out, "formulas.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(fhtab, file.path(out, "formula_herbs.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

s_out <- summ[, c("formula","n_herb","n_herb_matched","n_target","GC_hit","HT_hit")]
write.table(s_out, file.path(out, "formula_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\n=== 六方靶点统计 ===\n")
print(s_out, row.names = FALSE)

cat("\n=== GC 轴命中明细 ===\n")
for (i in seq_len(nrow(summ)))
  cat(sprintf("%-12s GC %2d/19: %s\n", summ$formula[i], summ$GC_hit[i], summ$GC_genes[i]))

cat("\n=== 未匹配药材 ===\n")
cat(if (length(unique(miss_all))) paste(unique(miss_all), collapse = ", ") else "无", "\n")
cat(sprintf("\n库内可用药材 %d 味；六方共用到 %d 味（匹配 %d）\n",
            length(lib_herbs), length(unique(unlist(FORMULAS))),
            length(unique(unlist(FORMULAS))) - length(unique(miss_all))))


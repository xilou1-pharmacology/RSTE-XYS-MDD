source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 05_specificity_score.R — 校正版特异性统计量
#
# 三个校正：
#   1. 频率过滤：剔除被 >UBIQ_FRAC 药材命中的"滥交靶点"
#   2. IDF 加权：罕见靶点权重高
#   3. 疾病条件化：只统计 ∩ 疾病 DEG 的基因
#
# 零分布：药材置换（同样味数 + 靶点数匹配）
# 输出：naive vs corrected 对比表 + 图
# =============================================================

root <- ROOT
UBIQ_FRAC <- 0.30      # 超过 30% 药材命中的基因视为滥交，剔除
NPERM <- 2000
set.seed(20261102)

ht <- read.delim(file.path(root, "data/processed/herb_targets_v2.tsv"), stringsAsFactors = FALSE)
herbs <- unique(ht$herb); genes_all <- unique(ht$gene)
byh <- split(ht$gene, ht$herb)
freq <- table(ht$gene) / length(herbs)
idf  <- -log(as.numeric(freq)); names(idf) <- names(freq)

GC <- c("NR3C1","FKBP4","HSP90AA1","HSP90AB1","NR3C2","HSD11B1","HSD11B2","CRH",
        "CRHR1","CRHR2","POMC","MC2R","SERPINA6","TSC22D3","SGK1","PER1","TDO2","TPH1","TPH2")
HT <- c("SLC6A4","TPH1","TPH2","HTR1A","HTR1B","HTR2A","HTR2C","HTR4","SLC6A3","MAOA","MAOB","DDC")

ubiq <- names(freq)[freq > UBIQ_FRAC]
cat(sprintf("滥交靶点(>%.0f%%药材命中): %d 个 -> %s\n", 100*UBIQ_FRAC,
            length(ubiq), paste(ubiq, collapse=", ")))
cat(sprintf("GC 轴剔除后剩: %s\n", paste(setdiff(GC, ubiq), collapse=", ")))

GC_c <- setdiff(GC, ubiq)
HT_c <- setdiff(HT, ubiq)

# 疾病集：真实数据到位前用模拟（1500 基因，含 6 个 GC 轴基因）
DEG <- NULL
deg_file <- file.path(root, "data/processed/deg_blood.tsv")
if (file.exists(deg_file)) {
  d <- read.delim(deg_file, stringsAsFactors = FALSE)
  DEG <- toupper(d[[1]])
  cat("使用真实 DEG:", length(DEG), "基因\n")
} else {
  set.seed(7)
  DEG <- unique(c(sample(genes_all, 1500), sample(GC_c, 6)))
  cat("【模拟 DEG】", length(DEG), "基因（含 6 个 GC 轴基因）\n")
}

FORMULAS <- list(
  "XYS"    = c("柴胡","当归","白芍","白术","茯苓","甘草","薄荷","生姜"),
  "CSS"    = c("柴胡","陈皮","川芎","香附","枳壳","白芍","甘草"),
  "YJF"    = c("香附","川芎","苍术","栀子","神曲"),
  "SZRT"   = c("酸枣仁","甘草","知母","茯苓","川芎"),
  "TWBXD"  = c("生地黄","人参","玄参","丹参","茯苓","远志","桔梗","五味子","当归","天冬","麦冬","柏子仁","酸枣仁"),
  "WD"     = c("半夏","竹茹","枳实","陈皮","甘草","茯苓"),
  "JT"     = c("黄连","肉桂"))

tgt <- function(hs) unique(unlist(byh[hs]))
score <- function(hs, axis, conditioned) {
  g <- tgt(hs)
  if (conditioned) g <- intersect(g, DEG)
  hit <- intersect(axis, g)
  list(n = length(hit), w = sum(idf[hit]))
}

run <- function(axis, conditioned, label) {
  out <- data.frame()
  for (nm in names(FORMULAS)) {
    hs <- FORMULAS[[nm]]; hs <- hs[hs %in% herbs]
    if (!length(hs)) next
    o <- score(hs, axis, conditioned)
    nn <- numeric(NPERM); ww <- numeric(NPERM)
    for (i in seq_len(NPERM)) {
      s <- sample(herbs, length(hs))
      r <- score(s, axis, conditioned)
      nn[i] <- r$n; ww[i] <- r$w
    }
    zn <- if (sd(nn) > 0) (o$n - mean(nn)) / sd(nn) else NA
    zw <- if (sd(ww) > 0) (o$w - mean(ww)) / sd(ww) else NA
    out <- rbind(out, data.frame(arm = label, formula = nm, n_herb = length(hs),
                                 obs_n = o$n, null_n = round(mean(nn), 2),
                                 z_count = round(zn, 2), z_idf = round(zw, 2)))
  }
  out
}

res <- rbind(
  run(GC, FALSE, "naive (all GC genes, no conditioning)"),
  run(GC_c, FALSE, "corrected A (ubiquitous removed)"),
  run(GC_c, TRUE,  "corrected B (A + disease-conditioned)"),
  run(HT_c, FALSE, "control axis (HT, corrected A)")
)

write.table(res, file.path(root, "results/specificity_corrected.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (a in unique(res$arm)) {
  cat("\n===", a, "===\n")
  s <- res[res$arm == a, ]; s <- s[order(-s$z_idf), ]
  print(s[, c("formula", "obs_n", "null_n", "z_count", "z_idf")], row.names = FALSE)
}

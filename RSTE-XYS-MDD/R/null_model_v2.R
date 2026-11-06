#!/usr/bin/env Rscript
# =============================================================
# null_model_v2.R — 四层置换框架（v2.0 核心）
#
#   L1 基因层 : 表达丰度五分位 × PPI 度数五分位 匹配
#   L2 GO 层  : 在同一 GO 条目内抽替代基因（控制基因家族偏倚）
#   L3 药材层 : 从 TCMSP 随机抽同样味数药材组"伪方"（主判定）
#   L4 参数层 : 四档 OB/DL 阈值重跑，看 z 稳定性
#
# 支持多轴：--axes gc_axis.gmt,ht_axis.gmt
# 输出：null_results.tsv / null_dist.tsv / F2_dual_axis.pdf
# =============================================================

CONFIG <- list(
  dir_in        = "data/",
  dir_out       = "results/",
  n_perm        = 1000,
  seed          = 20261102,
  n_stratum     = 5,

  formula_file  = "formulas.tsv",        # formula \t gene
  formula_herbs = "formula_herbs.tsv",   # formula \t herb
  herb_dir      = "herb_targets/",       # <herb>.tsv (col: gene)
  deg_file      = "deg_blood.tsv",
  module_file   = "wgcna_modules.tsv",
  expr_file     = "expr_blood.tsv",
  degree_file   = "ppi_degree.tsv",
  gene_go       = "gene_go.tsv",         # gene \t go_id
  axes          = c(gc = "gc_axis.gmt", ht = "ht_axis.gmt"),
  lfc_cut       = 0.5,
  padj_cut      = 0.05,
  min_module    = 30,
  layers        = c("L1", "L2", "L3")    # L4 由独立阈值文件驱动
)

set.seed(CONFIG$seed)
dir.create(CONFIG$dir_out, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) if (file.exists(f)) read.delim(f, stringsAsFactors = FALSE) else NULL

# ---------- 数据 ----------
formulas <- read.delim(file.path(CONFIG$dir_in, CONFIG$formula_file),
                       header = FALSE, col.names = c("formula", "gene"))
formulas$gene <- toupper(formulas$gene)

deg <- rd(file.path(CONFIG$dir_in, CONFIG$deg_file))
colnames(deg)[1:3] <- c("gene", "logFC", "adj.P")
deg$gene <- toupper(deg$gene)

mod <- rd(file.path(CONFIG$dir_in, CONFIG$module_file))
colnames(mod)[1:2] <- c("gene", "module"); mod$gene <- toupper(mod$gene)
tab <- table(mod$module)
keep <- names(tab)[tab >= CONFIG$min_module & names(tab) != "grey"]
mod_genes <- unique(mod$gene[mod$module %in% keep])

deg_sig <- deg$gene[abs(deg$logFC) >= CONFIG$lfc_cut & deg$adj.P < CONFIG$padj_cut]
disease <- intersect(deg_sig, mod_genes)

expr <- read.delim(file.path(CONFIG$dir_in, CONFIG$expr_file),
                   row.names = 1, check.names = FALSE)
expr_mean <- rowMeans(as.matrix(expr), na.rm = TRUE)

degree <- setNames(rep(0, length(rownames(expr))), toupper(rownames(expr)))
d <- rd(file.path(CONFIG$dir_in, CONFIG$degree_file))
if (!is.null(d)) { colnames(d)[1:2] <- c("gene","degree"); d$gene <- toupper(d$gene)
                   degree[d$gene] <- d$degree }

bg <- intersect(toupper(rownames(expr)), deg$gene)
message(sprintf("DEG=%d | 模块基因=%d | 疾病集=%d | 背景=%d",
                length(deg_sig), length(mod_genes), length(disease), length(bg)))

# ---------- 分层工具 ----------
qbin <- function(x, q = CONFIG$n_stratum) cut(rank(x, ties.method = "first"), q, labels = FALSE)
expr_bin <- qbin(expr_mean[bg]); deg_bin <- qbin(degree[bg])
stratum  <- setNames(paste(expr_bin, deg_bin, sep = "_"), bg)

# L2 GO 索引
go <- rd(file.path(CONFIG$dir_in, CONFIG$gene_go))
go_list <- NULL
if (!is.null(go)) {
  colnames(go)[1:2] <- c("gene","go"); go$gene <- toupper(go$gene)
  go <- go[go$gene %in% bg, ]
  go_list <- split(go$gene, go$go)
  gene2go <- split(go$go, go$gene)
  message(sprintf("GO 条目 %d，覆盖基因 %d", length(go_list), length(gene2go)))
}

# 药材靶点
herb_targets <- list()
if (dir.exists(file.path(CONFIG$dir_in, CONFIG$herb_dir))) {
  files <- list.files(file.path(CONFIG$dir_in, CONFIG$herb_dir), pattern = "\\.tsv$", full.names = TRUE)
  for (f in files) {
    nm <- sub("\\.tsv$", "", basename(f))
    t <- read.delim(f, stringsAsFactors = FALSE)
    herb_targets[[nm]] <- unique(toupper(t[[1]]))
  }
  message(sprintf("药材靶点 %d 味", length(herb_targets)))
}

formula_herbs <- read.delim(file.path(CONFIG$dir_in, CONFIG$formula_herbs),
                            stringsAsFactors = FALSE, header = FALSE)
colnames(formula_herbs) <- c("formula", "herb")

# ---------- 轴 ----------
read_gmt <- function(f) lapply(strsplit(readLines(f), "\t"), function(x) toupper(unique(x[-(1:2)])))
axes <- list()
for (nm in names(CONFIG$axes)) {
  p <- file.path(CONFIG$dir_in, CONFIG$axes[[nm]])
  if (file.exists(p)) axes[[nm]] <- unique(unlist(read_gmt(p)))
}

stat_of <- function(set) {
  list(overlap = length(intersect(set, disease)),
       axis = sapply(axes, function(a) length(intersect(intersect(set, disease), a))))
}

# ---------- L1：分层基因抽样 ----------
sample_L1 <- function(targets) {
  pool_by <- split(bg, stratum)
  st <- stratum[targets]; st[is.na(st)] <- names(pool_by)[1]
  out <- unlist(lapply(names(table(st)), function(s) {
    k <- table(st)[[s]]; p <- pool_by[[s]]
    if (is.null(p) || length(p) < k) p <- bg
    sample(p, k, replace = length(p) < k)
  }), use.names = FALSE)
  out
}

# ---------- L2：GO 内抽样 ----------
sample_L2 <- function(targets) {
  if (is.null(go_list)) return(NULL)
  out <- character(0)
  for (g in targets) {
    terms <- gene2go[[g]]
    pool <- if (is.null(terms)) bg else unique(unlist(go_list[terms]))
    pool <- setdiff(pool, out)
    out <- c(out, sample(pool, 1))
  }
  out
}

# ---------- L3：伪方抽样 ----------
sample_L3 <- function(k_herb, n_target_obs) {
  if (length(herb_targets) < k_herb) return(NULL)
  hs <- sample(names(herb_targets), k_herb)
  u <- unique(unlist(herb_targets[hs]))
  u <- intersect(u, bg)
  if (length(u) > n_target_obs) u <- sample(u, n_target_obs)
  u
}

# ---------- 主循环 ----------
rows <- list(); dists <- list()
for (f in unique(formulas$formula)) {
  tg <- intersect(formulas$gene[formulas$formula == f], bg)
  if (length(tg) < 10) { message("[skip] ", f); next }
  obs <- stat_of(tg)
  k_herb <- sum(formula_herbs$formula == f)
  message(sprintf("[%s] 靶点=%d 药材=%d overlap=%d", f, length(tg), k_herb, obs$overlap))

  for (lyr in CONFIG$layers) {
    sim <- switch(lyr,
      L1 = function() sample_L1(tg),
      L2 = function() sample_L2(tg),
      L3 = function() sample_L3(k_herb, length(tg)))
    if (lyr == "L2" && is.null(go_list)) next

    null_ov <- integer(CONFIG$n_perm); null_ax <- NULL
    ok <- 0
    for (i in seq_len(CONFIG$n_perm)) {
      s <- tryCatch(sim(), error = function(e) NULL)
      if (is.null(s) || length(s) < 10) next
      ok <- ok + 1
      st <- stat_of(s)
      null_ov[i] <- st$overlap
      if (is.null(null_ax)) null_ax <- matrix(0, 0, length(axes))
      null_ax <- rbind(null_ax, as.numeric(st$axis))
    }
    if (ok < 50) { message("  [", lyr, "] 有效置换不足"); next }
    null_ov <- null_ov[seq_len(ok)]; null_ax <- null_ax[seq_len(ok), , drop = FALSE]

    z <- (obs$overlap - mean(null_ov)) / sd(null_ov)
    p <- (1 + sum(null_ov >= obs$overlap)) / (1 + ok)

    for (j in seq_along(axes)) {
      axn <- names(axes)[j]
      zo <- (obs$axis[axn] - mean(null_ax[, j])) / sd(null_ax[, j])
      po <- (1 + sum(null_ax[, j] >= obs$axis[axn])) / (1 + ok)
      rows[[length(rows) + 1]] <- data.frame(
        formula = f, layer = lyr, axis = axn, n_target = length(tg),
        obs = as.integer(obs$axis[axn]), null_mean = round(mean(null_ax[, j]), 2),
        z = round(zo, 3), p = signif(po, 3), n_perm_ok = ok)
      dists[[length(dists) + 1]] <- data.frame(
        formula = f, layer = lyr, axis = axn, null_value = null_ax[, j],
        observed = obs$axis[axn])
    }
    # 总 overlap（不分轴）
    rows[[length(rows) + 1]] <- data.frame(
      formula = f, layer = lyr, axis = "overall", n_target = length(tg),
      obs = obs$overlap, null_mean = round(mean(null_ov), 2),
      z = round(z, 3), p = signif(p, 3), n_perm_ok = ok)
  }
}

res <- do.call(rbind, rows); dst <- do.call(rbind, dists)

# FDR：按层 × 轴分别校正
res$p_FDR <- NA
for (lyr in unique(res$layer)) for (ax in unique(res$axis)) {
  idx <- res$layer == lyr & res$axis == ax
  res$p_FDR[idx] <- signif(p.adjust(res$p[idx], "BH"), 3)
}

write.table(res, file.path(CONFIG$dir_out, "null_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(dst, file.path(CONFIG$dir_out, "null_dist.tsv"),    sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- F2：双轴 + 三层零分布 ----------
pdf(file.path(CONFIG$dir_out, "F2_dual_axis.pdf"), width = 11, height = 7)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (axn in names(axes)) {
  sub <- res[res$axis == axn & res$layer == "L3", ]
  sub <- sub[order(-sub$z), ]
  if (nrow(sub) == 0) next
  bp <- barplot(sub$z, names.arg = sub$formula, las = 2, col = "steelblue",
                main = paste0(axn, " axis | L3 herb-swap"), ylab = "z-score")
  abline(h = 2, lty = 2, col = "red")
  text(bp, sub$z, labels = sprintf("%.2f", sub$z), pos = 3, cex = 0.7)
}
for (lyr in c("L1", "L2", "L3")) {
  s <- dst[dst$layer == lyr & dst$axis == names(axes)[1], ]
  if (nrow(s) == 0) next
  hist(s$null_value, breaks = 40, col = "grey85", border = "grey60",
       main = paste0(names(axes)[1], " | null ", lyr), xlab = "axis gene count (null)")
  abline(v = unique(s$observed), col = "red", lwd = 2)
}
dev.off()

cat("\n=== 完成 ===\n")
print(head(res[order(res$p), ], 12))
cat("\n主判定看 layer=L3；L1/L2 为敏感性\n")

# 17_mr_sensitivity.R — cis-MR 敏感性分析
# IVW / 加权中位数 / MR-Egger(+截距) / Cochran Q / 留一法

suppressPackageStartupMessages({
  library(MendelianRandomization); library(data.table)
})

root <- "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
d <- fread(file.path(root, "results/mr_persnp.tsv"), data.table = FALSE)
cat("SNP:", nrow(d), "| 基因:", length(unique(d$gene)), "\n")

res <- list()
for (g in unique(d$gene)) {
  s <- d[d$gene == g, ]
  n <- nrow(s)
  if (n < 1) next
  bxy <- s$by; bxx <- s$bx; sxx <- s$sx; syy <- s$sy

  row <- data.frame(gene = g, nSNP = n, method = "Wald ratio",
                    b = NA, se = NA, P = NA, eg_icpt = NA, eg_icpt_P = NA,
                    Q = NA, Q_P = NA, loo_max = NA)

  if (n == 1) {
    row$b <- bxy / bxx
    row$se <- abs(syy / bxx)
    row$P <- 2 * pnorm(-abs(row$b / row$se))
    res[[g]] <- row
    next
  }

  mri <- mr_input(bx = bxx, bxse = sxx, by = bxy, byse = syy)
  out <- data.frame(gene = character(), nSNP = integer(), method = character(),
                    b = numeric(), se = numeric(), P = numeric(),
                    eg_icpt = numeric(), eg_icpt_P = numeric(),
                    Q = numeric(), Q_P = numeric(), loo_max = numeric())

  ivw <- mr_ivw(mri)
  out <- rbind(out, data.frame(gene = g, nSNP = n, method = "IVW",
                               b = ivw@Estimate, se = ivw@StdError, P = ivw@Pvalue,
                               eg_icpt = NA, eg_icpt_P = NA,
                               Q = ivw@Heter.Stat[1], Q_P = ivw@Heter.Stat[2], loo_max = NA))

  if (n >= 3) {
    wm <- mr_median(mri)
    out <- rbind(out, data.frame(gene = g, nSNP = n, method = "Weighted median",
                                 b = wm@Estimate, se = wm@StdError, P = wm@Pvalue,
                                 eg_icpt = NA, eg_icpt_P = NA, Q = NA, Q_P = NA, loo_max = NA))
    eg <- mr_egger(mri)
    out <- rbind(out, data.frame(gene = g, nSNP = n, method = "MR-Egger",
                                 b = eg@Estimate, se = eg@StdError.Est, P = eg@Pvalue.Est,
                                 eg_icpt = eg@Intercept, eg_icpt_P = eg@Pvalue.Int,
                                 Q = eg@Heter.Stat[1], Q_P = eg@Heter.Stat[2], loo_max = NA))
  }

  ## 留一法：剔除单 SNP 后 IVW 的最大相对变化
  if (n >= 3) {
    base <- ivw@Estimate
    deltas <- c()
    for (k in seq_len(n)) {
      mk <- mr_input(bx = bxx[-k], bxse = sxx[-k], by = bxy[-k], byse = syy[-k])
      ek <- tryCatch(mr_ivw(mk)@Estimate, error = function(e) NA)
      if (!is.na(ek) && base != 0) deltas <- c(deltas, abs((ek - base) / base))
    }
    out$loo_max[out$method == "IVW"] <- if (length(deltas)) max(deltas) else NA
  }

  res[[g]] <- out
}

R <- do.call(rbind, res)
R$OR <- round(exp(R$b), 3)
num <- c("b","se","eg_icpt","Q")
for (c in num) R[[c]] <- ifelse(is.na(R[[c]]), NA, round(R[[c]], 4))
R$P <- ifelse(is.na(R$P), NA, signif(R$P, 3))
R$eg_icpt_P <- ifelse(is.na(R$eg_icpt_P), NA, signif(R$eg_icpt_P, 3))
R$Q_P <- ifelse(is.na(R$Q_P), NA, signif(R$Q_P, 3))
R$loo_max <- ifelse(is.na(R$loo_max), NA, round(R$loo_max, 3))

R <- R[order(R$P), c("gene","nSNP","method","OR","b","se","P","eg_icpt","eg_icpt_P","Q","Q_P","loo_max")]
write.table(R, file.path(root, "results/mr_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

cat("\n=== 敏感性分析结果 ===\n")
print(R[R$gene %in% c("GSK3B","PPARG","NFKB1","JUN"), ], row.names = FALSE)
cat("\n完成 -> results/mr_sensitivity.tsv\n")

# 14_validate_GSE76826.R — 白细胞队列方向验证（组织匹配）
# GSE76826: leukocytes, GPL17077, 12 healthy vs 10 depression (排除 10 缓解期)

suppressPackageStartupMessages({ library(limma); library(data.table) })

root <- "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
rd   <- file.path(root, "data/raw/GSE76826")
outd <- file.path(root, "results")

lines <- readLines(gzfile(file.path(rd, "GSE76826_series_matrix.txt.gz")), warn = FALSE)
b <- grep("!series_matrix_table_begin", lines); e <- grep("!series_matrix_table_end", lines)
tab <- fread(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"), data.table = FALSE)
rownames(tab) <- tab[[1]]; tab <- tab[, -1]
cat("探针 × 样本:", nrow(tab), "×", ncol(tab), "\n")

diag <- state <- NULL
for (ln in grep("^!Sample_characteristics_ch1", lines, value = TRUE)) {
  v <- gsub('"', "", strsplit(ln, "\t")[[1]][-1])
  if (length(v) && grepl("^diagnosis:", v[1])) diag <- sub("diagnosis: ", "", v)
  if (length(v) && grepl("^state:", v[1])) state <- sub("state: ", "", v)
}
gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", lines, value = TRUE)[1], "\t")[[1]][-1])
keep <- (diag == "Healthy") | (diag == "Major depressive disorder" & state == "Depression")
cat("纳入:", sum(keep), "例 | 分组:\n"); print(table(diag[keep]))
tab <- tab[, gsm[keep]]
grp <- ifelse(diag[keep] == "Healthy", "HC", "MDD")

map <- fread(file.path(root, "data/processed/GPL17077_probe2gene.tsv"), data.table = FALSE, sep = "\t")
map <- map[map$probe %in% rownames(tab), ]
expr <- as.matrix(tab[map$probe, , drop = FALSE])
rownames(expr) <- toupper(map$symbol)
expr <- expr[order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE), ]
expr <- expr[!duplicated(rownames(expr)), ]
if (max(expr, na.rm = TRUE) > 50) { expr <- log2(expr + 1); cat("已 log2\n") }
expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), ]
cat("基因数:", nrow(expr), "\n")

design <- model.matrix(~ factor(grp, levels = c("HC", "MDD")))
fit <- lmFit(expr, design); fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
cat("FDR<0.05:", sum(tt$adj.P.Val < 0.05), "| P<0.05:", sum(tt$P.Value < 0.05), "\n")
write.table(tt[, c("gene","logFC","P.Value","adj.P.Val")],
            file.path(outd, "deg_GSE76826.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

disc <- fread(file.path(root, "data/processed/deg_128vs64.tsv"), data.table = FALSE)
disc$gene <- toupper(disc$gene)
mg <- merge(tt[, c("gene","logFC")], disc[, c("gene","logFC")], by = "gene",
            suffixes = c("_val","_disc"))
tgt <- toupper(trimws(fread(file.path(root, "data/processed/XYS_targets_OB30DL018.tsv"), header = FALSE)[[1]]))
mg$is_target <- mg$gene %in% tgt

sp <- function(a, b) {
  ct <- suppressWarnings(cor.test(a, b, method = "spearman"))
  c(rho = ct$estimate, p = ct$p.value)
}
r_all <- sp(mg$logFC_val, mg$logFC_disc)
sub <- mg[mg$is_target, ]; r_tgt <- sp(sub$logFC_val, sub$logFC_disc)
ns <- mg[!mg$is_target, ]; r_nt <- sp(ns$logFC_val, ns$logFC_disc)

cat("\n=== 方向一致性 (Spearman) ===\n")
cat(sprintf("  全部 %d 基因      rho = %+.3f, P = %.3g\n", nrow(mg), r_all["rho"], r_all["p"]))
cat(sprintf("  XYS 靶点 %d 基因  rho = %+.3f, P = %.3g\n", nrow(sub), r_tgt["rho"], r_tgt["p"]))
cat(sprintf("  非靶点 %d 基因    rho = %+.3f\n", nrow(ns), r_nt["rho"]))
cat(sprintf("  靶点同向比例: %.1f%% | 非靶点: %.1f%%\n",
            100*mean(sign(sub$logFC_val)==sign(sub$logFC_disc)),
            100*mean(sign(ns$logFC_val)==sign(ns$logFC_disc))))

res <- data.frame(
  comparison = c("all_genes","XYS_targets","non_targets"),
  n = c(nrow(mg), nrow(sub), nrow(ns)),
  spearman_rho = round(c(r_all["rho"], r_tgt["rho"], r_nt["rho"]), 4),
  P = signif(c(r_all["p"], r_tgt["p"], NA), 3),
  same_direction_pct = round(100*c(mean(sign(mg$logFC_val)==sign(mg$logFC_disc)),
                                   mean(sign(sub$logFC_val)==sign(sub$logFC_disc)),
                                   mean(sign(ns$logFC_val)==sign(ns$logFC_disc))), 1))
write.table(res, file.path(outd, "validation_direction_GSE76826.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

## marker 检查
chk <- c("CD8A","CD3D","CD14","LYZ","IL1B","TNF","S100A8","S100A9","GZMK","CD19")
cat("\n=== 关键基因方向 ===\n")
for (g in chk) {
  a <- sub[sub$gene == g, ]; b2 <- mg[mg$gene == g, ]
  if (nrow(b2)) cat(sprintf("  %-8s 发现 %+7.3f | 验证 %+7.3f\n", g, b2$logFC_disc[1], b2$logFC_val[1]))
}
cat("\n完成 -> results/validation_direction_GSE76826.tsv\n")

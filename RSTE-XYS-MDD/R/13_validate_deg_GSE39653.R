# 13_validate_deg_GSE39653.R — 基因层面方向验证
# 用 limma 在 GSE39653 (PBMC, 21 MDD / 24 HC) 做差异表达，
# 再与发现队列 (GSE98793, 全血, 128/64) 的 logFC 做相关。

suppressPackageStartupMessages({ library(limma); library(data.table) })

root <- "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
rd   <- file.path(root, "data/raw/GSE39653")
outd <- file.path(root, "results")

## ---- 表达矩阵 ----
lines <- readLines(gzfile(file.path(rd, "GSE39653_series_matrix.txt.gz")), warn = FALSE)
b <- grep("!series_matrix_table_begin", lines); e <- grep("!series_matrix_table_end", lines)
tab <- fread(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"), data.table = FALSE)
rownames(tab) <- tab[[1]]; tab <- tab[, -1]

## ---- 分组 ----
dis <- gsub('"', "", strsplit(grep("^!Sample_characteristics_ch1", lines, value = TRUE)[1], "\t")[[1]][-1])
dis <- sub("disease: ", "", dis)
gsm <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", lines, value = TRUE)[1], "\t")[[1]][-1])
keep <- dis %in% c("major depressive disorder", "healthy control")
tab <- tab[, gsm[keep]]; grp <- ifelse(dis[keep] == "major depressive disorder", "MDD", "HC")
cat("样本:", length(grp), "|", paste(names(table(grp)), table(grp), collapse = " / "), "\n")

## ---- 注释 ----
map <- fread(file.path(root, "data/processed/GPL10558_probe2gene.tsv"),
             data.table = FALSE, sep = "\t")
map <- map[map$probe %in% rownames(tab), ]
expr <- as.matrix(tab[map$probe, , drop = FALSE])
rownames(expr) <- toupper(map$symbol)
expr <- expr[order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE), ]
expr <- expr[!duplicated(rownames(expr)), ]
if (max(expr, na.rm = TRUE) > 50) expr <- log2(expr + 1)
expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), ]
cat("基因数:", nrow(expr), "\n")

## ---- limma ----
design <- model.matrix(~ factor(grp, levels = c("HC", "MDD")))
fit <- lmFit(expr, design); fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
cat("FDR<0.05:", sum(tt$adj.P.Val < 0.05), "| P<0.05:", sum(tt$P.Value < 0.05), "\n")

write.table(tt[, c("gene", "logFC", "P.Value", "adj.P.Val")],
            file.path(outd, "deg_GSE39653.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 与发现队列相关 ----
disc <- fread(file.path(root, "data/processed/deg_128vs64.tsv"), data.table = FALSE)
disc$gene <- toupper(disc$gene)
mg <- merge(tt[, c("gene", "logFC")], disc[, c("gene", "logFC")],
            by = "gene", suffixes = c("_val", "_disc"))
cat("两队列共有基因:", nrow(mg), "\n")

tgt <- toupper(trimws(fread(file.path(root, "data/processed/XYS_targets_OB30DL018.tsv"),
                            header = FALSE)[[1]]))
mg$is_target <- mg$gene %in% tgt

## 全部基因
r_all <- cor(mg$logFC_val, mg$logFC_disc, method = "spearman")
p_all <- cor.test(mg$logFC_val, mg$logFC_disc, method = "spearman")$p.value
## 靶点基因
sub <- mg[mg$is_target, ]
r_tgt <- cor(sub$logFC_val, sub$logFC_disc, method = "spearman")
p_tgt <- cor.test(sub$logFC_val, sub$logFC_disc, method = "spearman")$p.value
## 非靶点基因（对照）
ns <- mg[!mg$is_target, ]
r_nt <- cor(ns$logFC_val, ns$logFC_disc, method = "spearman")

cat(sprintf("\n=== 方向一致性 (Spearman) ===\n"))
cat(sprintf("  全部 %d 基因        rho = %+.3f, P = %.3g\n", nrow(mg), r_all, p_all))
cat(sprintf("  XYS 靶点 %d 基因    rho = %+.3f, P = %.3g\n", nrow(sub), r_tgt, p_tgt))
cat(sprintf("  非靶点 %d 基因      rho = %+.3f\n", nrow(ns), r_nt))

## 靶点集内方向一致比例
conc <- mean(sign(sub$logFC_val) == sign(sub$logFC_disc))
cat(sprintf("  靶点基因同向比例: %.1f%%\n", 100 * conc))
cat(sprintf("  非靶点同向比例:   %.1f%%\n", 100 * mean(sign(ns$logFC_val) == sign(ns$logFC_disc))))

res <- data.frame(
  comparison = c("all_genes", "XYS_targets", "non_targets"),
  n = c(nrow(mg), nrow(sub), nrow(ns)),
  spearman_rho = round(c(r_all, r_tgt, r_nt), 4),
  P = signif(c(p_all, p_tgt, NA), 3),
  same_direction_pct = round(100 * c(mean(sign(mg$logFC_val) == sign(mg$logFC_disc)), conc,
                                     mean(sign(ns$logFC_val) == sign(ns$logFC_disc))), 1)
)
write.table(res, file.path(outd, "validation_direction_GSE39653.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n完成 -> results/validation_direction_GSE39653.tsv\n")

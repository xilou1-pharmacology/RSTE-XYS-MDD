# 12_validate_GSE39653.R — 独立验证队列
# GSE39653: PBMC, 21 MDD / 24 HC (8 BD 排除), GPL10558
# 检验：XYS 靶点集分数、T 细胞、CD8、单核细胞 —— 方向是否与发现队列一致

suppressPackageStartupMessages({
  library(GSVA); library(GSEABase); library(data.table)
})

root <- "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
rd   <- file.path(root, "data/raw/GSE39653")
outd <- file.path(root, "results")

## ---------- 1. 表达矩阵 ----------
lines <- readLines(gzfile(file.path(rd, "GSE39653_series_matrix.txt.gz")), warn = FALSE)
b <- grep("!series_matrix_table_begin", lines); e <- grep("!series_matrix_table_end", lines)
tab <- fread(text = paste(lines[(b + 1):(e - 1)], collapse = "\n"), data.table = FALSE)
rownames(tab) <- tab[[1]]; tab <- tab[, -1]
cat("探针 × 样本:", nrow(tab), "×", ncol(tab), "\n")

## ---------- 2. 分组 ----------
meta <- grep("^!Sample_characteristics_ch1", lines, value = TRUE)[1]
dis  <- strsplit(meta, "\t")[[1]][-1]; dis <- gsub('"', "", dis); dis <- sub("disease: ", "", dis)
gsm  <- gsub('"', "", strsplit(grep("^!Sample_geo_accession", lines, value = TRUE)[1], "\t")[[1]][-1])
grp  <- data.frame(gsm = gsm, disease = dis, stringsAsFactors = FALSE)
print(table(grp$disease))
keep <- grp$disease %in% c("major depressive disorder", "healthy control")
grp  <- grp[keep, ]; tab <- tab[, grp$gsm, drop = FALSE]
grp$group <- ifelse(grp$disease == "major depressive disorder", "MDD", "HC")
cat("纳入分析:", nrow(grp), "例 |", table(grp$group), "\n")

## ---------- 3. 探针注释 ----------
map <- fread(file.path(root, "data/processed/GPL10558_probe2gene.tsv"),
             data.table = FALSE, sep = "\t", header = TRUE)
cat("注释列: probe | symbol\n")
map <- map[!is.na(map$symbol) & map$symbol != "" & map$symbol != "---", ]
map$symbol <- toupper(map$symbol)
map <- map[map$probe %in% rownames(tab), ]
cat("可映射探针:", nrow(map), "\n")

## ---------- 4. 折叠到基因（取平均表达最高的探针）----------
expr <- as.matrix(tab[map$probe, , drop = FALSE])
rownames(expr) <- map$symbol
ord <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
expr <- expr[ord, ]
expr <- expr[!duplicated(rownames(expr)), ]
cat("基因数:", nrow(expr), "\n")

## 是否需要 log2
rng <- range(expr, na.rm = TRUE)
cat("表达范围:", round(rng[1], 2), "-", round(rng[2], 2), "\n")
if (rng[2] > 50) { expr <- log2(expr + 1); cat("已做 log2 转换\n") }
expr <- expr[apply(expr, 1, function(x) all(is.finite(x))), ]

## ---------- 5. 基因集 ----------
tgt <- toupper(trimws(fread(file.path(root, "data/processed/XYS_targets_OB30DL018.tsv"),
                            header = FALSE)[[1]]))
markers <- list(
  XYS_targets = intersect(tgt, rownames(expr)),
  T_cell = intersect(c("CD3D","CD3E","CD3G","CD2","CD28","LCK","ITK","CD6","TRAC","TRBC2"), rownames(expr)),
  CD8_T  = intersect(c("CD8A","CD8B","GZMK","GZMA","PRF1","NKG7","CD3D","CD3E"), rownames(expr)),
  Monocyte = intersect(c("CD14","LYZ","FCGR3A","CSF1R","S100A8","S100A9","VCAN","FCN1","CD68"), rownames(expr))
)
cat("\n基因集大小:\n"); print(sapply(markers, length))

## ---------- 6. ssGSEA ----------
param <- ssgseaParam(expr, markers, normalize = TRUE)
sc <- gsva(param, verbose = FALSE)
cat("ssGSEA 完成:", dim(sc), "\n")

## ---------- 7. 统计 ----------
res <- data.frame()
for (s in rownames(sc)) {
  a <- sc[s, grp$group == "MDD"]; b2 <- sc[s, grp$group == "HC"]
  d <- mean(a) - mean(b2)
  p <- wilcox.test(a, b2)$p.value
  set.seed(1)
  bs <- replicate(2000, {
    ia <- sample(a, replace = TRUE); ib <- sample(b2, replace = TRUE)
    mean(ia) - mean(ib)
  })
  ci <- quantile(bs, c(0.025, 0.975))
  res <- rbind(res, data.frame(signature = s, n_gene = length(markers[[s]]),
                               mean_MDD = round(mean(a), 4), mean_HC = round(mean(b2), 4),
                               diff = round(d, 4), CI_lo = round(ci[1], 4), CI_hi = round(ci[2], 4),
                               P = signif(p, 3)))
}
res$FDR <- signif(p.adjust(res$P, "BH"), 3)
print(res)
write.table(res, file.path(outd, "validation_GSE39653.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(grp, file.path(outd, "validation_GSE39653_samples.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n完成 -> results/validation_GSE39653.tsv\n")

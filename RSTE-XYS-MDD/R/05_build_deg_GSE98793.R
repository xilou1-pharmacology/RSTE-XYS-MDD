source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 05_build_deg_GSE98793.R
# 从 series matrix 解析 → 注释 → 分组 → limma → DEG
#
# 分组：64 纯 MDD（anxiety=no） vs 64 健康对照
# 输出：
#   data/processed/expr_blood.tsv   (gene × sample, log2)
#   data/processed/deg_blood.tsv    (gene, logFC, AveExpr, t, P, adj.P)
#   data/processed/pheno_GSE98793.tsv
# =============================================================

suppressPackageStartupMessages(library(limma))
root <- ROOT
raw  <- file.path(root, "data/raw"); out <- file.path(root, "data/processed")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) cat(sprintf(...), "\n")

# ---------- 1. 读 series matrix 的 pheno ----------
con <- gzfile(file.path(raw, "GSE98793_series_matrix.txt.gz"), "rt")
lines <- readLines(con, warn = FALSE); close(con)
msg("series matrix 行数: %d", length(lines))

grab <- function(tag) {
  ln <- grep(paste0("^", tag), lines, value = TRUE)
  if (!length(ln)) return(NULL)
  lapply(ln, function(x) {
    v <- strsplit(x, "\t")[[1]][-1]
    gsub('^"|"$', "", v)
  })
}
gsm  <- grab("!Sample_geo_accession")[[1]]
ttl  <- grab("!Sample_title")[[1]]
chars <- grab("!Sample_characteristics_ch1")
msg("样本数: %d | characteristics 行数: %d", length(gsm), length(chars))

# 把 "key: value" 拆成字段
kv <- list()
for (row in chars) {
  key <- sub(":.*$", "", row[1])
  val <- sub("^[^:]*:\\s*", "", row)
  kv[[key]] <- val
}
msg("字段: %s", paste(names(kv), collapse = " | "))

# ---------- 2. 分组 ----------
group_raw <- kv[["subject group"]]
anxiety   <- kv[["anxiety"]]
if (is.null(group_raw)) group_raw <- kv[[grep("group", names(kv), value = TRUE)[1]]]
if (is.null(anxiety))   anxiety   <- kv[[grep("anxiety", names(kv), value = TRUE)[1]]]

pheno <- data.frame(gsm = gsm, title = ttl, group = group_raw,
                    anxiety = if (is.null(anxiety)) NA else anxiety,
                    stringsAsFactors = FALSE)
msg("\n分组分布:\n"); print(table(pheno$group, pheno$anxiety, useNA = "ifany"))

# 纯抑郁 = MDD 且 anxiety=no；对照 = CNTL 且 anxiety=no
is_mdd  <- grepl("MDD|case|depress", pheno$group, ignore.case = TRUE)
is_cntl <- grepl("CNTL|control|healthy", pheno$group, ignore.case = TRUE)
no_anx  <- is.na(pheno$anxiety) | grepl("^no$", pheno$anxiety, ignore.case = TRUE)

sel <- (is_mdd & no_anx) | (is_cntl & no_anx)
pheno$label <- NA
pheno$label[is_mdd & no_anx]  <- "MDD"
pheno$label[is_cntl & no_anx] <- "CTRL"
msg("入选: MDD=%d, CTRL=%d", sum(pheno$label=="MDD", na.rm=TRUE),
    sum(pheno$label=="CTRL", na.rm=TRUE))

# ---------- 3. 读表达矩阵 ----------
i1 <- grep("^!series_matrix_table_begin", lines); i2 <- grep("^!series_matrix_table_end", lines)
tab <- read.delim(text = lines[(i1 + 1):(i2 - 1)], header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE, quote = "")
msg("表达矩阵: %d 探针 × %d 样本", nrow(tab), ncol(tab) - 1)
expr <- as.matrix(tab[, -1]); rownames(expr) <- gsub('"', "", tab[[1]])
colnames(expr) <- gsub('^"|"$', "", colnames(expr))
if (max(expr, na.rm = TRUE) > 50) { expr <- log2(expr + 1); msg("已 log2 转换") }

# ---------- 4. 探针 → 基因 ----------
con2 <- gzfile(file.path(raw, "GPL570.annot.gz"), "rt")
al <- readLines(con2, warn = FALSE); close(con2)
j1 <- grep("^!platform_table_begin", al); j2 <- grep("^!platform_table_end", al)
ann <- read.delim(text = al[(j1 + 1):(j2 - 1)], header = TRUE,
                  check.names = FALSE, stringsAsFactors = FALSE, quote = "")
sym_col <- grep("Gene Symbol", colnames(ann), value = TRUE, ignore.case = TRUE)[1]
msg("注释表: %d 行 | 符号列: %s", nrow(ann), sym_col)
ann2 <- data.frame(ID = gsub('"', "", ann[["ID"]]), sym = ann[[sym_col]], stringsAsFactors = FALSE)
ann2$sym <- sub(" ?///.*$", "", ann2$sym)
ann2 <- ann2[!is.na(ann2$sym) & ann2$sym != "" & ann2$sym != "---", ]
ann2 <- ann2[!duplicated(ann2$ID), ]

common <- intersect(rownames(expr), ann2$ID)
msg("可注释探针: %d / %d", length(common), nrow(expr))
expr <- expr[common, , drop = FALSE]
sym  <- ann2$sym[match(rownames(expr), ann2$ID)]

# 多探针取均值最大者
ord <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
expr <- expr[ord, , drop = FALSE]; sym <- sym[ord]
keep <- !duplicated(sym)
expr <- expr[keep, , drop = FALSE]
rownames(expr) <- sym[keep]
msg("基因水平矩阵: %d 基因 × %d 样本", nrow(expr), ncol(expr))

# ---------- 5. limma ----------
sel_samples <- pheno$gsm[!is.na(pheno$label)]
expr <- expr[, intersect(colnames(expr), sel_samples), drop = FALSE]
lbl  <- pheno$label[match(colnames(expr), pheno$gsm)]
msg("用于分析: %d 样本（MDD=%d, CTRL=%d）", length(lbl), sum(lbl=="MDD"), sum(lbl=="CTRL"))

grp <- factor(lbl, levels = c("CTRL", "MDD"))
design <- model.matrix(~ grp)
fit <- lmFit(expr, design); fit <- eBayes(fit)
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "P")
tt$gene <- rownames(tt)

# ---------- 6. 输出 ----------
expr_out <- data.frame(gene = rownames(expr), expr, check.names = FALSE)
write.table(expr_out, file.path(out, "expr_blood.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tt[, c("gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
            file.path(out, "deg_blood.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pheno, file.path(out, "pheno_GSE98793.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- 7. 报告 ----------
sig05 <- sum(tt$adj.P.Val < 0.05)
sig05_lfc <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5)
sig05_lfc3 <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.3)
msg("\n=== DEG 结果 ===")
msg("adj.P<0.05: %d", sig05)
msg("adj.P<0.05 & |logFC|>=0.5: %d", sig05_lfc)
msg("adj.P<0.05 & |logFC|>=0.3: %d", sig05_lfc3)
msg("\nTop 10:")
print(head(tt[, c("gene","logFC","P.Value","adj.P.Val")], 10), row.names = FALSE)

# 关键基因
key <- c("NR3C1","FKBP5","SLC6A4","HTR2A","IL1B","TNF","BDNF","AKT1","GSK3B","CRP")
kk <- tt[tt$gene %in% key, c("gene","logFC","P.Value","adj.P.Val")]
if (nrow(kk)) { msg("\n关键基因:"); print(kk, row.names = FALSE) }


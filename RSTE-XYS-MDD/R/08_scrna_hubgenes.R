source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 08_scrna_hubgenes.R — 单细胞：hub 基因在哪类细胞表达
#
# 输入：data/raw/GSE336230_harmony.rds（GEO 提供的 Harmony 校正 Seurat 对象）
# 输出：
#   results/scrna_celltype_composition.tsv
#   results/scrna_hub_expression.tsv      （每细胞类型 × hub 基因 平均表达/阳性率）
#   results/scrna_target_score.tsv        （241 靶点的模块评分，按细胞类型）
#   figures/F_scRNA_dotplot.pdf
# =============================================================

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(Matrix)
})
root <- ROOT
out  <- file.path(root, "results"); fig <- file.path(root, "figures")
dir.create(out, showWarnings = FALSE); dir.create(fig, showWarnings = FALSE)

obj <- readRDS(file.path(root, "data/raw/GSE336230_harmony.rds"))
cat("对象类型:", paste(class(obj), collapse = "/"), "\n")
if (inherits(obj, "Seurat")) {
  cat("细胞数:", ncol(obj), "| 基因数:", nrow(obj), "\n")
  cat("Assays:", paste(Assays(obj), collapse = ", "), "\n")
  cat("meta 列:\n"); print(colnames(obj@meta.data))
} else {
  cat("不是 Seurat 对象，顶层元素:\n"); print(names(obj))
}

# ---- 找细胞类型列 ----
md <- obj@meta.data
cand <- grep("celltype|cell_type|cellType|annotation|ident|cluster|label", colnames(md),
             ignore.case = TRUE, value = TRUE)
cat("\n候选细胞类型列:", paste(cand, collapse = ", "), "\n")
ct_col <- if (length(cand)) cand[1] else "seurat_clusters"
if (!ct_col %in% colnames(md)) {
  if ("seurat_clusters" %in% colnames(md)) ct_col <- "seurat_clusters" else {
    obj <- tryCatch(FindClusters(obj, resolution = 0.5, verbose = FALSE), error = function(e) obj)
    ct_col <- "seurat_clusters"
  }
}
ct <- as.character(md[[ct_col]]); names(ct) <- rownames(md)
cat("使用列:", ct_col, "| 类型数:", length(unique(ct)), "\n")
print(sort(table(ct), decreasing = TRUE))

comp <- as.data.frame(table(`cell_type` = ct))
write.table(comp, file.path(out, "scrna_celltype_composition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- hub 基因表达 ----
hubs <- c("AKT1","TNF","IL1B","NFKB1","PTGS2","TP53","TGFB1","PIK3CG",
          "JUN","CASP3","MYC","IFNG","BCL2","MMP9","IL10")
obj <- tryCatch(JoinLayers(obj), error = function(e) obj)
DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else Assays(obj)[1]
expr <- tryCatch(LayerData(obj, assay = DefaultAssay(obj), layer = "data"),
                 error = function(e) GetAssayData(obj, slot = "data"))
hubs <- intersect(hubs, rownames(expr))
cat("\n可用 hub 基因:", length(hubs), "->", paste(hubs, collapse = ", "), "\n")

rows <- list()
for (g in hubs) {
  v <- as.numeric(expr[g, ])
  for (k in unique(ct)) {
    idx <- which(ct == k)
    rows[[length(rows) + 1]] <- data.frame(
      cell_type = k, gene = g, n_cell = length(idx),
      mean_expr = mean(v[idx]), pct_expr = 100 * mean(v[idx] > 0))
  }
}
hubtab <- do.call(rbind, rows)
write.table(hubtab, file.path(out, "scrna_hub_expression.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 241 靶点模块评分 ----
tg <- read.delim(file.path(root, "data/processed/XYS_targets_OB30DL018.tsv"))$gene
tg <- intersect(tg, rownames(expr))
cat("靶点可映射:", length(tg), "/241\n")
score <- Matrix::colMeans(expr[tg, , drop = FALSE])
scoretab <- data.frame(
  cell_type = names(tapply(score, ct, mean)),
  mean_target_score = as.numeric(tapply(score, ct, mean)),
  n_cell = as.numeric(tapply(score, ct, length)))
scoretab <- scoretab[order(-scoretab$mean_target_score), ]
write.table(scoretab, file.path(out, "scrna_target_score.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(scoretab, row.names = FALSE, digits = 4)

# ---- dotplot ----
hubtab$gene <- factor(hubtab$gene, levels = hubs)
p <- ggplot(hubtab, aes(gene, cell_type)) +
  geom_point(aes(size = pct_expr, colour = mean_expr)) +
  scale_colour_gradient(low = "grey90", high = "firebrick") +
  scale_size_area(max_size = 6) +
  theme_bw(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, colour = "Mean expr", size = "% expressing",
       title = "Hub gene expression across PBMC subsets (GSE336230)")
ggsave(file.path(fig, "F_scRNA_dotplot.pdf"), p, width = 7, height = 4)
cat("\n图已保存: figures/F_scRNA_dotplot.pdf\n")

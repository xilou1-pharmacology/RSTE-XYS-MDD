source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 09_scrna_final.R — 单细胞最终分析
#
# 数据：GSE336230（MDD 患者 PBMC，33,833 细胞）
#   condition = Low CRP / High CRP-Baseline / High CRP-Infliximab
#   cell_types = 自带注释（CD14+/CD16+ 单核、T/NK/B 等）
#
# 输出：
#   results/scrna_composition.tsv      各细胞类型在三组中的比例
#   results/scrna_hub_by_celltype.tsv  hub 基因在各细胞类型的表达
#   results/scrna_target_by_celltype.tsv  241 靶点评分
#   figures/F_scRNA_dotplot.pdf
#   figures/F_scRNA_composition.pdf
# =============================================================

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(Matrix)
})
root <- ROOT
out <- file.path(root,"results"); fig <- file.path(root,"figures")

obj <- readRDS(file.path(root,"data/raw/GSE336230_harmony.rds"))
md <- obj@meta.data
# 唯一标签：cluster + 注释
md$ct <- paste0("C", md$seurat_clusters, "_", md$cell_types)
obj@meta.data <- md
cat("细胞:", ncol(obj), "| 细胞类型:", length(unique(md$ct)), "\n")

# ---- 1. 组成（按 condition）----
tab <- as.data.frame(table(ct = md$ct, condition = md$condition))
tab$prop <- ave(tab$Freq, tab$condition, FUN = function(x) 100*x/sum(x))
write.table(tab, file.path(out,"scrna_composition.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ---- 2. hub 基因 × 细胞类型 ----
hubs <- c("TNF","IL1B","NFKB1","PTGS2","AKT1","TP53","TGFB1","PIK3CG","JUN","CASP3","MYC","IFNG","BCL2")
obj <- tryCatch(JoinLayers(obj), error=function(e) obj)
DefaultAssay(obj) <- "RNA"
expr <- tryCatch(LayerData(obj, assay="RNA", layer="data"),
                 error=function(e) GetAssayData(obj, slot="data"))
hubs <- intersect(hubs, rownames(expr))
cat("hub 基因:", length(hubs), "\n")

rows <- list()
for (g in hubs) {
  v <- as.numeric(expr[g,])
  for (k in unique(md$ct)) {
    i <- which(md$ct == k)
    rows[[length(rows)+1]] <- data.frame(cell_type=k, gene=g, n=length(i),
      mean_expr=mean(v[i]), pct=100*mean(v[i]>0))
  }
}
hubtab <- do.call(rbind, rows)
write.table(hubtab, file.path(out,"scrna_hub_by_celltype.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# ---- 3. 靶点评分 ----
tg <- read.delim(file.path(root,"data/processed/XYS_targets_OB30DL018.tsv"))$gene
tg <- intersect(tg, rownames(expr))
score <- Matrix::colMeans(expr[tg,,drop=FALSE])
st <- data.frame(cell_type=names(tapply(score, md$ct, mean)),
                 mean_target_score=as.numeric(tapply(score, md$ct, mean)),
                 n=as.numeric(tapply(score, md$ct, length)))
st <- st[order(-st$mean_target_score),]
write.table(st, file.path(out,"scrna_target_by_celltype.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
cat("\n靶点评分 Top10:\n"); print(head(st,10), row.names=FALSE, digits=4)

# ---- 4. 单核细胞比例（关键结果）----
mono <- tab[grepl("CD14|CD16", tab$ct),]
mono_sum <- aggregate(prop ~ condition, data=mono, FUN=sum)
cat("\n单核细胞占比（%）：\n"); print(mono_sum, row.names=FALSE, digits=3)

# ---- 5. 图 ----
hubtab$gene <- factor(hubtab$gene, levels=hubs)
p1 <- ggplot(hubtab, aes(gene, cell_type)) +
  geom_point(aes(size=pct, colour=mean_expr)) +
  scale_colour_gradient(low="grey90", high="firebrick") +
  scale_size_area(max_size=6) + theme_bw(base_size=8) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x=NULL,y=NULL,colour="Mean expr",size="% expressing",
       title="Hub genes across PBMC subsets (GSE336230)")
ggsave(file.path(fig,"F_scRNA_dotplot.pdf"), p1, width=7, height=5)

p2 <- ggplot(mono, aes(ct, prop, fill=condition)) +
  geom_col(position="dodge") + theme_bw(base_size=9) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x=NULL,y="% of cells",fill=NULL,title="Monocyte subsets by inflammatory status")
ggsave(file.path(fig,"F_scRNA_composition.pdf"), p2, width=6, height=4)
cat("\n图已保存\n")

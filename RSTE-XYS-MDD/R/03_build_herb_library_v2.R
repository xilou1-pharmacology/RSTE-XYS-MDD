source("config.R")
#!/usr/bin/env Rscript
# =============================================================
# 03_build_herb_library_v2.R
# 从平台内置三库合并构建 药材→靶点 表（TCMSP + HIT）
#
# 为什么合并：
#   TCMSP 经 OB>=30/DL>=0.18 过滤后，基因宇宙从 1442 缩到 446，
#   任何 300 靶点的方都会覆盖库的 2/3 → 汇聚检验饱和、失去分辨力。
#   HIT 补充 2195 个基因，并集 13060，GC 轴覆盖 19/20。
#
# 输出:
#   data/processed/herb_targets_v2.tsv    herb \t gene \t source
#   data/processed/herb_stats_v2.tsv
#   data/processed/axis_coverage.tsv      轴基因在各库的覆盖
# =============================================================

dat_dir <- Sys.getenv("RSTE_TCMS_DIR", file.path(ROOT, "data", "raw", "tcmsp"))
out_dir <- PROCESSED
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GC <- c("NR3C1","FKBP5","FKBP4","HSP90AA1","HSP90AB1","NR3C2","HSD11B1","HSD11B2","CRH",
        "CRHR1","CRHR2","POMC","MC2R","SERPINA6","TSC22D3","SGK1","PER1","TDO2","TPH1","TPH2")
HT <- c("SLC6A4","TPH1","TPH2","HTR1A","HTR1B","HTR2A","HTR2C","HTR4","SLC6A3","MAOA","MAOB","DDC")

message("=== 1. TCMSP ===")
tg  <- readRDS(file.path(dat_dir, "TCMSP_merged_targets.rds"))
ing <- readRDS(file.path(dat_dir, "TCMSP_merged_ingredients.rds"))
uni <- readRDS(file.path(dat_dir, "TCMSP-uniport.rds"))

map <- data.frame(tn = trimws(as.character(uni$target_name)),
                  g  = trimws(as.character(uni[["GENE name"]])), stringsAsFactors = FALSE)
map$g <- toupper(sub("[,; ].*$", "", map$g))
map <- map[!is.na(map$g) & map$g != "" & !is.na(map$tn), ]
map <- map[!duplicated(map$tn), ]

phys <- data.frame(MOL_ID = as.character(ing$MOL_ID),
                   ob = suppressWarnings(as.numeric(ing$ob)),
                   dl = suppressWarnings(as.numeric(ing$dl)), stringsAsFactors = FALSE)
phys <- phys[!duplicated(phys$MOL_ID), ]

tcm <- data.frame(herb = trimws(as.character(tg[["中药名"]])),
                  MOL_ID = as.character(tg$MOL_ID),
                  tn = trimws(as.character(tg$target_name)), stringsAsFactors = FALSE)
tcm <- merge(tcm, phys, by = "MOL_ID", all.x = TRUE)
tcm <- merge(tcm, map, by = "tn")
tcm_all <- unique(tcm[, c("herb", "g")]); names(tcm_all) <- c("herb", "gene")
tcm_all$source <- "TCMSP"
tcm_f   <- tcm[!is.na(tcm$ob) & tcm$ob >= 30 & !is.na(tcm$dl) & tcm$dl >= 0.18, ]
tcm_f   <- unique(tcm_f[, c("herb", "g")]); names(tcm_f) <- c("herb", "gene")
message(sprintf("  TCMSP 全量: %d 药 / %d 基因 %d 对", length(unique(tcm_all$herb)),
                length(unique(tcm_all$gene)), nrow(tcm_all)))
message(sprintf("  TCMSP OB/DL 过滤: %d 药 / %d 基因 %d 对", length(unique(tcm_f$herb)),
                length(unique(tcm_f$gene)), nrow(tcm_f)))

message("=== 2. HIT ===")
hit_t <- readRDS(file.path(dat_dir, "HIT_all_target_compound.rds"))
hit_h <- readRDS(file.path(dat_dir, "HIT_all_herb_compound.rds"))
ht <- data.frame(cid = as.character(hit_t[["Compound Id"]]),
                 gene = toupper(trimws(as.character(hit_t$Symbol))), stringsAsFactors = FALSE)
hh <- data.frame(cid = as.character(hit_h[["Compound Id"]]),
                 herb = trimws(as.character(hit_h[["Chinese Character"]])), stringsAsFactors = FALSE)
hh <- hh[!is.na(hh$herb) & hh$herb != "", ]
hh <- hh[!duplicated(hh$cid), ]
hit <- merge(ht, hh, by = "cid")
hit <- hit[!is.na(hit$gene) & hit$gene != "", c("herb", "gene")]
hit$source <- "HIT"
hit <- unique(hit)
message(sprintf("  HIT: %d 药 / %d 基因 %d 对", length(unique(hit$herb)),
                length(unique(hit$gene)), nrow(hit)))

message("=== 3. 合并 ===")
all <- rbind(tcm_all, hit)
ht2 <- unique(all[, c("herb", "gene")])
n_per <- table(ht2$herb)
keep <- names(n_per)[n_per >= 5 & n_per <= 1200]
ht2 <- ht2[ht2$herb %in% keep, ]
message(sprintf("  并集: %d 药 / %d 基因 %d 对（剔除 %d 味：靶点<5 或 >1200）",
                length(keep), length(unique(ht2$gene)), nrow(ht2), length(setdiff(names(n_per), keep))))

st <- data.frame(herb = names(n_per), n_target = as.integer(n_per))
st <- st[st$herb %in% keep, ]; st <- st[order(-st$n_target), ]

write.table(ht2, file.path(out_dir, "herb_targets_v2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(st, file.path(out_dir, "herb_stats_v2.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 轴覆盖 ----
U <- unique(ht2$gene)
cov <- data.frame(axis = c(rep("GC", length(GC)), rep("HT", length(HT))),
                  gene = c(GC, HT), in_library = c(GC %in% U, HT %in% U))
write.table(cov, file.path(out_dir, "axis_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nGC 轴覆盖 %d/%d\n", sum(GC %in% U), length(GC)))
cat("  缺失:", paste(setdiff(GC, U), collapse = ", "), "\n")
cat(sprintf("HT 轴覆盖 %d/%d\n", sum(HT %in% U), length(HT)))
cat("  缺失:", paste(setdiff(HT, U), collapse = ", "), "\n")

# ---- 伪方规模 ----
set.seed(1)
sim <- replicate(2000, length(unique(ht2$gene[ht2$herb %in% sample(keep, 8)])))
cat(sprintf("\n随机 8 味伪方靶点数：中位 %.0f（四分位 %.0f–%.0f）\n",
            median(sim), quantile(sim, .25), quantile(sim, .75)))
cat(sprintf("基因宇宙 %d；300 靶点覆盖 %.1f%%\n", length(U), 100 * 300 / length(U)))

cat("\n=== 完成 ===\n")


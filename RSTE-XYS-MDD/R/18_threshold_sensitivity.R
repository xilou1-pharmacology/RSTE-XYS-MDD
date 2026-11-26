# 18_threshold_sensitivity.R — OB/DL 阈值稳定性分析
# 对应 Allemailem 2026 五维评估中的 "stability"
# 对每个 (OB, DL) 组合重建 XYS 靶点集，跑同样的 GO BP 富集，比较结果

suppressPackageStartupMessages({ library(data.table) })

dat_dir <- "D:/LKFlow2.0_Windows/LKFlow2.0_Windows/LKFlow2/0.data"
root    <- "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
out_dir <- file.path(root, "results")
gmt_f   <- file.path(root, "data/raw/kegg/GO_Biological_Process_2021.gmt")

## ---- TCMSP ----
tg  <- readRDS(file.path(dat_dir, "TCMSP_merged_targets.rds"))
ing <- readRDS(file.path(dat_dir, "TCMSP_merged_ingredients.rds"))
uni <- readRDS(file.path(dat_dir, "TCMSP-uniport.rds"))
hit_t <- readRDS(file.path(dat_dir, "HIT_all_target_compound.rds"))

map <- data.frame(tn = trimws(as.character(uni$target_name)),
                  g  = trimws(as.character(uni[["GENE name"]])), stringsAsFactors = FALSE)
map$g <- toupper(sub("[,; ].*$", "", map$g))
map <- map[!is.na(map$g) & map$g != "" & !is.na(map$tn), ]
map <- map[!duplicated(map$tn), ]

phys <- data.frame(MOL_ID = as.character(ing$MOL_ID),
                   ob = suppressWarnings(as.numeric(ing$ob)),
                   dl = suppressWarnings(as.numeric(ing$dl)), stringsAsFactors = FALSE)
phys <- phys[!duplicated(phys$MOL_ID), ]

x <- data.frame(herb = trimws(as.character(tg[["中药名"]])),
                MOL_ID = as.character(tg$MOL_ID),
                tn = trimws(as.character(tg$target_name)), stringsAsFactors = FALSE)
x <- merge(x, phys, by = "MOL_ID", all.x = TRUE)
x <- merge(x, map, by = "tn")

XYS <- c("柴胡","当归","白芍","白术","茯苓","甘草","薄荷","生姜")
x8 <- x[x$herb %in% XYS, ]

cat("逍遥散 TCMSP 记录:", nrow(x8), "| 全量基因:", length(unique(x8$g)), "\n")

## ---- HIT（无 OB/DL，恒定加入）----
hit_h <- readRDS(file.path(dat_dir, "HIT_all_herb_compound.rds"))
ht <- data.frame(cid = as.character(hit_t[["Compound Id"]]),
                 gene = toupper(trimws(as.character(hit_t$Symbol))), stringsAsFactors = FALSE)
hh <- data.frame(cid = as.character(hit_h[["Compound Id"]]),
                 herb = trimws(as.character(hit_h[["Chinese Character"]])), stringsAsFactors = FALSE)
hh <- hh[!is.na(hh$herb) & hh$herb != "", ]; hh <- hh[!duplicated(hh$cid), ]
hitx <- merge(ht, hh, by = "cid")
hit_xys <- unique(hitx$gene[hitx$herb %in% XYS])
hit_xys <- hit_xys[!is.na(hit_xys) & hit_xys != ""]
cat("HIT 逍遥散基因:", length(hit_xys), "\n")

## ---- GMT ----
read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE)
  ln <- ln[nzchar(ln)]
  sets <- lapply(strsplit(ln, "\t"), function(p) unique(p[-c(1, 2)]))
  names(sets) <- vapply(strsplit(ln, "\t"), `[`, "", 1)
  sets
}
gmt <- read_gmt(gmt_f)
cat("GO BP 条目:", length(gmt), "\n")
universe <- unique(unlist(gmt))
N <- length(universe)
cat("注释宇宙 N =", N, "\n")

enrich <- function(genes, label) {
  q <- intersect(genes, universe); nq <- length(q)
  if (nq < 5) return(NULL)
  rows <- list()
  for (t in names(gmt)) {
    S <- intersect(gmt[[t]], universe); m <- length(S)
    if (m < 5 || m > N - 5) continue
    k <- length(intersect(q, S))
    if (k == 0) next
    p <- phyper(k - 1, m, N - m, nq, lower.tail = FALSE)
    rows[[t]] <- data.frame(term = t, hits = k, size = m, query = nq, P = p,
                            stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  if (is.null(d)) return(NULL)
  d$FDR <- p.adjust(d$P, "BH")
  d$threshold <- label
  d[order(d$P), ]
}

GRID <- list(
  c(30, 0.18), c(30, 0.15), c(20, 0.18), c(20, 0.15),
  c(15, 0.18), c(15, 0.15), c(0, 0)
)

all_res <- list(); summ <- list()
for (g in GRID) {
  ob <- g[1]; dl <- g[2]
  sel <- x8[!is.na(x8$ob) & x8$ob >= ob & !is.na(x8$dl) & x8$dl >= dl, ]
  tcm_genes <- unique(sel$g)
  genes <- unique(c(tcm_genes, hit_xys))
  lab <- sprintf("OB%g_DL%g", ob, dl)
  r <- enrich(genes, lab)
  if (is.null(r)) { cat(lab, ": 基因不足\n"); next }
  sig <- r$FDR < 0.05
  cat(sprintf("%-12s 基因 %4d | 显著条目 %5d\n", lab, length(genes), sum(sig)))
  summ[[lab]] <- data.frame(threshold = lab, OB = ob, DL = dl,
                            n_gene = length(genes), n_term = nrow(r), n_sig = sum(sig))
  all_res[[lab]] <- r
}

S <- do.call(rbind, summ)

## 与基线（OB30_DL0.18）比较显著条目重叠
base <- all_res[["OB30_DL0.18"]]
base_sig <- base$term[base$FDR < 0.05]
S$jaccard_vs_base <- NA
S$kept_pct <- NA
for (lab in S$threshold) {
  if (lab == "OB30_DL0.18") { S$jaccard_vs_base[S$threshold==lab] <- 1; S$kept_pct[S$threshold==lab] <- 100; next }
  s <- all_res[[lab]]$term[all_res[[lab]]$FDR < 0.05]
  inter <- length(intersect(s, base_sig)); uni <- length(union(s, base_sig))
  S$jaccard_vs_base[S$threshold==lab] <- round(inter/uni, 3)
  S$kept_pct[S$threshold==lab] <- round(100*length(intersect(s, base_sig))/length(base_sig), 1)
}

write.table(S, file.path(out_dir, "threshold_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== 汇总 ===\n"); print(S, row.names = FALSE)

## 关键通路在各阈值下的显著性
KEY <- c("citrate cycle", "steroid hormone biosynthetic", "asthma", "cytokine")
cat("\n=== 关键通路命中 ===\n")
for (lab in names(all_res)) {
  d <- all_res[[lab]]
  hits <- c()
  for (k in KEY) {
    i <- grep(k, d$term, ignore.case = TRUE)
    if (length(i)) hits <- c(hits, sprintf("%s(P=%.1e)", k, min(d$P[i]))) else hits <- c(hits, paste0(k, "(无)"))
  }
  cat(sprintf("%-12s %s\n", lab, paste(hits, collapse = " | ")))
}

cat("\n完成 -> results/threshold_sensitivity.tsv\n")

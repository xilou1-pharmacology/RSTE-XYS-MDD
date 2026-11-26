# 19_threshold_KEGG.R — OB/DL 阈值稳定性（KEGG，与主分析一致）
# 基线构造：TCMSP OB>=30 & DL>=0.18，逍遥散 8 味 → 241 基因（不含 HIT）
# 背景：KEGG_2021_Human 基因宇宙 N = 8,078

suppressPackageStartupMessages({ library(data.table) })

dat_dir <- "D:/LKFlow2.0_Windows/LKFlow2.0_Windows/LKFlow2/0.data"
root    <- "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
out_dir <- file.path(root, "results")
gmt_f   <- file.path(root, "data/raw/kegg/KEGG_2021_Human.gmt")

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

x <- data.frame(herb = trimws(as.character(tg[["中药名"]])),
                MOL_ID = as.character(tg$MOL_ID),
                tn = trimws(as.character(tg$target_name)), stringsAsFactors = FALSE)
x <- merge(x, phys, by = "MOL_ID")
x <- merge(x, map, by = "tn")

XYS <- c("柴胡","当归","白芍","白术","茯苓","甘草","薄荷","生姜")
x8 <- x[x$herb %in% XYS, ]

read_gmt <- function(f) {
  ln <- readLines(f, warn = FALSE); ln <- ln[nzchar(ln)]
  sets <- lapply(strsplit(ln, "\t"), function(p) unique(p[-c(1, 2)]))
  names(sets) <- vapply(strsplit(ln, "\t"), `[`, "", 1)
  sets
}
gmt <- read_gmt(gmt_f)
universe <- unique(unlist(gmt)); N <- length(universe)
cat("KEGG 条目:", length(gmt), "| 基因宇宙 N =", N, "\n")

enrich <- function(genes, label) {
  q <- intersect(genes, universe); nq <- length(q)
  rows <- list()
  for (t in names(gmt)) {
    S <- intersect(gmt[[t]], universe); m <- length(S)
    if (m < 5) next
    k <- length(intersect(q, S))
    if (k == 0) next
    p <- phyper(k - 1, m, N - m, nq, lower.tail = FALSE)
    rows[[t]] <- data.frame(term = t, hits = k, size = m, query = nq, P = p,
                            stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, rows)
  d$FDR <- p.adjust(d$P, "BH"); d$threshold <- label
  d[order(d$P), ]
}

GRID <- list(c(30, 0.18), c(30, 0.15), c(25, 0.18), c(20, 0.18),
             c(20, 0.15), c(15, 0.18), c(10, 0.18), c(0, 0))

all_res <- list(); summ <- list()
for (g in GRID) {
  sel <- x8[!is.na(x8$ob) & x8$ob >= g[1] & !is.na(x8$dl) & x8$dl >= g[2], ]
  genes <- unique(sel$g)
  lab <- sprintf("OB%g_DL%g", g[1], g[2])
  r <- enrich(genes, lab)
  sig <- r$FDR < 0.05
  cat(sprintf("%-12s 基因 %4d | 检验 %4d | 显著 %4d\n", lab, length(genes), nrow(r), sum(sig)))
  summ[[lab]] <- data.frame(threshold = lab, OB = g[1], DL = g[2], n_gene = length(genes),
                            n_tested = nrow(r), n_sig = sum(sig))
  all_res[[lab]] <- r
}

S <- do.call(rbind, summ)
base <- all_res[["OB30_DL0.18"]]; base_sig <- base$term[base$FDR < 0.05]
S$jaccard <- NA; S$retained_pct <- NA
for (lab in S$threshold) {
  s <- all_res[[lab]]$term[all_res[[lab]]$FDR < 0.05]
  S$jaccard[S$threshold == lab] <- round(length(intersect(s, base_sig)) /
                                         length(union(s, base_sig)), 3)
  S$retained_pct[S$threshold == lab] <- round(100 * length(intersect(s, base_sig)) /
                                              length(base_sig), 1)
}
write.table(S, file.path(out_dir, "threshold_sensitivity_KEGG.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== 汇总 ===\n"); print(S, row.names = FALSE)

cat("\n=== 基线 top 10 通路在各阈值下的 P ===\n")
top <- head(base$term, 10)
for (t in top) {
  vals <- sapply(names(all_res), function(l) {
    d <- all_res[[l]]; i <- which(d$term == t)
    if (!length(i)) NA else d$P[i]
  })
  cat(sprintf("%-52s %s\n", substr(t, 1, 52),
              paste(sprintf("%.1e", vals), collapse = " ")))
}
cat("\n完成 -> results/threshold_sensitivity_KEGG.tsv\n")

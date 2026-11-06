# config.R — 项目路径（从仓库根目录解析）
ROOT <- Sys.getenv("RSTE_ROOT", unset = normalizePath("."))
DATA <- file.path(ROOT, "data")
RAW <- file.path(DATA, "raw")
PROCESSED <- file.path(DATA, "processed")
RESULTS <- file.path(ROOT, "results")
FIGURES <- file.path(ROOT, "figures")
for (d in c(RAW, PROCESSED, RESULTS, FIGURES))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

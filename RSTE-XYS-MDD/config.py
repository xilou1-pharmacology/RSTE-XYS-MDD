"""config.py — 项目路径（从仓库根目录解析，不依赖任何绝对路径）"""
import os
from pathlib import Path

ROOT = Path(os.environ.get("RSTE_ROOT", Path(__file__).resolve().parent)).resolve()
DATA = ROOT / "data"
RAW = DATA / "raw"
PROCESSED = DATA / "processed"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
for _d in (RAW, PROCESSED, RESULTS, FIGURES):
    _d.mkdir(parents=True, exist_ok=True)

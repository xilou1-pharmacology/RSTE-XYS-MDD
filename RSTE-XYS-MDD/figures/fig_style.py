"""fig_style.py — 全图统一风格（Q1 期刊级）

用法:
    import fig_style as S
    S.setup()
    fig, ax = plt.subplots(...)
    S.panel(ax, "A")
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import matplotlib
matplotlib.use("Agg")
import os
os.environ.setdefault("MPLCONFIGDIR", __import__("tempfile").gettempdir())
import matplotlib.pyplot as plt
from matplotlib import rcParams

# ---- 配色（低饱和、色盲友好）----
PAL = {
    "monarch":   "#B2182B",
    "minister":  "#EF8A62",
    "assistant": "#67A9CF",
    "courier":   "#2166AC",
    "grey":      "#9E9E9E",
    "dark":      "#263238",
    "accent":    "#D6604D",
    "up":        "#C62828",
    "down":      "#1565C0",
    "neutral":   "#BDBDBD",
}
ROLE_COLOR = {"Monarch": PAL["monarch"], "Minister": PAL["minister"],
              "Assistant": PAL["assistant"], "Courier": PAL["courier"]}


def setup():
    rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": 8.5,
        "axes.titlesize": 9.5,
        "axes.labelsize": 8.5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7.5,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def panel(ax, letter, x=-0.16, y=1.07):
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="left")


def save(fig, path):
    fig.savefig(path + ".pdf")
    fig.savefig(path + ".png")
    print("saved", path + ".pdf/.png")

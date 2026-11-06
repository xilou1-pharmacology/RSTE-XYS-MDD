"""fig7_v2.py — Figure 7: compartment-resolved metabolomics (four panels)

A  plasma volcano (ST000062)
B  plasma versus CSF effect sizes for shared metabolites
C  metabolites most depleted in CSF
D  distribution of within-individual plasma-CSF correlations (n = 50)
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import sys, json

import fig_style as S
import numpy as np
import matplotlib.pyplot as plt

R = str(ROOT)
S.setup()

pl = [l.rstrip("\n").split("\t") for l in
      open(f"{R}/data/processed/metabolomics_ST000062.tsv", encoding="utf-8")][1:]
pl = [(r[0], float(r[1]), float(r[2]), float(r[3])) for r in pl if len(r) >= 4]
pc = [l.rstrip("\n").split("\t") for l in
      open(f"{R}/results/plasma_CSF_concordance.tsv", encoding="utf-8")][1:]
pc = [(r[0], float(r[1]), float(r[2]), float(r[3]), float(r[4])) for r in pc]
pairs = json.load(open(f"{R}/results/paired_plasma_csf.json", encoding="utf-8"))

KEY = {"Niacinamide": "Niacinamide", "Ketoleucine": "Ketoleucine", "Pyruvic acid": "Pyruvate",
       "Squalene": "Squalene", "Glutamic acid": "Glutamate", "Aspartic acid": "Aspartate",
       "2-Oxoglutaric acid": "2-OG"}

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.38, height_ratios=[1, 1])

# A
ax = fig.add_subplot(gs[0, 0])
x = np.array([r[1] for r in pl]); y = np.array([-np.log10(r[3]) for r in pl])
sig = y > -np.log10(0.05)
ax.scatter(x[~sig], y[~sig], s=10, c=S.PAL["neutral"], alpha=0.7, edgecolors="none")
ax.scatter(x[sig], y[sig], s=16, c=[S.PAL["up"] if v > 0 else S.PAL["down"] for v in x[sig]],
           alpha=0.9, edgecolors="white", linewidths=0.3)
ax.axhline(-np.log10(0.05), ls="--", lw=0.8, color=S.PAL["dark"])
ax.axvline(0, lw=0.5, color="k")
for i, r in enumerate(pl):
    if r[0] in KEY and y[i] > 0.9:
        ax.annotate(KEY[r[0]], (x[i], y[i]), fontsize=6.6, color=S.PAL["dark"],
                    textcoords="offset points", xytext=(3, 2))
ax.set_xlabel("log2 fold change (plasma)"); ax.set_ylabel("−log10 FDR")
ax.set_title("Plasma metabolome", loc="left")
S.panel(ax, "A")

# B
ax = fig.add_subplot(gs[0, 1])
px = np.array([r[1] for r in pc]); py = np.array([r[3] for r in pc])
ax.axhline(0, lw=0.6, color="grey"); ax.axvline(0, lw=0.6, color="grey")
ax.scatter(px, py, s=11, c=S.PAL["grey"], alpha=0.75, edgecolors="none")
for r in pc:
    if r[0] in KEY:
        ax.scatter(r[1], r[3], s=44, facecolors="none", edgecolors=S.PAL["accent"],
                   linewidths=1.0, zorder=4)
        off={"Aspartate":(-26,-9),"Glutamate":(6,-9),"Ketoleucine":(-10,9),
             "Niacinamide":(-34,-10),"Squalene":(-30,4),"Pyruvate":(-6,8),"2-OG":(-16,8)}.get(KEY[r[0]],(4,3))
        ax.annotate(KEY[r[0]], (r[1], r[3]), fontsize=6.0, color=S.PAL["accent"],
                    textcoords="offset points", xytext=off)
ax.set_xlabel("log2FC plasma"); ax.set_ylabel("log2FC CSF")
ax.text(0.03, 0.06, "ρ = −0.25, P = 0.009\nsign agreement 27.9%",
        transform=ax.transAxes, fontsize=7, color=S.PAL["dark"])
ax.set_title("Disease effects diverge", loc="left")
S.panel(ax, "B")

# C
ax = fig.add_subplot(gs[1, 0])
pc_sorted = sorted(pc, key=lambda r: r[3])[:10][::-1]
names = [r[0][:26] for r in pc_sorted]
vals = [r[3] for r in pc_sorted]
ax.barh(range(len(pc_sorted)), vals, color=S.PAL["down"], height=0.70)
ax.set_yticks(range(len(pc_sorted))); ax.set_yticklabels(names, fontsize=6.8)
ax.axvline(0, color="k", lw=0.6)
ax.set_xlabel("log2 fold change (CSF)")
ax.set_xlim(min(vals) * 1.28, 0.05)
ax.set_title("Most depleted in CSF", loc="left")
S.panel(ax, "C")

# D
ax = fig.add_subplot(gs[1, 1])
rv = np.array([r[1] for r in pairs])
ax.hist(rv, bins=np.arange(-0.5, 0.55, 0.05), color=S.ROLE_COLOR["Assistant"],
        edgecolor="white", linewidth=0.4)
ax.axvline(0, color="k", lw=0.7)
ax.axvline(float(np.median(rv)), color=S.PAL["accent"], lw=1.2, ls="--")
ax.text(0.03, 0.94, f"n = {len(rv)} metabolites\nmedian ρ = {np.median(rv):.2f}\n"
                    f"{100*np.mean(rv>0):.0f}% positive",
        transform=ax.transAxes, fontsize=7, va="top", color=S.PAL["dark"])
ax.set_xlabel("within-individual plasma–CSF Spearman ρ")
ax.set_ylabel("metabolites")
ax.set_title("Compartments are independent", loc="left")
S.panel(ax, "D")

S.save(fig, f"{R}/figures/F7_metabolomics")

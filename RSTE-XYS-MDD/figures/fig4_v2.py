"""fig4_v2.py — Figure 4: transcriptomic landscape and immune signatures

(a) DEG volcano (GSE98793, 128 cases vs 64 controls)
(b) ssGSEA immune signature differences
(c) distributions of the three informative signatures
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import fig_style as S
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = Path(str(ROOT))
S.setup()

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.44, wspace=0.46, height_ratios=[1, 0.85])

# ---------------- (a) volcano ----------------
ax = fig.add_subplot(gs[0, :])
deg = pd.read_csv(R / "data/processed/deg_128vs64.tsv", sep="\t")
fc = deg["logFC"].values
p = deg["P.Value"].values
y = -np.log10(np.clip(p, 1e-300, 1))
sig = (p < 0.05)
up, dn = sig & (fc > 0), sig & (fc < 0)
ax.scatter(fc[~sig], y[~sig], s=4, c="#CFD8DC", alpha=0.55, edgecolors="none", rasterized=True)
ax.scatter(fc[up], y[up], s=7, c=S.PAL["up"], alpha=0.8, edgecolors="none", rasterized=True)
ax.scatter(fc[dn], y[dn], s=7, c=S.PAL["down"], alpha=0.8, edgecolors="none", rasterized=True)
ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color="#78909C")
ax.axvline(0.5, ls=":", lw=0.6, color="#78909C")
ax.axvline(-0.5, ls=":", lw=0.6, color="#78909C")
lab = deg.loc[sig].copy(); lab["y"] = y[sig]
up_lab = lab[lab["logFC"] > 0].nlargest(2, "y")
dn_lab = lab[lab["logFC"] < 0].nlargest(2, "y")
for k, (_, r) in enumerate(pd.concat([up_lab, dn_lab]).iterrows()):
    dx = 4 if r["logFC"] > 0 else -4
    dy = 2 if k % 2 == 0 else -10
    ax.annotate(r["gene"], (r["logFC"], r["y"]), fontsize=6.2, color=S.PAL["dark"],
                textcoords="offset points", xytext=(dx, dy),
                ha="left" if r["logFC"] > 0 else "right")
ax.set_xlabel("log2 fold change (MDD / control)")
ax.set_ylabel("−log10 P")
ax.text(0.015, 0.96,
        f"{int(sig.sum())} genes at P < 0.05\n"
        f"{int((deg['adj.P.Val'] < 0.05).sum())} genes at FDR < 0.05",
        transform=ax.transAxes, fontsize=6.8, va="top", color=S.PAL["dark"])
ax.set_title("Differential expression", loc="left")
S.panel(ax, "a")

# ---------------- (b) immune bars ----------------
ax = fig.add_subplot(gs[1, 0])
imm = pd.read_csv(R / "results/immune_diff.tsv", sep="\t")
imm["short"] = (imm["signature"].str.replace("_et_al", "", regex=False)
                .str.replace("_", " ", regex=False))
imm = imm.sort_values("diff")
colors = [S.PAL["down"] if f < 0.05 else "#B0BEC5" for f in imm["FDR"]]
ax.barh(range(len(imm)), imm["diff"], color=colors, height=0.70)
ax.axvline(0, lw=0.6, color="k")
ax.set_yticks(range(len(imm)))
ax.set_yticklabels(imm["short"], fontsize=6.4)
for i, (d, f) in enumerate(zip(imm["diff"], imm["FDR"])):
    ax.text(d - 0.006 if d < 0 else d + 0.006, i, f"{f:.3f}", va="center",
            ha="right" if d < 0 else "left", fontsize=5.6,
            color=S.PAL["down"] if f < 0.05 else "#78909C")
ax.set_xlabel("ssGSEA difference (MDD − control)")
ax.set_xlim(imm["diff"].min() * 1.5, imm["diff"].max() * 1.5)
ax.set_title("Immune signatures", loc="left")
S.panel(ax, "b")

# ---------------- (c) boxplots ----------------
ax = fig.add_subplot(gs[1, 1])
ss = pd.read_csv(R / "data/processed/immune_ssgsea.tsv", sep="\t").set_index("signature")
ph = pd.read_csv(R / "data/processed/pheno_GSE98793.tsv", sep="\t")
ctrl = ph.loc[ph["group"].str.startswith("CNTL"), "gsm"].tolist()
case = ph.loc[ph["group"].str.startswith("CASE"), "gsm"].tolist()
pick = ["T_cells_Danaher_et_al", "CD8_Rooney_et_al", "Macrophages_Rooney_et_al"]
names = ["T cells", "CD8 T cells", "Macrophages"]
pos = np.arange(len(pick))
for k, (sig_name, nm) in enumerate(zip(pick, names)):
    a = pd.to_numeric(ss.loc[sig_name, [c for c in ctrl if c in ss.columns]], errors="coerce").dropna().values
    b = pd.to_numeric(ss.loc[sig_name, [c for c in case if c in ss.columns]], errors="coerce").dropna().values
    for off, dat, col in [(-0.17, a, "#B0BEC5"), (0.17, b, S.PAL["accent"])]:
        bp = ax.boxplot([dat], positions=[pos[k] + off], widths=0.27, patch_artist=True,
                        showfliers=False, medianprops=dict(color="k", lw=0.9))
        bp["boxes"][0].set(facecolor=col, edgecolor="#546E7A", linewidth=0.6)
    top = max(np.nanmax(a), np.nanmax(b))
    ax.text(pos[k], top + 0.07, f"{np.nanmean(b)-np.nanmean(a):+.2f}", ha="center",
            fontsize=6.8, color=S.PAL["dark"])
ax.set_xticks(pos); ax.set_xticklabels(names, fontsize=7)
ax.set_ylabel("ssGSEA score")
ax.legend(handles=[Line2D([], [], marker="s", ls="", ms=6, mfc="#B0BEC5", mec="#546E7A", label="control"),
                   Line2D([], [], marker="s", ls="", ms=6, mfc=S.PAL["accent"], mec="#546E7A", label="MDD")],
          frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.155), ncol=2, fontsize=6.4)
ax.set_title("Selected signatures", loc="left")
S.panel(ax, "c")

S.save(fig, str(R / "figures/F4_immune"))

import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
"""fig456.py — Figures 4 (immune), 5 (single cell), 6 (Mendelian randomization)"""
import sys, os, collections, csv

import fig_style as S
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = str(ROOT)
S.setup()

# =====================================================================
# Figure 4 — immune deconvolution
# =====================================================================
imm = [l.rstrip("\n").split("\t") for l in open(f"{R}/results/immune_diff.tsv", encoding="utf-8")][1:]
sig = [r[0].replace("_et_al", "").replace("_", " ") for r in imm]
diff = np.array([float(r[3]) for r in imm])
pfdr = np.array([float(r[5]) for r in imm])
order = np.argsort(diff)
sig = [sig[i] for i in order]; diff = diff[order]; pfdr = pfdr[order]

fig = plt.figure(figsize=(7.2, 3.4))
gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.36)

ax = fig.add_subplot(gs[0, 0])
cols = [S.PAL["down"] if (d < 0 and q < 0.05) else S.PAL["up"] if (d > 0 and q < 0.05)
        else ("#90A4AE" if q < 0.10 else S.PAL["neutral"]) for d, q in zip(diff, pfdr)]
ax.barh(range(len(sig)), diff, color=cols, height=0.7)
ax.axvline(0, color="k", lw=0.6)
ax.set_yticks(range(len(sig))); ax.set_yticklabels(sig, fontsize=6.2)
for i, (d, q) in enumerate(zip(diff, pfdr)):
    ax.text(d + (0.012 if d > 0 else -0.012), i, f"FDR {q:.3f}", va="center",
            ha="left" if d > 0 else "right", fontsize=5.2, color=S.PAL["dark"])
ax.set_xlabel("ssGSEA score difference (MDD − control)")
ax.set_xlim(min(diff) * 1.55, max(diff) * 1.75)
ax.set_title("Immune signatures in MDD blood (GSE98793)", loc="left")
S.panel(ax, "A")

# B: distributions
raw = [l.rstrip("\n").split("\t") for l in open(f"{R}/data/processed/immune_ssgsea.tsv", encoding="utf-8")]
hdr = raw[0][1:]
vals = {r[0]: np.array([float(x) for x in r[1:]]) for r in raw[1:]}
phen = [l.rstrip("\n").split("\t") for l in open(f"{R}/data/processed/pheno_GSE98793.tsv", encoding="utf-8")]
ph = phen[0]
ci = 2 if len(ph) > 2 else 1
lab = {r[0]: r[ci] for r in phen[1:]}
grp = np.array([1 if any(k in str(lab.get(s, "")).lower() for k in ["mdd", "case", "patient", "depress"]) else 0
                for s in hdr])

ax = fig.add_subplot(gs[0, 1])
picks = ["T_cells_Danaher_et_al", "CD8_Rooney_et_al", "Macrophages_Rooney_et_al"]
pos = np.arange(len(picks))
for i, k in enumerate(picks):
    v = vals[k]
    for g, off, c in [(0, -0.17, S.PAL["neutral"]), (1, 0.17, S.PAL["accent"])]:
        x = v[grp == g]
        if len(x) == 0: continue
        bp = ax.boxplot([x], positions=[pos[i] + off], widths=0.28, patch_artist=True,
                        showfliers=False, medianprops=dict(color="k", lw=0.7))
        for b in bp["boxes"]:
            b.set(facecolor=c, edgecolor="#546E7A", linewidth=0.5, alpha=0.85)
    d = v[grp == 1].mean() - v[grp == 0].mean()
    ax.text(pos[i], 0.98, f"{d:+.2f}", ha="center", fontsize=7.4, color=S.PAL["dark"],
            transform=ax.get_xaxis_transform())
ax.set_xticks(pos); ax.set_xticklabels(["T cells", "CD8 T", "Macrophage"], fontsize=6.5)
ax.set_ylabel("ssGSEA score")
ax.legend(handles=[Line2D([], [], marker="s", ls="", mfc=S.PAL["neutral"], mec="#546E7A", label="control"),
                   Line2D([], [], marker="s", ls="", mfc=S.PAL["accent"], mec="#546E7A", label="MDD")],
          frameon=False, loc="lower center", bbox_to_anchor=(0.5, -0.34), ncol=2)
ax.set_ylim(-1.10, 1.15)
ax.set_title("Selected signatures", loc="left")
S.panel(ax, "B")
plt.tight_layout()
S.save(fig, f"{R}/figures/F4_immune")

# =====================================================================
# Figure 5 — single cell
# =====================================================================
pat = [l.rstrip("\n").split("\t") for l in open(f"{R}/results/scrna_patient_level.tsv", encoding="utf-8")]
ph_ = pat[0]; pat = pat[1:]
tgtct = [l.rstrip("\n").split("\t") for l in open(f"{R}/results/scrna_target_by_celltype_merged.tsv", encoding="utf-8")][1:]
ctn = {r[0]: int(float(r[2])) for r in tgtct}
hub = [l.rstrip("\n").split("\t") for l in open(f"{R}/results/scrna_hub_by_celltype_merged.tsv", encoding="utf-8")][1:]

fig = plt.figure(figsize=(7.2, 6.2))
gs = fig.add_gridspec(2, 3, hspace=0.55, wspace=0.45, height_ratios=[1, 1.15])

# A composition
ax = fig.add_subplot(gs[0, 0])
items = sorted(ctn.items(), key=lambda x: -x[1])[:9]
ax.barh(range(len(items))[::-1], [v for _, v in items], color=S.PAL["assistant"], height=0.72)
ax.set_yticks(range(len(items))[::-1]); ax.set_yticklabels([k[:18] for k, _ in items], fontsize=5.6)
ax.set_xlabel("cells"); ax.set_title("PBMC composition", loc="left")
S.panel(ax, "A")

# B target score by cell type
ax = fig.add_subplot(gs[0, 1])
tt = sorted([(r[0], float(r[1])) for r in tgtct], key=lambda x: x[1])
ax.barh(range(len(tt)), [v for _, v in tt],
        color=[S.PAL["accent"] if "Monocyte" in k else S.PAL["neutral"] for k, _ in tt], height=0.72)
ax.set_yticks(range(len(tt))); ax.set_yticklabels([k[:18] for k, _ in tt], fontsize=5.2)
ax.set_xlabel("mean target-set score"); ax.set_title("XYS targets by subset", loc="left")
S.panel(ax, "B")

# C paired delta
ax = fig.add_subplot(gs[0, 2])
for r in pat:
    d = float(r[6])
    ax.plot([0, 1], [float(r[4]), float(r[5])], color="#B0BEC5", lw=0.7, zorder=1)
    ax.scatter([0], [float(r[4])], s=16, c=S.PAL["accent"], zorder=3, edgecolors="white", linewidths=0.3)
    ax.scatter([1], [float(r[5])], s=16, c=S.PAL["assistant"], zorder=3, edgecolors="white", linewidths=0.3)
ax.set_xticks([0, 1]); ax.set_xticklabels(["monocyte", "other"], fontsize=6.5)
ax.set_ylabel("target-set score"); ax.set_xlim(-0.3, 1.3)
ax.text(0.5, 0.96, "8/9 donors positive\npaired P = 0.008", transform=ax.transAxes,
        ha="center", va="top", fontsize=5.6)
ax.set_title("Within-donor contrast", loc="left")
S.panel(ax, "C")

# D hub dotplot
ax = fig.add_subplot(gs[1, :])
genes = sorted({r[1] for r in hub})
cts = sorted({r[0] for r in hub}, key=lambda k: -ctn.get(k, 0))[:10]
M = np.zeros((len(genes), len(cts))); P = np.zeros_like(M)
for r in hub:
    if r[0] in cts and r[1] in genes:
        i = genes.index(r[1]); j = cts.index(r[0])
        M[i, j] = float(r[3]); P[i, j] = float(r[4])
sc = ax.scatter(np.repeat(np.arange(len(cts)), len(genes)),
                np.tile(np.arange(len(genes)), len(cts)),
                s=(P.T.ravel()) * 1.6, c=M.T.ravel(), cmap="Reds", vmin=0,
                vmax=np.percentile(M, 99), edgecolors="#78909C", linewidths=0.3)
ax.set_xticks(range(len(cts))); ax.set_xticklabels([c[:16] for c in cts], rotation=35, ha="right", fontsize=5.8)
ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
cb = plt.colorbar(sc, ax=ax, fraction=0.025, pad=0.01); cb.set_label("mean expression", fontsize=6)
cb.ax.tick_params(labelsize=5.5)
ax.set_title("Hub-gene expression across PBMC subsets", loc="left")
S.panel(ax, "D")
plt.tight_layout()
S.save(fig, f"{R}/figures/F5_singlecell")

# =====================================================================
# Figure 6 — Mendelian randomization
# =====================================================================
mr = [l.rstrip("\n").split("\t") for l in open(f"{R}/results/mr_multiSNP_fixed.tsv", encoding="utf-8")][1:]
rows = [(r[0], float(r[4]), float(r[5]), float(r[7]), int(r[1]), float(r[2])) for r in mr]
rows.sort(key=lambda x: x[1])
fig = plt.figure(figsize=(7.2, 4.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1], wspace=0.34)

ax = fig.add_subplot(gs[0, 0])
for i, (g, b, se, p, ns, f) in enumerate(rows):
    lo, hi = b - 1.96 * se, b + 1.96 * se
    c = S.PAL["up"] if b > 0 else S.PAL["down"]
    sig = p < 0.0029
    ax.plot([np.exp(lo), np.exp(hi)], [i, i], color=c, lw=1.0 if sig else 0.7,
            alpha=1.0 if sig else 0.55)
    ax.scatter([np.exp(b)], [i], s=22 if sig else 13, c=c, edgecolors="white",
               linewidths=0.4, zorder=3, alpha=1.0 if sig else 0.75)
ax.axvline(1, color="k", lw=0.6, ls="--")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([f"{g} ({ns})" for g, _, _, _, ns, _ in rows], fontsize=6)
ax.set_xlabel("OR for MDD per SD increase in expression")
ax.set_xscale("log")
ax.set_xticks([0.6,0.8,1.0,1.25,1.6,2.0])
ax.set_xticklabels(["0.6","0.8","1.0","1.25","1.6","2.0"])
ax.minorticks_off()
ax.set_title("cis-eQTL Mendelian randomization", loc="left")
ax.text(0.98, 0.02, "Bonferroni P < 0.0029", transform=ax.transAxes, ha="right", fontsize=5.6)
S.panel(ax, "A")

ax = fig.add_subplot(gs[0, 1])
ax.scatter([r[5] for r in rows], [-np.log10(r[3]) for r in rows], s=16,
           c=[S.PAL["up"] if r[1] > 0 else S.PAL["down"] for r in rows],
           edgecolors="white", linewidths=0.3)
for r in rows:
    if r[3] < 0.01:
        ax.annotate(r[0], (r[5], -np.log10(r[3])), fontsize=5.4, textcoords="offset points",
                    xytext=(4, 2), color=S.PAL["dark"])
ax.axhline(-np.log10(0.0029), color=S.PAL["accent"], lw=0.7, ls="--")
ax.set_xlabel("minimum instrument F"); ax.set_ylabel("−log10 P")
ax.set_title("Instrument strength", loc="left")
S.panel(ax, "B")
plt.tight_layout()
S.save(fig, f"{R}/figures/F6_mr")
print("done")

"""fig3.py — Figure 3: specificity controls for conventional enrichment

A  top-10 pathway overlap between XYS and six comparator formulae
B  size-matched convergence z-scores against two disease gene sets
C  pathway hit frequency in random draws from the herbal vs protein-coding background
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

ov = json.load(open(f"{R}/results/fig3_overlap.json", encoding="utf-8"))
zz = json.load(open(f"{R}/results/fig3_z.json", encoding="utf-8"))
rf = json.load(open(f"{R}/results/fig3_random_freq.json", encoding="utf-8"))

NAME = {"逍遥散": "XYS", "柴胡疏肝散": "CSS", "越鞠丸": "YJW", "酸枣仁汤": "SZRT",
        "天王补心丹": "TWBXD", "温胆汤": "WD", "半夏白术天麻汤": "BZTM"}

fig = plt.figure(figsize=(7.2, 7.2))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=1.15, wspace=0.46)

# ---------------- A -----------------
ax = fig.add_subplot(gs[0, :])
forms = ov["formulas"]; pws = ov["pathways"]; M = np.array(ov["matrix"], dtype=float)
order = sorted(range(len(forms)), key=lambda i: (forms[i] != "逍遥散", -M[i].sum()))
M = M[order]; forms = [forms[i] for i in order]
xys_row = forms.index("逍遥散")
keep = [j for j in range(len(pws)) if M[xys_row, j] > 0]
M = M[:, keep]; pws = [pws[j] for j in keep]
nr, nc = M.shape
for i in range(nr):
    for j in range(nc):
        if M[i, j] > 0:
            ax.scatter(j, i, s=52, c=S.PAL["accent"] if forms[i] == "XYS" else "#546E7A",
                       edgecolors="white", linewidths=0.5, zorder=3)
        else:
            ax.scatter(j, i, s=20, facecolors="none", edgecolors="#CFD8DC",
                       linewidths=0.7, zorder=2)
for i in range(nr):
    ax.text(nc - 0.30, i, f"{int(M[i].sum())}", va="center", ha="left",
            fontsize=8, color=S.PAL["dark"], fontweight="bold")
ax.text(nc - 0.30, -0.95, "n", va="center", ha="left", fontsize=7.5, color=S.PAL["dark"])
ax.set_xticks(range(nc))
ax.set_xticklabels([p[:30] for p in pws], rotation=32, ha="right", fontsize=7.2)
ax.set_yticks(range(nr)); ax.set_yticklabels([NAME.get(f, f) for f in forms], fontsize=8)
ax.set_xlim(-0.6, nc + 0.35); ax.set_ylim(nr - 0.4, -1.2)
ax.tick_params(axis="y", length=0)
ax.spines["left"].set_visible(False)
ax.set_title("XYS top-10 pathways recovered by every comparator", loc="left")
S.panel(ax, "A", x=-0.10, y=1.06)

# ---------------- B -----------------
ax = fig.add_subplot(gs[1, 0])
lab = ["XYS", "CSS", "YJW", "SZRT", "TWBXD", "WD", "BZTM"]
key = ["逍遥散", "柴胡疏肝散", "越鞠丸", "酸枣仁汤", "天王补心丹", "温胆汤", "半夏白术天麻汤"]
pri = [zz["primary"].get(k, {}).get("z", 0) for k in key]
sen = [zz["sensitivity"].get(k, {}).get("z", 0) for k in key]
x = np.arange(len(key)); w = 0.38
ax.bar(x - w/2, pri, w, color="#B0BEC5", label="primary DEG set")
ax.bar(x + w/2, sen, w, color=S.PAL["accent"], label="sensitivity DEG set")
ax.set_xticks(x); ax.set_xticklabels(lab, fontsize=5.8, rotation=30, ha="right")
ax.set_ylabel("convergence z")
ax.axhline(0, color="k", lw=0.5)
ax.set_ylim(0,4.9)
ax.legend(frameon=False, fontsize=5.6, loc="upper center", ncol=2, columnspacing=0.8)
ax.set_title("Size-matched convergence", loc="left")
S.panel(ax, "B")

# ---------------- C -----------------
ax = fig.add_subplot(gs[1, 1])
rf2 = sorted(rf, key=lambda d: -(d["herbal"] - d["genome"]))[:10]
names = [d["pathway"][:30] for d in rf2]
h = [d["herbal"] for d in rf2]; g = [d["genome"] for d in rf2]
y = np.arange(len(rf2))
ax.barh(y - 0.19, h, 0.36, color=S.ROLE_COLOR["Minister"], label="herbal background")
ax.barh(y + 0.19, g, 0.36, color=S.ROLE_COLOR["Courier"], label="protein-coding background")
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=5.2)
ax.invert_yaxis()
ax.set_xlabel("hits in top-10 (of 200 random sets)")
ax.legend(frameon=False, fontsize=5.4, loc="lower right")
ax.set_title("Random draws reproduce the signal", loc="left")
S.panel(ax, "C")

plt.tight_layout()
S.save(fig, f"{R}/figures/F3_specificity")

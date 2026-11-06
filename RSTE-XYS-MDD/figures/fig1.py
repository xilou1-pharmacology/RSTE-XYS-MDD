"""fig1.py — Figure 1: Role-Stratified Target Enrichment (RSTE)

A  targets per herb: all versus unique, coloured by role
B  pathways reaching FDR < 0.05 under whole-genome versus herbal background
C  recovery of seven established herb-target associations
D  role-resolved enrichment: BH-FDR and null-model P
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import sys

import fig_style as S
import numpy as np
import matplotlib.pyplot as plt

R = str(ROOT)
S.setup()

HERBS = [("Chaihu", "Monarch", 761, 69), ("Danggui", "Minister", 183, 14),
         ("Baishao", "Minister", 570, 18), ("Baizhu", "Assistant", 231, 15),
         ("Fuling", "Assistant", 70, 3), ("Gancao", "Assistant", 287, 14),
         ("Bohe", "Courier", 260, 29), ("Shengjiang", "Courier", 880, 225)]

fig = plt.figure(figsize=(7.2, 6.0))
gs = fig.add_gridspec(2, 2, hspace=0.58, wspace=0.40, height_ratios=[1, 1])

# ---------- A ----------
ax = fig.add_subplot(gs[0, 0])
y = np.arange(len(HERBS))
cols = [S.ROLE_COLOR[r] for _, r, _, _ in HERBS]
ax.barh(y, [t for _, _, t, _ in HERBS], color=cols, alpha=0.30, height=0.66)
ax.barh(y, [u for _, _, _, u in HERBS], color=cols, alpha=1.0, height=0.66)
for i, (_, _, t, u) in enumerate(HERBS):
    ax.text(t + 16, i, f"{t} / {u}", va="center", fontsize=7.0, color=S.PAL["dark"])
ax.set_yticks(y); ax.set_yticklabels([h for h, _, _, _ in HERBS], fontsize=7.6)
ax.invert_yaxis()
ax.set_xlabel("number of targets  (all / unique)"); ax.set_xlim(0, 1010)

ax.set_title("Targets per herb", loc="left")
S.panel(ax, "A", x=-0.30, y=1.10)

# ---------- B ----------
ax = fig.add_subplot(gs[0, 1])
sig = [165, 7]; tot = [234, 26]
x = np.arange(2)
ax.bar(x, sig, color=[S.PAL["neutral"], S.PAL["accent"]], width=0.52)
for i, (s, t) in enumerate(zip(sig, tot)):
    ax.text(i, s + 6, f"{s} / {t}\n({100*s/t:.0f}%)", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels(["whole\ngenome", "herbal\nuniverse"], fontsize=7.6)
ax.set_ylabel("pathways with FDR < 0.05"); ax.set_ylim(0, 205)
ax.set_title("Background substitution", loc="left")
S.panel(ax, "B", x=-0.32, y=1.10)

# ---------- C ----------
ax = fig.add_subplot(gs[1, 0])
ax.set_axis_off()
rows = [("Gancao", "HSD11B2"), ("Gancao", "CYP11B2"), ("Bohe + Shengjiang", "TRPM8"),
        ("Shengjiang", "TRPA1"), ("Shengjiang", "MAOA"), ("Shengjiang", "MAOB"),
        ("Danggui", "citrate cycle")]
for i, (h, t) in enumerate(rows):
    yy = 0.90 - i * 0.125
    ax.scatter([0.025], [yy], s=26, marker="o", c=S.ROLE_COLOR["Courier"],
               transform=ax.transAxes, clip_on=False, zorder=5)
    ax.text(0.09, yy, h, fontsize=7.6, transform=ax.transAxes)
    ax.text(0.52, yy, "→ " + t, fontsize=7.6, transform=ax.transAxes,
            family="monospace", color=S.PAL["dark"])
ax.text(0.02, 1.04, "Seven established associations recovered (7/7)", fontsize=8,
        transform=ax.transAxes, fontweight="bold")
S.panel(ax, "C", x=-0.04, y=1.16)

# ---------- D ----------
ax = fig.add_subplot(gs[1, 1])
names = ["Citrate\ncycle\n(Danggui)", "Steroid\nbiosynthesis\n(Gancao)", "Cytokine\nsignalling\n(Bohe)"]
fdr = [0.016, 0.006, 0.016]; nullp = [0.065, 0.438, 0.015]
x = np.arange(3); w = 0.34
ax.bar(x - w/2, fdr, w, color=S.ROLE_COLOR["Minister"], label="BH-FDR")
ax.bar(x + w/2, nullp, w, color=S.PAL["neutral"], label="null-model P")
ax.axhline(0.05, ls="--", lw=0.9, color=S.ROLE_COLOR["Monarch"])
ax.text(-0.62, 0.068, "0.05", fontsize=7, color=S.ROLE_COLOR["Monarch"])
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=7.0)
ax.set_ylabel("P value"); ax.set_ylim(0, 0.52); ax.set_xlim(-0.7, 2.6)
ax.legend(frameon=False, fontsize=7, loc="upper left")
ax.set_title("Role-resolved enrichment", loc="left", pad=8)
S.panel(ax, "D", x=-0.30, y=1.16)

S.save(fig, f"{R}/figures/F1_RSTE_scheme")

"""fig8_v3.py — Figure 8: molecular docking (final, with 3D poses)

A  binding-energy matrix
B  ranked binding energies
C  contact-residue counts
D  three predicted binding poses (AKT1–paeoniflorin, PTGS2–quercetin, TNF–paeoniflorin)
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import sys, json, collections

import fig_style as S
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

R = str(ROOT)
PYMOL = str(ROOT)
S.setup()

dk = json.load(open(f"{R}/data/dock/docking_results.json", encoding="utf-8"))
rows = [(r, lg, v) for r, lg, v in dk if v is not None]
recs = ["TNF", "PTGS2", "AKT1"]
ligs = ["Luteolin", "Quercetin", "Kaempferol", "Paeoniflorin"]
M = np.full((len(recs), len(ligs)), np.nan)
for r, lg, v in rows:
    M[recs.index(r), ligs.index(lg)] = v
cols = {"TNF": S.ROLE_COLOR["Monarch"], "PTGS2": S.ROLE_COLOR["Minister"],
        "AKT1": S.ROLE_COLOR["Courier"]}

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.62], hspace=0.12, wspace=0.55,
                      top=0.97, bottom=0.03, left=0.13, right=0.99)

# ---------- A ----------
ax = fig.add_subplot(gs[0, 0])
im = ax.imshow(M, cmap="YlOrRd_r", vmin=-10.2, vmax=-7.0, aspect="auto")
ax.set_xticks(range(4)); ax.set_xticklabels(["Lut", "Que", "Kae", "Pae"], fontsize=7.5)
ax.set_yticks(range(3)); ax.set_yticklabels(recs, fontsize=8)
for i in range(3):
    for j in range(4):
        ax.text(j, i, f"{M[i,j]:.1f}", ha="center", va="center", fontsize=7,
                color="white" if M[i, j] < -9.2 else "black")
ax.set_title("Binding energy", loc="left")
S.panel(ax, "A", x=-0.34, y=1.13)

# ---------- B ----------
ax = fig.add_subplot(gs[0, 1])
rows_s = sorted(rows, key=lambda r: r[2])
y = np.arange(len(rows_s))
ax.barh(y, [r[2] for r in rows_s], color=[cols[r[0]] for r in rows_s], height=0.72)
ax.set_yticks(y)
ax.set_yticklabels([f"{r[0][:4]}–{r[1][:4]}" for r in rows_s], fontsize=6.4)
ax.axvline(-7.0, ls="--", lw=0.8, color=S.PAL["dark"])
ax.set_xlabel("ΔG (kcal/mol)"); ax.set_xlim(-10.8, 0)
ax.set_title("Ranked complexes", loc="left")
S.panel(ax, "B", x=-0.42, y=1.13)

# ---------- C ----------
ax = fig.add_subplot(gs[0, 2])
ct = [l.rstrip("\n").split("\t") for l in
      open(f"{R}/data/dock/contacts.tsv", encoding="utf-8")][1:]
cnt = collections.Counter()
for r in ct:
    cnt[f"{r[0]}–{r[1][:4].capitalize()}"] += 1
keys = [f"{r[0]}–{r[1][:4]}" for r in rows_s if f"{r[0]}–{r[1][:4]}" in cnt]
vals = [cnt[k] for k in keys]
yy = np.arange(len(keys))
ax.barh(yy, vals, color=[cols[k.split("–")[0]] for k in keys], height=0.72)
ax.set_yticks(yy); ax.set_yticklabels(keys, fontsize=6.4)
ax.set_xlabel("contact residues (<4 Å)")
ax.invert_yaxis()
ax.set_title("Buried interface", loc="left")
S.panel(ax, "C", x=-0.42, y=1.13)

# ---------- D: three poses ----------
POSES = [("AKT1_vs_paeoniflorin_pocket.png", "AKT1 – paeoniflorin", "−9.91"),
         ("PTGS2_vs_quercetin_pocket.png", "PTGS2 – quercetin", "−9.83"),
         ("TNF_vs_paeoniflorin_pocket.png", "TNF – paeoniflorin", "−8.32")]
for k, (fn, label, dg) in enumerate(POSES):
    ax = fig.add_subplot(gs[1, k])
    img = mpimg.imread(f"{PYMOL}/{fn}")
    h, w = img.shape[:2]
    cy, cx = int(h * 0.055), int(w * 0.055)
    ax.imshow(img[cy:h-cy, cx:w-cx])
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.6); sp.set_color("#B0BEC5")
    ax.set_xlabel(f"{label}\nΔG {dg} kcal/mol", fontsize=7.6, labelpad=3)
    if k == 0:
        S.panel(ax, "D", x=-0.10, y=1.10)

S.save(fig, f"{R}/figures/F8_docking")

"""fig2_v4.py — Figure 2: composition and connectivity (three clean panels)

A  most connected constituents (coloured by herb role)
B  degree distributions
C  compounds per target versus PPI degree
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import sys, collections

import fig_style as S
import numpy as np
import matplotlib.pyplot as plt

R = str(ROOT)
S.setup()

rows = [l.rstrip("\n").split("\t") for l in
        open(f"{R}/data/processed/compound_target_edges.tsv", encoding="utf-8")][1:]
comp, targ, mol2herb, mol2name = collections.defaultdict(set), collections.defaultdict(set), {}, {}
for mol, name, herb, g in rows:
    comp[mol].add(g); targ[g].add(mol); mol2herb[mol] = herb; mol2name[mol] = name

ROLE = {"柴胡": "Monarch", "当归": "Minister", "白芍": "Minister",
        "白术": "Assistant", "茯苓": "Assistant", "甘草": "Assistant",
        "薄荷": "Courier", "生姜": "Courier"}
ppi = {}
for i, l in enumerate(open(f"{R}/data/processed/ppi_degree.tsv", encoding="utf-8")):
    if i == 0: continue
    p = l.rstrip().split("\t")
    if len(p) >= 2: ppi[p[0]] = int(p[1])

fig = plt.figure(figsize=(7.2, 3.2))
gs = fig.add_gridspec(1, 3, wspace=0.52, width_ratios=[1.25, 1, 1])

# A
ax = fig.add_subplot(gs[0, 0])
top = sorted(((m, len(v)) for m, v in comp.items()), key=lambda x: -x[1])[:14][::-1]
ax.barh(range(len(top)), [v for _, v in top],
        color=[S.ROLE_COLOR[ROLE[mol2herb[m]]] for m, _ in top], height=0.72)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([mol2name[m][:20] for m, _ in top], fontsize=7.4)
for i, (m, v) in enumerate(top):
    ax.text(v + 2, i, str(v), va="center", fontsize=7)
ax.set_xlabel("number of targets"); ax.set_xlim(0, max(v for _, v in top) * 1.20)
for lab in ["Monarch", "Minister", "Assistant", "Courier"]:
    ax.scatter([], [], s=22, c=S.ROLE_COLOR[lab], label=lab, edgecolors="white", linewidths=0.3)
ax.legend(frameon=False, loc="lower right", fontsize=6.8, handletextpad=0.3)
ax.set_title("Most connected constituents", loc="left")
S.panel(ax, "A", x=-0.42, y=1.14)

# B
ax = fig.add_subplot(gs[0, 1])
cd = [len(comp[m]) for m in comp]; td = [len(targ[g]) for g in targ]
bins = np.arange(0, 150, 8)
ax.hist(cd, bins=bins, color=S.ROLE_COLOR["Minister"], alpha=0.9,
        edgecolor="white", linewidth=0.4, label="targets / compound")
ax.hist(td, bins=bins, color=S.ROLE_COLOR["Courier"], alpha=0.6,
        edgecolor="white", linewidth=0.4, label="compounds / target")
ax.set_xlabel("degree"); ax.set_ylabel("count")
ax.set_ylim(0, 245)
ax.legend(frameon=False, loc="upper right", fontsize=6.8)
ax.set_title("Degree distributions", loc="left")
S.panel(ax, "B", x=-0.36, y=1.14)

# C
ax = fig.add_subplot(gs[0, 2])
gx = np.array([len(targ[g]) for g in targ])
gy = np.array([ppi.get(g, 0) for g in targ])
ok = gy > 0
ax.scatter(gx[ok], gy[ok], s=13, c=S.PAL["neutral"], alpha=0.85, edgecolors="none")
for g, off in [("PTGS2", (-36, 2)), ("AKT1", (-38, -4)), ("NFKB1", (8, -13)),
               ("TNF", (8, 4)), ("NCOA7", (-38, 6))]:
    if g in targ and ppi.get(g):
        ax.scatter(len(targ[g]), ppi[g], s=38, facecolors="none",
                   edgecolors=S.PAL["accent"], linewidths=1.1, zorder=4)
        ax.annotate(g, (len(targ[g]), ppi[g]), fontsize=6.8, color=S.PAL["accent"],
                    textcoords="offset points", xytext=off)
ax.set_xlabel("compounds per target"); ax.set_ylabel("PPI degree")
ax.text(0.42, 0.20, "r = 0.13", transform=ax.transAxes, fontsize=7.6, va="top")
ax.set_title("Two notions of hub", loc="left")
S.panel(ax, "C", x=-0.36, y=1.14)

S.save(fig, f"{R}/figures/F2_compound_target")

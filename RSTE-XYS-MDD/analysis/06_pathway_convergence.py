#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
06_pathway_convergence.py (v2) — 通路层面收敛检验（Enrichr GMT 版）
"""
import os, math
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
from collections import defaultdict

ROOT = str(ROOT)
GMT_DIR = os.path.join(ROOT, "data/raw/kegg")
OUT = os.path.join(ROOT, "results")
os.makedirs(OUT, exist_ok=True)

def read_gmt(path):
    sets = {}
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3: continue
            term = p[0].strip()
            genes = {x.strip().upper() for x in p[2:] if x.strip()}
            if genes: sets[term] = genes
    return sets

KEGG = read_gmt(os.path.join(GMT_DIR, "KEGG_2021_Human.gmt"))
GOBP = read_gmt(os.path.join(GMT_DIR, "GO_Biological_Process_2021.gmt"))
print(f"KEGG 通路 {len(KEGG)} | GO BP {len(GOBP)}")

def read_col(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for i, l in enumerate(f):
            if i == 0: continue
            p = l.rstrip("\n").split("\t")
            if p and p[0].strip(): out.append(p[0].strip().upper())
    return out

XY = set(read_col(os.path.join(ROOT, "data/processed/XYS_targets_OB30DL018.tsv")))
print("逍遥散靶点:", len(XY))

rows = []
with open(os.path.join(ROOT, "data/processed/deg_blood.tsv"), encoding="utf-8") as f:
    f.readline()
    for l in f:
        p = l.rstrip("\n").split("\t")
        if len(p) >= 6:
            try: rows.append((p[0].upper(), float(p[1]), float(p[4]), float(p[5])))
            except ValueError: pass

UNIVERSE = {r[0] for r in rows}
print("检测基因:", len(UNIVERSE))

DEG05 = {r[0] for r in rows if r[2] < 0.05}
DEG_LFC = {r[0] for r in rows if r[2] < 0.05 and abs(r[1]) >= 0.25}
DEG_01 = {r[0] for r in rows if r[2] < 0.01}
print(f"DEG P<0.05: {len(DEG05)} | P<0.05&|logFC|>=0.25: {len(DEG_LFC)} | P<0.01: {len(DEG_01)}")

def hyper(k, K, n, N):
    if n == 0 or K == 0: return 1.0
    lo = max(0, n + K - N)
    if k <= lo: return 1.0
    return sum(math.comb(K, x) * math.comb(N - K, n - x)
               for x in range(k, min(n, K) + 1)) / math.comb(N, n)

def enrich(query, sets, universe, min_size=5, max_size=500):
    q = query & universe
    res = []
    N = len(universe); n = len(q)
    for term, genes in sets.items():
        mem = genes & universe
        if len(mem) < min_size or len(mem) > max_size: continue
        k = len(mem & q)
        if k < 2: continue
        res.append([term, k, len(mem), n, hyper(k, len(mem), n, N)])
    res.sort(key=lambda r: r[4])
    m = len(res); prev = 1.0
    for i in range(m - 1, -1, -1):
        v = min(prev, res[i][4] * m / (i + 1)); prev = v
        res[i].append(min(1.0, v))
    return res

def save(res, name):
    with open(os.path.join(OUT, name), "w", encoding="utf-8") as f:
        f.write("term\thits\tsize\tquery\tP\tFDR\n")
        for r in res: f.write("\t".join(map(str, r)) + "\n")
    sig = sum(1 for r in res if r[5] < 0.05)
    print(f"  {name:28s} {len(res):4d} 条, FDR<0.05: {sig}")

print("\n=== 富集 ===")
RX_K = enrich(XY, KEGG, UNIVERSE);  save(RX_K, "pw_XYS_KEGG.tsv")
RD_K = enrich(DEG05, KEGG, UNIVERSE); save(RD_K, "pw_DEG05_KEGG.tsv")
RD2_K = enrich(DEG_LFC, KEGG, UNIVERSE); save(RD2_K, "pw_DEGlfc_KEGG.tsv")
RX_G = enrich(XY, GOBP, UNIVERSE);  save(RX_G, "pw_XYS_GOBP.tsv")
RD_G = enrich(DEG05, GOBP, UNIVERSE); save(RD_G, "pw_DEG05_GOBP.tsv")

def top_terms(res, q=0.05, top=30):
    s = {r[0] for r in res if r[5] < q}
    return s if s else {r[0] for r in res[:top]}

for name, a, b in [("KEGG: XYS vs DEG(P<0.05)", RX_K, RD_K),
                   ("KEGG: XYS vs DEG(P<0.05,|FC|>=.25)", RX_K, RD2_K),
                   ("GOBP: XYS vs DEG(P<0.05)", RX_G, RD_G)]:
    sa, sb = top_terms(a), top_terms(b)
    sh = sa & sb
    print(f"\n[{name}] XYS={len(sa)} DEG={len(sb)} 共享={len(sh)}")
    for t in sorted(sh)[:15]: print("    -", t)

print("\n=== XYS 靶点 Top10 KEGG ===")
for r in RX_K[:10]: print(f"  {r[0][:55]:55s} hits={r[1]:3d} P={r[4]:.2g} FDR={r[5]:.2g}")

print("\n=== DEG Top10 KEGG ===")
for r in RD_K[:10]: print(f"  {r[0][:55]:55s} hits={r[1]:3d} P={r[4]:.2g} FDR={r[5]:.2g}")


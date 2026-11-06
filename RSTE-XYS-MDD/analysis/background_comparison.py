"""background_comparison.py — 计算配方级靶点集在两种背景下的 KEGG 富集

统一规则（与 Methods 一致）：
  可测通路 = 背景宇宙内基因数在 5–400 之间，且被查询集命中至少 1 个
  超几何检验 + Benjamini–Hochberg（全部可测通路内校正）
输出：results/background_comparison.tsv
"""
import sys, math, csv
from pathlib import Path
from collections import defaultdict

_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(_ROOT))
from config import ROOT, RESULTS, RAW, PROCESSED


def read_gmt(p):
    out = {}
    for line in open(p, encoding="utf-8"):
        f = line.rstrip("\n").split("\t")
        if len(f) > 2:
            out[f[0]] = {x.upper() for x in f[2:] if x}
    return out


def hypergeom_tail(Kc, N, n, x):
    return sum(math.comb(Kc, i) * math.comb(N - Kc, n - i)
               for i in range(x, min(Kc, n) + 1)) / math.comb(N, n)


K = read_gmt(RAW / "kegg" / "KEGG_2021_Human.gmt")
all_kegg = set().union(*K.values())

herb_t = set()
for line in open(PROCESSED / "herb_targets_v2.tsv", encoding="utf-8"):
    f = line.rstrip("\n").split("\t")
    if len(f) > 1:
        herb_t.add(f[1].upper())
OMEGA = herb_t & all_kegg

XYS = {l.split("\t")[0].strip().upper() for l in
       open(PROCESSED / "XYS_targets_OB30DL018.tsv", encoding="utf-8")} - {""}
XYS &= all_kegg


def enrich(universe, label):
    rows = []
    for pw, genes in K.items():
        Kc = len(genes & universe)
        if not (5 <= Kc <= 400):
            continue
        hit = len(XYS & genes)
        if hit == 0:
            continue
        rows.append([pw, hit, Kc, hypergeom_tail(Kc, len(universe), len(XYS), hit)])
    rows.sort(key=lambda r: r[3])
    m = len(rows)
    for i, r in enumerate(rows):
        r.append(min(1.0, r[3] * m / (i + 1)))
    sig = sum(1 for r in rows if r[4] < 0.05)
    print(f"{label}: 查询 {len(XYS)} | 宇宙 {len(universe)} | 可测 {m} | FDR<0.05 {sig} "
          f"({100*sig/m:.0f}%)")
    for r in rows[:8]:
        print(f"    {r[0][:44]:44s} hit={r[1]:3d} K={r[2]:4d} FDR={r[4]:.2e}")
    return rows, sig


print("=== 背景对比：同一 241 靶点集 ===")
g, sg = enrich(all_kegg, "全基因组")
print()
h, sh = enrich(OMEGA, "中药靶点宇宙")

with open(RESULTS / "background_comparison.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["pathway", "hits", "K_genome", "P_genome", "FDR_genome",
                "K_herbal", "P_herbal", "FDR_herbal"])
    gd = {r[0]: r for r in g}
    hd = {r[0]: r for r in h}
    for pw in sorted(set(gd) | set(hd)):
        gg, hh = gd.get(pw), hd.get(pw)
        w.writerow([pw,
                    (gg or hh)[1],
                    gg[2] if gg else "", f"{gg[3]:.3e}" if gg else "", f"{gg[4]:.3e}" if gg else "",
                    hh[2] if hh else "", f"{hh[3]:.3e}" if hh else "", f"{hh[4]:.3e}" if hh else ""])
print("\n已写出 results/background_comparison.tsv")
print(f"摘要用数字: 全基因组 {sg}/{len(g)} ({100*sg/len(g):.0f}%), "
      f"中药宇宙 {sh}/{len(h)} ({100*sh/len(h):.0f}%)")

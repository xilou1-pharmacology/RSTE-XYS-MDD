#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
10_rste.py — Role-Stratified Target Enrichment (RSTE)  v1.1

算法
----
1. 角色归并：R_role = ∪ T(h), h ∈ role
2. 独有靶点：U(h) = T(h) \ ∪_{h'≠h} T(h')      （单药独有）
             U(role) = R_role \ ∪_{role'≠role} R_role'  （角色独有）
3. 背景替换：Ω = 中药靶点宇宙（2,012 基因），而非全基因组（8,078）
4. 超几何富集 + BH-FDR（每味药内校正）
5. 正向对照：已知药理靶点应落在对应药/角色的独有集合内
6. 零模型：等规模随机基因集，比较「FDR<0.05 的通路数」与「最小 P 的秩」

输出
----
results/rste_enrichment.tsv
results/rste_positive_controls.tsv
results/rste_null.tsv
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT


import collections, math, random, csv, os

ROOT = str(ROOT)
LIB = f"{ROOT}/data/processed/herb_targets_v2.tsv"
KEGG = f"{ROOT}/data/raw/kegg/KEGG_2021_Human.gmt"
OUT = f"{ROOT}/results"
os.makedirs(OUT, exist_ok=True)

ROLES = {"柴胡": "Monarch", "当归": "Minister", "白芍": "Minister",
         "白术": "Assistant", "茯苓": "Assistant", "甘草": "Assistant",
         "薄荷": "Courier", "生姜": "Courier"}

herb2t = collections.defaultdict(set)
allt = set()
with open(LIB, encoding="utf-8") as f:
    next(f)
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 2:
            g = p[1].strip().upper()
            herb2t[p[0].strip()].add(g)
            allt.add(g)

def read_gmt(path):
    d = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        x = line.rstrip("\n").split("\t")
        if len(x) >= 3:
            g = {i.strip().upper() for i in x[2:] if i.strip()}
            if g:
                d[x[0].strip()] = g
    return d

K = read_gmt(KEGG)
OMEGA = allt & set().union(*K.values())
print(f"Ω（中药靶点宇宙 ∩ KEGG）= {len(OMEGA)}")

# ---------- 独有集合 ----------
uniq_herb = {}
for h in ROLES:
    others = set().union(*[herb2t[x] for x in ROLES if x != h])
    uniq_herb[h] = herb2t[h] - others

by_role = collections.defaultdict(set)
for h, r in ROLES.items():
    by_role[r] |= herb2t[h]
uniq_role = {}
for r in set(ROLES.values()):
    others = set().union(*[by_role[x] for x in by_role if x != r])
    uniq_role[r] = by_role[r] - others

for h, r in ROLES.items():
    print(f"  {r:9s} {h:3s} | T={len(herb2t[h]):4d} | 单药独有 {len(uniq_herb[h]):4d} | 角色独有 {len(uniq_role[r]):4d}")

# ---------- 富集 ----------
def hyper(k, Kc, n, N):
    return sum(math.comb(Kc, x) * math.comb(N - Kc, n - x)
               for x in range(k, min(n, Kc) + 1)) / math.comb(N, n)

def enrich(q):
    q = q & OMEGA
    n = len(q)
    res = []
    if n < 5:
        return res
    for term, genes in K.items():
        m = genes & OMEGA
        if len(m) < 5 or len(m) > 400:
            continue
        k = len(m & q)
        if k < 2:
            continue
        res.append([term, k, len(m), n, hyper(k, len(m), n, len(OMEGA))])
    res.sort(key=lambda r: r[4])
    prev = 1.0
    for i in range(len(res) - 1, -1, -1):
        v = min(prev, res[i][4] * len(res) / (i + 1))
        prev = v
        res[i].append(min(1.0, v))
    return res

rows = []
obs = {}
for h, r in ROLES.items():
    res = enrich(uniq_herb[h])
    obs[h] = res
    for x in res:
        rows.append([h, r, x[0], x[1], x[2], x[3], f"{x[4]:.3e}", f"{x[5]:.3e}"])
    print(f"  {h:3s} ({r:9s}) 可测 {len(res):3d} | FDR<0.05: {sum(1 for x in res if x[5]<0.05):2d}"
          + (f" | 首位 {res[0][0][:38]} (FDR={res[0][5]:.2g})" if res else ""))

with open(f"{OUT}/rste_enrichment.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["herb", "role", "pathway", "hits", "pathway_size", "query_size", "P", "FDR"])
    w.writerows(rows)

# ---------- 正向对照 ----------
# level: "herb" = 单药独有; "role" = 角色独有
CONTROLS = [
    ("甘草", "HSD11B2", "herb", "11β-HSD2; glycyrrhizin inhibits this enzyme causing pseudohyperaldosteronism"),
    ("甘草", "CYP11B2", "herb", "Aldosterone synthase"),
    ("薄荷", "TRPM8", "role", "Menthol cold receptor (Courier role: mint + ginger)"),
    ("生姜", "TRPA1", "herb", "Gingerol receptor"),
    ("生姜", "MAOA", "herb", "Monoamine oxidase A"),
    ("生姜", "MAOB", "herb", "Monoamine oxidase B"),
    ("当归", "Citrate cycle (TCA cycle)", "pathway", "TCA cycle enrichment"),
]
pc, hit = [], 0
for herb, target, level, note in CONTROLS:
    if level == "pathway":
        rec = any(x[0] == target and x[5] < 0.05 for x in obs.get(herb, []))
        where = "通路富集 FDR<0.05"
    elif level == "herb":
        rec = target in uniq_herb[herb]
        where = "单药独有"
    else:
        role = ROLES[herb]
        rec = target in uniq_role[role]
        holders = [h for h in ROLES if target in herb2t[h]]
        where = f"角色独有（{'+'.join(holders)}）"
    hit += rec
    pc.append([herb, target, note, "YES" if rec else "no", where])

with open(f"{OUT}/rste_positive_controls.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["herb", "target", "known_pharmacology", "recovered", "level"])
    w.writerows(pc)

print(f"\n=== 正向对照: {hit}/{len(CONTROLS)} ===")
for x in pc:
    print(f"  {'✔' if x[3]=='YES' else '✘'} {x[0]:3s} → {x[1]:26s} {x[4]:26s} {x[2][:48]}")

# ---------- 零模型：FDR<0.05 通路数 ----------
random.seed(42)
pool = sorted(OMEGA)
null = []
for h, r in ROLES.items():
    n = len(uniq_herb[h] & OMEGA)
    if n < 5:
        continue
    counts = []
    for _ in range(200):
        sub = set(random.sample(pool, n))
        res = enrich(sub)
        counts.append(sum(1 for x in res if x[5] < 0.05))
    so = sum(1 for x in obs[h] if x[5] < 0.05)
    mu = sum(counts) / len(counts)
    sd = (sum((c - mu) ** 2 for c in counts) / (len(counts) - 1)) ** 0.5
    p = (sum(1 for c in counts if c >= so) + 1) / (len(counts) + 1)
    null.append([h, r, n, so, f"{mu:.1f}", f"{sd:.1f}", f"{p:.3f}"])
    print(f"  零模型 {h:3s}: n={n:3d} 观测FDR<0.05={so:2d} | 随机 {mu:.1f}±{sd:.1f} | P={p:.3f}")

with open(f"{OUT}/rste_null.tsv", "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["herb", "role", "n_unique", "obs_sig", "null_mean", "null_sd", "P"])
    w.writerows(null)

print("\n完成 -> results/rste_enrichment.tsv, rste_positive_controls.tsv, rste_null.tsv")

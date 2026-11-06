#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
07_mr_cis.py — cis-eQTL 孟德尔随机化（本地文件，不联网）

暴露：eQTLGen cis-eQTL（FDR<0.05，按 P 排序）
结局：FinnGen R12 F5_DEPRESSIO
对象：PPI top hub 基因

逻辑：每个基因取最强 cis-eQTL 作工具变量（单 SNP → Wald ratio）
      F = (beta/se)^2，F>10 才保留
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT

import gzip, math, os, csv, json

ROOT = str(ROOT)
EQTL = os.path.join(ROOT, "data/raw/mr/eqtlgen_cis_eQTL_FDR05.txt.gz")
FINN = os.path.join(ROOT, "data/raw/mr/finngen_R12_F5_DEPRESSIO.gz")
DEGREE = os.path.join(ROOT, "data/processed/ppi_degree.tsv")
OUT = os.path.join(ROOT, "results")
os.makedirs(OUT, exist_ok=True)

# hub 基因：PPI 度数 top 20
hubs = []
with open(DEGREE, encoding="utf-8") as f:
    rd = csv.DictReader(f, delimiter="\t")
    for i, r in enumerate(rd):
        if i >= 20: break
        hubs.append(r["gene"])
hubs = set(hubs)
print("hub 基因:", sorted(hubs))

# ---- 1. eQTLGen：每个基因取最小 P ----
best = {}
n = 0
with gzip.open(EQTL, "rt", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    while True:
        try:
            line = f.readline()
        except EOFError:
            print("  [warn] eQTLGen 文件截断，继续使用已读部分"); break
        if not line: break
        n += 1
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        gs = p[idx["GeneSymbol"]]
        if gs not in hubs: continue
        try:
            pv = float(p[idx["Pvalue"]])
        except ValueError:
            continue
        if gs not in best or pv < best[gs]["p"]:
            best[gs] = dict(p=pv, snp=p[idx["SNP"]], ea=p[idx["AssessedAllele"]],
                            oa=p[idx["OtherAllele"]],
                            beta=float(p[idx["Zscore"]]) / 1.0,
                            z=float(p[idx["Zscore"]]), n=int(float(p[idx["NrSamples"]])))
print(f"eQTLGen 扫描 {n:,} 行，命中 {len(best)}/{len(hubs)} 个 hub 基因")

# eQTLGen 只给 Z，需 beta 与 se：se ≈ 1/sqrt(2*p*(1-p)*N)，用 MAF 近似不可得
# 采用标准做法：beta = Z / sqrt(2*MAF*(1-MAF)*(N + Z^2))，MAF 缺失时用 se = 1/sqrt(N) 近似
for g, d in best.items():
    N = d["n"] if d["n"] > 0 else 30000
    # 保守近似：se ≈ 1/sqrt(N)，beta ≈ Z/sqrt(N)
    d["se"] = 1.0 / math.sqrt(N)
    d["beta"] = d["z"] * d["se"]
    d["f"] = (d["beta"] / d["se"]) ** 2

keep = {g: d for g, d in best.items() if d["f"] > 10}
print(f"F>10 的工具变量: {len(keep)}/{len(best)}")

# ---- 2. FinnGen：按 rsid 提取 ----
want = {d["snp"] for d in keep.values()}
print("需要从 FinnGen 提取的 SNP:", len(want))
fin = {}
with gzip.open(FINN, "rt", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().lstrip("#").rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        rs = p[idx["rsids"]]
        if rs in want:
            fin[rs] = dict(beta=float(p[idx["beta"]]), se=float(p[idx["sebeta"]]),
                           ea=p[idx["alt"]], oa=p[idx["ref"]], p=float(p[idx["pval"]]),
                           gene=p[idx["nearest_genes"]])
print(f"FinnGen 命中 {len(fin)}/{len(want)} 个 SNP")

# ---- 3. Wald ratio ----
rows = []
for g, d in keep.items():
    rs = d["snp"]
    if rs not in fin: continue
    o = fin[rs]
    # 对齐等位基因
    if d["ea"] == o["ea"] and d["oa"] == o["oa"]:
        b_out = o["beta"]
    elif d["ea"] == o["oa"] and d["oa"] == o["ea"]:
        b_out = -o["beta"]
    else:
        continue  # 等位基因不一致（可能是链翻转），跳过
    b_exp = d["beta"]
    if b_exp == 0: continue
    ratio = b_out / b_exp
    se_ratio = abs(o["se"] / b_exp)
    z = ratio / se_ratio if se_ratio > 0 else 0
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    rows.append([g, rs, round(d["beta"], 4), round(d["se"], 4), round(d["f"], 1),
                 round(b_out, 4), round(o["se"], 4), round(ratio, 4), round(se_ratio, 4),
                 round(z, 2), f"{p:.3g}"])

rows.sort(key=lambda r: float(r[10]))
with open(os.path.join(OUT, "mr_results.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["gene", "SNP", "beta_exposure", "se_exposure", "F",
                "beta_outcome", "se_outcome", "Wald_ratio", "se_ratio", "z", "P"])
    w.writerows(rows)

print(f"\n=== cis-MR 结果（{len(rows)} 个基因）===")
print(f"{'基因':10s} {'SNP':14s} {'F':>7s} {'OR':>7s} {'P':>10s}")
for r in rows:
    print(f"{r[0]:10s} {r[1]:14s} {r[4]:7.1f} {math.exp(float(r[7])):7.3f} {r[10]:>10s}")

sig = [r for r in rows if float(r[10]) < 0.05]
print(f"\nP<0.05: {len(sig)} 个")
print("完成 ->", os.path.join(OUT, "mr_results.tsv"))


#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
16_mr_persnp.py — 重建逐 SNP 的 cis-MR 数据（供敏感性分析用）

规则（与稿件 Methods 一致）：
  - 17 个 hub 基因
  - eQTLGen cis-eQTL（FDR<0.05）
  - 100 kb 距离剪枝（每个 100 kb 窗口保留 P 最小者）
  - 与 FinnGen R12 F5_DEPRESSIO 按 rsid 匹配、等位基因对齐
  - beta_exp = Z/sqrt(N), se_exp = 1/sqrt(N)（eQTLGen 只给 Z）

输出：results/mr_persnp.tsv
"""
import gzip, math, os, csv

ROOT = "D:/LKFlow2.0_Windows/Ning data/XYS_MDD_project"
EQTL = os.path.join(ROOT, "data/raw/mr/eqtlgen_cis_eQTL_FDR05.txt.gz")
FINN = os.path.join(ROOT, "data/raw/mr/finngen_R12_F5_DEPRESSIO.gz")
OUT = os.path.join(ROOT, "results")
GENES = ["GSK3B","PPARG","NFKB1","JUN","MMP9","EGF","PIK3CG","IL1B","FOS","IFNG",
         "MYC","TNF","CASP3","IL10","TP53","PTGS2","AKT1"]
WINDOW = 100_000  # 100 kb

# ---- 1. eQTLGen: 收集全部 cis-eQTL ----
cand = {g: [] for g in GENES}
n = 0
with gzip.open(EQTL, "rt", encoding="utf-8", errors="replace") as f:
    hdr = f.readline().rstrip("\n").split("\t")
    idx = {c: i for i, c in enumerate(hdr)}
    while True:
        try:
            line = f.readline()
        except EOFError:
            print("  [warn] eQTLGen 文件截断，使用已读部分"); break
        if not line: break
        n += 1
        p = line.rstrip("\n").split("\t")
        if len(p) < len(hdr): continue
        gs = p[idx["GeneSymbol"]]
        if gs not in cand: continue
        try:
            pv = float(p[idx["Pvalue"]]); z = float(p[idx["Zscore"]])
            pos = int(float(p[idx["SNPPos"]])); nN = int(float(p[idx["NrSamples"]]))
        except (ValueError, IndexError):
            continue
        cand[gs].append(dict(snp=p[idx["SNP"]], pos=pos, z=z, p=pv, n=nN,
                             ea=p[idx["AssessedAllele"]], oa=p[idx["OtherAllele"]]))
print(f"eQTLGen 扫描 {n:,} 行")
print("每个基因的 cis-eQTL 数:", {g: len(v) for g, v in cand.items() if v})

# ---- 2. 100 kb 剪枝 ----
def prune(snps, window=WINDOW):
    snps = sorted(snps, key=lambda d: d["p"])
    kept = []
    for s in snps:
        if all(abs(s["pos"] - k["pos"]) >= window for k in kept):
            kept.append(s)
    return kept

sel = {}
for g, v in cand.items():
    if v: sel[g] = prune(v)
print("\n剪枝后每个基因的 SNP 数:", {g: len(v) for g, v in sel.items()})
print("合计:", sum(len(v) for v in sel.values()))

# ---- 3. FinnGen ----
want = {s["snp"] for v in sel.values() for s in v}
print(f"\n需要提取的 SNP: {len(want)}")
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
                           ea=p[idx["alt"]], oa=p[idx["ref"]], p=float(p[idx["pval"]]))
print(f"FinnGen 命中 {len(fin)}/{len(want)}")

# ---- 4. 组装 ----
rows, dropped = [], {"no_outcome": 0, "allele_mismatch": 0, "bad_f": 0}
for g, v in sel.items():
    for s in v:
        if s["snp"] not in fin:
            dropped["no_outcome"] += 1; continue
        o = fin[s["snp"]]
        if s["ea"] == o["ea"] and s["oa"] == o["oa"]:
            by = o["beta"]
        elif s["ea"] == o["oa"] and s["oa"] == o["ea"]:
            by = -o["beta"]
        else:
            dropped["allele_mismatch"] += 1; continue
        se_exp = 1.0 / math.sqrt(s["n"]) if s["n"] > 0 else 1.0 / math.sqrt(30000)
        bx = s["z"] * se_exp
        fstat = (bx / se_exp) ** 2
        if fstat <= 10:
            dropped["bad_f"] += 1; continue
        rows.append(dict(gene=g, snp=s["snp"], pos=s["pos"], bx=round(bx, 6), sx=round(se_exp, 6),
                         by=round(by, 6), sy=round(o["se"], 6), F=round(fstat, 1),
                         p_eqtl=f"{s['p']:.3g}"))

rows.sort(key=lambda r: (r["gene"], r["pos"]))
with open(os.path.join(OUT, "mr_persnp.tsv"), "w", encoding="utf-8", newline="") as f:
    w = csv.DictWriter(f, delimiter="\t", fieldnames=["gene","snp","pos","bx","sx","by","sy","F","p_eqtl"])
    w.writeheader(); w.writerows(rows)

from collections import Counter
cnt = Counter(r["gene"] for r in rows)
print(f"\n最终 SNP: {len(rows)} 条，覆盖 {len(cnt)} 个基因")
print("  每个基因:", dict(cnt))
print("  丢弃:", dropped)
print("  F 最小值: %.1f | 最大值: %.1f" % (min(r["F"] for r in rows), max(r["F"] for r in rows)))
print("完成 -> results/mr_persnp.tsv")

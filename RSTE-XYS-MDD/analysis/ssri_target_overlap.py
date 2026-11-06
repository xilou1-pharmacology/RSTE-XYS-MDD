#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ssri_target_overlap.py — 逍遥散靶点 vs 抗抑郁药靶点 重叠分析

用法:
    python ssri_target_overlap.py --targets targets_all.tsv --out results

数据源: ChEMBL REST API (https://www.ebi.ac.uk/chembl/api/data)
输出:
    ssri_targets.tsv        每个药的靶点表
    overlap_stats.tsv       每个药与逍遥散的重叠统计
    gc_axis_overlap.tsv     GC 轴基因是否被命中
"""
import argparse, json, os, sys, time
import urllib.request, urllib.parse
from math import comb

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
CACHE = ".chembl_cache"
os.makedirs(CACHE, exist_ok=True)

DRUGS = {
    "fluoxetine":   "Fluoxetine",
    "sertraline":   "Sertraline",
    "paroxetine":   "Paroxetine",
    "citalopram":   "Citalopram",
    "escitalopram": "Escitalopram",
    "venlafaxine":  "Venlafaxine",
    "duloxetine":   "Duloxetine",
    "mirtazapine":  "Mirtazapine",
}

GC_AXIS = ["NR3C1", "FKBP5", "FKBP4", "HSP90AA1", "HSP90AB1",
           "NR3C2", "HSD11B1", "HSD11B2", "CRH", "CRHR1",
           "CRHR2", "POMC", "MC2R", "SERPINA6", "TSC22D3",
           "SGK1", "PER1", "TDO2", "TPH1", "TPH2"]


def get(url, cache_key):
    path = os.path.join(CACHE, cache_key + ".json")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = json.loads(r.read().decode("utf-8"))
            with open(path, "w", encoding="utf-8") as f:
                json.dump(data, f)
            return data
        except Exception as e:
            if attempt == 2:
                print(f"  ! {url} -> {e}", file=sys.stderr)
                return {}
            time.sleep(2)


def chembl_id(drug):
    q = urllib.parse.quote(drug)
    d = get(f"{CHEMBL}/molecule?pref_name__iexact={q}&format=json", f"mol_{drug}")
    mols = d.get("molecules", [])
    return mols[0]["molecule_chembl_id"] if mols else None


def targets_of(mol_id, drug):
    d = get(f"{CHEMBL}/activity?molecule_chembl_id={mol_id}"
            f"&target_type=SINGLE%20PROTEIN&limit=1000&format=json", f"act_{drug}")
    tids = {a["target_chembl_id"] for a in d.get("activities", []) if a.get("target_chembl_id")}
    genes = set()
    for tid in tids:
        t = get(f"{CHEMBL}/target/{tid}?format=json", f"tgt_{tid}")
        for comp in t.get("target_components", []) or []:
            for syn in comp.get("target_component_synonyms", []) or []:
                if syn.get("syn_type") == "GENE_SYMBOL":
                    genes.add(syn["component_synonym"].upper())
    return genes


def fisher(a, b, c, d):
    """a=both, b=only herb, c=only drug, d=neither"""
    n = a + b + c + d
    def p(x):
        return comb(a + b, x) * comb(c + d, a + c - x) / comb(n, a + c)
    obs = p(a)
    return sum(p(x) for x in range(a, min(a + b, a + c) + 1) if p(x) <= obs + 1e-12)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="targets_all.tsv")
    ap.add_argument("--out", default="results")
    ap.add_argument("--gene-col", default="gene")
    ap.add_argument("--formula-col", default="formula")
    ap.add_argument("--formula", default="逍遥散")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    import csv
    herb = set()
    with open(args.targets, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row[args.formula_col] == args.formula:
                herb.add(row[args.gene_col].upper())
    print(f"{args.formula} targets: {len(herb)}")

    universe = 20000  # 蛋白编码基因近似值
    rows, all_t = [], {}
    for key, name in DRUGS.items():
        mid = chembl_id(name)
        if not mid:
            print(f"[skip] {name}: no ChEMBL id")
            continue
        tg = targets_of(mid, key)
        all_t[name] = tg
        both = herb & tg
        a, b = len(both), len(herb - tg)
        c = len(tg - herb)
        d = universe - a - b - c
        p = fisher(a, b, c, d) if a else 1.0
        rows.append([name, len(tg), a, len(both) / max(1, len(herb | tg)),
                     f"{p:.3g}", ",".join(sorted(both))[:200]])
        print(f"{name:14s} targets={len(tg):5d} overlap={a:3d} p={p:.3g}")

    with open(os.path.join(args.out, "overlap_stats.tsv"), "w", encoding="utf-8") as f:
        f.write("drug\tn_targets\toverlap\tjaccard\tp_fisher\tshared_genes\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    # GC 轴是否被各药命中
    with open(os.path.join(args.out, "gc_axis_overlap.tsv"), "w", encoding="utf-8") as f:
        f.write("gene\tin_xiaoyao\t" + "\t".join(DRUGS.values()) + "\n")
        for g in GC_AXIS:
            line = [g, "Y" if g in herb else "N"]
            for name in DRUGS.values():
                line.append("Y" if g in all_t.get(name, set()) else "N")
            f.write("\t".join(line) + "\n")

    print(f"\nDone -> {args.out}")


if __name__ == "__main__":
    main()

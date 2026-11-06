"""fig3_data.py — 计算 Figure 3 数据（numpy 向量化）

输出 results/fig3_z.json, fig3_overlap.json, fig3_draws.json
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import numpy as np, json, collections

R = str(ROOT)
def gmt(p):
    s = {}
    for l in open(p, encoding="utf-8", errors="replace"):
        x = l.rstrip("\n").split("\t")
        if len(x) >= 3:
            g = {i.strip().upper() for i in x[2:] if i.strip()}
            if g:
                s[x[0].strip()] = g
    return s


K = gmt(f"{R}/data/raw/kegg/KEGG_2021_Human.gmt")
BGS = set().union(*K.values())
BG = sorted(BGS)
gi = {g: i for i, g in enumerate(BG)}
N = len(BG)

paths = [(t, g & BGS) for t, g in K.items() if 5 <= len(g & BGS) <= 500]
P = np.zeros((len(paths), N), dtype=np.int16)
for i, (t, g) in enumerate(paths):
    for x in g:
        P[i, gi[x]] = 1
print(f"背景 {N} | 合格通路 {len(paths)}")


def vec(genes):
    v = np.zeros(N, dtype=np.int16)
    for g in genes:
        if g in gi:
            v[gi[g]] = 1
    return v


def shared_mat(Dv, Q):
    cq = Q @ P.T
    cd = Dv @ P.T
    return (((cq >= 5) & (cd >= 5)).sum(1))


form = collections.defaultdict(set)
for i, l in enumerate(open(f"{R}/data/processed/formulas.tsv", encoding="utf-8")):
    p = l.rstrip("\n").split("\t")
    if i == 0:
        hdr = p
        continue
    idx = {c: j for j, c in enumerate(hdr)}
    form[p[idx.get("formula", 0)]].add(p[idx.get("gene", 1)].strip().upper())


def degset(f):
    s = set()
    for i, l in enumerate(open(f"{R}/data/processed/{f}", encoding="utf-8")):
        if i == 0:
            continue
        p = l.rstrip("\n").split("\t")
        try:
            if float(p[4]) < 0.05 and abs(float(p[1])) >= 0.25:
                s.add(p[0].upper())
        except Exception:
            pass
    return s & BGS


def enrich_top(q, k=5):
    q = q & BGS
    v = vec(q)
    c = v @ P.T
    pv = []
    n = int(v.sum())
    for i, (t, g) in enumerate(paths):
        kk = int(c[i])
        if kk < 2:
            pv.append(1.0)
            continue
        from math import comb
        m = len(g)
        s = sum(comb(m, x) * comb(N - m, n - x) for x in range(kk, min(n, m) + 1)) / comb(N, n)
        pv.append(s)
    order = np.argsort(pv)[:k]
    return [paths[i][0] for i in order]


rng = np.random.default_rng(2024)
out = {}
for name, fn, size, reps in [("primary", "deg_128vs64.tsv", 156, 50),
                             ("sensitivity", "deg_blood.tsv", 209, 50)]:
    D = degset(fn)
    Dv = vec(D)
    row = {}
    for f, s in form.items():
        q = s & BGS
        if len(q) < size:
            continue
        zs, ps = [], []
        for _ in range(reps):
            sub = set(rng.choice(BG, size, replace=False).tolist())
            obs = shared_mat(Dv, vec(sub)[None, :])[0]
            Q = np.zeros((500, N), dtype=np.int16)
            for j in range(500):
                Q[j, rng.choice(N, size, replace=False)] = 1
            null = shared_mat(Dv, Q)
            mu, sd = null.mean(), null.std(ddof=1)
            zs.append((obs - mu) / sd if sd else 0)
            ps.append((np.sum(null >= obs) + 1) / 501)
        row[f] = dict(z=round(float(np.mean(zs)), 2), P=round(float(np.mean(ps)), 3))
    out[name] = row

json.dump(out, open(f"{R}/results/fig3_z.json", "w"), ensure_ascii=False, indent=1)
for k in out:
    print(f"\n=== {k} ===")
    for f, v in sorted(out[k].items(), key=lambda x: -x[1]["z"]):
        print(f"  {f:12s} z={v['z']:5.2f}  P={v['P']:.3f}")
    rk = sorted(out[k], key=lambda x: -out[k][x]["z"]).index("逍遥散") + 1
    print(f"  逍遥散排名 {rk}/{len(out[k])}")

# 随机抽样通路家族
rng2 = np.random.default_rng(3)
univ = set()
for i, l in enumerate(open(f"{R}/data/processed/herb_targets_v2.tsv", encoding="utf-8")):
    if i == 0:
        continue
    p = l.rstrip("\n").split("\t")
    if len(p) >= 2:
        univ.add(p[1].strip().upper())
univ &= BGS
herb = [enrich_top(set(rng2.choice(sorted(univ), 300, replace=False).tolist()), 5) for _ in range(3)]
geno = [enrich_top(set(rng2.choice(BG, 300, replace=False).tolist()), 5) for _ in range(3)]
json.dump({"herbal": herb, "genome": geno}, open(f"{R}/results/fig3_draws.json", "w"), ensure_ascii=False)
print("\nherbal:", herb[0])
print("genome:", geno[0])

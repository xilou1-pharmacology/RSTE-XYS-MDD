"""fig3_convergence.py — 通路收敛的规模匹配置换检验（numpy 向量化）

观测统计量：从方剂靶点集中随机抽 SIZE 个，计算与 DEG 集的共享通路数
零分布：从背景中随机抽 SIZE 个基因，重复
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
P = np.zeros((N, len(paths)), dtype=np.float32)
for i, (t, g) in enumerate(paths):
    for x in g:
        P[gi[x], i] = 1.0
print(f"背景 {N} | 通路 {len(paths)}", flush=True)


def vec(genes):
    v = np.zeros(N, dtype=np.float32)
    for g in genes:
        j = gi.get(g)
        if j is not None:
            v[j] = 1.0
    return v


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


rng = np.random.default_rng(2024)
out = {}
for name, fn, size, reps in [("primary", "deg_128vs64.tsv", 156, 30),
                             ("sensitivity", "deg_blood.tsv", 209, 30)]:
    Dv = vec(degset(fn))
    cd = Dv @ P
    keep = cd >= 5          # 只保留 DEG 命中 >=5 的通路
    row = {}
    for f, s in form.items():
        q = sorted((s & BGS))
        if len(q) < size:
            continue
        zs, ps = [], []
        for r in range(reps):
            sub = vec(set(rng.choice(q, size, replace=False).tolist()))
            obs = int((((sub @ P) >= 5) & keep).sum())
            Q = np.zeros((400, N), dtype=np.float32)
            idxs = rng.random((400, N)).argpartition(size, axis=1)[:, :size]
            np.put_along_axis(Q, idxs, 1.0, axis=1)
            null = (((Q @ P) >= 5) & keep).sum(1)
            mu, sd = null.mean(), null.std(ddof=1)
            zs.append((obs - mu) / sd if sd else 0)
            ps.append((int((null >= obs).sum()) + 1) / 401)
        row[f] = dict(z=round(float(np.mean(zs)), 2), P=round(float(np.mean(ps)), 3))
        print(f"  {name:12s} {f:12s} z={row[f]['z']:5.2f} P={row[f]['P']:.3f}", flush=True)
    out[name] = row

json.dump(out, open(f"{R}/results/fig3_z.json", "w"), ensure_ascii=False, indent=1)
for k in out:
    rk = sorted(out[k], key=lambda x: -out[k][x]["z"]).index("逍遥散") + 1
    print(f"{k}: 逍遥散排名 {rk}/{len(out[k])}", flush=True)
print("DONE", flush=True)

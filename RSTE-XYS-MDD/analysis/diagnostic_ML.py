"""diagnostic_ML.py — 诊断性能(ROC/AUC) + 机器学习驱动 hub 基因与靶点集的重叠

Requires data/processed/expr_blood_full.tsv, generated from GEO GSE98793 by
R/05_build_deg_GSE98793.R (48 MB; not redistributed here).

回答两个问题：
  1. 靶点集分数能否区分 MDD 与对照（AUC）？
  2. 用三种机器学习方法从数据里选出的 hub 基因，是否与 XYS 靶点集重叠？
"""
import sys
from pathlib import Path
_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import ROOT
import sys, warnings
from pathlib import Path
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegressionCV
from sklearn.svm import LinearSVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import RFE
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.preprocessing import StandardScaler
from scipy import stats

ROOT = Path(str(ROOT))
P = ROOT / "data/processed"

expr = pd.read_csv(P / "expr_blood_full.tsv", sep="\t", index_col=0)
pheno = pd.read_csv(P / "pheno_GSE98793.tsv", sep="\t")
print("表达矩阵:", expr.shape, "| 表型:", pheno.shape, "| 列:", list(pheno.columns)[:5])

# 对齐
pheno = pheno.set_index("gsm")
common = [s for s in expr.columns if s in pheno.index]
expr = expr[common]
y = pheno.loc[common, "group"].astype(str).str.upper().str.startswith("CASE").astype(int).values
print(f"样本 {len(common)} | 病例 {y.sum()} | 对照 {(y==0).sum()}")

# ---------- 1. 靶点集分数 ----------
tgt = {l.split("\t")[0].strip().upper() for l in
       open(P / "XYS_targets_OB30DL018.tsv", encoding="utf-8")} - {""}
present = [g for g in expr.index if str(g).upper() in tgt]
print(f"靶点集在矩阵中可检出: {len(present)}/{len(tgt)}")
z = expr.loc[present].apply(lambda r: (r - r.mean()) / (r.std() + 1e-9), axis=1)
score = z.mean(axis=0).values

auc = roc_auc_score(y, score)
# 置换检验
rng = np.random.default_rng(7)
null = [roc_auc_score(rng.permutation(y), score) for _ in range(2000)]
p_perm = (np.sum(np.array(null) >= auc) + 1) / 2001
# DeLong 近似用 bootstrap
boot = []
for _ in range(1000):
    idx = rng.integers(0, len(y), len(y))
    if len(np.unique(y[idx])) < 2:
        continue
    boot.append(roc_auc_score(y[idx], score[idx]))
ci = (np.percentile(boot, 2.5), np.percentile(boot, 97.5))
print(f"\n【靶点集分数诊断性能】AUC = {auc:.3f} (95% CI {ci[0]:.3f}–{ci[1]:.3f}), "
      f"置换 P = {p_perm:.4f}")

# ---------- 2. 机器学习 hub 基因 ----------
deg = pd.read_csv(P / "deg_128vs64.tsv", sep="\t")
gcol = "gene" if "gene" in deg.columns else deg.columns[0]
pcol = [c for c in deg.columns if c.lower() in ("p.value", "pvalue", "p")][0]
genes = deg.loc[deg[pcol] < 0.05, gcol].astype(str).str.upper().tolist()
genes = [g for g in genes if g in expr.index]
print(f"\nDEG (P<0.05) 用于机器学习: {len(genes)} 个基因")

X = expr.loc[genes].T.values
X = StandardScaler().fit_transform(X)
cv = StratifiedKFold(5, shuffle=True, random_state=7)

sel = {}
try:
    lasso = LogisticRegressionCV(Cs=20, penalty="l1", solver="liblinear", cv=5,
                                 scoring="roc_auc", max_iter=5000, random_state=7)
    lasso.fit(X, y)
    sel["LASSO"] = {genes[i] for i in np.where(np.abs(lasso.coef_[0]) > 0)[0]}
except Exception as e:
    print("LASSO 失败:", e)

try:
    svc = LinearSVC(penalty="l1", dual=False, C=0.05, max_iter=20000, random_state=7)
    svc.fit(X, y)
    sel["SVM-RFE"] = {genes[i] for i in np.argsort(-np.abs(svc.coef_[0]))[:100]}
except Exception as e:
    print("SVM 失败:", e)

try:
    rf = RandomForestClassifier(n_estimators=500, random_state=7, n_jobs=-1)
    rf.fit(X, y)
    imp = rf.feature_importances_
    sel["RandomForest"] = {genes[i] for i in np.argsort(-imp)[:100]}
except Exception as e:
    print("RF 失败:", e)

print("\n【三种方法选出的 hub 基因数】")
for k, v in sel.items():
    print(f"  {k:14s} {len(v):4d}")

if len(sel) >= 2:
    from functools import reduce
    cons = reduce(lambda a, b: a & b, sel.values())
    print(f"\n  ∩ 三法交集: {len(cons)} 个")
    for k, v in sel.items():
        ov = v & tgt
        exp = len(v) * len(tgt) / len(expr.index)
        print(f"  {k} 与 XYS 靶点集重叠: {len(ov)} (期望 {exp:.1f}) "
              f"| 富集倍数 {len(ov)/exp:.2f}" if exp > 0 else "")
    if cons:
        ovc = cons & tgt
        expc = len(cons) * len(tgt) / len(expr.index)
        print(f"  交集 hub 与 XYS 靶点重叠: {len(ovc)}/{len(cons)} (期望 {expc:.1f})")
        # 超几何
        p_h = stats.hypergeom.sf(len(ovc) - 1, len(expr.index), len(tgt), len(cons))
        print(f"  超几何 P = {p_h:.3g}")
        print("  重叠基因:", sorted(ovc)[:20])

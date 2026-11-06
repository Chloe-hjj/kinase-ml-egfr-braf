from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold
from features import fingerprints, scaffold

OUT = Path("results/figures")
OUT.mkdir(parents = True, exist_ok = True)

to_do = [("EGFR potency (wild type)", "data/processed/EGFR_wt_potency.csv", "pIC50"),
        ("BRAF potency (wild type)", "data/processed/BRAF_wt_potency.csv", "pIC50"),
        ("EGFR delta (L858R,T790M)", "data/processed/EGFR_L858R_T790M_pairs.csv", "delta"),
        ("BRAF delta (V600E)", "data/processed/BRAF_V600E_pairs.csv", "delta")]

fig, axes = plt.subplots(2, 2, figsize = (8, 8))
for ax, (title, path, target) in zip(axes.ravel(), to_do):
    df = pd.read_csv(path)
    X = fingerprints(df["smiles"].tolist())
    y = df[target].to_numpy()
    groups = df["smiles"].map(scaffold).to_numpy()
    pred = np.zeros_like(y, dtype = float)
    for tr, te in GroupKFold(n_splits = 5).split(X, y, groups):
        model = RandomForestRegressor(n_estimators = 300, n_jobs = -1, random_state = 0).fit(X[tr], y[tr])
        pred[te] = model.predict(X[te])
    lo, hi = min(y.min(), pred.min()), max(y.max(), pred.max())
    ax.scatter(y, pred, s = 6, alpha = 0.3, color = "#0072B2", edgecolors = "none")
    ax.plot([lo, hi], [lo, hi], color = "black", lw = 1)
    ax.set_title(f"{title}\nR² = {r2_score(y, pred):.2f}, Spearman = {spearmanr(y, pred)[0]:.2f}, n = {len(y)}", fontsize = 9)
    ax.set_xlabel("Measured")
    ax.set_ylabel("Predicted (out-of-fold)")
    print(title, "done")
fig.suptitle("Random forest, scaffold-grouped 5-fold cross-validation", fontsize = 11)
fig.tight_layout()
fig.savefig(OUT / "fig4_pred_vs_true.png", dpi = 150)

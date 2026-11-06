from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = Path("results/figures")
OUT.mkdir(parents = True, exist_ok = True)

TASKS = ["EGFR potency", "BRAF potency", "EGFR delta (L858R,T790M)", "BRAF delta (V600E)"]
LABELS = ["EGFR\npotency", "BRAF\npotency", "EGFR delta\n(L858R,T790M)", "BRAF delta\n(V600E)"]
COLORS = {"mean baseline": "#999999", "random forest": "#0072B2", "MLP": "#D55E00"}

res = pd.concat([pd.read_csv("results/rf_cv_folds.csv"),
                 pd.read_csv("results/mlp_cv_folds.csv")], ignore_index=True)
stats = res.groupby(["task", "split", "model"])["r2"].agg(["mean", "std"]).reset_index()

def get(task, split, model):
    row = stats[(stats.task == task) & (stats.split == split) & (stats.model == model)]
    return (row["mean"].iloc[0], row["std"].iloc[0]) if len(row) else (np.nan, np.nan)

# Figure 1: delta distributions
fig, axes = plt.subplots(1, 2, figsize = (9, 3.5), sharey = False)
for ax, (path, title) in zip(axes, [("data/processed/EGFR_L858R_T790M_pairs.csv", "EGFR: L858R,T790M vs wild type"),
                                    ("data/processed/BRAF_V600E_pairs.csv", "BRAF: V600E vs wild type")]):
    d = pd.read_csv(path)["delta"]
    ax.hist(d, bins = 30, color = "#0072B2", alpha = 0.85, edgecolor = "white")
    ax.axvline(0, color = "black", lw = 1)
    ax.axvline(d.median(), color = "#D55E00", lw = 1.5, ls = "--", label = f"median {d.median():+.2f}")
    ax.set_title(f"{title} (n = {len(d)})", fontsize = 10)
    ax.set_xlabel("delta pIC50 (mutant - wild type)")
    ax.legend(frameon = False)
axes[0].set_ylabel("Number of molecules")
fig.tight_layout()
fig.savefig(OUT / "fig1_delta_distribution.png", dpi = 150)
plt.close(fig)

# Figure 2: model comparison on scaffold-grouped splits
fig, ax = plt.subplots(figsize = (9, 4))
models = ["mean baseline", "random forest", "MLP"]
w = 0.26
for i, m in enumerate(models):
    vals = [get(t, "scaffold", m) for t in TASKS]
    ax.bar(np.arange(len(TASKS)) + (i - 1) * w, [v[0] for v in vals], w,
           yerr = [v[1] for v in vals], capsize = 3, color = COLORS[m], label = m)
ax.axhline(0, color = "black", lw = 0.8)
ax.set_xticks(np.arange(len(TASKS)))
ax.set_xticklabels(LABELS)
ax.set_ylabel("R² on held-out scaffolds")
ax.set_title("Model comparison, scaffold-grouped splits (mean ± std across splits)", fontsize = 10)
ax.legend(frameon = False)
fig.tight_layout()
fig.savefig(OUT / "fig2_model_comparison.png", dpi = 150)
plt.close(fig)

# Figure 3: optimism of random splits (random forest)
fig, ax = plt.subplots(figsize = (9, 4))
for i, (split, col) in enumerate([("random", "#56B4E9"), ("scaffold", "#0072B2")]):
    vals = [get(t, split, "random forest") for t in TASKS]
    ax.bar(np.arange(len(TASKS)) + (i - 0.5) * 0.35, [v[0] for v in vals], 0.35,
           yerr = [v[1] for v in vals], capsize = 3, color = col, label = f"{split} split")
ax.set_xticks(np.arange(len(TASKS)))
ax.set_xticklabels(LABELS)
ax.set_ylabel("R² (random forest)")
ax.set_title("Random splits overestimate performance on new chemistry", fontsize=10)
ax.legend(frameon = False)
fig.tight_layout()
fig.savefig(OUT / "fig3_random_vs_scaffold.png", dpi = 150)
plt.close(fig)
print("saved figures to", OUT)

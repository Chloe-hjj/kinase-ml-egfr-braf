import numpy as np, pandas as pd, torch
from torch import nn
from scipy.stats import spearmanr
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupShuffleSplit, ShuffleSplit

from features import fingerprints
from splits import get_fold

def make_model(n_in):
    return nn.Sequential(nn.Linear(n_in, 512), nn.ReLU(), nn.Dropout(0.2), nn.Linear(512, 128), nn.ReLU(), nn.Dropout(0.2), nn.Linear(128, 1))

def train_mlp(Xtr, ytr, Xva, yva, seed, max_epochs = 200, patience = 20, lr = 1e-3, wd = 1e-5, bs = 128):
    torch.manual_seed(seed)
    mu = float(ytr.mean())
    Xtr_t = torch.tensor(Xtr, dtype = torch.float32)
    ytr_t = torch.tensor(ytr - mu, dtype = torch.float32).unsqueeze(1)
    Xva_t = torch.tensor(Xva, dtype = torch.float32)
    yva_t = torch.tensor(yva - mu, dtype = torch.float32).unsqueeze(1)
    model = make_model(Xtr.shape[1])
    opt = torch.optim.Adam(model.parameters(), lr = lr, weight_decay = wd)
    loss_fn = nn.MSELoss()
    best, wait, best_state = float("inf"), 0, None
  
    for epoch in range(max_epochs):
        model.train()
        perm = torch.randperm(len(Xtr_t))
  
        for i in range(0, len(perm), bs):
            b = perm[i:i + bs]
            opt.zero_grad()
            loss_fn(model(Xtr_t[b]), ytr_t[b]).backward()
            opt.step()
        model.eval()
   
        with torch.no_grad():
            val = loss_fn(model(Xva_t), yva_t).item()
        
        if val < best - 1e-4:
            best, wait = val, 0
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
        else:
            wait += 1
            if wait >= patience:
                break

        model.load_state_dict(best_state)
        model.eval()

        return model, mu

def split_val(tr, groups_tr, seed):
    if groups_tr is not None:
        a, b = next(GroupShuffleSplit(n_splits = 1, test_size = 0.1, random_state = seed).split(tr, groups = groups_tr))
    else:
        a, b = next(ShuffleSplit(n_splits = 1, test_size = 0.1, random_state = seed).split(tr))
    return tr[a], tr[b]

to_do  = [("EGFR potency", "data/processed/EGFR_wt_potency.csv", "pIC50"),
        ("BRAF potency", "data/processed/BRAF_wt_potency.csv", "pIC50"),
        ("EGFR delta (L858R,T790M)", "data/processed/EGFR_L858R_T790M_pairs.csv", "delta"),
        ("BRAF delta (V600E)", "data/processed/BRAF_V600E_pairs.csv", "delta")]

if __name__ == "__main__":
    rows = []

    for label, path, target in to_do:
        df = pd.read_csv(path)
        X = fingerprints(df["smiles"].tolist())
        y = df[target].to_numpy()
      
        for split in ["scaffold", "random"]:
            for k, (tr, te, g) in enumerate(get_fold(df, split, seed = 0)):
                tr_fit, tr_val = split_val(tr, g, seed = k)
                model, mu = train_mlp(X[tr_fit], y[tr_fit], X[tr_val], y[tr_val], seed = k)
                with torch.no_grad():
                    pred = model(torch.tensor(X[te], dtype = torch.float32)).squeeze(1).numpy() + mu
                rows.append(dict(task = label, split = split, model = "MLP", fold = k,
                                 rmse = mean_squared_error(y[te], pred) ** 0.5,
                                 r2 = r2_score(y[te], pred),
                                 spearman = spearmanr(y[te], pred)[0]))

        print(label, split, "done")

    res = pd.DataFrame(rows)
    res.to_csv("results/mlp_cv_folds.csv", index = False)
    
    print(res.groupby(["task", "split", "model"])[["rmse", "r2", "spearman"]].agg(["mean", "std"]).round(2).to_string())



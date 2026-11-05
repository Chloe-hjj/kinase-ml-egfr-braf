import numpy as np, pandas as pd
from scipy.stats import spearmanr
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GroupKFold, KFold

from features import fingerprints, scaffold
from splits import get_fold

def eva(df, target, split, seed = 0):
    x = fingerprints(df["smiles"].tolist())
    y = df[target].to_numpy()
    rows = []

    for k, (tr, te, _) in enumerate(get_fold(df, split, seed)):
        models = [("mean baseline", DummyRegressor()),
                  ("random forest", RandomForestRegressor(n_estimators = 300, n_jobs=-1, random_state=seed))]

        for name, model in models:
            model.fit(x[tr], y[tr]); 
            pred = model.predict(x[te])
            rows.append(dict(split = split, model = name, fold = k, 
rmse = mean_squared_error(y[te], pred)**0.5, r2 = r2_score(y[te], pred), 
spearman = spearmanr(y[te], pred)[0] if name == "random forest" else np.nan))

    return pd.DataFrame(rows)

to_do = [("EGFR potency", "data/processed/EGFR_wt_potency.csv", "pIC50"),
        ("BRAF potency", "data/processed/BRAF_wt_potency.csv", "pIC50"),
        ("EGFR delta (L858R,T790M)", "data/processed/EGFR_L858R_T790M_pairs.csv", "delta"),
        ("BRAF delta (V600E)", "data/processed/BRAF_V600E_pairs.csv", "delta")]

if __name__ == "__main__":
    allres = []
    for label, path, target in to_do:
        df = pd.read_csv(path)
        for split in ["scaffold", "random"]:
            r = eva(df, target, split); r.insert(0, "task", label); allres.append(r)
            print(label, split, "done")

    res = pd.concat(allres)
    res.to_csv("results/rf_cv_folds.csv", index = False)
    
    print(res.groupby(["task", "split", "model"])[["rmse","r2", "spearman"]].agg(["mean", "std"]).round(2).to_string())

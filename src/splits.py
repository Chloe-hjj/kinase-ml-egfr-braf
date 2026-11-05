import numpy as np
from sklearn.model_selection import GroupShuffleSplit, KFold
from features import scaffold

def get_fold(df, split, seed = 0):
    idx = np.arange(len(df))
    if split == "scaffold":
        groups = df["smiles"].map(scaffold).to_numpy()
        gss = GroupShuffleSplit(n_splits = 10, test_size = 0.2, random_state = seed)
        return [(tr,te, groups[tr]) for tr, te in gss.split(idx, groups = groups)]

    kf = KFold(n_splits = 5, shuffle = True, random_state = seed)
    return [(tr, te, None) for tr, te in kf.split(idx)]


       

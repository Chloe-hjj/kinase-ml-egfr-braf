import numpy as np, pandas as pd
from make_pairs import std_smiles

def build_potency(name, assay_type = "B"):
    df = pd.read_csv(f"data/raw/{name}_ic50_raw.csv", low_memory=False)
    df = df[(df["standard_relation"] == "=") & (df["standard_value"] > 0) & (df["assay_type"] == assay_type) & (df["assay_variant_mutation"].isna())]
    df["pIC50"] = 9 - np.log10(df["standard_value"].astype(float))
    out = df.groupby("molecule_chembl_id").agg(pIC50 = ("pIC50", "median"), smiles =("canonical_smiles", "first")).reset_index()
    out["smiles"] = out["smiles"].map(std_smiles)
    
    return out.dropna(subset = ["smiles"]).drop_duplicates("smiles")


if __name__ == "__main__":
    for name in ["EGFR", "BRAF"]:
        d = build_potency(name)
        d.to_csv(f"data/processed/{name}_wt_potency.csv", index = False)
       
        print(name, len(d), "molecules")

   

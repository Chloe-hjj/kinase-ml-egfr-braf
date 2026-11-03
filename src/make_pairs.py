import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

chooser = rdMolStandardize.LargestFragmentChooser()

def std_smiles(smi):
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    return Chem.MolToSmiles(chooser.choose(mol)) if mol else None

def build_pairs(name, mutation, assay_type="B"):
    df = pd.read_csv(f"data/raw/{name}_ic50_raw.csv", low_memory=False)
    df = df[(df["standard_relation"] == "=") & (df["standard_value"] >0)
            & (df["assay_type"] == assay_type) & (df["assay_variant_mutation"] != "UNDEFINED MUTATION")].copy()

    df["pIC50"] = 9 - np.log10(df["standard_value"].astype(float))

    wt = df[df["assay_variant_mutation"].isna()]
    mut = df[df["assay_variant_mutation"] == mutation]
    key = ["molecule_chembl_id", "document_chembl_id"]
    wt_m = wt.groupby(key)["pIC50"].median().rename("pIC50_wt").reset_index()
    mut_m = mut.groupby(key)["pIC50"].median().rename("pIC50_mut").reset_index()

    pairs = wt_m.merge(mut_m, on=key)
    pairs = (pairs.groupby("molecule_chembl_id")[["pIC50_wt", "pIC50_mut"]]
                  .median().reset_index())          
    pairs["delta"] = pairs["pIC50_mut"] - pairs["pIC50_wt"]

    smi = df.drop_duplicates("molecule_chembl_id").set_index("molecule_chembl_id")["canonical_smiles"]
    pairs["smiles"] = pairs["molecule_chembl_id"].map(smi).map(std_smiles)
    pairs = pairs.dropna(subset = ["smiles"]).drop_duplicates("smiles")
    return pairs

if __name__ == "__main__":
    for name, mut, out in [("EGFR", "L858R,T790M", "EGFR_L858R_T790M"),
                           ("BRAF", "V600E", "BRAF_V600E")]:
        p = build_pairs(name, mut)
        p.to_csv(f"data/processed/{out}_pairs.csv", index = False)
        
        print(out, len(p), "molecules")
        print(p["delta"].describe().round(2), "\n")





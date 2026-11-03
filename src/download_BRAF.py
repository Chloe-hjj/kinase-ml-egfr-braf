# In case of being able to download only EGFR with download.py, use this for BRAF.  
import pandas as pd
from chembl_webresource_client.new_client import new_client

def show_candi(name):
    for t in new_client.target.search(name):
        if t["organism"] == "Homo sapiens" and t["target_type"] =="SINGLE PROTEIN":
            print(name, t["target_chembl_id"], t["pref_name"])

TARGETS = {"BRAF": "CHEMBL5145"}

def download(name, chembl_id):
    acts = new_client.activity.filter(
        target_chembl_id = chembl_id, standard_type = "IC50", standard_units = "nM")
    df = pd.DataFrame(list(acts))
    df.to_csv(f"data/raw/{name}_ic50_raw.csv", index=False)
    print(name, len(df), "rows")

if __name__ == "__main__":
    for name, cid in TARGETS.items():
        download(name, cid)

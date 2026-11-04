import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

_gen = rdFingerprintGenerator.GetMorganGenerator(radius = 2, fpSize = 2048)

def fp_array(smiles):
    mol = Chem.MolFromSmiles(smiles) 
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {smiles}")
    fp = _gen.GetFingerprint(mol)
    arr = np.zeros((2048), dtype = np.uint8)
    DataStructs.ConvertToNumpyArray(fp, arr)

    return arr

def fingerprints(smiles_list):
    return np.vstack([fp_array(s) for s in smiles_list])

def scaffold(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
 
   

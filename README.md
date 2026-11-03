# Kinase potency and resistance-risk prediction for EGFR and BRAF (RDKit + scikit-learn + PyTorch)

A small project that:
(1) builds a baseline model of kinase inhibitor potency from chemical structure, and
(2) extends it to a simple proxy for **resistance risk**: how much a compound loses potency
against a clinically relevant mutant compared with the wild-type kinase.

The emphasis is on honest evaluation, using scaffold-grouped cross-validation, and on stating
clear limits.

> **Status:** In progress. A learning project on cheminformatics and machine learning for
> kinase drug discovery. After reading about the ANA-KIN doctoral network and its work on
> kinase drug discovery beyond the ATP site, I wanted to apply what I learned in my master's
> in computational biology to a related problem.

## Questions
1. Can models trained on known inhibitors predict the potency of **chemically new**
   compounds for EGFR and BRAF, and how much does the evaluation split change the answer?
2. Can the loss of potency against a mutant (a resistance proxy) be predicted from structure?

## Data
- Source: ChEMBL (IC50 values in nM for human EGFR, CHEMBL203, and BRAF, CHEMBL5145).
- Filters: exact values only (`=` relation), **biochemical assays only** (ChEMBL assay type B),
  entries labelled "UNDEFINED MUTATION" removed.
- Cleaning: pIC50 = 9 - log10(IC50 in nM), SMILES standardised with RDKit (largest fragment),
  duplicates merged by median.
- **Wild type** means rows with no annotated mutation, so some may be unannotated mutants.
- **Resistance proxy:** for molecules measured on both wild-type EGFR and the
  L858R,T790M double mutant (a known resistance mutant) **in the same paper**,
  `delta_pIC50 = pIC50(mutant) - pIC50(wild type)`. Negative values mean a loss of potency
  on the mutant.
- **BRAF** is analysed in the same way with V600E, but V600E is the drug target mutation,
  so this is a measure of **mutant selectivity**, not resistance.

| Target | Raw activities | Exact (`=`) | Paired molecules (same paper, biochemical) | Median delta (IQR) |
|---|---|---|---|---|
| EGFR (WT vs L858R,T790M) | 25,105 | 19,646 | 621 | +0.62 (-0.28 to +1.29) |
| BRAF (WT vs V600E) | 10,831 | 9,105 | 302 | +0.31 (+0.05 to +0.68) |

Unique molecules for the wild-type potency task: *to be updated*.

## Methods
- **Features:** Morgan fingerprints (radius 2, 2048 bits), RDKit.
- **Evaluation:** 5-fold cross-validation where molecules sharing a Bemis-Murcko scaffold are
  kept in the same fold (scaffold-grouped), compared with a random 5-fold split.
- **Models:**
  - Mean-value baseline.
  - Classical machine learning: random forest (scikit-learn).
  - Neural network: small MLP (PyTorch, Adam, MSE loss, early stopping), *planned*.
- **Metrics:** RMSE, R², Spearman correlation, reported as mean and standard deviation
  across folds.

## Results
*To be updated.*

### Potency baseline (wild-type pIC50)
| Target | Model | Split | RMSE | R² | Spearman |
|---|---|---|---|---|---|
| EGFR | Random forest | scaffold | | | |
| EGFR | MLP | scaffold | | | |
| BRAF | Random forest | scaffold | | | |
| BRAF | MLP | scaffold | | | |

### Resistance proxy (delta pIC50)
| Task | Model | Split | RMSE | R² | Spearman |
|---|---|---|---|---|---|
| EGFR delta | Random forest | scaffold | | | |
| BRAF delta | Random forest | scaffold | | | |

## Limitations
- Public ChEMBL data are noisy: IC50 values come from different assays and labs.
- Most data describe ATP-site inhibitors, so the models mostly learn ATP-site chemistry
  and do not predict allosteric or interface binding.
- The resistance target is a proxy, not a validated resistance measure, and it depends on
  compounds measured in both forms.
- The EGFR set is biased towards mutant-selective inhibitors, because many compounds were
  designed against T790M.
- Wild-type rows are those with no annotated mutation, so some may be unannotated mutants.
- The paired sets are small (621 and 302 molecules), so results will be noisy.
- No prospective experimental validation was done.

## Scheme
- [x] Download ChEMBL data for EGFR and BRAF
- [x] Explore data and count wild-type/mutant pairs (go/no-go)
- [x] Build paired datasets
- [ ] Wild-type potency dataset
- [ ] Fingerprints and scaffold-grouped cross-validation
- [ ] Random forest baseline
- [ ] PyTorch MLP
- [ ] Interpretation (substructure highlighting)
- [ ] Reproducibility check from a fresh clone

## How to reproduce
*To be updated once the scripts are finished and tested from a fresh clone.*

```bash
conda env create -f environment.yml
conda activate kinase-ml
python src/download.py       # ChEMBL download (can fail when the server is busy; retry or download the kinases separately!)
python src/make_pairs.py
```

## About
Hyunjeong (Chloe) Jung, M.Sc. Computational Biology (UPM).
https://www.linkedin.com/in/hjungchloe/

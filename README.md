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

Unique molecules for the wild-type potency task:
- EGFR 9341 molecules
- BRAF 2969 molecules

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

## Results (random forest)

### Potency baseline (wild-type pIC50)
| Target | Model | Split | RMSE | R² | Spearman |
|---|---|---|---|---|---|
| EGFR | Mean baseline | scaffold | 1.39 ± 0.10 | -0.01 ± 0.02 | n/a |
| EGFR | Random forest | scaffold | 0.87 ± 0.05 | 0.60 ± 0.06 | 0.79 ± 0.04 |
| EGFR | Random forest | random | 0.74 ± 0.04 | 0.72 ± 0.02 | 0.85 ± 0.01 |
| BRAF | Mean baseline | scaffold | 1.33 ± 0.05 | -0.03 ± 0.03 | n/a |
| BRAF | Random forest | scaffold | 0.70 ± 0.07 | 0.71 ± 0.07 | 0.84 ± 0.03 |
| BRAF | Random forest | random | 0.60 ± 0.05 | 0.79 ± 0.03 | 0.89 ± 0.02 |

### Resistance proxy (delta pIC50)
| Task | Model | Split | RMSE | R² | Spearman |
|---|---|---|---|---|---|
| EGFR delta | Mean baseline | scaffold | 1.20 ± 0.14 | -0.07 ± 0.11 | n/a |
| EGFR delta | Random forest | scaffold | 0.83 ± 0.16 | 0.46 ± 0.22 | 0.70 ± 0.13 |
| EGFR delta | Random forest | random | 0.71 ± 0.11 | 0.64 ± 0.10 | 0.80 ± 0.06 |
| BRAF delta | Mean baseline | scaffold | 0.61 ± 0.16 | -0.11 ± 0.10 | n/a |
| BRAF delta | Random forest | scaffold | 0.52 ± 0.12 | 0.17 ± 0.24 | 0.49 ± 0.21 |
| BRAF delta | Random forest | random | 0.44 ± 0.06 | 0.46 ± 0.06 | 0.61 ± 0.04 |

MLP results: *to be updated*.

### What the resuts show??
- Potency can be predicted from structure, even for new scaffolds. The random forest clearly
  is better than the mean baseline for both kinases.
- Random splits are pretty optimistic. Moving from a random to a scaffold-grouped split lowers
  R² from 0.72 to 0.60 for EGFR and from 0.79 to 0.71 for BRAF, and raises the error.
- The EGFR resistance proxy has learnable signal on new scaffolds (R² 0.46, Spearman 0.70),
  with large variation between splits.
- The BRAF V600E selectivity signal is weak on new scaffolds (R² 0.17 ± 0.24). The delta
  values vary little (SD about 0.6), so they are close to the experimental noise level, and the
  model generalises poorly to new chemistry.

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
- [x] Wild-type potency dataset
- [x] Fingerprints and scaffold-grouped cross-validation
- [x] Random forest baseline
- [ ] PyTorch MLP
- [ ] Interpretation (substructure highlighting)
- [ ] Reproducibility check from a fresh clone

## How to reproduce
```bash
conda env create -f environment.yml
conda activate kinase-ml
python src/download.py (src/download_BRAF.py)      # ChEMBL download (can fail when the server is busy; retry or download the kinases separately!)
python src/make_pairs.py
python src/make_potency.py
python src/cv_rf.py
```

## About
Hyunjeong (Chloe) Jung, M.Sc. Computational Biology (UPM).
https://www.linkedin.com/in/hjungchloe/

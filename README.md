# Kinase potency and resistance-risk prediction for EGFR and BRAF (RDKit + scikit-learn + PyTorch)

A small project that:
(1) builds a baseline model of kinase inhibitor potency from chemical structure, and
(2) extends it to a simple proxy for **resistance risk**: how much a compound loses potency
against a clinically relevant mutant compared with the wild-type kinase.

The emphasis is on honest evaluation, using scaffold splits, and on stating clear limits.

> **Status:** In progress. A learning project on cheminformatics and machine learning for
> kinase drug discovery. After reading about the ANA-KIN doctoral network and its work on
> kinase drug discovery beyond the ATP site, I wanted to apply what I learned in my master's
> in computational biology to a related problem.

## Questions
1. Can models trained on known inhibitors predict the potency of **chemically new**
   compounds for EGFR and BRAF, and how much does the evaluation split change the answer?
2. Can the loss of potency against a mutant (a resistance proxy) be predicted from structure?

## Data
- Source: ChEMBL, IC50 measurements (nM), binding/functional assays.
- Cleaning: exact values only (`=` relation), pIC50 = 9 - log10(IC50 in nM), SMILES
  standardised with RDKit, duplicates merged by median.
- **Resistance proxy:** for compounds measured on both the wild-type kinase and a mutant
  (EGFR mutation(s) / BRAF V600E), the target is
  `delta_pIC50 = pIC50(mutant) - pIC50(wild type)`. Mutant and wild-type assays will be
  identified from assay annotations and descriptions.

| Target | Activities | Unique molecules | Paired WT/mutant compounds |
|---|---|---|---|
| EGFR | [ ] | [ ] | [ ] |
| BRAF | [ ] | [ ] | [ ] |

## Methods
- **Features:** Morgan fingerprints (radius 2, 2048 bits), RDKit.
- **Split:** Bemis-Murcko scaffold split (80/10/10), with a random split for comparison.
- **Models:**
  - Classical machine learning: random forest (scikit-learn).
  - Neural network: small MLP (PyTorch, Adam, MSE loss, early stopping).
- **Metrics:** RMSE, R², Spearman correlation on the held-out test set.

## Results
*To be updated.*

### Potency baseline
| Target | Model | Split | RMSE | R² | Spearman |
|---|---|---|---|---|---|
| EGFR | Random forest | scaffold | [ ] | [ ] | [ ] |
| EGFR | MLP | scaffold | [ ] | [ ] | [ ] |
| BRAF | Random forest | scaffold | [ ] | [ ] | [ ] |
| BRAF | MLP | scaffold | [ ] | [ ] | [ ] |

### Resistance proxy (delta pIC50)
| Target | Model | Split | RMSE | R² | Spearman |
|---|---|---|---|---|---|
| [ ] | [ ] | scaffold | [ ] | [ ] | [ ] |

## Limitations
- Public ChEMBL data are noisy: IC50 values come from different assays and labs.
- Most data describe ATP-site inhibitors, so the models mostly learn ATP-site chemistry
  and do not predict allosteric or interface binding.
- The resistance target is a proxy, not a validated resistance measure, and it depends on
  compounds measured in both forms.
- No prospective experimental validation was done.

## How to reproduce
*To be added once the scripts are finished and tested from a fresh clone.*

## About
Hyunjeong (Chloe) Jung, M.Sc. Computational Biology (UPM).
https://www.linkedin.com/in/hjungchloe/

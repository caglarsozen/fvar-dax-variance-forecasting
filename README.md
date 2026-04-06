# Leak-safe functional forecasting of daily variance for DAX

This repository contains the locked manuscript pipeline for one-step-ahead forecasting of DAX Rogers--Satchell variance using leak-safe functional predictors and benchmark models.

## Main configuration

- Target: DAX (`^GDAXI`)
- Main paper: `PRUNE095 + ref20 + DAX_CLOSE`
- Appendix robustness: `PRUNE095 + ref20 + GLOBAL_CLOSE`

## Repository structure

```text
.
├─ R/
│  └─ fvar_variance_pipeline_locked.R
├─ manuscript/
│  └─ paper.tex
├─ README.md
├─ LICENSE
├─ .gitignore
└─ CITATION.cff

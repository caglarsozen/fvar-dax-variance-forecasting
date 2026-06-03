# Leak-safe functional forecasting of daily DAX variance

This repository contains replication materials for the manuscript:

**Leak-safe functional forecasting of daily variance: window-to-curve predictors for DAX Rogers--Satchell variance**

The main empirical pipeline is provided in:


R/01_run_full_pipeline.R


## Main manuscript configuration

The locked manuscript configuration in the script is:


OUT_BASE   <- "outfvar_DAX_FINAL_PRUNE095_LOCKED_MINORREV_ENET"
TARGET_RAW <- "^GDAXI"
TARGET_NAME <- "DAX"

CFG_MAIN   <- "PRUNE095"
REFIT_MAIN <- 20L

AVAIL_MAIN <- "DAX_CLOSE"
AVAIL_ROB  <- "GLOBAL_CLOSE"


The main paper reports the conservative **DAX-Close** convention. The **Global-Close** convention is used as supplementary timing robustness.

## What the pipeline does

The script implements:

- Yahoo Finance data retrieval;
- Rogers--Satchell DAX variance construction from OHLC data;
- strict intersection-calendar alignment;
- no imputation or forward filling;
- DAX-Close and Global-Close information-availability protocols;
- DEV-only smoothing, scaling, screening, and design choices;
- Elastic Net lag-stack benchmark;
- minimal ablation comparing raw lag-stack regularization with the smoothed fVAR-X/PLS representation;
- recursive expanding-window forecasting;
- COMMON-aligned TEST evaluation;
- QLIKE, MAE, and RMSE tables;
- Diebold--Mariano tests with Newey--West HAC correction;
- moving-block bootstrap confidence intervals;
- residual diagnostics;
- manuscript-ready CSV, TeX, and figure outputs.

## How to run

From the repository root:


source("R/00_install_packages.R")
source("R/01_run_full_pipeline.R")
source("R/99_session_info.R")


or:


source("R/run_all.R")

## Expected output folder

Running the pipeline creates:


outfvar_DAX_FINAL_PRUNE095_LOCKED_MINORREV_ENET/


with subfolders for data, tables, TeX tables, plots, and debug files.

## Expected headline interpretation

Under the baseline DAX-Close convention and COMMON-aligned TEST evaluation:

1. HAR-X(logVIX) has the lowest average QLIKE.
2. HAR ranks second.
3. Elastic Net lag stack ranks third.
4. fVAR-X(+logv_t) is the strongest functional specification.
5. The functional family is informative and competitive, but not dominant.
6. HAR and HAR-X(logVIX) have the clearest multiplicity-robust inferential support relative to GARCH(1,1).

## Reproducibility notes

The credibility of the empirical design depends on the exact timing, filtering, alignment, and preprocessing rules. The repository is therefore organized so that readers can audit:

- raw data retrieval;
- strict intersection calendar;
- DAX-Close and Global-Close timing maps;
- DEV/TEST split;
- DEV-only smoothing and Elastic Net screening;
- recursive fitting and component selection;
- COMMON-aligned TEST construction;
- QLIKE, MAE, RMSE, DM, bootstrap, and diagnostic outputs.

## Data

All market data are publicly retrieved from Yahoo Finance. No proprietary data are required.

## License

Please choose the final license before making the repository public.

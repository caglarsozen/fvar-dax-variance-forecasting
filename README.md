# Leak-safe functional forecasting of daily DAX variance

This repository contains replication materials for the manuscript:

**Leak-safe functional forecasting of daily variance: window-to-curve predictors for DAX Rogers--Satchell variance**

The main empirical replication script is:


fvar_variance_pipeline_locked_v2.R


## Main manuscript configuration

The locked manuscript configuration in the script is:


OUT_BASE    <- "outfvar_DAX_FINAL_PRUNE095_LOCKED_MINORREV_ENET"

TARGET_RAW  <- "^GDAXI"
TARGET_NAME <- "DAX"

CFG_MAIN    <- "PRUNE095"
REFIT_MAIN  <- 20L

AVAIL_MAIN  <- "DAX_CLOSE"
AVAIL_ROB   <- "GLOBAL_CLOSE"
```

The main paper reports the conservative **DAX-Close** information-availability convention. The more permissive **Global-Close** convention is retained as a supplementary timing-robustness comparison.

## What the pipeline does

The script implements the full empirical design used in the manuscript:

* Yahoo Finance data retrieval;
* Rogers--Satchell DAX variance construction from OHLC data;
* strict intersection-calendar alignment;
* no imputation or forward filling;
* DAX-Close and Global-Close information-availability protocols;
* DEV-only smoothing, scaling, screening, and design choices;
* Elastic Net lag-stack benchmark;
* minimal ablation comparing a raw regularized lag stack with the smoothed fVAR-X/PLS representation;
* recursive expanding-window forecasting;
* COMMON-aligned TEST evaluation;
* QLIKE, MAE, and RMSE tables;
* Diebold--Mariano tests with Newey--West HAC correction;
* moving-block bootstrap confidence intervals;
* residual diagnostics;
* manuscript-ready CSV, TeX, and figure outputs.

## How to run

From the repository root, run:


source("fvar_variance_pipeline_locked_v2.R")


The script downloads the required public market data and creates the manuscript output folder automatically.

## Expected output folder

Running the pipeline creates:


outfvar_DAX_FINAL_PRUNE095_LOCKED_MINORREV_ENET/


This folder contains the generated data objects, tables, TeX files, figures, and diagnostic outputs used in the manuscript.

## Expected headline interpretation

Under the baseline DAX-Close convention and COMMON-aligned TEST evaluation:

1. HAR-X(logVIX) has the lowest average QLIKE.
2. HAR ranks second.
3. Elastic Net lag stack ranks third.
4. fVAR-X(+logv_t) is the strongest functional specification.
5. The functional family is informative and competitive, but not dominant.
6. HAR and HAR-X(logVIX) have the clearest multiplicity-robust inferential support relative to GARCH(1,1).

## Reproducibility notes

The credibility of the empirical design depends on the exact timing, filtering, alignment, and preprocessing rules. The repository is organized so that readers can audit:

* raw data retrieval;
* strict intersection-calendar construction;
* DAX-Close and Global-Close timing maps;
* DEV/TEST split;
* DEV-only smoothing and Elastic Net screening;
* recursive fitting and component selection;
* COMMON-aligned TEST construction;
* QLIKE, MAE, RMSE, Diebold--Mariano, bootstrap, and diagnostic outputs.

## Data

All market data are publicly retrieved from Yahoo Finance. No proprietary data are required.

## License

This repository is distributed under the license provided in the `LICENSE` file.

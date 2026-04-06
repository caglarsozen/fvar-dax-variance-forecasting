###############################################################################
# fVAR Variance Forecasting 
# Target: DAX
#
# MAIN PAPER:        PRUNE095 + ref20 + DAX_CLOSE
# APPENDIX ROBUST.:  PRUNE095 + ref20 + GLOBAL_CLOSE
#
# Included:
# (1) Leak-safe preprocessing: DEV-only smoothing / EN / scaling / component selection
# (2) Information availability protocol:
#       - DAX_CLOSE    = MAIN
#       - GLOBAL_CLOSE = APPENDIX ROBUSTNESS
# (3) Stronger benchmarks: HAR, HAR-X(logVIX), GARCH, EGARCH, GJR-GARCH
# (4) Inference rigor: DM one-sided + two-sided + Holm/BH; NW-lag sensitivity
# (5) COMMON-aligned evaluation: strict intersection calendar; no LOCF
# (6) Protocol comparisons:
#       - selected functional vs GARCH
#       - full model set, own common sets
#       - relative-to-GARCH
# (7) Diagnostics:
#       - multi-model residual diagnostics (ARCH-LM, Ljung-Box on z_t^2)
# (8) Locked manuscript configuration:
#       - no design sweep selection in final reporting
#       - baseline fixed at PRUNE095 / ref20 / DAX_CLOSE
#
# Corrections in this version:
#   * fixed two-block table rbind issue
#   * fixed retained-variable table truncation
#   * fixed bootstrap B=10 consistency between main text and appendix sensitivity
###############################################################################

options(stringsAsFactors = FALSE)
set.seed(20250105)

# ============================= USER SETTINGS ==================================
OUT_BASE    <- "outfvar_DAX_FINAL_PRUNE095_LOCKED"

TARGET_RAW  <- "^GDAXI"
TARGET_NAME <- "DAX"

start_date  <- as.Date("2019-01-02")
end_date    <- as.Date("2025-12-31")

W_CURVE     <- 20L
GRID_EVAL   <- 30L

CFG_MAIN    <- "PRUNE095"
REFIT_MAIN  <- 20L

AVAIL_MAIN  <- "DAX_CLOSE"
AVAIL_ROB   <- "GLOBAL_CLOSE"

BOOT_B_GRID   <- c(5L, 10L, 20L)
DM_LAG_GRID   <- c(3L, 5L, 10L)

TEST_RATIO    <- 0.25
VAL_RATIO_DEV <- 0.20

EN_ALPHA  <- 0.5
DM_NW_LAG <- 5L

USE_BSPLINE_SMOOTH    <- TRUE
NBASIS_GRID           <- c(10L, 12L, 15L, 18L, 20L)
LAMBDA_GRID           <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)
GCV_SUBSAMPLE_WINDOWS <- 60L
REP_SERIES_GCV        <- "GC"

PLS_MAX_COMP <- 8L
PLS_NSEG     <- 5L

PLOT_DPI <- 600

# ============================= LABEL HELPERS ==================================
CFG_LABELS <- c(
  "CORE"     = "Core candidate construction",
  "PRUNE095" = "PRUNE095 candidate construction",
  "FULL"     = "Full candidate construction"
)

cfg_pretty <- function(x){
  if(x %in% names(CFG_LABELS)) unname(CFG_LABELS[x]) else x
}

# ============================= PACKAGES =======================================
REQUIRED_PACKAGES <- c(
  "data.table", "quantmod", "xts", "zoo", "glmnet",
  "pls", "forecast", "ggplot2", "fda", "rugarch"
)

check_required_packages <- function(pkgs = REQUIRED_PACKAGES){
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing_pkgs) > 0L){
    stop(
      paste0(
        "Missing required packages: ",
        paste(missing_pkgs, collapse = ", "),
        "\nInstall them first, for example via:\n",
        "install.packages(c(",
        paste(sprintf('\"%s\"', missing_pkgs), collapse = ", "),
        "), repos = \"https://cloud.r-project.org\")"
      )
    )
  }
}

check_required_packages()

suppressPackageStartupMessages({
  library(data.table)
  library(quantmod)
  library(xts)
  library(zoo)
  library(glmnet)
  library(pls)
  library(forecast)
  library(ggplot2)
  library(fda)
  library(rugarch)
})

# ============================= DIR HELPERS ====================================
mk_run_dirs <- function(OUT){
  for(d in c("data","tables","tex","plots","debug")){
    dir.create(file.path(OUT, d), showWarnings = FALSE, recursive = TRUE)
  }
}

dir.create(OUT_BASE, showWarnings = FALSE, recursive = TRUE)
mk_run_dirs(OUT_BASE)

write_csv <- function(df, OUT, name){
  fn <- file.path(OUT, "tables", paste0(name, ".csv"))
  dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(as.data.table(df), fn)
  fn
}

# ============================= TEX HELPERS ====================================
tex_escape <- function(x){
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("([%&#])", "\\\\\\1", x, perl = TRUE)
  x
}

pad_even_dt <- function(DT){
  DT <- as.data.table(copy(DT))
  if(nrow(DT) %% 2 == 1L){
    blank <- as.list(rep("", ncol(DT)))
    names(blank) <- names(DT)
    DT <- rbind(DT, as.data.table(blank), fill = TRUE)
  }
  DT
}

write_longtable_2block <- function(df, OUT, name, caption, label, left_cols, right_cols){
  stopifnot(length(left_cols) == 3, length(right_cols) == 3)
  
  df2 <- as.data.frame(df, stringsAsFactors = FALSE)
  if(nrow(df2) == 0){
    df2 <- data.frame(
      V1 = "", V2 = "", V3 = "",
      stringsAsFactors = FALSE
    )
    names(df2) <- left_cols
  }
  
  if(nrow(df2) %% 2 == 1){
    blank_row <- as.list(rep("", ncol(df2)))
    names(blank_row) <- names(df2)
    df2 <- rbind(df2, blank_row)
  }
  
  half <- nrow(df2) / 2
  LHS  <- df2[1:half, , drop = FALSE]
  RHS  <- df2[(half + 1):nrow(df2), , drop = FALSE]
  
  tex <- c(
    "\\begin{longtable}{L M P L M P}",
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{3pt}",
    "\\renewcommand{\\arraystretch}{1.03}",
    sprintf("\\caption{%s}\\label{%s}\\\\", caption, label),
    "\\toprule",
    sprintf("\\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s}\\\\",
            left_cols[1], left_cols[2], left_cols[3],
            right_cols[1], right_cols[2], right_cols[3]),
    "\\midrule",
    "\\endfirsthead",
    sprintf("\\caption[]{%s (continued).}\\\\", caption),
    "\\toprule",
    sprintf("\\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s} & \\textbf{%s}\\\\",
            left_cols[1], left_cols[2], left_cols[3],
            right_cols[1], right_cols[2], right_cols[3]),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{6}{r}{\\textit{Continued on next page}}\\\\",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot"
  )
  
  for(i in seq_len(half)){
    row <- c(
      tex_escape(LHS[i, left_cols[1]]),  tex_escape(LHS[i, left_cols[2]]),  tex_escape(LHS[i, left_cols[3]]),
      tex_escape(RHS[i, right_cols[1]]), tex_escape(RHS[i, right_cols[2]]), tex_escape(RHS[i, right_cols[3]])
    )
    tex <- c(tex, sprintf("%s & %s & %s & %s & %s & %s\\\\",
                          row[1], row[2], row[3], row[4], row[5], row[6]))
  }
  
  tex <- c(tex, "\\end{longtable}")
  
  fn <- file.path(OUT, "tex", paste0(name, ".tex"))
  dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
  writeLines(tex, fn)
  fn
}

write_longtable_simple <- function(df, OUT, name, caption, label, align = NULL){
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  
  if(nrow(df) == 0){
    df <- as.data.frame(matrix("", nrow = 1, ncol = ncol(df)))
    names(df) <- colnames(df)
  }
  
  if(is.null(align)){
    if(ncol(df) == 1L){
      align <- "l"
    } else {
      align <- paste(c("l", rep("c", ncol(df) - 1L)), collapse = "")
    }
  }
  
  hdr <- paste(sprintf("\\textbf{%s}", names(df)), collapse = " & ")
  
  tex <- c(
    sprintf("\\begin{longtable}{%s}", align),
    "\\footnotesize",
    "\\setlength{\\tabcolsep}{4pt}",
    "\\renewcommand{\\arraystretch}{1.03}",
    sprintf("\\caption{%s}\\label{%s}\\\\", caption, label),
    "\\toprule",
    paste0(hdr, "\\\\"),
    "\\midrule",
    "\\endfirsthead",
    sprintf("\\caption[]{%s (continued).}\\\\", caption),
    "\\toprule",
    paste0(hdr, "\\\\"),
    "\\midrule",
    "\\endhead",
    "\\midrule",
    sprintf("\\multicolumn{%d}{r}{\\textit{Continued on next page}}\\\\", ncol(df)),
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot"
  )
  
  for(i in seq_len(nrow(df))){
    row <- tex_escape(unlist(df[i, , drop = FALSE], use.names = FALSE))
    tex <- c(tex, paste0(paste(row, collapse = " & "), "\\\\"))
  }
  
  tex <- c(tex, "\\end{longtable}")
  
  fn <- file.path(OUT, "tex", paste0(name, ".tex"))
  dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
  writeLines(tex, fn)
  fn
}

# ============================= DISPLAY HELPERS ================================
display_model <- function(x){
  out <- x
  out[x == "Persistence(v_t)"] <- "Persistence($v_t$)"
  out[x == "AR(1)-logv"] <- "AR(1)-logv"
  out[x == "HAR"] <- "HAR"
  out[x == "HAR-X(logVIX)"] <- "HAR-X(logVIX)"
  out[x == "GARCH(1,1)"] <- "GARCH(1,1)"
  out[x == "EGARCH(1,1)"] <- "EGARCH(1,1)"
  out[x == "GJR-GARCH(1,1)"] <- "GJR-GARCH(1,1)"
  out[x == "fVAR"] <- "fVAR"
  out[x == "fVAR-X(+logv_t)"] <- "fVAR-X(+logv$_t$)"
  out[x == "fVAR-Lagged"] <- "fVAR-Lagged"
  out[x == "fVAR-Lagged-X"] <- "fVAR-Lagged-X"
  out
}

fmt4 <- function(x) sprintf("%.4f", x)
fmt6 <- function(x) sprintf("%.6f", x)
fmt_p_diag <- function(x){
  ifelse(!is.finite(x), "",
         ifelse(x < 1e-4, "$<0.0001$", formatC(x, format = "f", digits = 4)))
}
fmt_p_std <- function(x){
  ifelse(!is.finite(x), "", formatC(x, format = "f", digits = 4))
}
fmt_sci_or_fixed <- function(x){
  ifelse(!is.finite(x), "",
         ifelse(abs(x) < 1e-4, format(x, scientific = TRUE), formatC(x, format = "f", digits = 4)))
}

avail_mode_latex <- function(x){
  if(identical(x, "DAX_CLOSE")) "\\DAXCLOSE" else "\\GLOBALCLOSE"
}

# ============================= SYMBOL MAPS ====================================
SYMS_PRICE <- c(
  "SPY", "^FTSE", "^GDAXI", "^N225",
  "GC=F", "BZ=F", "SI=F", "HG=F",
  "DX-Y.NYB", "^VIX"
)
SYMS_RATE <- c("^TNX", "^IRX")
ALL_SYMS  <- unique(c(SYMS_PRICE, SYMS_RATE))

clean_name <- function(sym){
  nm <- sym
  nm <- gsub("\\^", "", nm)
  nm <- gsub("-Y\\.NYB", "", nm)
  nm <- gsub("=F", "", nm)
  nm <- gsub("\\.", "_", nm)
  nm
}

NAME_MAP <- setNames(sapply(ALL_SYMS, clean_name), ALL_SYMS)
NAME_MAP["^GDAXI"]   <- "DAX"
NAME_MAP["^N225"]    <- "Nikkei"
NAME_MAP["^VIX"]     <- "VIX"
NAME_MAP["DX-Y.NYB"] <- "DXY"
NAME_MAP["^TNX"]     <- "TNX"
NAME_MAP["^IRX"]     <- "IRX"

is_rate_sym <- function(sym) sym %in% SYMS_RATE

AFTER_DAX_CLOSE_RAW <- c("SPY","^VIX","DX-Y.NYB","^TNX","^IRX","GC=F","SI=F","HG=F","BZ=F")
is_after_dax_close_raw <- function(sym) sym %in% AFTER_DAX_CLOSE_RAW

IDX_GROUP  <- c("SPY","^FTSE","^GDAXI","^N225")
CORE_GROUP <- c("SPY","^FTSE","^GDAXI","^N225","GC=F","BZ=F","SI=F","HG=F","DX-Y.NYB","^VIX","^TNX","^IRX")

make_predictors_raw <- function(target_raw, cfg = c("FULL","PRUNE095","CORE")){
  cfg <- match.arg(cfg)
  if(cfg == "CORE")     return(unique(c(CORE_GROUP, target_raw)))
  if(cfg == "FULL")     return(unique(c(ALL_SYMS, target_raw)))
  if(cfg == "PRUNE095") return(unique(c(ALL_SYMS, target_raw)))
  unique(c(ALL_SYMS, target_raw))
}

# ============================= DOWNLOAD ONCE ==================================
get_ohlc <- function(sym, from, to){
  tryCatch(
    suppressWarnings(getSymbols(sym, src = "yahoo", from = from, to = to,
                                auto.assign = FALSE, warnings = FALSE)),
    error = function(e) NULL
  )
}

get_close <- function(ohlc){
  if(is.null(ohlc)) return(NULL)
  if(any(grepl("Adjusted", colnames(ohlc), ignore.case = TRUE))) Ad(ohlc) else Cl(ohlc)
}

message("=== Downloading OHLC/Close data once ===")
ohlc_list  <- list()
close_list <- list()

for(sym in ALL_SYMS){
  message("Downloading: ", sym)
  ohlc_list[[sym]]  <- get_ohlc(sym, start_date, end_date)
  close_list[[sym]] <- get_close(ohlc_list[[sym]])
}

ok_syms <- names(close_list)[sapply(close_list, function(x) !is.null(x) && NROW(x) > 50)]
if(!(TARGET_RAW %in% ok_syms)) stop("Target symbol download failed: ", TARGET_RAW)

# ============================= CORE HELPERS ===================================
lag_xts_k <- function(x, k = 1L){
  if(k <= 0L) return(x)
  core <- coredata(x)
  if(is.null(dim(core))) core <- matrix(core, ncol = 1L)
  if(nrow(core) <= k){
    out <- xts(matrix(NA_real_, nrow = nrow(core), ncol = ncol(core)), order.by = index(x))
    colnames(out) <- colnames(x)
    return(out)
  }
  lagged <- rbind(
    matrix(NA_real_, nrow = k, ncol = ncol(core)),
    core[seq_len(nrow(core) - k), , drop = FALSE]
  )
  out <- xts(lagged, order.by = index(x))
  colnames(out) <- colnames(x)
  out
}

rv_roger_satchell <- function(ohlc){
  if(is.null(ohlc)) return(NULL)
  O <- Op(ohlc); H <- Hi(ohlc); L <- Lo(ohlc); C <- Cl(ohlc)
  rs <- log(H/O) * log(H/C) + log(L/O) * log(L/C)
  rs <- pmax(as.numeric(rs), 1e-12)
  xts(rs, order.by = index(ohlc))
}

desc_stats <- function(x){
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if(length(x) < 5){
    return(c(SD = NA, MAD = NA, Min = NA, Max = NA, Skew = NA, Kurtosis = NA))
  }
  s  <- sd(x)
  m  <- mean(x)
  sk <- mean((x - m)^3) / (s^3 + 1e-12)
  ku <- mean((x - m)^4) / (s^4 + 1e-12)
  c(SD = s, MAD = mad(x), Min = min(x), Max = max(x), Skew = sk, Kurtosis = ku)
}

qlike_vec <- function(v, vhat){
  eps  <- 1e-12
  v    <- pmax(as.numeric(v), eps)
  vhat <- pmax(as.numeric(vhat), eps)
  log(vhat) + (v / vhat)
}

metric_var <- function(v, vhat){
  v    <- as.numeric(v)
  vhat <- as.numeric(vhat)
  data.frame(
    MAE_var  = mean(abs(v - vhat), na.rm = TRUE),
    RMSE_var = sqrt(mean((v - vhat)^2, na.rm = TRUE)),
    QLIKE    = mean(qlike_vec(v, vhat), na.rm = TRUE)
  )
}

arch_lm <- function(z, m = 10L){
  z2 <- as.numeric(z)^2
  z2 <- z2[is.finite(z2)]
  if(length(z2) <= m + 10L) return(list(stat = NA, p = NA, n = length(z2)))
  
  y <- z2[(m + 1):length(z2)]
  X <- sapply(1:m, function(k) z2[(m + 1 - k):(length(z2) - k)])
  fit <- lm(y ~ X)
  R2  <- summary(fit)$r.squared
  LM  <- length(y) * R2
  p   <- 1 - pchisq(LM, df = m)
  list(stat = as.numeric(LM), p = as.numeric(p), n = length(y))
}

nw_var <- function(d, L){
  d <- d[is.finite(d)]
  n <- length(d)
  if(n < (L + 5L)) return(var(d, na.rm = TRUE))
  d <- d - mean(d)
  gamma0 <- sum(d * d) / n
  v <- gamma0
  for(l in 1:L){
    gam <- sum(d[(l + 1):n] * d[1:(n - l)]) / n
    w   <- 1 - l / (L + 1)
    v   <- v + 2 * w * gam
  }
  v
}

dm_loss <- function(loss1, loss2, alternative = c("less","greater","two.sided"), L = 5L){
  alternative <- match.arg(alternative)
  ok <- is.finite(loss1) & is.finite(loss2)
  d  <- (loss1 - loss2)[ok]
  n  <- length(d)
  if(n < 40L) return(list(stat = NA, p = NA, n = n))
  
  v  <- nw_var(d, L = L)
  se <- sqrt(v / n)
  if(!is.finite(se) || se <= 0) return(list(stat = NA, p = NA, n = n))
  
  stat <- mean(d) / se
  df   <- n - 1L
  
  if(alternative == "less")      p <- pt(stat, df = df)
  if(alternative == "greater")   p <- 1 - pt(stat, df = df)
  if(alternative == "two.sided") p <- 2 * min(pt(stat, df = df), 1 - pt(stat, df = df))
  
  list(stat = as.numeric(stat), p = as.numeric(p), n = n)
}

mbb_mean_diff <- function(diff_series, B = 10L, R = 2000L){
  d <- as.numeric(diff_series)
  d <- d[is.finite(d)]
  n <- length(d)
  if(n < 60L) return(list(mean = NA, p_one = NA, ci = c(NA, NA), B = B, n = n))
  
  nb     <- ceiling(n / B)
  starts <- 1:(n - B + 1L)
  blocks <- lapply(starts, function(s) d[s:(s + B - 1L)])
  
  means <- numeric(R)
  for(r in seq_len(R)){
    idx <- sample.int(length(blocks), size = nb, replace = TRUE)
    x   <- unlist(blocks[idx], use.names = FALSE)
    x   <- x[1:n]
    means[r] <- mean(x)
  }
  
  m0    <- mean(d)
  p_one <- mean(means >= 0)
  ci    <- as.numeric(quantile(means, probs = c(0.025, 0.975), na.rm = TRUE))
  
  list(mean = m0, p_one = p_one, ci = ci, B = B, n = n)
}

standardize_by_idx <- function(X, idx_fit){
  X   <- as.matrix(X)
  mu  <- colMeans(X[idx_fit, , drop = FALSE], na.rm = TRUE)
  sdv <- apply(X[idx_fit, , drop = FALSE], 2, sd, na.rm = TRUE)
  sdv[!is.finite(sdv) | sdv < 1e-12] <- 1
  Xs <- sweep(X, 2, mu, "-")
  Xs <- sweep(Xs, 2, sdv, "/")
  list(Xs = Xs, center = mu, scale = sdv)
}

make_blocked_foldid <- function(n, K = 5L){
  brks <- floor(seq(0, n, length.out = K + 1L))
  foldid <- integer(n)
  for(k in 1:K){
    i1 <- brks[k] + 1L
    i2 <- brks[k + 1L]
    if(i1 <= i2) foldid[i1:i2] <- k
  }
  foldid
}

# ============================= MODEL MAP / FAMILY HELPERS =====================
MODEL_MAP <- c(
  "Persistence(v_t)" = "Persistence",
  "AR(1)-logv"       = "AR1_logv",
  "HAR"              = "HAR",
  "HAR-X(logVIX)"    = "HAR_X",
  "GARCH(1,1)"       = "GARCH11",
  "EGARCH(1,1)"      = "EGARCH11",
  "GJR-GARCH(1,1)"   = "GJR11",
  "fVAR"             = "fVAR",
  "fVAR-X(+logv_t)"  = "fVAR_X",
  "fVAR-Lagged"      = "fVAR_L",
  "fVAR-Lagged-X"    = "fVAR_LX"
)

FUNCTIONAL_LABELS <- c("fVAR", "fVAR-X(+logv_t)", "fVAR-Lagged", "fVAR-Lagged-X")

collect_metric_table_from_M <- function(M){
  V_true <- as.numeric(M[, "Actual"])
  
  out <- rbindlist(lapply(names(MODEL_MAP), function(lbl){
    col_nm <- MODEL_MAP[[lbl]]
    data.table(
      Model = lbl,
      metric_var(V_true, as.numeric(M[, col_nm]))
    )
  }), fill = TRUE)
  
  q_garch <- out[Model == "GARCH(1,1)", QLIKE]
  out[, Delta_QLIKE_vs_GARCH := QLIKE - q_garch]
  out[, Rank_QLIKE := frank(QLIKE, ties.method = "first")]
  out[]
}

choose_best_functional <- function(tab_metrics){
  cand <- copy(tab_metrics[Model %in% FUNCTIONAL_LABELS])
  cand[which.min(QLIKE)]
}

make_protocol_modelset_compare <- function(tab_main, tab_rob,
                                           main_name = "DAX_CLOSE",
                                           rob_name  = "GLOBAL_CLOSE"){
  x <- merge(
    tab_main[, .(Model,
                 QLIKE_MAIN = QLIKE,
                 Rank_MAIN  = Rank_QLIKE)],
    tab_rob[,  .(Model,
                 QLIKE_ROB  = QLIKE,
                 Rank_ROB   = Rank_QLIKE)],
    by = "Model", all = TRUE
  )
  
  setnames(x,
           old = c("QLIKE_MAIN","Rank_MAIN","QLIKE_ROB","Rank_ROB"),
           new = c(paste0("QLIKE_", main_name),
                   paste0("Rank_",  main_name),
                   paste0("QLIKE_", rob_name),
                   paste0("Rank_",  rob_name)))
  
  x[, Delta_MAIN_minus_ROB := get(paste0("QLIKE_", main_name)) - get(paste0("QLIKE_", rob_name))]
  setorderv(x, cols = c(paste0("Rank_", main_name), paste0("Rank_", rob_name)))
  x[]
}

make_protocol_vs_garch_compare <- function(tab_main, tab_rob,
                                           main_name = "DAX_CLOSE",
                                           rob_name  = "GLOBAL_CLOSE"){
  x <- merge(
    tab_main[, .(Model, Delta_vs_GARCH_MAIN = Delta_QLIKE_vs_GARCH)],
    tab_rob[,  .(Model, Delta_vs_GARCH_ROB  = Delta_QLIKE_vs_GARCH)],
    by = "Model", all = TRUE
  )
  
  setnames(x,
           old = c("Delta_vs_GARCH_MAIN","Delta_vs_GARCH_ROB"),
           new = c(paste0("Delta_vs_GARCH_", main_name),
                   paste0("Delta_vs_GARCH_", rob_name)))
  
  x[, DeltaGap_MAIN_minus_ROB :=
      get(paste0("Delta_vs_GARCH_", main_name)) -
      get(paste0("Delta_vs_GARCH_", rob_name))]
  setorderv(x, cols = paste0("Delta_vs_GARCH_", main_name))
  x[]
}

make_crossprotocol_common_compare <- function(M_main, M_rob,
                                              main_name = "DAX_CLOSE",
                                              rob_name  = "GLOBAL_CLOSE"){
  common_dates <- intersect(index(M_main), index(M_rob))
  if(length(common_dates) < 50L) return(NULL)
  
  M_main_cc <- M_main[common_dates]
  M_rob_cc  <- M_rob[common_dates]
  
  tab_main_cc <- collect_metric_table_from_M(M_main_cc)
  tab_rob_cc  <- collect_metric_table_from_M(M_rob_cc)
  
  full_cc <- make_protocol_modelset_compare(tab_main_cc, tab_rob_cc,
                                            main_name = main_name, rob_name = rob_name)
  rel_cc  <- make_protocol_vs_garch_compare(tab_main_cc, tab_rob_cc,
                                            main_name = main_name, rob_name = rob_name)
  
  list(
    dates       = common_dates,
    tab_main_cc = tab_main_cc,
    tab_rob_cc  = tab_rob_cc,
    full_cc     = full_cc,
    rel_cc      = rel_cc
  )
}

# ============================= PANEL BUILD ====================================
prune_by_corr <- function(ret_panel, sym_names, thr = 0.95){
  keep <- sym_names
  if(length(keep) <= 2L) return(keep)
  
  repeat {
    C <- suppressWarnings(cor(ret_panel[, keep, drop = FALSE], use = "pairwise.complete.obs"))
    C[lower.tri(C, diag = TRUE)] <- NA_real_
    mx <- max(abs(C), na.rm = TRUE)
    if(!is.finite(mx) || mx <= thr) break
    
    ij <- which(abs(C) == mx, arr.ind = TRUE)[1, ]
    a  <- keep[ij[1]]
    b  <- keep[ij[2]]
    ma <- mean(abs(C[a, ]), na.rm = TRUE)
    mb <- mean(abs(C[b, ]), na.rm = TRUE)
    drop <- ifelse(ma >= mb, a, b)
    
    keep <- setdiff(keep, drop)
    if(length(keep) <= 2L) break
  }
  
  keep
}

build_panel <- function(selected_raw_syms, target_raw, cfg_name,
                        avail_mode = c("GLOBAL_CLOSE","DAX_CLOSE")){
  avail_mode <- match.arg(avail_mode)
  
  closes <- lapply(selected_raw_syms, function(s) close_list[[s]])
  names(closes) <- selected_raw_syms
  
  close_panel <- do.call(merge, c(closes, all = FALSE))
  colnames(close_panel) <- sapply(selected_raw_syms, function(s) NAME_MAP[[s]])
  close_panel <- na.omit(close_panel)
  
  ret_list <- list()
  for(s in selected_raw_syms){
    nm   <- NAME_MAP[[s]]
    x    <- close_panel[, nm]
    xnum <- as.numeric(x)
    if(is_rate_sym(s)){
      r <- diff(xnum)
      r <- xts(r, order.by = index(x)[-1])
    } else {
      xnum <- pmax(xnum, 1e-12)
      r <- diff(log(xnum))
      r <- xts(r, order.by = index(x)[-1])
    }
    ret_list[[nm]] <- r
  }
  
  ret_panel <- do.call(merge, c(ret_list, all = FALSE))
  ret_panel <- na.omit(ret_panel)
  
  close_panel <- close_panel[index(ret_panel), , drop = FALSE]
  
  if(avail_mode == "DAX_CLOSE"){
    lag_raw   <- intersect(selected_raw_syms, selected_raw_syms[sapply(selected_raw_syms, is_after_dax_close_raw)])
    lag_names <- sapply(lag_raw, function(s) NAME_MAP[[s]])
    lag_names <- intersect(lag_names, colnames(ret_panel))
    
    targ_nm <- NAME_MAP[[target_raw]]
    lag_names <- setdiff(lag_names, targ_nm)
    
    if(length(lag_names) > 0){
      ret_panel[, lag_names]   <- lag_xts_k(ret_panel[, lag_names, drop = FALSE], k = 1L)
      close_panel[, lag_names] <- lag_xts_k(close_panel[, lag_names, drop = FALSE], k = 1L)
    }
  }
  
  ok_ret   <- complete.cases(ret_panel)
  ok_close <- complete.cases(close_panel)
  idx_keep <- index(ret_panel)[ok_ret & ok_close]
  ret_panel   <- ret_panel[idx_keep, , drop = FALSE]
  close_panel <- close_panel[idx_keep, , drop = FALSE]
  
  if(cfg_name == "PRUNE095"){
    idx_raw   <- intersect(IDX_GROUP, selected_raw_syms)
    idx_names <- sapply(idx_raw, function(s) NAME_MAP[[s]])
    idx_names <- intersect(idx_names, colnames(ret_panel))
    non_idx   <- setdiff(colnames(ret_panel), idx_names)
    
    if(length(idx_names) >= 3L){
      keep_idx  <- prune_by_corr(ret_panel, idx_names, thr = 0.95)
      keep_cols <- c(non_idx, keep_idx)
      ret_panel   <- ret_panel[, keep_cols, drop = FALSE]
      close_panel <- close_panel[, keep_cols, drop = FALSE]
    }
  }
  
  target_nm <- NAME_MAP[[target_raw]]
  if(!(target_nm %in% colnames(ret_panel))){
    stop("Target missing after panel build: ", target_nm)
  }
  
  list(close_panel = close_panel, ret_panel = ret_panel, target_nm = target_nm)
}

# ============================= MODEL HELPERS ==================================
rollmat <- function(x, W){
  x <- as.numeric(x)
  n <- length(x)
  if(n < W) stop("Series shorter than W.")
  N <- n - W + 1L
  M <- matrix(NA_real_, nrow = N, ncol = W)
  for(i in seq_len(N)) M[i, ] <- x[i:(i + W - 1L)]
  M
}

win_features <- function(M){
  W <- ncol(M)
  t <- 1:W
  mu   <- rowMeans(M)
  sdv  <- apply(M, 1, sd)
  last <- M[, W]
  slope <- apply(M, 1, function(y) coef(lm(y ~ t))[2])
  cbind(mu = mu, sd = sdv, last = last, slope = slope)
}

oos_forecast_expanding <- function(y, X = NULL, idx_eval, refit_every, fit_fun, pred_fun){
  yhat <- rep(NA_real_, length(idx_eval))
  fit  <- NULL
  
  for(j in seq_along(idx_eval)){
    i <- idx_eval[j]
    if(j == 1L || ((j - 1L) %% refit_every) == 0L){
      tr  <- 1:(i - 1L)
      fit <- tryCatch(fit_fun(y[tr], if(!is.null(X)) X[tr, , drop = FALSE] else NULL),
                      error = function(e) NULL)
    }
    yhat[j] <- tryCatch(pred_fun(fit, if(!is.null(X)) X[i, , drop = FALSE] else NULL),
                        error = function(e) NA_real_)
  }
  
  yhat
}

fit_ar1 <- function(ytr, Xtr = NULL){
  forecast::Arima(ytr, order = c(1, 0, 0), include.mean = TRUE)
}

pred_ar1_1step <- function(fit, Xnew = NULL){
  if(is.null(fit)) return(NA_real_)
  as.numeric(forecast::forecast(fit, h = 1)$mean[1])
}

pls_fit_best_blocked <- function(X, y, max_comp = 8L, nseg = 5L){
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  max_comp <- min(max_comp, p, n - 5L)
  if(max_comp < 1L) stop("Not enough data for PLS.")
  
  brks   <- floor(seq(0, n, length.out = nseg + 1L))
  cv_err <- rep(0, max_comp)
  cv_n   <- rep(0, max_comp)
  
  for(s in 2:nseg){
    tr_end <- brks[s - 1L]
    va_st  <- brks[s - 1L] + 1L
    va_en  <- brks[s]
    if(tr_end < 30L || va_st > va_en) next
    
    Xtr <- X[1:tr_end, , drop = FALSE]
    ytr <- y[1:tr_end]
    Xva <- X[va_st:va_en, , drop = FALSE]
    yva <- y[va_st:va_en]
    
    df_tr <- data.frame(y = ytr, Xtr, check.names = FALSE)
    fit   <- try(plsr(y ~ ., data = df_tr, ncomp = max_comp,
                      method = "oscorespls", scale = FALSE),
                 silent = TRUE)
    if(inherits(fit, "try-error")) next
    
    df_va <- data.frame(Xva, check.names = FALSE)
    for(k in 1:max_comp){
      ph <- try(as.numeric(predict(fit, newdata = df_va, ncomp = k)), silent = TRUE)
      if(inherits(ph, "try-error")) next
      mse <- mean((yva - ph)^2, na.rm = TRUE)
      if(is.finite(mse)){
        cv_err[k] <- cv_err[k] + mse
        cv_n[k]   <- cv_n[k] + 1L
      }
    }
  }
  
  cv_mse <- cv_err / pmax(cv_n, 1L)
  best_k <- which.min(cv_mse)
  
  df_full <- data.frame(y = y, X, check.names = FALSE)
  fit_full <- plsr(y ~ ., data = df_full, ncomp = best_k,
                   method = "oscorespls", scale = FALSE)
  
  list(fit = fit_full, ncomp = best_k, xnames = colnames(X), cv_mse = cv_mse)
}

fit_pls <- function(ytr, Xtr){
  pls_fit_best_blocked(Xtr, ytr, max_comp = PLS_MAX_COMP, nseg = PLS_NSEG)
}

pred_pls_1step <- function(obj, Xnew){
  if(is.null(obj)) return(NA_real_)
  Xdf <- as.data.frame(Xnew, check.names = FALSE)
  colnames(Xdf) <- obj$xnames
  as.numeric(predict(obj$fit, newdata = Xdf, ncomp = obj$ncomp))
}

fit_lm <- function(ytr, Xtr){
  df <- data.frame(y = ytr, Xtr, check.names = FALSE)
  lm(y ~ ., data = df)
}

pred_lm_1step <- function(fit, Xnew){
  if(is.null(fit)) return(NA_real_)
  df <- as.data.frame(Xnew, check.names = FALSE)
  as.numeric(predict(fit, newdata = df))
}

har_design <- function(logv, t_end_idx, Xextra = NULL){
  t0 <- t_end_idx
  lag1   <- logv[t0]
  mean5  <- sapply(t0, function(tt) mean(logv[(tt - 4L):tt],  na.rm = TRUE))
  mean22 <- sapply(t0, function(tt) mean(logv[(tt - 21L):tt], na.rm = TRUE))
  
  X <- cbind(lag1 = lag1, mean5 = mean5, mean22 = mean22)
  if(!is.null(Xextra)) X <- cbind(X, Xextra)
  X
}

garch_forecast_oos_map <- function(ret_target, t_end_idx, idx_eval, refit_every, spec){
  r    <- as.numeric(ret_target)
  yhat <- rep(NA_real_, length(idx_eval))
  fit  <- NULL
  
  for(j in seq_along(idx_eval)){
    i        <- idx_eval[j]
    t_origin <- t_end_idx[i]
    
    if(j == 1L || ((j - 1L) %% refit_every) == 0L){
      tr   <- 1:t_origin
      r_tr <- r[tr]
      fit  <- try(rugarch::ugarchfit(spec, r_tr, solver = "hybrid",
                                     fit.control = list(scale = 1)),
                  silent = TRUE)
      if(inherits(fit, "try-error")) fit <- NULL
    }
    
    if(is.null(fit)){
      yhat[j] <- NA_real_
      next
    }
    
    fc <- try(rugarch::ugarchforecast(fit, n.ahead = 1), silent = TRUE)
    if(inherits(fc, "try-error")){
      yhat[j] <- NA_real_
      next
    }
    
    sig <- as.numeric(rugarch::sigma(fc))
    yhat[j] <- sig^2
  }
  
  yhat
}

z_from_vhat <- function(r, vhat){
  as.numeric(r) / sqrt(pmax(as.numeric(vhat), 1e-12))
}

# ============================= ONE RUN ========================================
run_one <- function(cfg_name, refit_every, avail_mode = c("GLOBAL_CLOSE","DAX_CLOSE")){
  avail_mode <- match.arg(avail_mode)
  
  run_tag <- paste0("run_", TARGET_NAME, "__", cfg_name, "__ref", refit_every, "__", avail_mode)
  OUT <- file.path(OUT_BASE, run_tag)
  mk_run_dirs(OUT)
  
  message("\n=== RUN: target=DAX cfg=", cfg_name,
          " refit=", refit_every,
          " avail=", avail_mode, " ===")
  
  sel_raw <- make_predictors_raw(TARGET_RAW, cfg = cfg_name)
  sel_raw <- intersect(sel_raw, ok_syms)
  sel_raw <- unique(c(sel_raw, TARGET_RAW))
  
  pan <- build_panel(sel_raw, TARGET_RAW, cfg_name, avail_mode = avail_mode)
  close_panel <- pan$close_panel
  ret_panel   <- pan$ret_panel
  target_nm   <- pan$target_nm
  
  write.csv(data.frame(Date = index(ret_panel), coredata(ret_panel)),
            file.path(OUT, "data", "panel_returns_common_calendar.csv"),
            row.names = FALSE)
  
  rv <- rv_roger_satchell(ohlc_list[[TARGET_RAW]])
  rv <- rv[index(ret_panel)]
  rv <- na.omit(rv)
  
  ret_panel   <- ret_panel[index(rv), , drop = FALSE]
  close_panel <- close_panel[index(rv), , drop = FALSE]
  logv        <- log(as.numeric(rv) + 1e-12)
  
  write.csv(data.frame(Date = index(rv), RV = as.numeric(rv), logRV = logv),
            file.path(OUT, "data", "target_RS_variance.csv"),
            row.names = FALSE)
  
  # ----- Descriptive stats ----------------------------------------------------
  DESC <- t(apply(ret_panel, 2, desc_stats))
  DESC <- data.frame(Series = rownames(DESC), round(DESC, 6), row.names = NULL)
  write_csv(DESC, OUT, "A1_desc_returns_raw")
  
  desc_tex <- DESC[, c("Series","SD","MAD","Min","Max")]
  write_longtable_2block(
    desc_tex,
    OUT, "A1_desc_returns",
    caption = sprintf("Descriptive statistics of transformed daily predictors under the strict intersection calendar (no imputation), reported for the locked baseline design under %s.", avail_mode_latex(avail_mode)),
    label   = sprintf("tab:desc_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    left_cols  = c("Series","SD","MAD"),
    right_cols = c("Series","Min","Max")
  )
  
  # ----- Functional dataset ---------------------------------------------------
  retX <- ret_panel
  Tn   <- nrow(retX)
  
  W_START <- max(W_CURVE, 22L)
  t_end   <- W_START:(Tn - 1L)
  N       <- length(t_end)
  if(N < 250L) stop("Too few observations after intersection calendar.")
  
  PRED_VARS <- colnames(retX)
  
  Xwin_list <- lapply(PRED_VARS, function(v) rollmat(retX[, v], W_CURVE))
  names(Xwin_list) <- PRED_VARS
  
  k_idx <- t_end - (W_CURVE - 1L)
  Xwin_aligned <- lapply(Xwin_list, function(M) M[k_idx, , drop = FALSE])
  
  y_logv <- logv[t_end + 1L]
  y_rv   <- as.numeric(rv)[t_end + 1L]
  
  # ----- Splits ---------------------------------------------------------------
  n_all     <- N
  n_test    <- floor(TEST_RATIO * n_all)
  n_dev     <- n_all - n_test
  n_val     <- floor(VAL_RATIO_DEV * n_dev)
  train_end <- n_dev - n_val
  
  idx_val    <- (train_end + 1L):n_dev
  idx_test   <- (n_dev + 1L):n_all
  idx_dev    <- 1:n_dev
  
  curve_grid <- 1:W_CURVE
  eval_grid  <- seq(1, W_CURVE, length.out = GRID_EVAL)
  
  smooth_linear <- function(wvec){
    approx(x = curve_grid, y = wvec, xout = eval_grid, rule = 2)$y
  }
  
  # ----- Representative series for DEV-only GCV -------------------------------
  rep_var <- if(REP_SERIES_GCV %in% names(Xwin_aligned)) REP_SERIES_GCV else names(Xwin_aligned)[1]
  
  DESIGN_META <- data.table(
    Target                       = TARGET_NAME,
    CandidatePredictorSet        = cfg_name,
    PredictorSetConstruction     = cfg_pretty(cfg_name),
    AvailabilityProtocol         = avail_mode,
    RepresentativeSeries_GCV     = rep_var,
    GCV_Nbasis_Grid              = paste(NBASIS_GRID, collapse = ","),
    GCV_Lambda_Grid              = paste(format(LAMBDA_GRID, scientific = TRUE), collapse = ","),
    PLS_Component_Grid           = paste(seq_len(PLS_MAX_COMP), collapse = ","),
    PLS_NSegments                = PLS_NSEG,
    PLS_BlockedValidation        = "Forward contiguous blocked validation: segment s validates on block s; training uses blocks 1:(s-1)",
    DEV_Block_Split              = sprintf("After holding out TEST, DEV has %d observations; internal validation uses the last %d observations of DEV (%.0f%%).",
                                           n_dev, n_val, 100 * VAL_RATIO_DEV)
  )
  write_csv(DESIGN_META, OUT, "A0_design_meta")
  
  # ----- DEV-only GCV for smoothing ------------------------------------------
  OPT_NB <- NBASIS_GRID[1]
  OPT_LA <- LAMBDA_GRID[3]
  
  if(USE_BSPLINE_SMOOTH){
    Mrep    <- Xwin_aligned[[rep_var]]
    
    dev_rows <- idx_dev
    n_sub    <- min(GCV_SUBSAMPLE_WINDOWS, length(dev_rows))
    sub_idx  <- sort(sample(dev_rows, size = n_sub, replace = FALSE))
    
    gcv_grid <- CJ(nbasis = NBASIS_GRID, lambda = LAMBDA_GRID)
    gcv_grid[, gcv := NA_real_]
    
    for(ii in 1:nrow(gcv_grid)){
      nb  <- gcv_grid$nbasis[ii]
      lam <- gcv_grid$lambda[ii]
      
      basis_tmp <- create.bspline.basis(rangeval = c(1, W_CURVE), nbasis = nb)
      fdParobj  <- fdPar(basis_tmp, Lfdobj = int2Lfd(2), lambda = lam)
      
      gcv_vals <- rep(NA_real_, n_sub)
      for(k in seq_along(sub_idx)){
        yk <- matrix(Mrep[sub_idx[k], ], ncol = 1)
        sm <- try(smooth.basis(argvals = curve_grid, y = yk, fdParobj = fdParobj), silent = TRUE)
        gcv_vals[k] <- if(inherits(sm, "try-error")) NA_real_ else as.numeric(sm$gcv)
      }
      gcv_grid$gcv[ii] <- mean(gcv_vals, na.rm = TRUE)
    }
    
    best_row <- gcv_grid[which.min(gcv), ]
    OPT_NB <- best_row$nbasis
    OPT_LA <- best_row$lambda
    
    write_csv(gcv_grid, OUT, "A2_gcv_grid_raw")
    
    gcv_tex <- copy(gcv_grid)
    gcv_tex[, `Basis size` := nbasis]
    gcv_tex[, `\\lambda` := fifelse(
      lambda == 1e-6, "$10^{-6}$",
      fifelse(lambda == 1e-5, "$10^{-5}$",
              fifelse(lambda == 1e-4, "$10^{-4}$",
                      fifelse(lambda == 1e-3, "$10^{-3}$", "$10^{-2}$")))
    )]
    gcv_tex[, GCV := formatC(gcv, format = "f", digits = 6)]
    gcv_tex <- gcv_tex[, .(`Basis size`, `\\lambda`, GCV)]
    
    write_longtable_2block(
      gcv_tex,
      OUT, "A2_gcv_grid",
      caption = sprintf("Representative-window GCV grid for B-spline smoothing on DEV ($W=%d$), reported for the locked baseline design under %s. The representative predictor series is %s.", W_CURVE, avail_mode_latex(avail_mode), rep_var),
      label = sprintf("tab:gcv_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
      left_cols  = c("Basis size","\\lambda","GCV"),
      right_cols = c("Basis size","\\lambda","GCV")
    )
    
    gcv_plot_dt <- copy(gcv_grid)
    gcv_plot_dt[, lambda_s := format(lambda, scientific = TRUE)]
    
    p_gcv <- ggplot(gcv_plot_dt, aes(x = lambda_s, y = factor(nbasis), fill = gcv)) +
      geom_tile() +
      labs(x = "lambda", y = "nbasis",
           title = sprintf("GCV heatmap (DEV-only; rep=%s; W=%d)", rep_var, W_CURVE)) +
      theme_minimal(base_size = 14)
    
    ggsave(file.path(OUT, "plots", "F4_gcv_heatmap.png"),
           p_gcv, width = 11, height = 5.5, dpi = PLOT_DPI)
    
    smooth_bspline <- function(wvec){
      basis_use <- create.bspline.basis(rangeval = c(1, W_CURVE), nbasis = OPT_NB)
      fdParobj  <- fdPar(basis_use, Lfdobj = int2Lfd(2), lambda = OPT_LA)
      yk <- matrix(wvec, ncol = 1)
      sm <- try(smooth.basis(argvals = curve_grid, y = yk, fdParobj = fdParobj), silent = TRUE)
      if(inherits(sm, "try-error")) return(smooth_linear(wvec))
      as.numeric(eval.fd(eval_grid, sm$fd))
    }
    
    WINDOW_TO_GRID <- smooth_bspline
  } else {
    WINDOW_TO_GRID <- smooth_linear
  }
  
  # ----- Construct functional/grid features ----------------------------------
  Xflat <- NULL
  for(v in names(Xwin_aligned)){
    M  <- Xwin_aligned[[v]]
    M2 <- t(apply(M, 1, WINDOW_TO_GRID))
    colnames(M2) <- paste0(v, "_g", seq_len(GRID_EVAL))
    Xflat <- cbind(Xflat, M2)
  }
  
  Xsum <- NULL
  feat_names <- c("mu","sd","last","slope")
  for(v in names(Xwin_aligned)){
    F <- win_features(Xwin_aligned[[v]])
    colnames(F) <- paste0(v, "_", feat_names)
    Xsum <- cbind(Xsum, F)
  }
  
  # ----- DEV-only scaling -----------------------------------------------------
  Xflat_std <- standardize_by_idx(Xflat, idx_fit = idx_dev)$Xs
  Xsum_std  <- standardize_by_idx(Xsum,  idx_fit = idx_dev)$Xs
  
  # ----- DEV-only blocked Elastic Net ----------------------------------------
  foldid_dev <- make_blocked_foldid(length(idx_dev), K = 5L)
  en_fit <- cv.glmnet(
    Xsum_std[idx_dev, , drop = FALSE],
    y_logv[idx_dev],
    alpha   = EN_ALPHA,
    foldid  = foldid_dev,
    grouped = FALSE
  )
  
  co <- as.matrix(coef(en_fit, s = "lambda.min"))
  sel_feat <- rownames(co)[abs(co[, 1]) > 1e-10]
  sel_feat <- setdiff(sel_feat, "(Intercept)")
  sel_vars <- unique(sub("_(mu|sd|last|slope)$", "", sel_feat))
  if(length(sel_vars) == 0L) sel_vars <- names(Xwin_aligned)
  
  png(file.path(OUT, "plots", "F3_elasticnet_cv.png"),
      width = 11, height = 6, units = "in", res = PLOT_DPI)
  par(mar = c(4,4,2,1) + 0.2)
  plot(en_fit, main = "Elastic Net CV (DEV-only; blocked folds)")
  dev.off()
  
  SELTAB <- data.table(`Retained variable` = sel_vars)
  SELTAB[, blank1 := ""]
  SELTAB[, blank2 := ""]
  SELTAB <- SELTAB[, .(`Retained variable`, blank1, blank2)]
  SELTAB <- pad_even_dt(SELTAB)
  write_csv(SELTAB, OUT, "A3_selected_vars_raw")
  
  write_longtable_2block(
    SELTAB,
    OUT, "A3_selected_vars",
    caption = "Elastic Net screening on DEV with blocked cross-validation. Reported entries are the variables retained for fVAR feature construction. In the locked baseline design, all twelve predictors are retained.",
    label = sprintf("tab:sel_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    left_cols  = c("Retained variable","blank1","blank2"),
    right_cols = c("Retained variable","blank1","blank2")
  )
  
  keep_cols  <- unlist(lapply(sel_vars, function(vv) grep(paste0("^", vv, "_g"), colnames(Xflat_std), value = TRUE)))
  Xflat_sel  <- Xflat_std[, keep_cols, drop = FALSE]
  
  # ----- Lagged functional designs -------------------------------------------
  Xflat_lag <- cbind(Xflat_sel[-1, , drop = FALSE], Xflat_sel[-N, , drop = FALSE])
  colnames(Xflat_lag) <- c(
    paste0(colnames(Xflat_sel), "__t"),
    paste0(colnames(Xflat_sel), "__tminus1")
  )
  
  y_logv_lag <- y_logv[-1]
  y_rv_lag   <- y_rv[-1]
  
  logv_t     <- logv[t_end]
  logv_t_std <- standardize_by_idx(matrix(logv_t, ncol = 1), idx_fit = idx_dev)$Xs[, 1]
  
  Xflat_X     <- cbind(logv_t = as.numeric(logv_t_std), Xflat_sel)
  Xflat_lag_X <- cbind(logv_t = as.numeric(logv_t_std[-1]), Xflat_lag)
  
  n_all_L     <- length(y_logv_lag)
  n_test_L    <- floor(TEST_RATIO * n_all_L)
  n_dev_L     <- n_all_L - n_test_L
  n_val_L     <- floor(VAL_RATIO_DEV * n_dev_L)
  train_end_L <- n_dev_L - n_val_L
  idx_val_L   <- (train_end_L + 1L):n_dev_L
  idx_test_L  <- (n_dev_L + 1L):n_all_L
  
  # ----- HAR / HAR-X(logVIX) --------------------------------------------------
  Xhar <- har_design(logv, t_end_idx = t_end, Xextra = NULL)
  
  x_vix <- NULL
  if("VIX" %in% colnames(close_panel)){
    x_vix <- log(pmax(as.numeric(close_panel[t_end, "VIX"]), 1e-12))
    x_vix <- standardize_by_idx(matrix(x_vix, ncol = 1), idx_fit = idx_dev)$Xs[, 1]
  }
  
  XharX <- if(!is.null(x_vix)) {
    har_design(logv, t_end_idx = t_end, Xextra = cbind(logVIX = x_vix))
  } else {
    Xhar
  }
  
  # ----- Validation forecasts -------------------------------------------------
  pred_ar1_val    <- oos_forecast_expanding(y_logv, NULL,            idx_val,   refit_every, fit_ar1, pred_ar1_1step)
  pred_har_val    <- oos_forecast_expanding(y_logv, Xhar,            idx_val,   refit_every, fit_lm,  pred_lm_1step)
  pred_harX_val   <- oos_forecast_expanding(y_logv, XharX,           idx_val,   refit_every, fit_lm,  pred_lm_1step)
  
  pred_fvar_val   <- oos_forecast_expanding(y_logv, Xflat_sel,       idx_val,   refit_every, fit_pls, pred_pls_1step)
  pred_fvarX_val  <- oos_forecast_expanding(y_logv, Xflat_X,         idx_val,   refit_every, fit_pls, pred_pls_1step)
  pred_fvarL_val  <- oos_forecast_expanding(y_logv_lag, Xflat_lag,   idx_val_L, refit_every, fit_pls, pred_pls_1step)
  pred_fvarLX_val <- oos_forecast_expanding(y_logv_lag, Xflat_lag_X, idx_val_L, refit_every, fit_pls, pred_pls_1step)
  
  vhat_ar1_val    <- exp(pred_ar1_val)
  vhat_har_val    <- exp(pred_har_val)
  vhat_harX_val   <- exp(pred_harX_val)
  vhat_fvar_val   <- exp(pred_fvar_val)
  vhat_fvarX_val  <- exp(pred_fvarX_val)
  vhat_fvarL_val  <- exp(pred_fvarL_val)
  vhat_fvarLX_val <- exp(pred_fvarLX_val)
  
  v_t           <- as.numeric(rv)[t_end]
  vhat_pers_val <- v_t[idx_val]
  
  # ----- GARCH family ---------------------------------------------------------
  spec_sgarch <- rugarch::ugarchspec(
    variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  spec_egarch <- rugarch::ugarchspec(
    variance.model = list(model = "eGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  spec_gjr <- rugarch::ugarchspec(
    variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
    mean.model     = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"
  )
  
  vhat_garch_val  <- garch_forecast_oos_map(retX[, target_nm], t_end, idx_val, refit_every, spec_sgarch)
  vhat_egarch_val <- garch_forecast_oos_map(retX[, target_nm], t_end, idx_val, refit_every, spec_egarch)
  vhat_gjr_val    <- garch_forecast_oos_map(retX[, target_nm], t_end, idx_val, refit_every, spec_gjr)
  
  # ----- Test forecasts -------------------------------------------------------
  pred_ar1_test    <- oos_forecast_expanding(y_logv, NULL,            idx_test,   refit_every, fit_ar1, pred_ar1_1step)
  pred_har_test    <- oos_forecast_expanding(y_logv, Xhar,            idx_test,   refit_every, fit_lm,  pred_lm_1step)
  pred_harX_test   <- oos_forecast_expanding(y_logv, XharX,           idx_test,   refit_every, fit_lm,  pred_lm_1step)
  
  pred_fvar_test   <- oos_forecast_expanding(y_logv, Xflat_sel,       idx_test,   refit_every, fit_pls, pred_pls_1step)
  pred_fvarX_test  <- oos_forecast_expanding(y_logv, Xflat_X,         idx_test,   refit_every, fit_pls, pred_pls_1step)
  pred_fvarL_test  <- oos_forecast_expanding(y_logv_lag, Xflat_lag,   idx_test_L, refit_every, fit_pls, pred_pls_1step)
  pred_fvarLX_test <- oos_forecast_expanding(y_logv_lag, Xflat_lag_X, idx_test_L, refit_every, fit_pls, pred_pls_1step)
  
  vhat_ar1_test    <- exp(pred_ar1_test)
  vhat_har_test    <- exp(pred_har_test)
  vhat_harX_test   <- exp(pred_harX_test)
  vhat_fvar_test   <- exp(pred_fvar_test)
  vhat_fvarX_test  <- exp(pred_fvarX_test)
  vhat_fvarL_test  <- exp(pred_fvarL_test)
  vhat_fvarLX_test <- exp(pred_fvarLX_test)
  
  vhat_pers_test   <- v_t[idx_test]
  
  vhat_garch_test  <- garch_forecast_oos_map(retX[, target_nm], t_end, idx_test, refit_every, spec_sgarch)
  vhat_egarch_test <- garch_forecast_oos_map(retX[, target_nm], t_end, idx_test, refit_every, spec_egarch)
  vhat_gjr_test    <- garch_forecast_oos_map(retX[, target_nm], t_end, idx_test, refit_every, spec_gjr)
  
  # ----- COMMON-aligned evaluation -------------------------------------------
  dates_oos     <- index(rv)[t_end + 1L]
  dates_oos_lag <- dates_oos[-1]
  
  V_val  <- xts(y_rv[idx_val],  order.by = dates_oos[idx_val])
  V_test <- xts(y_rv[idx_test], order.by = dates_oos[idx_test])
  
  S_val <- list(
    Persistence = xts(vhat_pers_val,   order.by = dates_oos[idx_val]),
    AR1_logv    = xts(vhat_ar1_val,    order.by = dates_oos[idx_val]),
    HAR         = xts(vhat_har_val,    order.by = dates_oos[idx_val]),
    HAR_X       = xts(vhat_harX_val,   order.by = dates_oos[idx_val]),
    GARCH11     = xts(vhat_garch_val,  order.by = dates_oos[idx_val]),
    EGARCH11    = xts(vhat_egarch_val, order.by = dates_oos[idx_val]),
    GJR11       = xts(vhat_gjr_val,    order.by = dates_oos[idx_val]),
    fVAR        = xts(vhat_fvar_val,   order.by = dates_oos[idx_val]),
    fVAR_X      = xts(vhat_fvarX_val,  order.by = dates_oos[idx_val]),
    fVAR_L      = xts(vhat_fvarL_val,  order.by = dates_oos_lag[idx_val_L]),
    fVAR_LX     = xts(vhat_fvarLX_val, order.by = dates_oos_lag[idx_val_L])
  )
  
  S_test <- list(
    Persistence = xts(vhat_pers_test,   order.by = dates_oos[idx_test]),
    AR1_logv    = xts(vhat_ar1_test,    order.by = dates_oos[idx_test]),
    HAR         = xts(vhat_har_test,    order.by = dates_oos[idx_test]),
    HAR_X       = xts(vhat_harX_test,   order.by = dates_oos[idx_test]),
    GARCH11     = xts(vhat_garch_test,  order.by = dates_oos[idx_test]),
    EGARCH11    = xts(vhat_egarch_test, order.by = dates_oos[idx_test]),
    GJR11       = xts(vhat_gjr_test,    order.by = dates_oos[idx_test]),
    fVAR        = xts(vhat_fvar_test,   order.by = dates_oos[idx_test]),
    fVAR_X      = xts(vhat_fvarX_test,  order.by = dates_oos[idx_test]),
    fVAR_L      = xts(vhat_fvarL_test,  order.by = dates_oos_lag[idx_test_L]),
    fVAR_LX     = xts(vhat_fvarLX_test, order.by = dates_oos_lag[idx_test_L])
  )
  
  M_val  <- na.omit(do.call(merge, c(list(Actual = V_val),  S_val,  all = FALSE)))
  M_test <- na.omit(do.call(merge, c(list(Actual = V_test), S_test, all = FALSE)))
  
  write.csv(data.frame(Date = index(M_test), coredata(M_test)),
            file.path(OUT, "data", "paths_TEST_common_aligned.csv"),
            row.names = FALSE)
  
  # ----- Metrics --------------------------------------------------------------
  tab_val  <- collect_metric_table_from_M(M_val)
  tab_test <- collect_metric_table_from_M(M_test)
  
  best_fun_val  <- choose_best_functional(tab_val)
  best_fun_test <- choose_best_functional(tab_test)
  
  best_fun_val_label  <- best_fun_val$Model[1]
  best_fun_val_col    <- MODEL_MAP[[best_fun_val_label]]
  best_fun_test_label <- best_fun_test$Model[1]
  best_fun_test_col   <- MODEL_MAP[[best_fun_test_label]]
  
  write_csv(tab_test, OUT, "T2_var_oos_metrics_TEST_raw")
  
  # manuscript-ready QLIKE table
  tab_qlike_tex <- copy(tab_test[, .(Model, QLIKE, Delta_QLIKE_vs_GARCH)])
  setorder(tab_qlike_tex, QLIKE)
  tab_qlike_tex[, Model := display_model(Model)]
  tab_qlike_tex[, QLIKE := fmt4(QLIKE)]
  tab_qlike_tex[, `\\Delta QLIKE vs GARCH` := fmt4(Delta_QLIKE_vs_GARCH)]
  tab_qlike_tex <- tab_qlike_tex[, .(Model, QLIKE, `\\Delta QLIKE vs GARCH`)]
  tab_qlike_tex <- pad_even_dt(tab_qlike_tex)
  
  write_longtable_2block(
    tab_qlike_tex,
    OUT, "T1_var_oos_qlike",
    caption = "Out-of-sample variance forecasting performance under QLIKE (lower is better) on the \\texttt{COMMON-aligned TEST} set. Results are reported under the baseline \\DAXCLOSE\\ information-availability convention. The final column reports differences relative to GARCH(1,1), computed as model QLIKE minus GARCH QLIKE.",
    label = sprintf("tab:var_oos_qlike_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    left_cols  = c("Model","QLIKE","\\Delta QLIKE vs GARCH"),
    right_cols = c("Model","QLIKE","\\Delta QLIKE vs GARCH")
  )
  
  # manuscript-ready MAE/RMSE
  tab_mrmse_tex <- copy(tab_test[, .(Model, MAE_var, RMSE_var)])
  tab_mrmse_tex[, Model := display_model(Model)]
  tab_mrmse_tex[, MAE := sprintf("%.7f", MAE_var)]
  tab_mrmse_tex[, RMSE := sprintf("%.7f", RMSE_var)]
  tab_mrmse_tex <- tab_mrmse_tex[, .(Model, MAE, RMSE)]
  tab_mrmse_tex <- pad_even_dt(tab_mrmse_tex)
  
  write_longtable_2block(
    tab_mrmse_tex,
    OUT, "T2_var_oos_mae_rmse",
    caption = "Out-of-sample variance forecasting performance under MAE and RMSE on the variance scale, evaluated on the \\texttt{COMMON-aligned TEST} set. Results are reported under the baseline \\DAXCLOSE\\ information-availability convention.",
    label = sprintf("tab:var_oos_mae_rmse_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    left_cols  = c("Model","MAE","RMSE"),
    right_cols = c("Model","MAE","RMSE")
  )
  
  # ----- DM tests vs GARCH ----------------------------------------------------
  V_true <- as.numeric(M_test[, "Actual"])
  loss_g <- qlike_vec(V_true, as.numeric(M_test[, "GARCH11"]))
  
  dm_one <- function(model_col, model_name){
    loss_m <- qlike_vec(V_true, as.numeric(M_test[, model_col]))
    dmL <- dm_loss(loss_m, loss_g, alternative = "less",      L = DM_NW_LAG)
    dm2 <- dm_loss(loss_m, loss_g, alternative = "two.sided", L = DM_NW_LAG)
    data.table(
      Comparison = paste(model_name, "vs GARCH(1,1)"),
      DM         = dmL$stat,
      P_one      = dmL$p,
      P_two      = dm2$p,
      N          = dmL$n
    )
  }
  
  dm_tab <- rbindlist(list(
    dm_one("AR1_logv", "AR(1)-logv"),
    dm_one("HAR",      "HAR"),
    dm_one("HAR_X",    "HAR-X(logVIX)"),
    dm_one("EGARCH11", "EGARCH(1,1)"),
    dm_one("GJR11",    "GJR-GARCH(1,1)"),
    dm_one("fVAR",     "fVAR"),
    dm_one("fVAR_X",   "fVAR-X(+logv_t)"),
    dm_one("fVAR_L",   "fVAR-Lagged"),
    dm_one("fVAR_LX",  "fVAR-Lagged-X")
  ), fill = TRUE)
  
  dm_tab[, Holm := p.adjust(P_one, method = "holm")]
  dm_tab[, BH   := p.adjust(P_one, method = "BH")]
  
  write_csv(dm_tab, OUT, "T3_dm_qlike_TEST_raw")
  
  dm_tab_tex <- copy(dm_tab)
  dm_tab_tex[, Comparison := gsub("fVAR-X\\(\\+logv_t\\)", "fVAR-X(+logv$_t$)", Comparison)]
  dm_tab_tex[, DM := fmt4(DM)]
  dm_tab_tex[, `$p_{\\text{one}}$` := fmt4(P_one)]
  dm_tab_tex[, `$p_{\\text{two}}$` := fmt4(P_two)]
  dm_tab_tex[, Holm := fmt4(Holm)]
  dm_tab_tex[, BH := fmt4(BH)]
  dm_tab_tex <- dm_tab_tex[, .(Comparison, DM, `$p_{\\text{one}}$`, `$p_{\\text{two}}$`, Holm, BH, N)]
  
  write_longtable_simple(
    dm_tab_tex,
    OUT, "T3_dm_qlike",
    caption = "Diebold--Mariano tests on QLIKE loss differentials relative to GARCH(1,1), evaluated on the \\texttt{COMMON-aligned TEST} set under the baseline \\DAXCLOSE\\ convention. One-sided and two-sided $p$-values are reported together with Holm- and Benjamini--Hochberg-adjusted one-sided values.",
    label = sprintf("tab:dm_qlike_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    align = "lcccccc"
  )
  
  # ----- Bootstrap CI + sensitivity (computed once, reused everywhere) -------
  loss_bestf <- qlike_vec(V_true, as.numeric(M_test[, best_fun_val_col]))
  diff_bestf <- loss_bestf - loss_g
  
  boot_sens_raw <- rbindlist(lapply(BOOT_B_GRID, function(BB){
    bt <- mbb_mean_diff(diff_bestf, B = BB, R = 2000L)
    data.table(
      B        = BB,
      MeanDiff = bt$mean,
      CI_L     = bt$ci[1],
      CI_U     = bt$ci[2]
    )
  }), fill = TRUE)
  
  boot10_row <- boot_sens_raw[B == 10L]
  boot_tab_ci <- data.table(
    Comparison = paste(best_fun_val_label, "minus GARCH(1,1)"),
    MeanDiff   = boot10_row$MeanDiff,
    CI_L       = boot10_row$CI_L,
    CI_U       = boot10_row$CI_U
  )
  write_csv(boot_tab_ci, OUT, "T4_boot_ci_B10_raw")
  
  boot_tab_ci_tex <- copy(boot_tab_ci)
  boot_tab_ci_tex[, Comparison := gsub("fVAR-X\\(\\+logv_t\\)", "fVAR-X(+logv$_t$)", Comparison)]
  boot_tab_ci_tex[, `Mean diff.` := fmt4(MeanDiff)]
  boot_tab_ci_tex[, `CI lower` := fmt4(CI_L)]
  boot_tab_ci_tex[, `CI upper` := fmt4(CI_U)]
  boot_tab_ci_tex <- boot_tab_ci_tex[, .(Comparison, `Mean diff.`, `CI lower`, `CI upper`)]
  
  write_longtable_simple(
    boot_tab_ci_tex,
    OUT, "T4_boot_ci",
    caption = "Moving-block bootstrap 95\\% confidence interval for the mean QLIKE difference on the \\texttt{COMMON-aligned TEST} set (validation-selected functional representative minus GARCH). Results are reported under the baseline \\DAXCLOSE\\ convention with block length $B=10$.",
    label = sprintf("tab:boot_ci_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    align = "lccc"
  )
  
  # ----- Diagnostics ----------------------------------------------------------
  r_next        <- as.numeric(retX[, target_nm])[t_end + 1L]
  r_test_xts    <- xts(r_next, order.by = dates_oos)
  r_test_common <- r_test_xts[index(M_test)]
  
  diag_models <- list(
    "GARCH(1,1)"      = "GARCH11",
    "EGARCH(1,1)"     = "EGARCH11",
    "HAR"             = "HAR",
    "HAR-X(logVIX)"   = "HAR_X",
    "fVAR-X(+logv_t)" = "fVAR_X",
    "fVAR-Lagged-X"   = "fVAR_LX"
  )
  
  diag_rows <- lapply(names(diag_models), function(lbl){
    col_nm <- diag_models[[lbl]]
    vhat_i <- as.numeric(M_test[, col_nm])
    z_i    <- z_from_vhat(as.numeric(r_test_common), vhat_i)
    
    a10 <- arch_lm(z_i, m = 10L)
    lb  <- try(Box.test(as.numeric(z_i)^2, lag = 10L, type = "Ljung-Box"), silent = TRUE)
    lbp <- if(inherits(lb, "try-error")) NA_real_ else as.numeric(lb$p.value)
    
    data.table(
      Model      = lbl,
      ARCH_LM_10 = a10$stat,
      P_10       = a10$p,
      LB_p_10    = lbp,
      N          = a10$n
    )
  })
  
  diag_tab <- rbindlist(diag_rows, fill = TRUE)
  write_csv(diag_tab, OUT, "A4_diag_TEST_raw")
  
  diag_tab_tex <- copy(diag_tab)
  diag_tab_tex[, Model := display_model(Model)]
  diag_tab_tex[, `ARCH--LM(10)` := fmt4(ARCH_LM_10)]
  diag_tab_tex[, `$p$-value` := fmt_p_diag(P_10)]
  diag_tab_tex[, `LB $p$-value` := fmt_p_diag(LB_p_10)]
  diag_tab_tex <- diag_tab_tex[, .(Model, `ARCH--LM(10)`, `$p$-value`, `LB $p$-value`, N)]
  
  write_longtable_simple(
    diag_tab_tex,
    OUT, "A4_diag_TEST",
    caption = "Residual diagnostics based on standardized return residuals $z_t=r_t/\\sqrt{\\hat v_t}$: ARCH--LM(10) and Ljung--Box tests on $z_t^2$ (TEST; COMMON-aligned). Reported as forecasting-oriented sanity checks rather than adequacy tests.",
    label = sprintf("tab:diag_%s_%s_ref%d_%s", TARGET_NAME, cfg_name, refit_every, avail_mode),
    align = "lcccc"
  )
  
  # ----- Figures --------------------------------------------------------------
  df_path <- data.frame(
    Date         = index(M_test),
    Actual       = as.numeric(M_test[, "Actual"]),
    GARCH        = as.numeric(M_test[, "GARCH11"]),
    HAR          = as.numeric(M_test[, "HAR"]),
    EGARCH       = as.numeric(M_test[, "EGARCH11"]),
    BestFunction = as.numeric(M_test[, best_fun_val_col])
  )
  
  p1 <- ggplot(df_path, aes(x = Date)) +
    geom_line(aes(y = Actual, color = "Actual"), linewidth = 1.2) +
    geom_line(aes(y = GARCH, color = "GARCH(1,1)"), linewidth = 1.0, linetype = "dashed") +
    geom_line(aes(y = EGARCH, color = "EGARCH(1,1)"), linewidth = 0.9, linetype = "dotdash") +
    geom_line(aes(y = HAR, color = "HAR"), linewidth = 0.9, linetype = "dotted") +
    geom_line(aes(y = BestFunction, color = display_model(best_fun_val_label)), linewidth = 1.0) +
    labs(
      y = "RS variance", color = NULL,
      title = sprintf("Variance paths (TEST, common aligned) ??? Target=%s | cfg=%s | refit=%d | avail=%s",
                      TARGET_NAME, cfg_name, refit_every, avail_mode)
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  
  ggsave(file.path(OUT, "plots", "F1_paths_TEST.png"),
         p1, width = 11, height = 5.5, dpi = PLOT_DPI)
  
  loss_dt <- data.table(
    Date           = index(M_test),
    GARCH          = loss_g,
    HAR            = qlike_vec(V_true, as.numeric(M_test[, "HAR"])),
    EGARCH         = qlike_vec(V_true, as.numeric(M_test[, "EGARCH11"])),
    BestFunctional = qlike_vec(V_true, as.numeric(M_test[, best_fun_val_col]))
  )
  
  setnames(loss_dt, "BestFunctional", display_model(best_fun_val_label))
  
  loss_long <- melt(loss_dt, id.vars = "Date", variable.name = "Model", value.name = "QLIKE")
  x_lo <- as.numeric(quantile(loss_long$QLIKE, 0.005, na.rm = TRUE))
  x_hi <- as.numeric(quantile(loss_long$QLIKE, 0.995, na.rm = TRUE))
  
  p2 <- ggplot(loss_long, aes(x = QLIKE, fill = Model)) +
    geom_density(alpha = 0.30) +
    coord_cartesian(xlim = c(x_lo, x_hi)) +
    labs(x = "QLIKE (lower is better)", y = "Density",
         title = "QLIKE distribution (TEST, common aligned; central 99%)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "top")
  
  ggsave(file.path(OUT, "plots", "F2_qlike_density_TEST.png"),
         p2, width = 11, height = 5.0, dpi = PLOT_DPI)
  
  # ----- Validation summary ---------------------------------------------------
  qlike_garch_val <- tab_val[Model == "GARCH(1,1)", QLIKE]
  qlike_bestf_val <- best_fun_val$QLIKE[1]
  
  loss_g_val     <- qlike_vec(as.numeric(M_val[, "Actual"]), as.numeric(M_val[, "GARCH11"]))
  loss_bestf_val <- qlike_vec(as.numeric(M_val[, "Actual"]), as.numeric(M_val[, best_fun_val_col]))
  
  dm_val_less <- dm_loss(loss_bestf_val, loss_g_val, alternative = "less",      L = DM_NW_LAG)
  dm_val_two  <- dm_loss(loss_bestf_val, loss_g_val, alternative = "two.sided", L = DM_NW_LAG)
  
  run_sum <- data.table(
    Target                   = TARGET_NAME,
    Cfg                      = cfg_name,
    PredictorSetConstruction = cfg_pretty(cfg_name),
    Refit                    = refit_every,
    Avail                    = avail_mode,
    RepresentativeSeries_GCV = rep_var,
    ValidationBlockRule      = "Last 20% of DEV (chronological internal validation block)",
    BestFunctional_VAL       = best_fun_val_label,
    QLIKE_BestFunctional_VAL = qlike_bestf_val,
    QLIKE_GARCH_VAL          = qlike_garch_val,
    MeanDiffBestF_VAL        = qlike_bestf_val - qlike_garch_val,
    DM_statBestF_VAL         = dm_val_less$stat,
    DM_p_lessBestF_VAL       = dm_val_less$p,
    DM_p_twoBestF_VAL        = dm_val_two$p,
    DM_N_VAL                 = dm_val_less$n,
    BestFunctional_TEST      = best_fun_test_label,
    N_total                  = N,
    N_dev                    = n_dev,
    N_val_block              = n_val,
    N_test_block             = n_test,
    N_val_common             = NROW(M_val),
    N_test_common            = NROW(M_test)
  )
  
  write_csv(run_sum, OUT, "run_summary")
  
  sink(file.path(OUT, "debug", "sessionInfo.txt"))
  print(sessionInfo())
  sink()
  
  list(
    OUT                = OUT,
    summary            = run_sum,
    M_val              = M_val,
    M_test             = M_test,
    tab_val            = tab_val,
    tab_test           = tab_test,
    best_fun_val_label = best_fun_val_label,
    best_fun_val_col   = best_fun_val_col,
    diff_bestf         = diff_bestf,
    boot_sens_raw      = boot_sens_raw
  )
}

# ============================= LOCKED MAIN RUN ================================
MAIN_RUN <- try(run_one(CFG_MAIN, REFIT_MAIN, avail_mode = AVAIL_MAIN), silent = TRUE)
if(inherits(MAIN_RUN, "try-error")){
  stop("MAIN run failed: ", MAIN_RUN)
}

MAIN_OUT <- MAIN_RUN$OUT

message("\n================ LOCKED MAIN RUN ================")
print(MAIN_RUN$summary)
message("MAIN_OUT = ", MAIN_OUT)

write_csv(copy(MAIN_RUN$summary), OUT_BASE, "MAIN_run_summary")

# ============================= LOCKED ROBUSTNESS RUN ==========================
ROB_RUN <- try(run_one(CFG_MAIN, REFIT_MAIN, avail_mode = AVAIL_ROB), silent = TRUE)
if(inherits(ROB_RUN, "try-error")){
  stop("ROBUSTNESS run failed: ", ROB_RUN)
}

ROB_OUT <- ROB_RUN$OUT

message("\n================ LOCKED ROBUSTNESS RUN ================")
print(ROB_RUN$summary)
message("ROB_OUT = ", ROB_OUT)

write_csv(copy(ROB_RUN$summary), OUT_BASE, "ROBUSTNESS_run_summary")

# ============================= LOCKED VALIDATION SUMMARY ======================
LOCKED_SUMMARY <- data.table(
  `Predictor-set construction` = "PRUNE095 candidate construction",
  Refit                        = REFIT_MAIN,
  `Best functional model`      = display_model(MAIN_RUN$summary$BestFunctional_VAL),
  `Mean diff. vs GARCH`        = fmt4(MAIN_RUN$summary$MeanDiffBestF_VAL),
  `DM $p$-value`               = fmt_sci_or_fixed(MAIN_RUN$summary$DM_p_lessBestF_VAL),
  `$N_{val,\\mathrm{common}}$`  = MAIN_RUN$summary$N_val_common,
  `$N_{test,\\mathrm{common}}$` = MAIN_RUN$summary$N_test_common
)
write_csv(LOCKED_SUMMARY, OUT_BASE, "A1_sweep_summary_VAL_locked")
write_longtable_simple(
  LOCKED_SUMMARY,
  OUT_BASE, "A1_sweep_summary_VAL_locked",
  caption = "Internal-validation summary under the locked baseline manuscript design and the MAIN availability protocol (\\DAXCLOSE). $N_{val,\\mathrm{common}}$ denotes the common-aligned count within the DEV-internal validation block; $N_{test,\\mathrm{common}}$ is shown separately only for reference.",
  label = "tab:sweep_summary_VAL",
  align = "lcccccc"
)

# ============================= PROTOCOL ROBUSTNESS TABLES =====================
tab_main <- copy(MAIN_RUN$tab_test)
tab_rob  <- copy(ROB_RUN$tab_test)
M_main   <- MAIN_RUN$M_test
M_rob    <- ROB_RUN$M_test

# A4a) selected best functional vs GARCH
q_main_g  <- tab_main[Model == "GARCH(1,1)", QLIKE]
q_main_bf <- tab_main[Model == MAIN_RUN$best_fun_val_label, QLIKE]

q_rob_g   <- tab_rob[Model == "GARCH(1,1)", QLIKE]
q_rob_bf  <- tab_rob[Model == ROB_RUN$best_fun_val_label, QLIKE]

AVTAB_BF <- data.table(
  Protocol                    = c("\\DAXCLOSE", "\\GLOBALCLOSE"),
  `Selected functional model` = c(display_model(MAIN_RUN$best_fun_val_label),
                                  display_model(ROB_RUN$best_fun_val_label)),
  `QLIKE (GARCH)`             = c(fmt4(q_main_g), fmt4(q_rob_g)),
  `QLIKE (selected)`          = c(fmt4(q_main_bf), fmt4(q_rob_bf)),
  `Mean diff.`                = c(fmt4(q_main_bf - q_main_g), fmt4(q_rob_bf - q_rob_g)),
  Notes                       = c("Baseline reported timing convention",
                                  "Supplementary timing convention")
)
write_csv(AVTAB_BF, OUT_BASE, "A4a_bestfunctional_protocol_compare")
write_longtable_simple(
  AVTAB_BF,
  OUT_BASE, "A4a_bestfunctional_protocol_compare",
  caption = "Protocol comparison for GARCH and the validation-selected functional representative. Baseline convention = \\DAXCLOSE; supplementary convention = \\GLOBALCLOSE.",
  label = "tab:availability_protocol_compare_bestfunctional",
  align = "lccccc"
)

# A4b) full model-set own-common comparison
PROT_FULL <- make_protocol_modelset_compare(
  tab_main, tab_rob,
  main_name = AVAIL_MAIN,
  rob_name  = AVAIL_ROB
)

PROT_FULL_TEX <- copy(PROT_FULL)
PROT_FULL_TEX[, Model := display_model(Model)]
PROT_FULL_TEX[, `QLIKE (\\DAXCLOSE)` := fmt4(get("QLIKE_DAX_CLOSE"))]
PROT_FULL_TEX[, `Rank (\\DAXCLOSE)`  := get("Rank_DAX_CLOSE")]
PROT_FULL_TEX[, `QLIKE (\\GLOBALCLOSE)` := fmt4(get("QLIKE_GLOBAL_CLOSE"))]
PROT_FULL_TEX[, `Rank (\\GLOBALCLOSE)`  := get("Rank_GLOBAL_CLOSE")]
PROT_FULL_TEX[, `\\Delta QLIKE` := fmt4(Delta_MAIN_minus_ROB)]
PROT_FULL_TEX <- PROT_FULL_TEX[, .(Model, `QLIKE (\\DAXCLOSE)`, `Rank (\\DAXCLOSE)`,
                                   `QLIKE (\\GLOBALCLOSE)`, `Rank (\\GLOBALCLOSE)`, `\\Delta QLIKE`)]
write_csv(PROT_FULL_TEX, OUT_BASE, "A4b_protocol_modelset_owncommon")
write_longtable_simple(
  PROT_FULL_TEX,
  OUT_BASE, "A4b_protocol_modelset_owncommon",
  caption = "Protocol comparison for the full model set using each protocol's own COMMON-aligned TEST set. Lower QLIKE is better; ranks are computed within protocol. Here, $\\Delta$QLIKE = QLIKE(\\DAXCLOSE) - QLIKE(\\GLOBALCLOSE).",
  label = "tab:protocol_modelset_owncommon",
  align = "lccccc"
)

# A4c) relative-to-GARCH protocol comparison
PROT_REL <- make_protocol_vs_garch_compare(
  tab_main, tab_rob,
  main_name = AVAIL_MAIN,
  rob_name  = AVAIL_ROB
)

PROT_REL_TEX <- copy(PROT_REL)
PROT_REL_TEX[, Model := display_model(Model)]
PROT_REL_TEX[, `$\\Delta$QLIKE vs GARCH (\\DAXCLOSE)` := fmt4(get("Delta_vs_GARCH_DAX_CLOSE"))]
PROT_REL_TEX[, `$\\Delta$QLIKE vs GARCH (\\GLOBALCLOSE)` := fmt4(get("Delta_vs_GARCH_GLOBAL_CLOSE"))]
PROT_REL_TEX[, Gap := fmt4(DeltaGap_MAIN_minus_ROB)]
PROT_REL_TEX <- PROT_REL_TEX[, .(Model, `$\\Delta$QLIKE vs GARCH (\\DAXCLOSE)`,
                                 `$\\Delta$QLIKE vs GARCH (\\GLOBALCLOSE)`, Gap)]
write_csv(PROT_REL_TEX, OUT_BASE, "A4c_protocol_vs_garch_owncommon")
write_longtable_simple(
  PROT_REL_TEX,
  OUT_BASE, "A4c_protocol_vs_garch_owncommon",
  caption = "Protocol comparison of relative performance versus GARCH(1,1), using each protocol's own COMMON-aligned TEST set. More negative values indicate larger gains over GARCH. Here, Gap = [$\\Delta$QLIKE vs GARCH under \\DAXCLOSE] - [$\\Delta$QLIKE vs GARCH under \\GLOBALCLOSE].",
  label = "tab:protocol_vs_garch_owncommon",
  align = "lccc"
)

# Optional verification only
CROSS_OBJ <- make_crossprotocol_common_compare(
  M_main, M_rob,
  main_name = AVAIL_MAIN,
  rob_name  = AVAIL_ROB
)
if(!is.null(CROSS_OBJ)){
  write_csv(CROSS_OBJ$full_cc, OUT_BASE, "A4d_protocol_modelset_crosscommon_raw")
  write_csv(CROSS_OBJ$rel_cc,  OUT_BASE, "A4e_protocol_vs_garch_crosscommon_raw")
}

# ============================= APPENDIX SENSITIVITY ===========================
# A2) Bootstrap block-size sensitivity (reused from MAIN_RUN; no recomputation)
BOOT_SENS_RAW <- copy(MAIN_RUN$boot_sens_raw)
BOOT_SENS <- BOOT_SENS_RAW[, .(
  `Block length $B$` = B,
  `Mean diff.`       = fmt4(MeanDiff),
  `95\\% CI`         = paste0("[", fmt4(CI_L), ", ", fmt4(CI_U), "]")
)]
write_csv(BOOT_SENS, OUT_BASE, "A2_bootB_sensitivity")
write_longtable_simple(
  BOOT_SENS,
  OUT_BASE, "A2_bootB_sensitivity",
  caption = "Moving-block bootstrap sensitivity over block length $B$ on TEST for the validation-selected functional representative relative to GARCH.",
  label = "tab:bootB_sensitivity",
  align = "ccc"
)

# A5) DM NW-lag sensitivity
M_best      <- MAIN_RUN$M_test
V_true_best <- as.numeric(M_best[, "Actual"])
loss_g_best <- qlike_vec(V_true_best, as.numeric(M_best[, "GARCH11"]))
loss_f_best <- qlike_vec(V_true_best, as.numeric(M_best[, MAIN_RUN$best_fun_val_col]))

DM_SENS <- rbindlist(lapply(DM_LAG_GRID, function(L){
  d1 <- dm_loss(loss_f_best, loss_g_best, alternative = "less",      L = L)
  d2 <- dm_loss(loss_f_best, loss_g_best, alternative = "two.sided", L = L)
  data.table(
    `NW lag $L$`          = L,
    DM                    = fmt4(d1$stat),
    `$p_{\\mathrm{one}}$` = fmt4(d1$p),
    `$p_{\\mathrm{two}}$` = fmt4(d2$p),
    `$N$`                 = d1$n
  )
}), fill = TRUE)

write_csv(DM_SENS, OUT_BASE, "A5_dm_nwlag_sensitivity")
write_longtable_simple(
  DM_SENS,
  OUT_BASE, "A5_dm_nwlag_sensitivity",
  caption = "Diebold--Mariano sensitivity over the Newey--West lag $L$ on TEST for the validation-selected functional representative relative to GARCH.",
  label = "tab:dm_nwlag_sensitivity",
  align = "ccccc"
)

# ============================= MANIFEST =======================================
to_tex_path <- function(x) gsub("\\\\", "/", x)

write_manifest_locked <- function(main_out, out_base){
  tex_main <- c(
    "% ===================== LOCKED MANUSCRIPT MANIFEST =====================",
    "% MAIN TABLES (PRUNE095, ref20, DAX_CLOSE)",
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "T1_var_oos_qlike.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "T2_var_oos_mae_rmse.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "T3_dm_qlike.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "T4_boot_ci.tex"))),
    "",
    "% MAIN FIGURES",
    sprintf("%% F1: %s", to_tex_path(file.path(main_out, "plots", "F1_paths_TEST.png"))),
    sprintf("%% F2: %s", to_tex_path(file.path(main_out, "plots", "F2_qlike_density_TEST.png"))),
    "",
    "% APPENDIX TABLES",
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "A1_desc_returns.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "A2_gcv_grid.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "A3_selected_vars.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(main_out, "tex", "A4_diag_TEST.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(out_base, "tex", "A1_sweep_summary_VAL_locked.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(out_base, "tex", "A2_bootB_sensitivity.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(out_base, "tex", "A4a_bestfunctional_protocol_compare.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(out_base, "tex", "A4b_protocol_modelset_owncommon.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(out_base, "tex", "A4c_protocol_vs_garch_owncommon.tex"))),
    sprintf("\\input{%s}", to_tex_path(file.path(out_base, "tex", "A5_dm_nwlag_sensitivity.tex"))),
    "",
    "% OPTIONAL verification CSV outputs",
    sprintf("%% %s", to_tex_path(file.path(out_base, "tables", "A4d_protocol_modelset_crosscommon_raw.csv"))),
    sprintf("%% %s", to_tex_path(file.path(out_base, "tables", "A4e_protocol_vs_garch_crosscommon_raw.csv"))),
    "",
    "% APPENDIX FIGURES",
    sprintf("%% F3: %s", to_tex_path(file.path(main_out, "plots", "F3_elasticnet_cv.png"))),
    sprintf("%% F4: %s", to_tex_path(file.path(main_out, "plots", "F4_gcv_heatmap.png"))),
    "% =================================================================="
  )
  
  fn <- file.path(out_base, "tex", "manifest_locked_prune095.tex")
  dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
  writeLines(tex_main, fn)
  fn
}

manifest_fn <- write_manifest_locked(MAIN_OUT, OUT_BASE)

message("\nLocked manifest written: ", normalizePath(manifest_fn))
message("\nDONE ??? Locked outputs under: ", normalizePath(OUT_BASE))

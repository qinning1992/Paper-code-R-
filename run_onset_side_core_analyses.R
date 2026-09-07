#!/usr/bin/env Rscript
# =============================================================================
# Onset-side-stratified re-analysis of the core analyses
# Manuscript: Differential Structural Profiles of the Hippocampus and Thalamus
#             and Their Cognitive Relevance in Unilateral MTLE (PNP-S-26-00652)
#
# Written in response to Reviewer #1 (Comments 2, 4, 6) and Reviewer #2
# (Comment 1), who asked that the key analyses be repeated separately in
# left-onset and right-onset MTLE rather than only in the pooled cohort.
#
# Four core analyses are repeated by seizure-onset side:
#   A. Additive models of hippocampal and thalamic volume (M0-M3)
#   B. HC-referenced ipsilateral-versus-contralateral atrophy z scores
#   C. Cognitive performance, left-onset versus right-onset
#   D. Fine-grained structural abnormality mapping versus healthy controls
#   E. Supporting predictor x onset-side interaction tests (full sample)
#
# Statistical choices are kept identical to the original project scripts:
# 5,000 Freedman-Lane residual permutations, the same covariate sets, and
# family-wise Benjamini-Hochberg FDR correction.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(purrr)
  library(readr); library(stringr); library(ggplot2); library(scales); library(tibble)
})

options(stringsAsFactors = FALSE, scipen = 999)
# Ensure Delta / multiplication / circle glyphs render in figure text.
invisible(suppressWarnings(Sys.setlocale("LC_CTYPE",
  Sys.getenv("R_LOCALE_CTYPE", unset = "C.UTF-8"))))
set.seed(2026)

N_PERM   <- 5000L
N_BOOT   <- 2000L
MIN_N    <- 15L
OUTDIR    <- "out"
TABLE_DIR <- file.path(OUTDIR, "01_tables")
FIG_DIR   <- file.path(OUTDIR, "02_figures")
QC_DIR    <- file.path(OUTDIR, "03_qc")
for (p in c(TABLE_DIR, FIG_DIR, QC_DIR)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUTDIR, "run_log.txt")
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)
log_msg <- function(...) {
  z <- paste0(...)
  cat(z, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", z, "\n", sep = "",
      file = LOG_FILE, append = TRUE)
}
write_table <- function(x, filename) readr::write_excel_csv(x, file.path(TABLE_DIR, filename), na = "")
save_plot <- function(p, stem, width, height) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p, width = width, height = height,
         units = "in", dpi = 600, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf, bg = "white")
}
star_from_q <- function(q) case_when(
  is.finite(q) & q < 0.001 ~ "***", is.finite(q) & q < 0.01 ~ "**",
  is.finite(q) & q < 0.05 ~ "*",    is.finite(q) & q < 0.10 ~ "○", TRUE ~ "")

# ---------------------------------------------------------------- helpers ----
normalize_side_lr <- function(x) {
  out <- rep(NA_character_, length(x))
  out[is.finite(x) & x == 1] <- "L"; out[is.finite(x) & x == 0] <- "R"; out
}
hedges_g <- function(x, y) {
  nx <- length(x); ny <- length(y)
  if (nx < 2 || ny < 2) return(NA_real_)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  d <- (mean(x) - mean(y)) / sp
  d * (1 - 3 / (4 * (nx + ny) - 9))
}
hedges_g_se <- function(g, nx, ny) sqrt((nx + ny) / (nx * ny) + g^2 / (2 * (nx + ny)))
effective_covars <- function(d, covars) {
  keep <- character(0)
  for (cv in intersect(covars, names(d))) {
    x <- d[[cv]]
    if (is.factor(x) || is.character(x)) {
      if (length(unique(x[!is.na(x)])) >= 2) keep <- c(keep, cv)
    } else {
      xn <- suppressWarnings(as.numeric(x))
      if (sum(is.finite(xn)) >= MIN_N && is.finite(sd(xn, na.rm = TRUE)) &&
          sd(xn, na.rm = TRUE) > 0) keep <- c(keep, cv)
    }
  }
  keep
}
residualize_one <- function(y, d, covars) {
  y <- as.numeric(y); covars <- effective_covars(d, covars)
  if (length(covars) == 0) return(y - mean(y, na.rm = TRUE))
  X <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse = " + "))), data = d)
  y - as.numeric(X %*% lm.fit(X, y)$coefficients)
}
row_sum_complete <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (length(vars) == 0) return(rep(NA_real_, nrow(df)))
  x <- as.matrix(df[, vars, drop = FALSE]); storage.mode(x) <- "double"
  out <- rowSums(x, na.rm = FALSE); out[rowSums(is.finite(x)) < length(vars)] <- NA_real_; out
}
make_ipsi_contra_pairs <- function(df, side_col, left_cols, right_cols,
                                   left_prefix, right_prefix, out_prefix) {
  lsuf <- sub(paste0("^", left_prefix), "", left_cols)
  ipsi <- contra <- character(0)
  for (suf in lsuf) {
    lc <- paste0(left_prefix, suf); rc <- paste0(right_prefix, suf)
    ic <- paste0("ipsi_", out_prefix, "_", suf); cc <- paste0("contra_", out_prefix, "_", suf)
    df[[ic]] <- ifelse(df[[side_col]] == "L", df[[lc]], df[[rc]])
    df[[cc]] <- ifelse(df[[side_col]] == "L", df[[rc]], df[[lc]])
    ipsi <- c(ipsi, ic); contra <- c(contra, cc)
  }
  list(df = df, ipsi_cols = ipsi, contra_cols = contra)
}

# Vectorised within-stratum permutation of residuals (B columns at once).
perm_within_blocks <- function(res0, blocks, B) {
  n <- length(res0)
  out <- matrix(res0, nrow = n, ncol = B)
  for (ii in blocks) {
    m <- length(ii)
    if (m > 1) {
      ord <- apply(matrix(runif(m * B), nrow = m), 2, order)
      out[ii, ] <- matrix(res0[ii][ord], nrow = m)
    }
  }
  out
}

# Freedman-Lane test for a single coefficient (identical to the original script).
freedman_lane_term <- function(y, X0, X1, term_index, strata = NULL, B = N_PERM, seed = 1L) {
  n <- length(y)
  if (n < MIN_N || qr(X0)$rank < ncol(X0) || qr(X1)$rank < ncol(X1))
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p_perm = NA_real_))
  fit0 <- lm.fit(X0, y); fit1 <- lm.fit(X1, y)
  yhat0 <- as.numeric(X0 %*% fit0$coefficients); res0 <- y - yhat0
  XtX_inv <- solve(crossprod(X1)); df1 <- n - ncol(X1)
  res1 <- y - as.numeric(X1 %*% fit1$coefficients); s2 <- sum(res1^2) / df1
  se_obs <- sqrt(s2 * XtX_inv[term_index, term_index])
  t_obs <- fit1$coefficients[term_index] / se_obs
  if (is.null(strata)) strata <- rep("all", n)
  blocks <- split(seq_len(n), as.character(strata))
  set.seed(seed)
  res_perm <- perm_within_blocks(res0, blocks, B)
  Yp <- yhat0 + res_perm
  coef_p <- XtX_inv %*% crossprod(X1, Yp)
  resid_p <- Yp - X1 %*% coef_p
  s2_p <- colSums(resid_p^2) / df1
  t_p <- coef_p[term_index, ] / sqrt(s2_p * XtX_inv[term_index, term_index])
  p_perm <- (1 + sum(abs(t_p[is.finite(t_p)]) >= abs(t_obs))) / (1 + sum(is.finite(t_p)))
  c(beta = unname(fit1$coefficients[term_index]), se = unname(se_obs),
    t = unname(t_obs), p_perm = unname(p_perm))
}

# Freedman-Lane test for the incremental R-squared of a nested model comparison.
freedman_lane_dR2 <- function(y, X0, X1, strata = NULL, B = N_PERM, seed = 1L) {
  n <- length(y)
  if (n < MIN_N || qr(X0)$rank < ncol(X0) || qr(X1)$rank < ncol(X1))
    return(c(dR2 = NA_real_, F = NA_real_, df1 = NA_real_, df2 = NA_real_, p_perm = NA_real_))
  rss <- function(Y, X) {
    B_ <- solve(crossprod(X), crossprod(X, Y)); colSums((Y - X %*% B_)^2)
  }
  Y <- matrix(y, ncol = 1)
  tss <- sum((y - mean(y))^2)
  rss0 <- rss(Y, X0); rss1 <- rss(Y, X1)
  dR2_obs <- as.numeric((rss0 - rss1) / tss)
  k <- ncol(X1) - ncol(X0); df2 <- n - ncol(X1)
  F_obs <- as.numeric(((rss0 - rss1) / k) / (rss1 / df2))
  fit0 <- lm.fit(X0, y); yhat0 <- as.numeric(X0 %*% fit0$coefficients); res0 <- y - yhat0
  if (is.null(strata)) strata <- rep("all", n)
  blocks <- split(seq_len(n), as.character(strata))
  set.seed(seed)
  res_perm <- perm_within_blocks(res0, blocks, B)
  Yp <- yhat0 + res_perm
  r0 <- rss(Yp, X0); r1 <- rss(Yp, X1)
  F_p <- ((r0 - r1) / k) / (r1 / df2)
  p_perm <- (1 + sum(F_p[is.finite(F_p)] >= F_obs)) / (1 + sum(is.finite(F_p)))
  c(dR2 = dR2_obs, F = F_obs, df1 = k, df2 = df2, p_perm = p_perm)
}

# Nonparametric bootstrap CI for incremental R-squared.
boot_dR2_ci <- function(d, y_var, vars0, vars1, B = N_BOOT, seed = 1L) {
  set.seed(seed)
  n <- nrow(d)
  y <- as.numeric(d[[y_var]])
  mm <- function(vs) model.matrix(as.formula(paste0("~ ", paste(sprintf("`%s`", vs), collapse = " + "))), data = d)
  X0 <- mm(vars0); X1 <- mm(vars1)
  r2q <- function(yy, XX) {
    q <- qr(XX)
    if (q$rank < ncol(XX)) return(NA_real_)
    tss <- sum((yy - mean(yy))^2)
    if (!is.finite(tss) || tss == 0) return(NA_real_)
    1 - sum((yy - qr.fitted(q, yy))^2) / tss
  }
  out <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    a <- r2q(y[idx], X1[idx, , drop = FALSE])
    c0 <- r2q(y[idx], X0[idx, , drop = FALSE])
    out[b] <- a - c0
  }
  q <- quantile(out, c(.025, .975), na.rm = TRUE)
  c(lo = unname(q[1]), hi = unname(q[2]), n_ok = sum(is.finite(out)))
}

# =============================================================================
# 1. Data preparation
# =============================================================================
log_msg("Loading data")
pat_raw <- read_excel(file.path("code and data", "df_raw_clean.xlsx"))
hc_raw  <- read_excel(file.path("code and data", "df_hc_clean.xlsx"))

pat <- pat_raw %>% mutate(
  side_lr      = normalize_side_lr(.data[["side"]]),
  Onset        = factor(ifelse(side_lr == "L", "Left-onset", "Right-onset"),
                        levels = c("Left-onset", "Right-onset")),
  sex          = factor(as.character(sex)),
  Age          = as.numeric(Age),
  Duration_Dis = as.numeric(Duration_Dis),
  Edu_year     = as.numeric(Edu_year),
  eTIV         = as.numeric(eTIV),
  HS_status    = ifelse(suppressWarnings(as.numeric(Pathology_Group)) == 1, "HS", "non-HS"))
hc <- hc_raw %>% mutate(sex = factor(as.character(sex)),
                        Age = as.numeric(Age), eTIV = as.numeric(eTIV))
stopifnot(!any(is.na(pat$side_lr)))
log_msg("Patients: left-onset=", sum(pat$side_lr == "L"),
        ", right-onset=", sum(pat$side_lr == "R"), "; HC=", nrow(hc))
log_msg("HS status: HS=", sum(pat$HS_status == "HS"), ", non-HS=", sum(pat$HS_status == "non-HS"))

## ---- whole hippocampal and summed thalamic volumes (additive-model predictors)
thal_suffix <- c("AV","VA","VAmc","VLa","VLp","VPL","VM","CL","CeM","CM","Pf",
                 "MDm","MDl","LD","LP","MGN","LGN","L_Sg","MV_Re","Pc","Pt","Pu_Total")
pat$Thal_Left_SummedVolume  <- row_sum_complete(pat, paste0("Left_",  thal_suffix))
pat$Thal_Right_SummedVolume <- row_sum_complete(pat, paste0("Right_", thal_suffix))
pat$Hippo_Left_Volume  <- as.numeric(pat[["Left-Subcort-Hippocampus"]])
pat$Hippo_Right_Volume <- as.numeric(pat[["Right-Subcort-Hippocampus"]])

pat <- pat %>% mutate(
  Hippo_Ipsi_Volume        = ifelse(side_lr == "L", Hippo_Left_Volume,  Hippo_Right_Volume),
  Hippo_Contra_Volume      = ifelse(side_lr == "L", Hippo_Right_Volume, Hippo_Left_Volume),
  Thal_Ipsi_SummedVolume   = ifelse(side_lr == "L", Thal_Left_SummedVolume,  Thal_Right_SummedVolume),
  Thal_Contra_SummedVolume = ifelse(side_lr == "L", Thal_Right_SummedVolume, Thal_Left_SummedVolume),
  Thal_Global_SummedVolume = rowMeans(cbind(Thal_Left_SummedVolume, Thal_Right_SummedVolume), na.rm = FALSE))

## ---- cognitive outcome definitions
scale_meta <- tribble(
  ~Scale, ~Scale_Display, ~Is_Normed, ~Domain,
  "WAIS-\u603b\u667a\u5546","WAIS FSIQ",TRUE,"Global/other",
  "\u8a00\u8bed\u7406\u89e3\u6307\u6570VCI","WAIS VCI",TRUE,"Verbal",
  "\u77e5\u89c9\u63a8\u7406\u6307\u6570PRI","WAIS PRI",TRUE,"Visuospatial",
  "\u5de5\u4f5c\u8bb0\u5fc6\u6307\u6570WMI","WAIS WMI",TRUE,"Global/other",
  "\u52a0\u5de5\u901f\u5ea6\u6307\u6570PSI","WAIS PSI",TRUE,"Global/other",
  "\u603b\u8bb0\u5fc6\u5546","WMS FSMQ",TRUE,"Global/other",
  "\u542c\u89c9\u8bb0\u5fc6\u6307\u6570","WMS AMI",TRUE,"Verbal",
  "\u89c6\u89c9\u8bb0\u5fc6\u6307\u6570","WMS VMI",TRUE,"Visuospatial",
  "\u5373\u523b\u8bb0\u5fc6\u6307\u6570","WMS IMI",TRUE,"Global/other",
  "\u5ef6\u8fdf\u8bb0\u5fc6\u6307\u6570","WMS DMI",TRUE,"Global/other",
  "BNT-\u547d\u540d","BNT Naming",FALSE,"Verbal",
  "VA-immediate","VA Immediate",FALSE,"Verbal",
  "VA-recall","VA Delay",FALSE,"Verbal",
  "VTrials I-V Total","VTrials I-V Total",FALSE,"Verbal",
  "CAVLT-slope","CAVLT Slope",FALSE,"Verbal",
  "FA-immediate","FA Immediate",FALSE,"Visuospatial",
  "FA-recall","FA Delay",FALSE,"Visuospatial",
  "FTrials I-V Total","FTrials I-V Total",FALSE,"Visuospatial",
  "AFLT-slope","AFLT Slope",FALSE,"Visuospatial")

# Outcome family 1: the age-standardised Wechsler indices used in the original
# additive models (n = 43-44; low power once split by onset side).
OUTCOMES_NORMED <- c("\u5de5\u4f5c\u8bb0\u5fc6\u6307\u6570WMI","\u603b\u8bb0\u5fc6\u5546","\u542c\u89c9\u8bb0\u5fc6\u6307\u6570",
                     "\u89c6\u89c9\u8bb0\u5fc6\u6307\u6570","\u5373\u523b\u8bb0\u5fc6\u6307\u6570","\u5ef6\u8fdf\u8bb0\u5fc6\u6307\u6570")
# Outcome family 2: raw memory/naming scores available in all 71 patients, so
# the onset-side strata retain 38/33 patients. Age is added as a covariate
# because these are raw rather than age-standardised scores.
OUTCOMES_RAW <- c("BNT-\u547d\u540d","VTrials I-V Total","VA-recall",
                  "FTrials I-V Total","FA-recall")

COVARS_NORMED <- c("Edu_year","Duration_Dis","sex","eTIV")
COVARS_RAW    <- c("Edu_year","Duration_Dis","sex","eTIV","Age")

model_sets <- tribble(
  ~Set, ~Set_Display, ~Hippo, ~Thal,
  "Primary_Ipsi_Global","Ipsilateral hippocampus + global thalamus",
     "Hippo_Ipsi_Volume","Thal_Global_SummedVolume",
  "Primary_Contra_Global","Contralateral hippocampus + global thalamus",
     "Hippo_Contra_Volume","Thal_Global_SummedVolume",
  "Sensitivity_Ipsi_Ipsi","Ipsilateral hippocampus + ipsilateral thalamus",
     "Hippo_Ipsi_Volume","Thal_Ipsi_SummedVolume",
  "Sensitivity_Contra_Ipsi","Contralateral hippocampus + ipsilateral thalamus",
     "Hippo_Contra_Volume","Thal_Ipsi_SummedVolume")

# =============================================================================
# A. Additive models, stratified by seizure-onset side
# =============================================================================
log_msg("A. Additive models stratified by onset side")

run_additive_one <- function(d, y_var, covars, hippo_var, thal_var, seed_base) {
  vars <- unique(c(y_var, covars, hippo_var, thal_var))
  d <- d[, vars, drop = FALSE] %>% drop_na()
  cv <- effective_covars(d, covars)
  n <- nrow(d)
  na_row <- tibble(N = n, dR2_H_alone = NA_real_, dR2_T_alone = NA_real_,
                   dR2_H_beyond_T = NA_real_, dR2_T_beyond_H = NA_real_,
                   p_H_alone = NA_real_, p_T_alone = NA_real_,
                   p_H_beyond_T = NA_real_, p_T_beyond_H = NA_real_,
                   R2_M0 = NA_real_, R2_M1 = NA_real_, R2_M2 = NA_real_, R2_M3 = NA_real_,
                   df2_M3 = NA_real_,
                   CI_H_beyond_T_lo = NA_real_, CI_H_beyond_T_hi = NA_real_,
                   CI_H_alone_lo = NA_real_, CI_H_alone_hi = NA_real_,
                   CI_T_beyond_H_lo = NA_real_, CI_T_beyond_H_hi = NA_real_,
                   CI_T_alone_lo = NA_real_, CI_T_alone_hi = NA_real_)
  if (n < MIN_N) return(na_row)
  y <- as.numeric(d[[y_var]])
  mm <- function(vs) model.matrix(as.formula(paste0("~ ", paste(sprintf("`%s`", vs), collapse = " + "))), data = d)
  X0 <- mm(cv); X1 <- mm(c(cv, hippo_var)); X2 <- mm(c(cv, thal_var)); X3 <- mm(c(cv, hippo_var, thal_var))
  r2 <- function(X) { f <- lm.fit(X, y); 1 - sum((y - X %*% f$coefficients)^2) / sum((y - mean(y))^2) }
  t_H_alone   <- freedman_lane_dR2(y, X0, X1, strata = d$sex, seed = seed_base + 1L)
  t_T_alone   <- freedman_lane_dR2(y, X0, X2, strata = d$sex, seed = seed_base + 2L)
  t_H_beyondT <- freedman_lane_dR2(y, X2, X3, strata = d$sex, seed = seed_base + 3L)
  t_T_beyondH <- freedman_lane_dR2(y, X1, X3, strata = d$sex, seed = seed_base + 4L)
  ci_hbt <- boot_dR2_ci(d, y_var, c(cv, thal_var),  c(cv, hippo_var, thal_var), seed = seed_base + 5L)
  ci_ha  <- boot_dR2_ci(d, y_var, cv,               c(cv, hippo_var),            seed = seed_base + 6L)
  ci_tbh <- boot_dR2_ci(d, y_var, c(cv, hippo_var), c(cv, hippo_var, thal_var), seed = seed_base + 7L)
  ci_ta  <- boot_dR2_ci(d, y_var, cv,               c(cv, thal_var),             seed = seed_base + 8L)
  tibble(N = n,
         dR2_H_alone = t_H_alone[["dR2"]], dR2_T_alone = t_T_alone[["dR2"]],
         dR2_H_beyond_T = t_H_beyondT[["dR2"]], dR2_T_beyond_H = t_T_beyondH[["dR2"]],
         p_H_alone = t_H_alone[["p_perm"]], p_T_alone = t_T_alone[["p_perm"]],
         p_H_beyond_T = t_H_beyondT[["p_perm"]], p_T_beyond_H = t_T_beyondH[["p_perm"]],
         R2_M0 = r2(X0), R2_M1 = r2(X1), R2_M2 = r2(X2), R2_M3 = r2(X3),
         df2_M3 = t_H_beyondT[["df2"]],
         CI_H_beyond_T_lo = ci_hbt[["lo"]], CI_H_beyond_T_hi = ci_hbt[["hi"]],
         CI_H_alone_lo = ci_ha[["lo"]], CI_H_alone_hi = ci_ha[["hi"]],
         CI_T_beyond_H_lo = ci_tbh[["lo"]], CI_T_beyond_H_hi = ci_tbh[["hi"]],
         CI_T_alone_lo = ci_ta[["lo"]], CI_T_alone_hi = ci_ta[["hi"]])
}

strata_def <- list(`All patients` = NULL, `Left-onset` = "Left-onset", `Right-onset` = "Right-onset")
fam_def <- list(
  list(name = "Age-standardised Wechsler indices", outcomes = OUTCOMES_NORMED, covars = COVARS_NORMED),
  list(name = "Raw memory and naming scores",      outcomes = OUTCOMES_RAW,    covars = COVARS_RAW))

additive_rows <- list(); ctr <- 0L
for (fam in fam_def) for (si in seq_along(strata_def)) for (mi in seq_len(nrow(model_sets)))
  for (y_var in fam$outcomes) {
    ctr <- ctr + 1L
    sname <- names(strata_def)[si]; sval <- strata_def[[si]]
    d <- if (is.null(sval)) pat else pat %>% filter(Onset == sval)
    res <- run_additive_one(d, y_var, fam$covars, model_sets$Hippo[mi], model_sets$Thal[mi],
                            seed_base = 500000L + ctr * 10L)
    additive_rows[[length(additive_rows) + 1L]] <- bind_cols(
      tibble(Outcome_Family = fam$name, Stratum = sname,
             Set = model_sets$Set[mi], Set_Display = model_sets$Set_Display[mi],
             Scale = y_var,
             Scale_Display = scale_meta$Scale_Display[match(y_var, scale_meta$Scale)],
             Domain = scale_meta$Domain[match(y_var, scale_meta$Scale)]), res)
  }
additive_raw <- bind_rows(additive_rows)
log_msg("  additive model rows: ", nrow(additive_raw))

additive_res <- additive_raw %>%
  pivot_longer(c(dR2_H_alone, dR2_T_alone, dR2_H_beyond_T, dR2_T_beyond_H),
               names_to = "Contrast", values_to = "Delta_R2") %>%
  mutate(Contrast = recode(Contrast,
           dR2_H_alone = "Hippocampus alone (M1-M0)", dR2_T_alone = "Thalamus alone (M2-M0)",
           dR2_H_beyond_T = "Hippocampus beyond thalamus (M3-M2)",
           dR2_T_beyond_H = "Thalamus beyond hippocampus (M3-M1)"),
         P_perm = case_when(
           Contrast == "Hippocampus alone (M1-M0)" ~ p_H_alone,
           Contrast == "Thalamus alone (M2-M0)" ~ p_T_alone,
           Contrast == "Hippocampus beyond thalamus (M3-M2)" ~ p_H_beyond_T,
           TRUE ~ p_T_beyond_H),
         CI_lo = case_when(Contrast == "Hippocampus beyond thalamus (M3-M2)" ~ CI_H_beyond_T_lo,
                           Contrast == "Hippocampus alone (M1-M0)" ~ CI_H_alone_lo,
                           Contrast == "Thalamus beyond hippocampus (M3-M1)" ~ CI_T_beyond_H_lo,
                           TRUE ~ CI_T_alone_lo),
         CI_hi = case_when(Contrast == "Hippocampus beyond thalamus (M3-M2)" ~ CI_H_beyond_T_hi,
                           Contrast == "Hippocampus alone (M1-M0)" ~ CI_H_alone_hi,
                           Contrast == "Thalamus beyond hippocampus (M3-M1)" ~ CI_T_beyond_H_hi,
                           TRUE ~ CI_T_alone_hi)) %>%
  group_by(Outcome_Family, Stratum, Set, Contrast) %>%
  mutate(q_BH = p.adjust(P_perm, method = "BH")) %>% ungroup() %>%
  mutate(Signif = star_from_q(q_BH)) %>%
  select(Outcome_Family, Stratum, Set, Set_Display, Scale, Scale_Display, Domain,
         N, df2_M3, Contrast, Delta_R2, CI_lo, CI_hi, P_perm, q_BH, Signif,
         R2_M0, R2_M1, R2_M2, R2_M3)
write_table(additive_res, "TableL_A1_Additive_models_by_onset_side.csv")

# =============================================================================
# E. Supporting predictor x onset-side interaction tests (full sample)
# =============================================================================
# Stratified models cannot by themselves establish that the two onset-side
# strata differ; a formal interaction term is required. These full-sample tests
# retain n = 43-71 and are reported alongside the stratified estimates.
log_msg("E. Predictor x onset-side interaction tests (full sample)")

run_interaction_one <- function(y_var, covars, pred_var, seed) {
  vars <- unique(c(y_var, covars, pred_var, "Onset"))
  d <- pat[, vars, drop = FALSE] %>% drop_na()
  cv <- effective_covars(d, covars)
  if (nrow(d) < MIN_N) return(tibble(N = nrow(d), Beta_Int = NA_real_, SE = NA_real_,
                                     t = NA_real_, P_perm = NA_real_))
  d$Pz <- as.numeric(scale(d[[pred_var]])); d$Yz <- as.numeric(scale(d[[y_var]]))
  rhs0 <- paste(c(sprintf("`%s`", cv), "Pz", "Onset"), collapse = " + ")
  X0 <- model.matrix(as.formula(paste0("~ ", rhs0)), data = d)
  X1 <- model.matrix(as.formula(paste0("~ ", rhs0, " + Pz:Onset")), data = d)
  idx <- grep(":", colnames(X1))
  if (length(idx) != 1) return(tibble(N = nrow(d), Beta_Int = NA_real_, SE = NA_real_,
                                      t = NA_real_, P_perm = NA_real_))
  tst <- freedman_lane_term(d$Yz, X0, X1, idx, strata = d$sex, seed = seed)
  tibble(N = nrow(d), Beta_Int = tst[["beta"]], SE = tst[["se"]],
         t = tst[["t"]], P_perm = tst[["p_perm"]])
}

int_grid <- bind_rows(
  expand_grid(Outcome_Family = "Age-standardised Wechsler indices",
              Scale = OUTCOMES_NORMED,
              Predictor = c("Hippo_Ipsi_Volume","Thal_Global_SummedVolume")),
  expand_grid(Outcome_Family = "Raw memory and naming scores",
              Scale = OUTCOMES_RAW,
              Predictor = c("Hippo_Ipsi_Volume","Thal_Global_SummedVolume")))

interaction_res <- map_dfr(seq_len(nrow(int_grid)), function(i) {
  fam <- int_grid$Outcome_Family[i]
  cv <- if (fam == "Raw memory and naming scores") COVARS_RAW else COVARS_NORMED
  bind_cols(int_grid[i, ],
            run_interaction_one(int_grid$Scale[i], cv, int_grid$Predictor[i],
                                seed = 700000L + i * 7L))
}) %>%
  mutate(Scale_Display = scale_meta$Scale_Display[match(Scale, scale_meta$Scale)],
         Predictor_Display = recode(Predictor,
           Hippo_Ipsi_Volume = "Ipsilateral hippocampal volume",
           Thal_Global_SummedVolume = "Global thalamic summed volume"),
         CI_lo = Beta_Int - 1.96 * SE, CI_hi = Beta_Int + 1.96 * SE) %>%
  group_by(Outcome_Family, Predictor) %>%
  mutate(q_BH = p.adjust(P_perm, method = "BH")) %>% ungroup() %>%
  mutate(Signif = star_from_q(q_BH)) %>%
  select(Outcome_Family, Scale, Scale_Display, Predictor, Predictor_Display,
         N, Beta_Int, SE, CI_lo, CI_hi, t, P_perm, q_BH, Signif)
write_table(interaction_res, "TableL_A2_Predictor_by_onset_side_interactions.csv")

# =============================================================================
# B. HC-referenced ipsilateral-versus-contralateral atrophy z scores, by onset side
# =============================================================================
log_msg("B. Ipsilateral-versus-contralateral atrophy z scores by onset side")

hip_fine <- c("Hippocampal_tail","subiculum_comb","CA1_comb","CA3_comb","CA4_comb",
              "GC_ML_DG_comb","molecular_layer_HP_comb","HATA","fimbria",
              "presubiculum_comb","parasubiculum")
thal_fine <- thal_suffix
pair_tbl <- bind_rows(
  tibble(Family = "Hippocampal subfields", L = paste0("lh_", hip_fine),
         R = paste0("rh_", hip_fine), Region = hip_fine),
  tibble(Family = "Thalamic nuclei", L = paste0("Left_", thal_fine),
         R = paste0("Right_", thal_fine), Region = thal_fine)) %>%
  filter(L %in% names(pat), R %in% names(pat), L %in% names(hc), R %in% names(hc))

# HC-referenced atrophy z: standardise each ROI against a control-based model
# adjusted for age, sex and eTIV, then sign-reverse so that positive = atrophy.
atrophy_z <- function(roi) {
  dh <- hc[, c("Age","sex","eTIV",roi)] %>% drop_na()
  if (nrow(dh) < 8) return(rep(NA_real_, nrow(pat)))
  fit <- lm(as.formula(paste0("`", roi, "` ~ Age + sex + eTIV")), data = dh)
  s <- summary(fit)$sigma
  dp <- pat[, c("Age","sex","eTIV",roi)]
  pred <- suppressWarnings(predict(fit, newdata = dp))
  -1 * (as.numeric(dp[[roi]]) - pred) / s
}
for (i in seq_len(nrow(pair_tbl))) {
  pat[[paste0("z_", pair_tbl$L[i])]] <- atrophy_z(pair_tbl$L[i])
  pat[[paste0("z_", pair_tbl$R[i])]] <- atrophy_z(pair_tbl$R[i])
}

ipsi_contra_z <- map_dfr(seq_len(nrow(pair_tbl)), function(i) {
  zl <- pat[[paste0("z_", pair_tbl$L[i])]]; zr <- pat[[paste0("z_", pair_tbl$R[i])]]
  zi <- ifelse(pat$side_lr == "L", zl, zr); zc <- ifelse(pat$side_lr == "L", zr, zl)
  map_dfr(c("Left-onset","Right-onset"), function(g) {
    k <- pat$Onset == g & is.finite(zi) & is.finite(zc)
    a <- zi[k]; b <- zc[k]
    if (sum(k) < 8) return(tibble())
    tt <- t.test(a, b, paired = TRUE)
    dz <- mean(a - b) / sd(a - b)
    tibble(Stratum = g, Family = pair_tbl$Family[i],
           Region = str_replace_all(str_remove(pair_tbl$Region[i], "_comb$"), "_", " "),
           N = sum(k), Mean_Ipsi_z = mean(a), Mean_Contra_z = mean(b),
           Mean_Diff_z = mean(a - b), SE_Diff = sd(a - b) / sqrt(sum(k)),
           Cohen_dz = dz, t = unname(tt$statistic), df = unname(tt$parameter),
           p_value = tt$p.value)
  })
}) %>%
  group_by(Stratum, Family) %>% mutate(q_BH = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>% mutate(Signif = star_from_q(q_BH))
write_table(ipsi_contra_z, "TableL_B1_Ipsi_vs_contra_atrophy_z_by_onset_side.csv")
log_msg("  ipsi/contra rows: ", nrow(ipsi_contra_z),
        "; corrected: ", sum(ipsi_contra_z$q_BH < 0.05, na.rm = TRUE))

## ---- B2. Is the ipsilateral-contralateral asymmetry itself different by onset
## side once HS composition is accounted for?  Right-onset patients are more
## often HS (22/33) than left-onset patients (17/38), and HS status is the main
## driver of ipsilateral hippocampal atrophy, so the unadjusted B1 comparison is
## confounded.  Here the within-patient asymmetry (ipsi minus contra atrophy z)
## is modelled as a function of onset side with HS status held constant.
log_msg("B2. Onset-side difference in ipsi-contra asymmetry, adjusted for HS status")
hs_by_onset <- pat %>% count(Onset, HS_status) %>% pivot_wider(names_from = HS_status, values_from = n, values_fill = 0)
write_table(hs_by_onset, "TableL_B0_HS_composition_by_onset_side.csv")
log_msg("  HS composition by onset side: ", paste(capture.output(print(hs_by_onset)), collapse = " | "))

asym_res <- map_dfr(seq_len(nrow(pair_tbl)), function(i) {
  zl <- pat[[paste0("z_", pair_tbl$L[i])]]; zr <- pat[[paste0("z_", pair_tbl$R[i])]]
  zi <- ifelse(pat$side_lr == "L", zl, zr); zc <- ifelse(pat$side_lr == "L", zr, zl)
  d <- tibble(Asym = zi - zc, Onset = pat$Onset, HS = factor(pat$HS_status), sex = pat$sex) %>% drop_na()
  if (nrow(d) < MIN_N) return(tibble())
  X0 <- model.matrix(~ HS, data = d)
  X1 <- model.matrix(~ HS + Onset, data = d)
  idx <- grep("^Onset", colnames(X1))
  tst <- freedman_lane_term(d$Asym, X0, X1, idx, strata = d$sex, seed = 950000L + i)
  mL <- mean(d$Asym[d$Onset == "Left-onset"]); mR <- mean(d$Asym[d$Onset == "Right-onset"])
  tibble(Family = pair_tbl$Family[i],
         Region = str_replace_all(str_remove(pair_tbl$Region[i], "_comb$"), "_", " "),
         N = nrow(d), Mean_Asym_Left = mL, Mean_Asym_Right = mR,
         Unadjusted_Diff_R_minus_L = mR - mL,
         Beta_Onset_HSadj = tst[["beta"]], SE = tst[["se"]],
         t = tst[["t"]], P_perm = tst[["p_perm"]])
}) %>% group_by(Family) %>% mutate(q_BH = p.adjust(P_perm, method = "BH")) %>%
  ungroup() %>% mutate(Signif = star_from_q(q_BH))
write_table(asym_res, "TableL_B2_Asymmetry_onset_side_effect_HSadjusted.csv")
log_msg("  asymmetry ROIs with onset-side effect surviving FDR after HS adjustment: ",
        sum(asym_res$q_BH < 0.05, na.rm = TRUE), " / ", nrow(asym_res))

# =============================================================================
# C. Cognitive performance, left-onset versus right-onset
# =============================================================================
log_msg("C. Cognitive performance left-onset versus right-onset")
cog_lr <- map_dfr(seq_len(nrow(scale_meta)), function(i) {
  v <- scale_meta$Scale[i]
  d <- pat[, c("Onset", v)] %>% drop_na()
  x <- as.numeric(d[[v]][d$Onset == "Left-onset"]); y <- as.numeric(d[[v]][d$Onset == "Right-onset"])
  if (length(x) < 5 || length(y) < 5) return(tibble())
  tt <- t.test(x, y, var.equal = FALSE); g <- hedges_g(x, y)
  se <- hedges_g_se(g, length(x), length(y))
  tibble(Scale = v, Scale_Display = scale_meta$Scale_Display[i], Domain = scale_meta$Domain[i],
         N_Left = length(x), Mean_Left = mean(x), SD_Left = sd(x),
         N_Right = length(y), Mean_Right = mean(y), SD_Right = sd(y),
         Hedges_g = g, CI_lo = g - 1.96 * se, CI_hi = g + 1.96 * se,
         Welch_t = unname(tt$statistic), df = unname(tt$parameter), p_value = tt$p.value)
}) %>% mutate(q_BH = p.adjust(p_value, method = "BH"), Signif = star_from_q(q_BH))
write_table(cog_lr, "TableL_C1_Cognition_left_vs_right_onset.csv")
log_msg("  cognitive scales tested: ", nrow(cog_lr),
        "; surviving FDR: ", sum(cog_lr$q_BH < 0.05, na.rm = TRUE))

# =============================================================================
# D. Fine-grained structural abnormality mapping by onset side versus HC
# =============================================================================
log_msg("D. Structural abnormality mapping by onset side versus HC")
struct_tbl <- bind_rows(
  tibble(Family = "Hippocampal subfields",
         ROI = c(paste0("lh_", hip_fine), paste0("rh_", hip_fine))),
  tibble(Family = "Thalamic nuclei",
         ROI = c(paste0("Left_", thal_fine), paste0("Right_", thal_fine)))) %>%
  filter(ROI %in% names(pat), ROI %in% names(hc)) %>%
  mutate(Hemisphere = ifelse(str_detect(ROI, "^(lh_|Left_)"), "Left", "Right"),
         Region = str_replace_all(str_remove(str_remove(ROI, "^(lh_|rh_|Left_|Right_)"), "_comb$"), "_", " "),
         ROI_Display = paste(Hemisphere, Region))

struct_all <- bind_rows(
  hc %>% mutate(Group = "HC"),
  pat %>% mutate(Group = as.character(Onset))) %>%
  mutate(Group = factor(Group, levels = c("HC","Left-onset","Right-onset")))

struct_res <- map_dfr(c("Left-onset","Right-onset"), function(g) {
  map_dfr(seq_len(nrow(struct_tbl)), function(i) {
    roi <- struct_tbl$ROI[i]
    d <- struct_all %>% filter(Group %in% c("HC", g)) %>%
      select(Group, Age, sex, eTIV, all_of(roi)) %>% drop_na()
    d$Group <- relevel(droplevels(d$Group), ref = "HC")
    X0 <- model.matrix(~ eTIV + Age + sex, data = d)
    X1 <- model.matrix(~ eTIV + Age + sex + Group, data = d)
    idx <- grep("^Group", colnames(X1))
    tst <- freedman_lane_term(as.numeric(d[[roi]]), X0, X1, idx, strata = d$sex,
                              seed = 900000L + which(c("Left-onset","Right-onset") == g) * 1000L + i)
    yres <- residualize_one(d[[roi]], d, c("eTIV","Age","sex"))
    tibble(Contrast = paste0(g, " vs HC"), Family = struct_tbl$Family[i],
           ROI = roi, ROI_Display = struct_tbl$ROI_Display[i],
           Hemisphere = struct_tbl$Hemisphere[i], N_Patient = sum(d$Group == g),
           N_HC = sum(d$Group == "HC"), Beta = tst[["beta"]], SE = tst[["se"]],
           t = tst[["t"]], P_perm = tst[["p_perm"]],
           Hedges_g = hedges_g(yres[d$Group == g], yres[d$Group == "HC"]))
  })
}) %>% group_by(Contrast, Family) %>% mutate(q_BH = p.adjust(P_perm, method = "BH")) %>%
  ungroup() %>% mutate(Signif = star_from_q(q_BH))
write_table(struct_res, "TableL_D1_Structural_mapping_by_onset_side_vs_HC.csv")
log_msg("  structural rows: ", nrow(struct_res))

# =============================================================================
# Summary counts
# =============================================================================
summary_counts <- bind_rows(
  additive_res %>% filter(!is.na(q_BH)) %>%
    group_by(Analysis = "Additive models", Outcome_Family, Stratum, Set = Set_Display, Contrast) %>%
    summarise(N_tests = n(), N_q05 = sum(q_BH < 0.05), N_q10 = sum(q_BH < 0.10),
              Min_q = min(q_BH), .groups = "drop"),
  ipsi_contra_z %>% group_by(Analysis = "Ipsi vs contra atrophy z", Outcome_Family = NA_character_,
                             Stratum, Set = Family, Contrast = "Ipsi - contra") %>%
    summarise(N_tests = n(), N_q05 = sum(q_BH < 0.05), N_q10 = sum(q_BH < 0.10),
              Min_q = min(q_BH), .groups = "drop"),
  cog_lr %>% group_by(Analysis = "Cognition left vs right", Outcome_Family = NA_character_,
                      Stratum = "All patients", Set = "Cognitive battery",
                      Contrast = "Left - right") %>%
    summarise(N_tests = n(), N_q05 = sum(q_BH < 0.05), N_q10 = sum(q_BH < 0.10),
              Min_q = min(q_BH), .groups = "drop"),
  struct_res %>% group_by(Analysis = "Structural mapping vs HC", Outcome_Family = NA_character_,
                          Stratum = Contrast, Set = Family, Contrast = "Patient - HC") %>%
    summarise(N_tests = n(), N_q05 = sum(q_BH < 0.05, na.rm = TRUE),
              N_q10 = sum(q_BH < 0.10, na.rm = TRUE), Min_q = min(q_BH, na.rm = TRUE), .groups = "drop"))
write_table(summary_counts, "TableL_E1_Summary_counts_all_analyses.csv")

# =============================================================================
# Figures
# =============================================================================
log_msg("Generating figures")
PAL_ONSET <- c(`Left-onset` = "#B2182B", `Right-onset` = "#2166AC", `All patients` = "#4D4D4D")

## Figure L-A: incremental R2 by onset side (primary model set, both families)
figA_dat <- additive_res %>%
  filter(Set == "Primary_Ipsi_Global") %>%
  mutate(Stratum = factor(Stratum, levels = c("All patients","Left-onset","Right-onset")),
         Contrast = factor(Contrast, levels = c("Hippocampus alone (M1-M0)",
                                                "Hippocampus beyond thalamus (M3-M2)",
                                                "Thalamus alone (M2-M0)",
                                                "Thalamus beyond hippocampus (M3-M1)")),
         Scale_Display = factor(Scale_Display, levels = rev(unique(scale_meta$Scale_Display))))
pA <- ggplot(figA_dat, aes(x = Delta_R2, y = Scale_Display, fill = Stratum)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "grey55", linewidth = 0.35) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi), position = position_dodge(width = 0.78),
                 height = 0.22, linewidth = 0.35, colour = "grey25", na.rm = TRUE) +
  geom_text(aes(label = Signif, x = pmax(Delta_R2, 0) + 0.012),
            position = position_dodge(width = 0.78), size = 3.1, hjust = 0, colour = "black") +
  facet_grid(Outcome_Family ~ Contrast, scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 26)) +
  scale_fill_manual(values = PAL_ONSET) +
  scale_x_continuous(limits = c(-0.06, 0.45), breaks = seq(0, 0.4, 0.1)) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 9),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Additive models of hippocampal and thalamic volume, by seizure-onset side",
       subtitle = paste0("Incremental explained variance with bootstrap 95% CI; Freedman-Lane ",
                         "permutation, BH-FDR within stratum, model set and contrast"),
       x = "Incremental R\u00b2 (\u0394R\u00b2)", y = NULL, fill = "Stratum")
save_plot(pA, "Figure_LA_Additive_models_by_onset_side", 13.5, 9.5)

## Figure L-B: ipsilateral versus contralateral atrophy z by onset side
figB_dat <- ipsi_contra_z %>%
  pivot_longer(c(Mean_Ipsi_z, Mean_Contra_z), names_to = "Side", values_to = "Z") %>%
  mutate(Side = recode(Side, Mean_Ipsi_z = "Ipsilateral", Mean_Contra_z = "Contralateral"),
         Side = factor(Side, levels = c("Ipsilateral","Contralateral")))
figB_lab <- ipsi_contra_z %>% group_by(Stratum, Family, Region) %>%
  summarise(y = max(Mean_Ipsi_z, Mean_Contra_z) + 0.12, lab = Signif[1], .groups = "drop")
pB <- ggplot(figB_dat, aes(x = Region, y = Z, fill = Side)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(data = figB_lab, aes(x = Region, y = y, label = lab), inherit.aes = FALSE,
            size = 3.2, colour = "black") +
  facet_grid(Stratum ~ Family, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = c(Ipsilateral = "#B2182B", Contralateral = "#4393C3")) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, colour = "black", size = 8),
        axis.text.y = element_text(colour = "black"), legend.position = "bottom",
        strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "HC-referenced atrophy z scores, ipsilateral versus contralateral, by seizure-onset side",
       subtitle = "Positive values indicate smaller-than-expected volume (greater atrophy); paired t tests, BH-FDR within family and stratum",
       x = NULL, y = "HC-referenced atrophy z score", fill = NULL)
save_plot(pB, "Figure_LB_Ipsi_contra_atrophy_z_by_onset_side", 14, 8.5)

## Figure L-C: cognition left versus right onset
figC_dat <- cog_lr %>%
  mutate(Scale_Display = factor(Scale_Display, levels = rev(scale_meta$Scale_Display)),
         FDR = ifelse(q_BH < 0.05, "q < .05", "q ≥ .05"))
pC <- ggplot(figC_dat, aes(x = Hedges_g, y = Scale_Display)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbarh(aes(xmin = CI_lo, xmax = CI_hi, colour = Domain), height = 0.18, linewidth = 0.55) +
  geom_point(aes(colour = Domain, shape = FDR), size = 2.4) +
  scale_colour_manual(values = c(Verbal = "#B2182B", Visuospatial = "#2166AC",
                                 `Global/other` = "#555555")) +
  scale_shape_manual(values = c(`q < .05` = 17, `q ≥ .05` = 16)) +
  theme_minimal(base_size = 12) +
  theme(axis.text = element_text(colour = "black"), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Cognitive performance: left-onset versus right-onset MTLE",
       subtitle = "Hedges' g (left minus right); positive values indicate higher scores in left-onset MTLE",
       x = "Hedges' g (95% CI)", y = NULL, colour = "Domain", shape = "FDR")
save_plot(pC, "Figure_LC_Cognition_left_vs_right_onset", 9.5, 8.5)

## Figure L-D: structural abnormality mapping by onset side
figD_dat <- struct_res %>%
  mutate(Contrast = factor(Contrast, levels = c("Left-onset vs HC","Right-onset vs HC")),
         ROI_Display = factor(ROI_Display, levels = rev(sort(unique(ROI_Display)))))
pD <- ggplot(figD_dat, aes(x = Contrast, y = ROI_Display, fill = Hedges_g)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = Signif), size = 3, colour = "black") +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-1.8, 1.8), oob = scales::squish) +
  theme_minimal(base_size = 10) +
  theme(axis.text = element_text(colour = "black"), axis.text.y = element_text(size = 7),
        strip.text = element_text(face = "bold"), legend.position = "right",
        plot.title = element_text(face = "bold", size = 14)) +
  labs(title = "Fine-grained structural abnormality mapping by seizure-onset side",
       subtitle = "Covariate-adjusted Hedges' g versus healthy controls; negative (blue) = smaller volume in patients",
       x = NULL, y = NULL, fill = "Hedges' g")
save_plot(pD, "Figure_LD_Structural_mapping_by_onset_side", 8.5, 12)

# =============================================================================
# QC
# =============================================================================
qc <- tibble(
  Item = c("R version","Permutations","Bootstrap resamples","Patients","Left-onset",
           "Right-onset","Healthy controls","HS","non-HS",
           "n WAIS WMI","n WMS FSMQ","n BNT/CAVLT/AFLT"),
  Value = c(R.version.string, N_PERM, N_BOOT, nrow(pat), sum(pat$side_lr == "L"),
            sum(pat$side_lr == "R"), nrow(hc), sum(pat$HS_status == "HS"),
            sum(pat$HS_status == "non-HS"),
            sum(is.finite(as.numeric(pat[["\u5de5\u4f5c\u8bb0\u5fc6\u6307\u6570WMI"]]))),
            sum(is.finite(as.numeric(pat[["\u603b\u8bb0\u5fc6\u5546"]]))),
            sum(is.finite(as.numeric(pat[["BNT-\u547d\u540d"]])))))
write_table(qc, "TableL_QC_run_metadata.csv")
writeLines(capture.output(sessionInfo()), file.path(QC_DIR, "sessionInfo.txt"))
log_msg("Done.")

#!/usr/bin/env Rscript

# Laterality sensitivity analyses requested during peer review.
# This script extends the original project scripts while preserving their core
# statistical choices: 5,000 permutations, covariate adjustment, family-wise
# BH-FDR, and the blue-white-red figure palette.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(scales)
})

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(2026)

N_PERM <- 5000L
MIN_N <- 15L
BASE_DIR <- normalizePath(".", winslash = "/", mustWork = TRUE)
OUTDIR <- file.path("analysis_outputs", "laterality_full_sensitivity")
TABLE_DIR <- file.path(OUTDIR, "01_tables")
FIG_DIR <- file.path(OUTDIR, "02_figures")
QC_DIR <- file.path(OUTDIR, "03_qc")
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(QC_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUTDIR, "run_log.txt")
if (file.exists(LOG_FILE)) file.remove(LOG_FILE)
log_msg <- function(...) {
  z <- paste0(...)
  cat(z, "\n", sep = "")
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", z, "\n",
      sep = "", file = LOG_FILE, append = TRUE)
}

write_table <- function(x, filename) {
  readr::write_excel_csv(x, file.path(TABLE_DIR, filename), na = "")
}

save_plot <- function(p, stem, width, height) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p,
         width = width, height = height, units = "in", dpi = 600,
         bg = "white")
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p,
         width = width, height = height, units = "in",
         device = grDevices::cairo_pdf, bg = "white")
}

star_from_q <- function(q) {
  case_when(
    is.finite(q) & q < 0.001 ~ "***",
    is.finite(q) & q < 0.01  ~ "**",
    is.finite(q) & q < 0.05  ~ "*",
    is.finite(q) & q < 0.10  ~ "○",
    TRUE ~ ""
  )
}

normalize_side_lr <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    xl <- tolower(trimws(x))
    out <- rep(NA_character_, length(xl))
    out[xl %in% c("1", "l", "left", "左", "左侧")] <- "L"
    out[xl %in% c("0", "r", "right", "右", "右侧")] <- "R"
    return(out)
  }
  out <- rep(NA_character_, length(x))
  out[is.finite(x) & x == 1] <- "L"
  out[is.finite(x) & x == 0] <- "R"
  out
}

hedges_g <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  nx <- length(x)
  ny <- length(y)
  if (nx < 5 || ny < 5) return(NA_real_)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  if (!is.finite(sp) || sp < 1e-12) return(NA_real_)
  d <- (mean(x) - mean(y)) / sp
  J <- 1 - 3 / (4 * (nx + ny) - 9)
  J * d
}

hedges_g_se <- function(g, nx, ny) {
  if (!is.finite(g) || nx < 5 || ny < 5) return(NA_real_)
  sqrt((nx + ny) / (nx * ny) + g^2 / (2 * (nx + ny - 2)))
}

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
  y <- as.numeric(y)
  covars <- effective_covars(d, covars)
  if (length(covars) == 0) return(y - mean(y, na.rm = TRUE))
  X <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse = " + "))), data = d)
  fit <- lm.fit(X, y)
  y - as.numeric(X %*% fit$coefficients)
}

perm_corr <- function(x, y, B = N_PERM, seed = 1L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < MIN_N || sd(x) == 0 || sd(y) == 0) return(c(r = NA_real_, p = NA_real_))
  xz <- as.numeric(scale(x))
  yz <- as.numeric(scale(y))
  obs <- sum(xz * yz) / (n - 1)
  idx <- replicate(B, sample.int(n, n, replace = FALSE))
  perm_r <- as.numeric(crossprod(xz, matrix(yz[idx], nrow = n))) / (n - 1)
  # Kept identical to the original structure-cognition script so that the
  # reviewer package reproduces the already reported Table S6/S7 results.
  p <- mean(abs(perm_r) >= abs(obs))
  c(r = obs, p = p)
}

# Vectorized Freedman-Lane test for one model coefficient.
freedman_lane_term <- function(y, X0, X1, term_index, strata = NULL,
                               B = N_PERM, seed = 1L) {
  n <- length(y)
  if (n < MIN_N || qr(X0)$rank < ncol(X0) || qr(X1)$rank < ncol(X1)) {
    return(c(beta = NA_real_, se = NA_real_, t = NA_real_, p_perm = NA_real_))
  }
  fit0 <- lm.fit(X0, y)
  fit1 <- lm.fit(X1, y)
  yhat0 <- as.numeric(X0 %*% fit0$coefficients)
  res0 <- y - yhat0
  XtX_inv <- solve(crossprod(X1))
  df1 <- n - ncol(X1)
  res1 <- y - as.numeric(X1 %*% fit1$coefficients)
  s2 <- sum(res1^2) / df1
  se_obs <- sqrt(s2 * XtX_inv[term_index, term_index])
  t_obs <- fit1$coefficients[term_index] / se_obs

  if (is.null(strata)) strata <- rep("all", n)
  strata <- as.character(strata)
  blocks <- split(seq_len(n), strata)
  set.seed(seed)
  res_perm <- matrix(NA_real_, nrow = n, ncol = B)
  for (b in seq_len(B)) {
    rp <- res0
    for (ii in blocks) {
      if (length(ii) > 1) rp[ii] <- res0[sample(ii, length(ii), replace = FALSE)]
    }
    res_perm[, b] <- rp
  }
  Yp <- yhat0 + res_perm
  coef_p <- XtX_inv %*% crossprod(X1, Yp)
  resid_p <- Yp - X1 %*% coef_p
  s2_p <- colSums(resid_p^2) / df1
  t_p <- coef_p[term_index, ] / sqrt(s2_p * XtX_inv[term_index, term_index])
  p_perm <- (1 + sum(abs(t_p[is.finite(t_p)]) >= abs(t_obs))) /
    (1 + sum(is.finite(t_p)))
  c(beta = unname(fit1$coefficients[term_index]),
    se = unname(se_obs), t = unname(t_obs), p_perm = unname(p_perm))
}

row_sum_complete <- function(df, vars) {
  vars <- intersect(vars, names(df))
  if (length(vars) == 0) return(rep(NA_real_, nrow(df)))
  x <- as.matrix(df[, vars, drop = FALSE])
  storage.mode(x) <- "double"
  out <- rowSums(x, na.rm = FALSE)
  out[rowSums(is.finite(x)) < length(vars)] <- NA_real_
  out
}

make_ipsi_contra_pairs <- function(df, side_col, left_cols, right_cols,
                                   left_prefix, right_prefix, out_prefix) {
  lsuf <- sub(paste0("^", left_prefix), "", left_cols)
  rsuf <- sub(paste0("^", right_prefix), "", right_cols)
  if (!setequal(lsuf, rsuf)) stop("Left/right ROI suffixes differ.")
  ipsi <- contra <- character(0)
  for (suf in lsuf) {
    lc <- paste0(left_prefix, suf)
    rc <- paste0(right_prefix, suf)
    ic <- paste0("ipsi_", out_prefix, "_", suf)
    cc <- paste0("contra_", out_prefix, "_", suf)
    df[[ic]] <- ifelse(df[[side_col]] == "L", df[[lc]], df[[rc]])
    df[[cc]] <- ifelse(df[[side_col]] == "L", df[[rc]], df[[lc]])
    ipsi <- c(ipsi, ic)
    contra <- c(contra, cc)
  }
  list(df = df, ipsi_cols = ipsi, contra_cols = contra)
}

log_msg("Starting laterality sensitivity package")
log_msg("Base directory: ", BASE_DIR)
log_msg("N permutations: ", N_PERM)

pat_raw <- read_excel(file.path("code and data", "df_raw_clean.xlsx"))
hc_raw <- read_excel(file.path("code and data", "df_hc_clean.xlsx"))

scale_meta <- tibble::tribble(
  ~Scale, ~Scale_Display, ~Is_Normed, ~Domain,
  "WAIS-总智商", "WAIS FSIQ", TRUE, "Global/other",
  "言语理解指数VCI", "WAIS VCI", TRUE, "Verbal",
  "知觉推理指数PRI", "WAIS PRI", TRUE, "Visuospatial",
  "工作记忆指数WMI", "WAIS WMI", TRUE, "Global/other",
  "加工速度指数PSI", "WAIS PSI", TRUE, "Global/other",
  "总记忆商", "WMS FSMQ", TRUE, "Global/other",
  "听觉记忆指数", "WMS AMI", TRUE, "Verbal",
  "视觉记忆指数", "WMS VMI", TRUE, "Visuospatial",
  "即刻记忆指数", "WMS IMI", TRUE, "Global/other",
  "延迟记忆指数", "WMS DMI", TRUE, "Global/other",
  "BNT-命名", "BNT Naming", FALSE, "Verbal",
  "VA-immediate", "VA Immediate", FALSE, "Verbal",
  "VA-recall", "VA Delay", FALSE, "Verbal",
  "VTrials I-V Total", "VTrials I-V Total", FALSE, "Verbal",
  "CAVLT-slope", "CAVLT Slope", FALSE, "Verbal",
  "FA-immediate", "FA Immediate", FALSE, "Visuospatial",
  "FA-recall", "FA Delay", FALSE, "Visuospatial",
  "FTrials I-V Total", "FTrials I-V Total", FALSE, "Visuospatial",
  "AFLT-slope", "AFLT Slope", FALSE, "Visuospatial"
)

missing_scales <- setdiff(scale_meta$Scale, names(pat_raw))
if (length(missing_scales) > 0) {
  stop("Missing cognitive scale columns: ", paste(missing_scales, collapse = ", "))
}

pat <- pat_raw %>%
  mutate(
    side_lr = normalize_side_lr(.data[["side"]]),
    side_lr = factor(side_lr, levels = c("R", "L")),
    sex = factor(as.character(sex)),
    Age = as.numeric(Age),
    Duration_Dis = as.numeric(Duration_Dis),
    Edu_year = as.numeric(Edu_year),
    eTIV = as.numeric(eTIV),
    HS_status = ifelse(suppressWarnings(as.numeric(Pathology_Group)) == 1,
                       "HS", "non-HS")
  )
hc <- hc_raw %>%
  mutate(
    sex = factor(as.character(sex)),
    Age = as.numeric(Age),
    eTIV = as.numeric(eTIV)
  )

if (any(is.na(pat$side_lr))) stop("Patient side could not be fully coded as L/R.")
log_msg("Patient counts: Left-onset=", sum(pat$side_lr == "L"),
        ", Right-onset=", sum(pat$side_lr == "R"),
        ", HC=", nrow(hc))

# Low-dimensional native ROIs used for cognition analyses.
hip_segment_suffix <- c("Whole_hippocampal_head", "Whole_hippocampal_body", "Hippocampal_tail")
hip_segment_display <- c(
  Whole_hippocampal_head = "Hippocampal head",
  Whole_hippocampal_body = "Hippocampal body",
  Hippocampal_tail = "Hippocampal tail"
)
thal_group_def <- list(
  AnteriorDorsalGroup = c("AV", "LD"),
  MediodorsalNuclei = c("MDm", "MDl"),
  IntralaminarNuclei = c("CL", "CM", "CeM", "Pf", "Pc"),
  PulvinarLPComplex = c("Pu_Total", "L_Sg", "LP"),
  VentralNuclei = c("VA", "VAmc", "VLa", "VLp", "VPL", "VM"),
  MidlineNuclei = c("MV_Re", "Pt"),
  GeniculateNuclei = c("LGN", "MGN")
)
thal_group_display <- c(
  AnteriorDorsalGroup = "Anterior-dorsal group",
  MediodorsalNuclei = "Mediodorsal nuclei",
  IntralaminarNuclei = "Intralaminar nuclei",
  PulvinarLPComplex = "Pulvinar-LP complex",
  VentralNuclei = "Ventral nuclei",
  MidlineNuclei = "Midline nuclei",
  GeniculateNuclei = "Geniculate nuclei"
)

for (g in names(thal_group_def)) {
  pat[[paste0("LeftGrp_", g)]] <- row_sum_complete(pat, paste0("Left_", thal_group_def[[g]]))
  pat[[paste0("RightGrp_", g)]] <- row_sum_complete(pat, paste0("Right_", thal_group_def[[g]]))
}

hip_l <- paste0("lh_", hip_segment_suffix)
hip_r <- paste0("rh_", hip_segment_suffix)
thal_l <- paste0("LeftGrp_", names(thal_group_def))
thal_r <- paste0("RightGrp_", names(thal_group_def))

hp_pair <- make_ipsi_contra_pairs(pat, "side_lr", hip_l, hip_r, "lh_", "rh_", "hp")
pat <- hp_pair$df
th_pair <- make_ipsi_contra_pairs(pat, "side_lr", thal_l, thal_r,
                                  "LeftGrp_", "RightGrp_", "thalgrp")
pat <- th_pair$df

native_roi_tbl <- bind_rows(
  tibble(Family = "Hippocampal segments", ROI = c(hip_l, hip_r)),
  tibble(Family = "Thalamic composite groups", ROI = c(thal_l, thal_r))
) %>%
  mutate(
    Hemisphere = ifelse(str_detect(ROI, "^(lh_|LeftGrp_)"), "Left", "Right"),
    Suffix = str_remove(ROI, "^(lh_|rh_|LeftGrp_|RightGrp_)"),
    ROI_Display = case_when(
      Family == "Hippocampal segments" ~ paste(Hemisphere, unname(hip_segment_display[Suffix])),
      TRUE ~ paste(Hemisphere, unname(thal_group_display[Suffix]))
    )
  )

ipsi_roi_tbl <- bind_rows(
  tibble(Family = "Hippocampal segments", ROI = c(hp_pair$ipsi_cols, hp_pair$contra_cols)),
  tibble(Family = "Thalamic composite groups", ROI = c(th_pair$ipsi_cols, th_pair$contra_cols))
) %>%
  mutate(
    Orientation = ifelse(str_starts(ROI, "ipsi_"), "Ipsi", "Contra"),
    Suffix = str_remove(ROI, "^(ipsi|contra)_(hp|thalgrp)_"),
    ROI_Display = case_when(
      Family == "Hippocampal segments" ~ paste(Orientation, unname(hip_segment_display[Suffix])),
      TRUE ~ paste(Orientation, unname(thal_group_display[Suffix]))
    )
  )

# -----------------------------------------------------------------------------
# 1) Left-onset versus right-onset clinical and cognitive characteristics.
# -----------------------------------------------------------------------------
log_msg("Running left-vs-right clinical/cognitive comparisons")

continuous_meta <- bind_rows(
  tibble(Variable = c("Age", "Duration_Dis", "Edu_year", "eTIV"),
         Display = c("Age", "Disease duration", "Education", "eTIV"),
         Type = "Clinical", Domain = "Clinical"),
  scale_meta %>% transmute(Variable = Scale, Display = Scale_Display,
                           Type = "Cognitive", Domain = Domain)
)

continuous_results <- map_dfr(seq_len(nrow(continuous_meta)), function(i) {
  v <- continuous_meta$Variable[i]
  d <- pat %>% select(side_lr, all_of(v)) %>% drop_na()
  x <- as.numeric(d[[v]][d$side_lr == "L"])
  y <- as.numeric(d[[v]][d$side_lr == "R"])
  tt <- tryCatch(t.test(x, y, var.equal = FALSE), error = function(e) NULL)
  g <- hedges_g(x, y)
  se_g <- hedges_g_se(g, length(x), length(y))
  tibble(
    Type = continuous_meta$Type[i], Domain = continuous_meta$Domain[i],
    Variable = v, Variable_Display = continuous_meta$Display[i],
    N_Left = length(x), Mean_Left = mean(x), SD_Left = sd(x),
    N_Right = length(y), Mean_Right = mean(y), SD_Right = sd(y),
    Mean_Difference_Left_minus_Right = mean(x) - mean(y),
    Hedges_g_Left_minus_Right = g,
    Hedges_g_CI_Lower = g - 1.96 * se_g,
    Hedges_g_CI_Upper = g + 1.96 * se_g,
    Welch_t = if (is.null(tt)) NA_real_ else unname(tt$statistic),
    Welch_df = if (is.null(tt)) NA_real_ else unname(tt$parameter),
    p_value = if (is.null(tt)) NA_real_ else tt$p.value
  )
}) %>%
  group_by(Type) %>% mutate(q_BH_within_type = p.adjust(p_value, method = "BH")) %>%
  ungroup()

categorical_results <- bind_rows(
  {
    z <- table(pat$side_lr, pat$sex)
    tibble(Variable = "sex", Levels = paste(colnames(z), collapse = " / "),
           Left_counts = paste(z["L", ], collapse = " / "),
           Right_counts = paste(z["R", ], collapse = " / "),
           Fisher_p = fisher.test(z)$p.value)
  },
  {
    z <- table(pat$side_lr, pat$HS_status)
    tibble(Variable = "HS_status", Levels = paste(colnames(z), collapse = " / "),
           Left_counts = paste(z["L", ], collapse = " / "),
           Right_counts = paste(z["R", ], collapse = " / "),
           Fisher_p = fisher.test(z)$p.value)
  }
)

write_table(continuous_results, "Table_L1_Left_vs_Right_continuous_characteristics.csv")
write_table(categorical_results, "Table_L1b_Left_vs_Right_categorical_characteristics.csv")

clinical_fig_dat <- continuous_results %>%
  filter(Type == "Cognitive") %>%
  mutate(Variable_Display = factor(Variable_Display, levels = rev(scale_meta$Scale_Display)),
         FDR = ifelse(q_BH_within_type < 0.05, "q < .05", "q ≥ .05"))

p_clin <- ggplot(clinical_fig_dat,
                 aes(x = Hedges_g_Left_minus_Right, y = Variable_Display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(xmin = Hedges_g_CI_Lower, xmax = Hedges_g_CI_Upper,
                    color = Domain), orientation = "y", width = 0.18,
                linewidth = 0.55) +
  geom_point(aes(color = Domain, shape = FDR), size = 2.4) +
  scale_color_manual(values = c(Verbal = "#B2182B", Visuospatial = "#2166AC",
                                `Global/other` = "#555555")) +
  scale_shape_manual(values = c(`q < .05` = 17, `q ≥ .05` = 16)) +
  theme_minimal(base_size = 12) +
  theme(axis.text.y = element_text(color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "Cognitive performance: left-onset versus right-onset TLE",
       subtitle = "Hedges' g (Left - Right); positive values indicate higher scores in left-onset TLE",
       x = "Hedges' g (95% approximate CI)", y = "", color = "Domain",
       shape = "Cognitive FDR")
save_plot(p_clin, "Figure_L1_Cognitive_Left_vs_Right_forest", 9.5, 8.5)

# -----------------------------------------------------------------------------
# 2) Native left/right fine-grained structural effects by onset side versus HC.
# -----------------------------------------------------------------------------
log_msg("Running native structural effects: left-onset vs HC and right-onset vs HC")

hip_fine_suffix <- c("Hippocampal_tail", "subiculum_comb", "CA1_comb", "CA3_comb",
                     "CA4_comb", "GC_ML_DG_comb", "molecular_layer_HP_comb",
                     "HATA", "fimbria", "presubiculum_comb", "parasubiculum")
thal_fine_suffix <- c("AV", "VA", "VAmc", "VLa", "VLp", "VPL", "VM", "CL",
                      "CeM", "CM", "Pf", "MDm", "MDl", "LD", "LP", "Pu_Total",
                      "MGN", "LGN", "L_Sg", "MV_Re", "Pc", "Pt")
struct_roi_tbl <- bind_rows(
  tibble(Family = "Hippocampus subfields",
         ROI = c(paste0("lh_", hip_fine_suffix), paste0("rh_", hip_fine_suffix))),
  tibble(Family = "Thalamus nuclei",
         ROI = c(paste0("Left_", thal_fine_suffix), paste0("Right_", thal_fine_suffix)))
) %>%
  filter(ROI %in% names(pat), ROI %in% names(hc)) %>%
  mutate(
    Hemisphere = ifelse(str_detect(ROI, "^(lh_|Left_)"), "Left", "Right"),
    Region = str_remove(ROI, "^(lh_|rh_|Left_|Right_)"),
    Region = str_replace_all(Region, "_comb$", ""),
    Region = str_replace_all(Region, "_", " "),
    ROI_Display = paste(Hemisphere, Region)
  )

struct_all <- bind_rows(
  hc %>% mutate(Group = "HC", side_lr = NA),
  pat %>% mutate(Group = ifelse(side_lr == "L", "Left-onset", "Right-onset"))
) %>% mutate(Group = factor(Group, levels = c("HC", "Left-onset", "Right-onset")))

run_structural_contrast <- function(group_level, contrast_index) {
  map_dfr(seq_len(nrow(struct_roi_tbl)), function(i) {
    roi <- struct_roi_tbl$ROI[i]
    d <- struct_all %>%
      filter(Group %in% c("HC", group_level)) %>%
      select(Group, Age, sex, eTIV, all_of(roi)) %>% drop_na()
    d$Group <- relevel(droplevels(d$Group), ref = "HC")
    X0 <- model.matrix(~ eTIV + Age + sex, data = d)
    X1 <- model.matrix(~ eTIV + Age + sex + Group, data = d)
    idx <- grep("^Group", colnames(X1))
    tst <- if (length(idx) == 1) {
      freedman_lane_term(as.numeric(d[[roi]]), X0, X1, idx,
                         strata = d$sex, B = N_PERM,
                         seed = 300000L + contrast_index * 1000L + i)
    } else c(beta = NA_real_, se = NA_real_, t = NA_real_, p_perm = NA_real_)
    yres <- residualize_one(d[[roi]], d, c("eTIV", "Age", "sex"))
    g <- hedges_g(yres[d$Group == group_level], yres[d$Group == "HC"])
    tibble(
      Contrast = paste0(group_level, " vs HC"),
      Family = struct_roi_tbl$Family[i], ROI = roi,
      ROI_Display = struct_roi_tbl$ROI_Display[i],
      Hemisphere = struct_roi_tbl$Hemisphere[i],
      N_Total = nrow(d), N_HC = sum(d$Group == "HC"),
      N_Patient = sum(d$Group == group_level),
      beta = tst["beta"], se = tst["se"], t = tst["t"],
      p_perm = tst["p_perm"], Hedges_g_adjusted = g
    )
  }) %>%
    group_by(Contrast, Family) %>%
    mutate(q_perm_family = p.adjust(p_perm, method = "BH")) %>%
    ungroup() %>% mutate(star = star_from_q(q_perm_family))
}

struct_results <- bind_rows(
  run_structural_contrast("Left-onset", 1L),
  run_structural_contrast("Right-onset", 2L)
)
write_table(struct_results, "Table_L2_Native_structural_effects_by_onset_side_vs_HC.csv")

struct_order <- struct_roi_tbl$ROI_Display
lim_g <- quantile(abs(struct_results$Hedges_g_adjusted), 0.98, na.rm = TRUE)
if (!is.finite(lim_g) || lim_g < 0.5) lim_g <- 0.5
hm_struct <- struct_results %>%
  mutate(Contrast = factor(Contrast, levels = c("Left-onset vs HC", "Right-onset vs HC")),
         ROI_Display = factor(ROI_Display, levels = rev(struct_order)))
p_struct <- ggplot(hm_struct, aes(x = Contrast, y = ROI_Display,
                                  fill = Hedges_g_adjusted)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = star), size = 3.4, fontface = "bold") +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-lim_g, lim_g), oob = squish, na.value = "grey90") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(size = 12, face = "bold", angle = 20, hjust = 1,
                                   color = "black"),
        axis.text.y = element_text(size = 8.2, color = "black"),
        strip.text.y = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 15, face = "bold"),
        legend.position = "right") +
  labs(title = "Native left/right structural effects by seizure-onset side",
       subtitle = paste0("Adjusted Hedges' g; BH-FDR within Family × Contrast; Nperm = ",
                         N_PERM, ".\nNegative values indicate smaller volume than HC."),
       x = "", y = "Native-hemisphere ROI", fill = "Adjusted\nHedges' g")
save_plot(p_struct, "Figure_L2_Native_structural_effects_by_onset_side", 9.5, 15)

# -----------------------------------------------------------------------------
# 3) Onset-side-stratified ipsi/contra correlations and native L/R correlations.
# -----------------------------------------------------------------------------
get_covars <- function(is_normed, include_side = FALSE) {
  z <- c("sex", "Duration_Dis", "Edu_year", "eTIV")
  if (include_side) z <- c(z, "side_lr")
  if (isTRUE(is_normed)) z else c(z, "Age")
}

run_corr_table <- function(dat, roi_tbl, analysis, subgroup, include_side = FALSE,
                           seed_offset = 0L) {
  out <- vector("list", nrow(roi_tbl) * nrow(scale_meta))
  k <- 1L
  for (i in seq_len(nrow(roi_tbl))) {
    roi <- roi_tbl$ROI[i]
    for (j in seq_len(nrow(scale_meta))) {
      sc <- scale_meta$Scale[j]
      covars <- get_covars(scale_meta$Is_Normed[j], include_side)
      d <- dat %>% select(all_of(unique(c(roi, sc, covars)))) %>% drop_na()
      covars <- effective_covars(d, covars)
      if (nrow(d) >= MIN_N) {
        xr <- residualize_one(d[[roi]], d, covars)
        yr <- residualize_one(d[[sc]], d, covars)
        z <- perm_corr(xr, yr, B = N_PERM,
                       seed = seed_offset + i * 100L + j)
      } else z <- c(r = NA_real_, p = NA_real_)
      out[[k]] <- tibble(
        Analysis = analysis, Subgroup = subgroup,
        Family = roi_tbl$Family[i], ROI = roi,
        ROI_Display = roi_tbl$ROI_Display[i],
        Scale = sc, Scale_Display = scale_meta$Scale_Display[j],
        Domain = scale_meta$Domain[j], Is_Normed = scale_meta$Is_Normed[j],
        CovarsUsed = paste(covars, collapse = " + "), n = nrow(d),
        r = z["r"], perm_p = z["p"]
      )
      k <- k + 1L
    }
  }
  bind_rows(out) %>%
    group_by(Analysis, Subgroup, Scale, Family) %>%
    mutate(q_perm = p.adjust(perm_p, method = "BH")) %>%
    ungroup() %>% mutate(star = star_from_q(q_perm))
}

log_msg("Running onset-side-stratified ipsi/contra correlations")
set.seed(2026)
side_corr <- bind_rows(
  run_corr_table(filter(pat, side_lr == "L"), ipsi_roi_tbl,
                 "Onset-side-stratified ipsi/contra", "Left-onset",
                 include_side = FALSE, seed_offset = 400000L),
  run_corr_table(filter(pat, side_lr == "R"), ipsi_roi_tbl,
                 "Onset-side-stratified ipsi/contra", "Right-onset",
                 include_side = FALSE, seed_offset = 500000L)
)
write_table(side_corr, "Table_L3_Onset_side_stratified_ipsi_contra_correlations.csv")

log_msg("Running native left/right correlations with onset side as a covariate")
native_corr <- run_corr_table(pat, native_roi_tbl,
                              "Native left/right hemisphere", "All patients",
                              include_side = TRUE, seed_offset = 600000L)
write_table(native_corr, "Table_L4_Native_left_right_correlations.csv")

scale_order <- scale_meta$Scale_Display
ipsi_order <- ipsi_roi_tbl$ROI_Display
lim_r_side <- max(0.2, quantile(abs(side_corr$r), 0.98, na.rm = TRUE))
hm_side <- side_corr %>%
  mutate(Subgroup = factor(Subgroup, levels = c("Left-onset", "Right-onset")),
         ROI_Display = factor(ROI_Display, levels = rev(ipsi_order)),
         Scale_Display = factor(Scale_Display, levels = scale_order))
p_side <- ggplot(hm_side, aes(x = Scale_Display, y = ROI_Display, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = star), size = 3.0, fontface = "bold") +
  facet_grid(Family ~ Subgroup, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-lim_r_side, lim_r_side), oob = squish,
                       na.value = "grey90") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8.5,
                                   face = "bold", color = "black"),
        axis.text.y = element_text(size = 8.4, color = "black"),
        strip.text = element_text(size = 10.5, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "Onset-side-stratified ipsi/contra structure-cognition associations",
       subtitle = paste0("Partial r; BH-FDR within Scale × Family; Nperm = ", N_PERM,
                         ". ○ q<.10, * q<.05, ** q<.01, *** q<.001"),
       x = "", y = "ROI", fill = "partial r")
save_plot(p_side, "Figure_L3_Onset_stratified_ipsi_contra_correlations", 16, 10)

native_order <- native_roi_tbl$ROI_Display
lim_r_native <- max(0.2, quantile(abs(native_corr$r), 0.98, na.rm = TRUE))
hm_native <- native_corr %>%
  mutate(ROI_Display = factor(ROI_Display, levels = rev(native_order)),
         Scale_Display = factor(Scale_Display, levels = scale_order))
p_native <- ggplot(hm_native, aes(x = Scale_Display, y = ROI_Display, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = star), size = 3.1, fontface = "bold") +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-lim_r_native, lim_r_native), oob = squish,
                       na.value = "grey90") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8.5,
                                   face = "bold", color = "black"),
        axis.text.y = element_text(size = 8.5, color = "black"),
        strip.text.y = element_text(size = 10.5, face = "bold"),
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "Native left/right structure-cognition associations",
       subtitle = paste0("Partial r adjusted for onset side and clinical covariates; ",
                         "BH-FDR within Scale × Family; Nperm = ", N_PERM),
       x = "", y = "Native-hemisphere ROI", fill = "partial r")
save_plot(p_native, "Figure_L4_Native_left_right_correlations", 13, 9)

# -----------------------------------------------------------------------------
# 4) Focused ROI × onset-side interaction models for verbal/visuospatial tests.
# -----------------------------------------------------------------------------
log_msg("Running native ROI by onset-side interaction models")
interaction_scales <- scale_meta %>% filter(Domain %in% c("Verbal", "Visuospatial"))

run_interaction <- function(roi_index, scale_index) {
  roi <- native_roi_tbl$ROI[roi_index]
  sc <- interaction_scales$Scale[scale_index]
  covars <- get_covars(interaction_scales$Is_Normed[scale_index], include_side = FALSE)
  d <- pat %>%
    select(side_lr, all_of(unique(c(roi, sc, covars)))) %>% drop_na()
  covars <- effective_covars(d, covars)
  if (nrow(d) < MIN_N || length(unique(d$side_lr)) < 2) {
    return(tibble(Family = native_roi_tbl$Family[roi_index], ROI = roi,
                  ROI_Display = native_roi_tbl$ROI_Display[roi_index], Scale = sc,
                  Scale_Display = interaction_scales$Scale_Display[scale_index],
                  Domain = interaction_scales$Domain[scale_index], n = nrow(d),
                  N_Left = sum(d$side_lr == "L"), N_Right = sum(d$side_lr == "R"),
                  beta_interaction = NA_real_, se_interaction = NA_real_,
                  CI_lower = NA_real_, CI_upper = NA_real_, t_interaction = NA_real_,
                  p_parametric = NA_real_, p_perm = NA_real_,
                  slope_Right = NA_real_, slope_Left = NA_real_))
  }
  d$side_lr <- factor(d$side_lr, levels = c("R", "L"))
  d$ROI_z <- as.numeric(scale(as.numeric(d[[roi]])))
  d$Y_z <- as.numeric(scale(as.numeric(d[[sc]])))
  reduced_terms <- c("ROI_z", "side_lr", covars)
  full_terms <- c("ROI_z * side_lr", covars)
  X0 <- model.matrix(as.formula(paste0("~ ", paste(reduced_terms, collapse = " + "))), data = d)
  X1 <- model.matrix(as.formula(paste0("~ ", paste(full_terms, collapse = " + "))), data = d)
  idx_int <- grep("ROI_z:side_lrL|side_lrL:ROI_z", colnames(X1))
  if (length(idx_int) != 1 || qr(X1)$rank < ncol(X1)) {
    tst <- c(beta = NA_real_, se = NA_real_, t = NA_real_, p_perm = NA_real_)
    slope_r <- slope_l <- p_param <- NA_real_
  } else {
    seed <- 700000L + roi_index * 100L + scale_index
    tst <- freedman_lane_term(d$Y_z, X0, X1, idx_int, strata = d$sex,
                              B = N_PERM, seed = seed)
    fit <- lm.fit(X1, d$Y_z)
    r <- d$Y_z - as.numeric(X1 %*% fit$coefficients)
    df_res <- nrow(X1) - ncol(X1)
    V <- sum(r^2) / df_res * solve(crossprod(X1))
    idx_roi <- match("ROI_z", colnames(X1))
    slope_r <- fit$coefficients[idx_roi]
    slope_l <- fit$coefficients[idx_roi] + fit$coefficients[idx_int]
    p_param <- 2 * pt(-abs(tst["t"]), df = df_res)
  }
  tibble(
    Family = native_roi_tbl$Family[roi_index], ROI = roi,
    ROI_Display = native_roi_tbl$ROI_Display[roi_index], Scale = sc,
    Scale_Display = interaction_scales$Scale_Display[scale_index],
    Domain = interaction_scales$Domain[scale_index],
    CovarsUsed = paste(covars, collapse = " + "), n = nrow(d),
    N_Left = sum(d$side_lr == "L"), N_Right = sum(d$side_lr == "R"),
    beta_interaction = tst["beta"], se_interaction = tst["se"],
    CI_lower = tst["beta"] - 1.96 * tst["se"],
    CI_upper = tst["beta"] + 1.96 * tst["se"],
    t_interaction = tst["t"], p_parametric = p_param,
    p_perm = tst["p_perm"], slope_Right = slope_r, slope_Left = slope_l
  )
}

interaction_results <- map_dfr(seq_len(nrow(native_roi_tbl)), function(i) {
  map_dfr(seq_len(nrow(interaction_scales)), function(j) run_interaction(i, j))
}) %>%
  group_by(Scale, Family) %>%
  mutate(q_perm = p.adjust(p_perm, method = "BH")) %>%
  ungroup() %>% mutate(star = star_from_q(q_perm))
write_table(interaction_results, "Table_L5_Native_ROI_by_onset_side_interactions.csv")

forest_dat <- interaction_results %>% filter(is.finite(q_perm)) %>% arrange(q_perm, p_perm)
if (sum(forest_dat$q_perm < 0.10, na.rm = TRUE) >= 8) {
  forest_dat <- forest_dat %>% filter(q_perm < 0.10)
} else {
  forest_dat <- forest_dat %>% slice_head(n = 15)
}
forest_dat <- forest_dat %>%
  mutate(Label = paste0(ROI_Display, " | ", Scale_Display),
         Label = factor(Label, levels = rev(unique(Label))),
         q_group = ifelse(q_perm < 0.05, "q < .05", "q ≥ .05"))
p_int <- ggplot(forest_dat, aes(x = beta_interaction, y = Label, color = Family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper), orientation = "y",
                width = 0.16, linewidth = 0.55) +
  geom_point(aes(shape = q_group), size = 2.5) +
  scale_color_manual(values = c(`Hippocampal segments` = "#B2182B",
                                `Thalamic composite groups` = "#2166AC")) +
  scale_shape_manual(values = c(`q < .05` = 17, `q ≥ .05` = 16)) +
  theme_minimal(base_size = 11.5) +
  theme(axis.text.y = element_text(size = 8.8, color = "black"),
        axis.text.x = element_text(color = "black"),
        legend.position = "bottom",
        plot.title = element_text(size = 15, face = "bold")) +
  labs(title = "Native ROI × seizure-onset side interactions",
       subtitle = "Standardized interaction coefficient (Left-onset slope minus Right-onset slope); top FDR-ranked results",
       x = "Standardized interaction beta (95% CI)", y = "", color = "ROI family",
       shape = "Permutation FDR")
save_plot(p_int, "Figure_L5_ROI_by_onset_side_interaction_forest", 11.5,
          max(6.5, 0.34 * nrow(forest_dat) + 2.2))

# -----------------------------------------------------------------------------
# QC and concise result summaries.
# -----------------------------------------------------------------------------
summary_counts <- bind_rows(
  struct_results %>% group_by(Analysis = "Native structure vs HC", Subgroup = Contrast,
                              Family) %>%
    summarise(N_tests = sum(is.finite(p_perm)), N_q05 = sum(q_perm_family < 0.05, na.rm = TRUE),
              N_q10 = sum(q_perm_family < 0.10, na.rm = TRUE),
              Min_q = min(q_perm_family, na.rm = TRUE), .groups = "drop"),
  side_corr %>% group_by(Analysis, Subgroup, Family) %>%
    summarise(N_tests = sum(is.finite(perm_p)), N_q05 = sum(q_perm < 0.05, na.rm = TRUE),
              N_q10 = sum(q_perm < 0.10, na.rm = TRUE),
              Min_q = min(q_perm, na.rm = TRUE), .groups = "drop"),
  native_corr %>% group_by(Analysis, Subgroup, Family) %>%
    summarise(N_tests = sum(is.finite(perm_p)), N_q05 = sum(q_perm < 0.05, na.rm = TRUE),
              N_q10 = sum(q_perm < 0.10, na.rm = TRUE),
              Min_q = min(q_perm, na.rm = TRUE), .groups = "drop"),
  interaction_results %>% group_by(Analysis = "ROI by onset-side interaction",
                                   Subgroup = "All patients", Family) %>%
    summarise(N_tests = sum(is.finite(p_perm)), N_q05 = sum(q_perm < 0.05, na.rm = TRUE),
              N_q10 = sum(q_perm < 0.10, na.rm = TRUE),
              Min_q = min(q_perm, na.rm = TRUE), .groups = "drop")
)
write_table(summary_counts, "Table_L6_Analysis_summary_counts.csv")

sample_qc <- tibble(
  Metric = c("N patients", "N left-onset", "N right-onset", "N HC",
             "N HS left-onset", "N HS right-onset", "N non-HS left-onset",
             "N non-HS right-onset"),
  Value = c(nrow(pat), sum(pat$side_lr == "L"), sum(pat$side_lr == "R"), nrow(hc),
            sum(pat$side_lr == "L" & pat$HS_status == "HS"),
            sum(pat$side_lr == "R" & pat$HS_status == "HS"),
            sum(pat$side_lr == "L" & pat$HS_status == "non-HS"),
            sum(pat$side_lr == "R" & pat$HS_status == "non-HS"))
)
write_table(sample_qc, "QC_sample_counts.csv")

availability_qc <- bind_rows(
  tibble(Data = "Patients", Variable = c(struct_roi_tbl$ROI, scale_meta$Scale),
         N_nonmissing = sapply(c(struct_roi_tbl$ROI, scale_meta$Scale),
                               function(v) sum(is.finite(suppressWarnings(as.numeric(pat[[v]]))))),
         N_total = nrow(pat)),
  tibble(Data = "HC", Variable = struct_roi_tbl$ROI,
         N_nonmissing = sapply(struct_roi_tbl$ROI,
                               function(v) sum(is.finite(suppressWarnings(as.numeric(hc[[v]]))))),
         N_total = nrow(hc))
)
readr::write_excel_csv(availability_qc, file.path(QC_DIR, "QC_data_availability.csv"), na = "")
writeLines(capture.output(sessionInfo()), file.path(QC_DIR, "sessionInfo.txt"), useBytes = TRUE)

manifest <- tibble(
  Component = c("Clinical/cognitive left-right comparison", "Native structure vs HC",
                "Onset-stratified ipsi/contra correlations", "Native L/R correlations",
                "ROI by onset-side interactions"),
  Main_table = c("Table_L1_Left_vs_Right_continuous_characteristics.csv",
                 "Table_L2_Native_structural_effects_by_onset_side_vs_HC.csv",
                 "Table_L3_Onset_side_stratified_ipsi_contra_correlations.csv",
                 "Table_L4_Native_left_right_correlations.csv",
                 "Table_L5_Native_ROI_by_onset_side_interactions.csv"),
  Main_figure = c("Figure_L1_Cognitive_Left_vs_Right_forest",
                  "Figure_L2_Native_structural_effects_by_onset_side",
                  "Figure_L3_Onset_stratified_ipsi_contra_correlations",
                  "Figure_L4_Native_left_right_correlations",
                  "Figure_L5_ROI_by_onset_side_interaction_forest")
)
write_table(manifest, "Analysis_manifest.csv")

log_msg("Completed all analyses")
log_msg("Tables: ", normalizePath(TABLE_DIR, winslash = "/"))
log_msg("Figures: ", normalizePath(FIG_DIR, winslash = "/"))
log_msg("QC: ", normalizePath(QC_DIR, winslash = "/"))
print(summary_counts, n = Inf)

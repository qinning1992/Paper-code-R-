#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

set.seed(20260824)

patient_file <- file.path("code and data", "df_raw_clean.xlsx")
hc_file <- file.path("code and data", "df_hc_clean.xlsx")
out_dir <- "PLS"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

n_perm <- 5000
n_boot <- 2000

to_num <- function(x) suppressWarnings(as.numeric(x))

clean_roi_name <- function(x) {
  hemi <- dplyr::case_when(
    grepl("^lh_", x) ~ "left",
    grepl("^rh_", x) ~ "right",
    grepl("^ipsi_", x) ~ "ipsilateral",
    grepl("^contra_", x) ~ "contralateral",
    TRUE ~ "unknown"
  )
  region <- x
  region <- sub("^lh_lh_", "", region)
  region <- sub("^rh_rh_", "", region)
  region <- sub("^ipsi_", "", region)
  region <- sub("^contra_", "", region)
  region <- sub("_volume$", "", region)
  data.frame(hemi = hemi, region = region, stringsAsFactors = FALSE)
}

make_roi_label <- function(hemi, region) {
  prefix <- dplyr::case_when(
    hemi == "left" ~ "L ",
    hemi == "right" ~ "R ",
    hemi == "ipsilateral" ~ "Ipsi ",
    hemi == "contralateral" ~ "Contra ",
    TRUE ~ ""
  )
  paste0(prefix, region)
}

prepare_data <- function(patient_df, hc_df) {
  cortical_cols <- grep("^(lh_lh_|rh_rh_).+_volume$", names(patient_df), value = TRUE)
  cortical_cols <- intersect(cortical_cols, names(hc_df))
  if (length(cortical_cols) != 68) {
    stop("Expected 68 shared cortical volume variables, found ", length(cortical_cols), ".")
  }

  pat <- patient_df %>%
    mutate(
      group = "MTLE",
      group_num = 1,
      side_num = to_num(.data[["side"]]),
      Age = to_num(.data[["Age"]]),
      sex = factor(.data[["sex"]]),
      eTIV = to_num(.data[["eTIV"]])
    )

  hc <- hc_df %>%
    mutate(
      group = "HC",
      group_num = 0,
      side_num = NA_real_,
      Age = to_num(.data[["Age"]]),
      sex = factor(.data[["sex"]]),
      eTIV = to_num(.data[["eTIV"]])
    )

  keep_cols <- c("group", "group_num", "side_num", "Age", "sex", "eTIV", cortical_cols)
  bind_rows(pat[, keep_cols], hc[, keep_cols]) %>%
    mutate(row_id = row_number())
}

residualize_and_scale <- function(dat, roi_cols) {
  needed <- c("group", "group_num", "Age", "sex", "eTIV", roi_cols)
  dat <- dat[complete.cases(dat[, needed]), needed]
  mm <- model.matrix(~ Age + sex + eTIV, data = dat)
  y_mat <- as.matrix(dat[, roi_cols])
  storage.mode(y_mat) <- "double"

  resid_mat <- matrix(NA_real_, nrow = nrow(y_mat), ncol = ncol(y_mat))
  colnames(resid_mat) <- roi_cols
  for (j in seq_along(roi_cols)) {
    fit <- lm.fit(mm, y_mat[, j])
    resid_mat[, j] <- fit$residuals
  }
  scaled_mat <- scale(resid_mat, center = TRUE, scale = TRUE)
  list(dat = dat, x = as.matrix(scaled_mat))
}

make_contrast <- function(group_num) {
  n_patient <- sum(group_num == 1)
  n_hc <- sum(group_num == 0)
  y <- ifelse(group_num == 1, 1 / n_patient, -1 / n_hc)
  y / sqrt(sum(y^2))
}

fit_pls <- function(x, group_num) {
  y <- make_contrast(group_num)
  cov_vec <- as.vector(crossprod(y, x))
  sv <- sqrt(sum(cov_vec^2))
  sal <- cov_vec / sv
  scores <- as.vector(x %*% sal)

  if (mean(scores[group_num == 1]) > mean(scores[group_num == 0])) {
    sal <- -sal
    scores <- -scores
    cov_vec <- -cov_vec
  }

  list(
    singular_value = sv,
    salience = sal,
    scores = scores,
    covariance_vector = cov_vec
  )
}

cohen_d <- function(x, g) {
  x1 <- x[g == 1]
  x0 <- x[g == 0]
  n1 <- length(x1)
  n0 <- length(x0)
  pooled_sd <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / (n1 + n0 - 2))
  (mean(x1) - mean(x0)) / pooled_sd
}

permutation_test <- function(x, group_num, n_iter) {
  obs <- fit_pls(x, group_num)$singular_value
  perm_sv <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    perm_sv[i] <- fit_pls(x, sample(group_num, length(group_num), replace = FALSE))$singular_value
  }
  p <- (sum(perm_sv >= obs) + 1) / (n_iter + 1)
  list(observed = obs, permuted = perm_sv, p_value = p)
}

bootstrap_salience <- function(dat, roi_cols, original_salience, n_iter) {
  idx_patient <- which(dat$group_num == 1)
  idx_hc <- which(dat$group_num == 0)
  boot_mat <- matrix(NA_real_, nrow = n_iter, ncol = length(roi_cols))
  colnames(boot_mat) <- roi_cols

  for (i in seq_len(n_iter)) {
    boot_idx <- c(
      sample(idx_patient, length(idx_patient), replace = TRUE),
      sample(idx_hc, length(idx_hc), replace = TRUE)
    )
    boot_dat <- dat[boot_idx, , drop = FALSE]
    prep <- residualize_and_scale(boot_dat, roi_cols)
    fit <- fit_pls(prep$x, prep$dat$group_num)
    sal <- fit$salience
    if (sum(sal * original_salience) < 0) sal <- -sal
    boot_mat[i, ] <- sal
  }

  boot_se <- apply(boot_mat, 2, sd, na.rm = TRUE)
  boot_low <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
  boot_high <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)
  list(se = boot_se, low = boot_low, high = boot_high, samples = boot_mat)
}

run_analysis <- function(all_dat, roi_cols, analysis_name, patient_filter = NULL, run_seed = NULL) {
  if (!is.null(run_seed)) set.seed(run_seed)

  if (is.null(patient_filter)) {
    dat <- all_dat
  } else {
    dat <- all_dat %>% filter(group == "HC" | patient_filter(.))
  }

  prep <- residualize_and_scale(dat, roi_cols)
  dat2 <- prep$dat
  x <- prep$x
  fit <- fit_pls(x, dat2$group_num)
  perm <- permutation_test(x, dat2$group_num, n_perm)
  boot <- bootstrap_salience(dat2, roi_cols, fit$salience, n_boot)

  score_dat <- data.frame(
    analysis = analysis_name,
    group = dat2$group,
    group_num = dat2$group_num,
    score = fit$scores
  )

  roi_meta <- clean_roi_name(roi_cols)
  salience_dat <- data.frame(
    analysis = analysis_name,
    roi = roi_cols,
    hemi = roi_meta$hemi,
    region = roi_meta$region,
    salience = fit$salience,
    bootstrap_se = boot$se,
    bootstrap_ratio = fit$salience / boot$se,
    bootstrap_ci_low = boot$low,
    bootstrap_ci_high = boot$high,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      abs_bootstrap_ratio = abs(bootstrap_ratio),
      rank_abs_bsr = rank(-abs_bootstrap_ratio, ties.method = "first"),
      direction_in_patients = ifelse(salience > 0, "lower in patients", "higher in patients")
    ) %>%
    arrange(rank_abs_bsr)

  global_dat <- data.frame(
    analysis = analysis_name,
    n_patient = sum(dat2$group_num == 1),
    n_hc = sum(dat2$group_num == 0),
    n_total = nrow(dat2),
    n_roi = length(roi_cols),
    n_perm = n_perm,
    n_bootstrap = n_boot,
    singular_value = fit$singular_value,
    permutation_p = perm$p_value,
    patient_score_mean = mean(fit$scores[dat2$group_num == 1]),
    hc_score_mean = mean(fit$scores[dat2$group_num == 0]),
    patient_minus_hc_score = mean(fit$scores[dat2$group_num == 1]) - mean(fit$scores[dat2$group_num == 0]),
    score_cohen_d = cohen_d(fit$scores, dat2$group_num)
  )

  perm_dat <- data.frame(
    analysis = analysis_name,
    permuted_singular_value = perm$permuted,
    observed_singular_value = perm$observed
  )

  list(global = global_dat, salience = salience_dat, scores = score_dat, perm = perm_dat)
}

make_aligned_cortical_data <- function(all_dat, roi_cols, target_side) {
  if (!target_side %in% c(0, 1)) {
    stop("target_side must be 1 for left-onset or 0 for right-onset.")
  }

  left_cols <- roi_cols[grepl("^lh_lh_", roi_cols)]
  regions <- sub("^lh_lh_", "", left_cols)
  regions <- sub("_volume$", "", regions)
  right_cols <- paste0("rh_rh_", regions, "_volume")
  missing_right <- setdiff(right_cols, roi_cols)
  if (length(missing_right) > 0) {
    stop("Missing right-hemisphere homologous cortical variables: ", paste(missing_right, collapse = ", "))
  }

  dat <- all_dat %>% filter(group == "HC" | side_num == target_side)
  if (target_side == 1) {
    ipsi_source <- left_cols
    contra_source <- right_cols
  } else {
    ipsi_source <- right_cols
    contra_source <- left_cols
  }

  aligned <- dat[, c("group", "group_num", "side_num", "Age", "sex", "eTIV")]
  ipsi_cols <- paste0("ipsi_", regions, "_volume")
  contra_cols <- paste0("contra_", regions, "_volume")
  aligned[, ipsi_cols] <- dat[, ipsi_source]
  aligned[, contra_cols] <- dat[, contra_source]

  list(dat = aligned, roi_cols = c(ipsi_cols, contra_cols))
}

patients <- read_excel(patient_file)
hcs <- read_excel(hc_file)
all_dat <- prepare_data(patients, hcs)
roi_cols <- grep("^(lh_lh_|rh_rh_).+_volume$", names(patients), value = TRUE)
roi_cols <- intersect(roi_cols, names(hcs))

diagnostics <- data.frame(
  item = c(
    "raw_patient_rows", "raw_hc_rows", "valid_patient_age_etiv",
    "valid_hc_age_etiv", "left_mtle_side_1", "right_mtle_side_0",
    "shared_cortical_volume_rois", "covariates"
  ),
  value = c(
    nrow(patients),
    nrow(hcs),
    sum(!is.na(to_num(patients[["Age"]])) & !is.na(to_num(patients[["eTIV"]]))),
    sum(!is.na(to_num(hcs[["Age"]])) & !is.na(to_num(hcs[["eTIV"]]))),
    sum(to_num(patients[["side"]]) == 1, na.rm = TRUE),
    sum(to_num(patients[["side"]]) == 0, na.rm = TRUE),
    length(roi_cols),
    "Age, sex, eTIV"
  )
)
write_csv(diagnostics, file.path(out_dir, "PLS_CorticalSensitivity_diagnostics.csv"))

results <- list(
  run_analysis(all_dat, roi_cols, "All MTLE vs HC", run_seed = 1001),
  run_analysis(
    all_dat, roi_cols, "Left MTLE vs HC",
    patient_filter = function(d) d$side_num == 1,
    run_seed = 1002
  ),
  run_analysis(
    all_dat, roi_cols, "Right MTLE vs HC",
    patient_filter = function(d) d$side_num == 0,
    run_seed = 1003
  )
)

aligned_left <- make_aligned_cortical_data(all_dat, roi_cols, target_side = 1)
aligned_right <- make_aligned_cortical_data(all_dat, roi_cols, target_side = 0)

aligned_results <- list(
  run_analysis(
    aligned_left$dat, aligned_left$roi_cols,
    "Left MTLE vs HC, ipsi/contra aligned",
    run_seed = 1002
  ),
  run_analysis(
    aligned_right$dat, aligned_right$roi_cols,
    "Right MTLE vs HC, ipsi/contra aligned",
    run_seed = 1003
  )
)

aligned_mapping <- data.frame(
  analysis = c(
    "Left MTLE vs HC, ipsi/contra aligned",
    "Right MTLE vs HC, ipsi/contra aligned"
  ),
  patient_side_definition = c("side = 1, left-onset MTLE", "side = 0, right-onset MTLE"),
  patient_ipsilateral_source = c("left cortical parcels", "right cortical parcels"),
  patient_contralateral_source = c("right cortical parcels", "left cortical parcels"),
  hc_reference_ipsi_label = c(
    "left cortical parcels, used as the left-onset reference orientation",
    "right cortical parcels, used as the right-onset reference orientation"
  ),
  hc_reference_contra_label = c(
    "right cortical parcels, used as the left-onset reference orientation",
    "left cortical parcels, used as the right-onset reference orientation"
  ),
  note = "Healthy controls have no true epileptogenic side; ipsi/contra labels in HCs are orientation references for the corresponding patient subgroup.",
  stringsAsFactors = FALSE
)

global_results <- bind_rows(lapply(results, `[[`, "global"))
salience_results <- bind_rows(lapply(results, `[[`, "salience"))
score_results <- bind_rows(lapply(results, `[[`, "scores"))
perm_results <- bind_rows(lapply(results, `[[`, "perm"))

aligned_global_results <- bind_rows(lapply(aligned_results, `[[`, "global"))
aligned_salience_results <- bind_rows(lapply(aligned_results, `[[`, "salience"))
aligned_score_results <- bind_rows(lapply(aligned_results, `[[`, "scores"))
aligned_perm_results <- bind_rows(lapply(aligned_results, `[[`, "perm"))

write_csv(global_results, file.path(out_dir, "PLS_CorticalSensitivity_global_results.csv"))
write_csv(salience_results, file.path(out_dir, "PLS_CorticalSensitivity_bootstrap_saliences.csv"))
write_csv(
  salience_results %>% group_by(analysis) %>% slice_min(rank_abs_bsr, n = 20) %>% ungroup(),
  file.path(out_dir, "PLS_CorticalSensitivity_top20_saliences.csv")
)
write_csv(score_results, file.path(out_dir, "PLS_CorticalSensitivity_scores.csv"))
write_csv(perm_results, file.path(out_dir, "PLS_CorticalSensitivity_permutation_distribution.csv"))

write_csv(aligned_global_results, file.path(out_dir, "PLS_CorticalSensitivity_aligned_global_results.csv"))
write_csv(aligned_salience_results, file.path(out_dir, "PLS_CorticalSensitivity_aligned_bootstrap_saliences.csv"))
write_csv(
  aligned_salience_results %>% group_by(analysis) %>% slice_min(rank_abs_bsr, n = 20) %>% ungroup(),
  file.path(out_dir, "PLS_CorticalSensitivity_aligned_top20_saliences.csv")
)
write_csv(aligned_score_results, file.path(out_dir, "PLS_CorticalSensitivity_aligned_scores.csv"))
write_csv(aligned_perm_results, file.path(out_dir, "PLS_CorticalSensitivity_aligned_permutation_distribution.csv"))
write_csv(aligned_mapping, file.path(out_dir, "PLS_CorticalSensitivity_aligned_mapping.csv"))

p_perm <- ggplot(perm_results, aes(x = permuted_singular_value)) +
  geom_histogram(bins = 45, fill = "grey80", color = "white") +
  geom_vline(aes(xintercept = observed_singular_value), color = "#B2182B", linewidth = 0.8) +
  facet_wrap(~ analysis, scales = "free") +
  labs(x = "Permuted singular value", y = "Permutation count") +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_permutation.png"), p_perm, width = 9, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_permutation.pdf"), p_perm, width = 9, height = 4.8)

p_scores <- ggplot(score_results, aes(x = group, y = score, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
  facet_wrap(~ analysis, scales = "free_x") +
  scale_fill_manual(values = c(HC = "#4D9221", MTLE = "#B2182B")) +
  labs(x = NULL, y = "Cortical-volume LV score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_scores.png"), p_scores, width = 8.5, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_scores.pdf"), p_scores, width = 8.5, height = 4.8)

top_for_plot <- salience_results %>%
  group_by(analysis) %>%
  slice_min(rank_abs_bsr, n = 20) %>%
  ungroup() %>%
  mutate(
    roi_label = make_roi_label(hemi, region),
    roi_label = reorder(roi_label, bootstrap_ratio)
  )

p_sal <- ggplot(top_for_plot, aes(x = roi_label, y = bootstrap_ratio, fill = bootstrap_ratio > 0)) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(~ analysis, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B")) +
  labs(x = NULL, y = "Bootstrap ratio") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_top20_saliences.png"), p_sal, width = 10, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_top20_saliences.pdf"), p_sal, width = 10, height = 7.5)

p_aligned_perm <- ggplot(aligned_perm_results, aes(x = permuted_singular_value)) +
  geom_histogram(bins = 45, fill = "grey80", color = "white") +
  geom_vline(aes(xintercept = observed_singular_value), color = "#B2182B", linewidth = 0.8) +
  facet_wrap(~ analysis, scales = "free") +
  labs(x = "Permuted singular value", y = "Permutation count") +
  theme_classic(base_size = 11)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_aligned_permutation.png"), p_aligned_perm, width = 8.5, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_aligned_permutation.pdf"), p_aligned_perm, width = 8.5, height = 4.8)

p_aligned_scores <- ggplot(aligned_score_results, aes(x = group, y = score, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.7) +
  facet_wrap(~ analysis, scales = "free_x") +
  scale_fill_manual(values = c(HC = "#4D9221", MTLE = "#B2182B")) +
  labs(x = NULL, y = "Aligned cortical-volume LV score") +
  theme_classic(base_size = 11) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_aligned_scores.png"), p_aligned_scores, width = 8.5, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_aligned_scores.pdf"), p_aligned_scores, width = 8.5, height = 4.8)

aligned_top_for_plot <- aligned_salience_results %>%
  group_by(analysis) %>%
  slice_min(rank_abs_bsr, n = 20) %>%
  ungroup() %>%
  mutate(
    roi_label = make_roi_label(hemi, region),
    roi_label = reorder(roi_label, bootstrap_ratio)
  )

p_aligned_sal <- ggplot(aligned_top_for_plot, aes(x = roi_label, y = bootstrap_ratio, fill = bootstrap_ratio > 0)) +
  geom_col(width = 0.75) +
  coord_flip() +
  facet_wrap(~ analysis, scales = "free_y") +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B")) +
  labs(x = NULL, y = "Bootstrap ratio") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_aligned_top20_saliences.png"), p_aligned_sal, width = 10, height = 6.5, dpi = 300)
ggsave(file.path(out_dir, "PLS_CorticalSensitivity_aligned_top20_saliences.pdf"), p_aligned_sal, width = 10, height = 6.5)

sink(file.path(out_dir, "PLS_CorticalSensitivity_run_log.txt"))
cat("Cortical-volume PLS sensitivity analysis\n")
cat("Run date:", as.character(Sys.time()), "\n")
cat("Input files:\n")
cat(" -", patient_file, "\n")
cat(" -", hc_file, "\n")
cat("Analyses: All MTLE vs HC; Left MTLE vs HC; Right MTLE vs HC\n")
cat("Aligned analyses: Left MTLE vs HC with left cortex treated as ipsilateral; Right MTLE vs HC with right cortex treated as ipsilateral\n")
cat("Covariates residualized from cortical volumes: Age, sex, eTIV\n")
cat("Permutations:", n_perm, "\n")
cat("Stratified bootstraps:", n_boot, "\n\n")
print(diagnostics)
cat("\nGlobal results:\n")
print(global_results)
cat("\nAligned global results:\n")
print(aligned_global_results)
cat("\nAligned side mapping:\n")
print(aligned_mapping)
sink()

print(global_results)
print(aligned_global_results)
message("PLS cortical sensitivity analysis completed. Outputs written to: ", normalizePath(out_dir))

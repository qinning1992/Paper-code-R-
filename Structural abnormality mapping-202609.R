# ==============================================================
# High-resolution structural abnormality mapping analysis
# FINAL CLEAN VERSION
# Outputs ONLY:
#   P3: All PAT / HS / Other vs HC  (3 files)
#   P4: MTLE / Left MTLE / Right MTLE vs HC (3 files)
# ==============================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

# ==========================================================
# 0) Required input objects
# ==========================================================
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found (patients)")
if(!exists("df_hc_clean"))  stop("❌ df_hc_clean not found (healthy controls)")

df_pat <- df_raw_clean
df_hc0 <- df_hc_clean

# IMPORTANT ASSUMPTION FOR P4:
# All patient rows in df_raw_clean are MTLE patients.
# side coding: Left = 1, Right = 0.

# ==========================================================
# 1) Settings
# ==========================================================
.p3_base_output_dir <- "C:/Users/86150/Documents/HIPP_THAL"
output_dir <- file.path(.p3_base_output_dir, "structural_abnormality")

# Final six files only
.allowed_outputs <- c(
  # P3
  "P3_familyFDR_AllPAT_HS_Other_vs_HC.csv",
  "FigP3_heatmap_AllPAT_HS_Other_familyFDR_HedgesG_star.png",
  "P3_familyFDR_summary_counts.csv",
  # P4
  "P4_familyFDR_MTLE_LeftMTLE_RightMTLE_vs_HC.csv",
  "FigP4_heatmap_MTLE_LeftMTLE_RightMTLE_familyFDR_HedgesG_star.png",
  "P4_familyFDR_summary_counts.csv"
)

# Clean output folder first, so the final folder contains ONLY the six files above
if(dir.exists(output_dir)){
  unlink(
    list.files(output_dir, full.names=TRUE, all.files=TRUE, no..=TRUE),
    recursive=TRUE,
    force=TRUE
  )
} else {
  dir.create(output_dir, recursive=TRUE)
}

cat("\n[Output folder cleaned]\n", output_dir, "\n")
cat("[Expected final outputs]\n - ", paste(.allowed_outputs, collapse="\n - "), "\n", sep="")

set.seed(2026)

# Column names
col_id   <- "ID"
col_age  <- "Age"
col_sex  <- "sex"
col_etiv <- "eTIV"
col_side <- "side"
col_path_group <- "Pathology_Group"

# Permutation
n_perm <- 5000
perm_seed <- 999

# Minimum subgroup size
min_n_group <- 8
print_every_roi <- 10

# Use the SAME covariates for P3 and P4.
# Do NOT include side as a covariate because HC has no meaningful seizure side,
# and side itself defines Left/Right MTLE subgroups.
covars_main <- c(col_etiv, col_age, col_sex)

# Pathology coding: Pathology_Group == 1 -> HS; all others -> nonHS/Other
make_pathology_class <- function(x){
  out <- dplyr::case_when(
    x == 1 ~ "G1_HS",
    TRUE   ~ "G2_nonHS"
  )
  factor(out, levels=c("G1_HS","G2_nonHS"))
}

# ==========================================================
# 2) ROI list
# ==========================================================
roi_list_user <- c(
  # ---------------- Cortex N=68 ----------------
  "lh_lh_cuneus_volume","lh_lh_entorhinal_volume","lh_lh_fusiform_volume", "lh_lh_caudalmiddlefrontal_volume", "lh_lh_inferiorparietal_volume",
  "lh_lh_inferiortemporal_volume","lh_lh_isthmuscingulate_volume","lh_lh_lateraloccipital_volume",
  "lh_lh_lateralorbitofrontal_volume","lh_lh_lingual_volume","lh_lh_medialorbitofrontal_volume",
  "lh_lh_middletemporal_volume","lh_lh_parahippocampal_volume","lh_lh_paracentral_volume",
  "lh_lh_parsopercularis_volume","lh_lh_parsorbitalis_volume","lh_lh_parstriangularis_volume",
  "lh_lh_pericalcarine_volume","lh_lh_postcentral_volume","lh_lh_posteriorcingulate_volume",
  "lh_lh_precentral_volume","lh_lh_precuneus_volume","lh_lh_rostralanteriorcingulate_volume",
  "lh_lh_rostralmiddlefrontal_volume","lh_lh_superiorfrontal_volume","lh_lh_superiorparietal_volume",
  "lh_lh_superiortemporal_volume","lh_lh_supramarginal_volume","lh_lh_frontalpole_volume",
  "lh_lh_temporalpole_volume","lh_lh_transversetemporal_volume","lh_lh_insula_volume",
  "lh_lh_bankssts_volume", "lh_lh_caudalanteriorcingulate_volume",
  
  "rh_rh_bankssts_volume","rh_rh_caudalanteriorcingulate_volume","rh_rh_caudalmiddlefrontal_volume",
  "rh_rh_cuneus_volume","rh_rh_entorhinal_volume","rh_rh_fusiform_volume","rh_rh_inferiorparietal_volume",
  "rh_rh_inferiortemporal_volume","rh_rh_isthmuscingulate_volume","rh_rh_lateraloccipital_volume",
  "rh_rh_lateralorbitofrontal_volume","rh_rh_lingual_volume","rh_rh_medialorbitofrontal_volume",
  "rh_rh_middletemporal_volume","rh_rh_parahippocampal_volume","rh_rh_paracentral_volume",
  "rh_rh_parsopercularis_volume","rh_rh_parsorbitalis_volume","rh_rh_parstriangularis_volume",
  "rh_rh_pericalcarine_volume","rh_rh_postcentral_volume","rh_rh_posteriorcingulate_volume",
  "rh_rh_precentral_volume","rh_rh_precuneus_volume","rh_rh_rostralanteriorcingulate_volume",
  "rh_rh_rostralmiddlefrontal_volume","rh_rh_superiorfrontal_volume","rh_rh_superiorparietal_volume",
  "rh_rh_superiortemporal_volume","rh_rh_supramarginal_volume","rh_rh_frontalpole_volume",
  "rh_rh_temporalpole_volume","rh_rh_transversetemporal_volume","rh_rh_insula_volume",
  
  
  # ---------------- Hippocampus subfields N= 22----------------
  "Left_Hippocampal_tail","Left_subiculum_comb","Left_CA1_comb","Left_CA2/3_comb","Left_CA4_comb",
  "Left_GC_ML_DG_comb","Left_molecular_layer_HP_comb","Left_HATA","Left_fimbria","Left_presubiculum_comb",
  "Left_parasubiculum",
  
  "Right_Hippocampal_tail","Right_subiculum_comb","Right_CA1_comb","Right_CA2/3_comb","Right_CA4_comb",
  "Right_GC_ML_DG_comb","Right_molecular_layer_HP_comb","Right_HATA","Right_fimbria","Right_presubiculum_comb",
  "Right_parasubiculum",
  
  # ---------------- Thalamus nuclei N= 44 ----------------
  "Left_AV","Left_VA","Left_VAmc","Left_VLa","Left_VLp","Left_VPL","Left_VM","Left_CL","Left_CeM","Left_CM","Left_Pf",
  "Left_MDm","Left_MDl","Left_LD","Left_LP","Left_Pu_Total","Left_MGN","Left_LGN","Left_L_Sg",
  "Left_MV_Re","Left_Pc","Left_Pt",
  
  "Right_AV","Right_VA","Right_VAmc","Right_VLa","Right_VLp","Right_VPL","Right_VM","Right_CL","Right_CeM","Right_CM","Right_Pf",
  "Right_MDm","Right_MDl","Right_LD","Right_LP","Right_Pu_Total","Right_MGN","Right_LGN","Right_L_Sg",
  "Right_MV_Re","Right_Pc","Right_Pt",
  
  # ===================== NEW: Subcortical N=16 (aseg-style summary) =====================
  "Left-Subcort-Thalamus",
  "Left-Subcort-Caudate",
  "Left-Subcort-Putamen",
  "Left-Subcort-Pallidum",
  "Left-Subcort-Hippocampus",
  "Left-Subcort-Amygdala",
  "Left-Subcort-Accumbens-area",
  "Left-Subcort-VentralDC",
  "Right-Subcort-Thalamus",
  "Right-Subcort-Caudate",
  "Right-Subcort-Putamen",
  "Right-Subcort-Pallidum",
  "Right-Subcort-Hippocampus",
  "Right-Subcort-Amygdala",
  "Right-Subcort-Accumbens-area",
  "Right-Subcort-VentralDC"
)

# 必须剔除：hippocampal fissure（CSF，不是脑组织） + hemi-eTIV（避免重复混入ROI）
roi_list_user <- roi_list_user[!grepl("hippocampal_fissure", roi_list_user, ignore.case = TRUE)]
roi_list_user <- roi_list_user[!(roi_list_user %in% c("lh_eTIV","rh_eTIV"))]

# ==========================================================
# 3) Prepare merged dataset
# ==========================================================
if(!(col_id %in% colnames(df_pat))) df_pat[[col_id]] <- paste0("PAT_", seq_len(nrow(df_pat)))
if(!(col_id %in% colnames(df_hc0))) df_hc0[[col_id]] <- paste0("HC_", seq_len(nrow(df_hc0)))

# Required covariates
for(cc in covars_main){
  if(!(cc %in% colnames(df_pat))) stop("❌ Missing patient covariate: ", cc)
  if(!(cc %in% colnames(df_hc0))) stop("❌ Missing HC covariate: ", cc)
}

# Factorize sex
if(col_sex %in% colnames(df_pat)) df_pat[[col_sex]] <- as.factor(df_pat[[col_sex]])
if(col_sex %in% colnames(df_hc0)) df_hc0[[col_sex]] <- as.factor(df_hc0[[col_sex]])

# Need side in patients for Left/Right MTLE
if(!(col_side %in% colnames(df_pat))){
  stop("❌ Patients df has no 'side' column; cannot run Left/Right MTLE analysis.")
}
if(!(col_side %in% colnames(df_hc0))){
  df_hc0[[col_side]] <- NA
}

# Pathology group
if(col_path_group %in% colnames(df_pat)){
  df_pat$Group4 <- as.character(make_pathology_class(df_pat[[col_path_group]]))
} else {
  stop("❌ Pathology_Group not found in df_raw_clean; P3 HS/Other analysis cannot be created.")
}

df_hc0$Group4 <- "HC"
df_pat$Group2 <- "PAT"
df_hc0$Group2 <- "HC"

to_int_safe <- function(x){
  if(is.factor(x)) x <- as.character(x)
  suppressWarnings(as.integer(x))
}

df_all <- bind_rows(df_hc0, df_pat) %>%
  mutate(
    Group2 = factor(Group2, levels=c("HC","PAT")),
    Group4 = factor(Group4, levels=c("HC","G1_HS","G2_nonHS")),
    side_num = to_int_safe(.data[[col_side]]),
    Laterality = case_when(
      Group2 == "HC" ~ "HC",
      is.finite(side_num) & side_num == 1 ~ "LeftOnset",
      is.finite(side_num) & side_num == 0 ~ "RightOnset",
      TRUE ~ NA_character_
    ),
    Laterality = factor(Laterality, levels=c("HC","LeftOnset","RightOnset"))
  )

cat("\n[Sanity check] Pathology groups:\n")
print(table(df_all$Group4, useNA="ifany"))
cat("\n[Sanity check] MTLE laterality groups:\n")
print(table(df_all$Laterality, useNA="ifany"))

if(sum(df_all$Laterality=="HC", na.rm=TRUE) < 5){
  stop("❌ Too few HC after merge.")
}
if(sum(df_all$Laterality=="LeftOnset", na.rm=TRUE) < min_n_group){
  stop("❌ Too few Left MTLE patients. Check side coding: Left should be 1.")
}
if(sum(df_all$Laterality=="RightOnset", na.rm=TRUE) < min_n_group){
  stop("❌ Too few Right MTLE patients. Check side coding: Right should be 0.")
}

# ROI columns available in BOTH patient and HC datasets
roi_cols <- intersect(roi_list_user, intersect(colnames(df_pat), colnames(df_hc0)))
if(length(roi_cols) < 10) stop("❌ ROI matched too few. Check ROI column names.")
cat("\nROI used:", length(roi_cols), "\n")

# ==========================================================
# 4) Statistical helper functions
# ==========================================================
hedges_g <- function(x, y){
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  nx <- length(x); ny <- length(y)
  if(nx < 5 || ny < 5) return(NA_real_)
  mx <- mean(x); my <- mean(y)
  sx <- var(x); sy <- var(y)
  sp <- sqrt(((nx-1)*sx + (ny-1)*sy) / (nx+ny-2))
  if(!is.finite(sp) || sp < 1e-12) return(NA_real_)
  d <- (mx - my)/sp
  J <- 1 - (3/(4*(nx+ny)-9))
  J*d
}

freedman_lane_fast <- function(d, y_col, covars, group_var, g_level, ref_level="HC",
                               B=2000, seed=999, strata_var=NULL){
  set.seed(seed)
  
  need <- unique(c(y_col, covars, group_var, strata_var))
  need <- need[need %in% colnames(d)]
  d2 <- d %>% select(all_of(need)) %>% drop_na()
  if(nrow(d2) < 20) return(list(p_perm=NA_real_, t_obs=NA_real_, n=nrow(d2)))
  
  d2[[group_var]] <- factor(d2[[group_var]])
  if(!(ref_level %in% levels(d2[[group_var]]))) return(list(p_perm=NA_real_, t_obs=NA_real_, n=nrow(d2)))
  if(!(g_level %in% levels(d2[[group_var]]))) return(list(p_perm=NA_real_, t_obs=NA_real_, n=nrow(d2)))
  
  d2 <- d2 %>% filter(.data[[group_var]] %in% c(ref_level, g_level))
  d2[[group_var]] <- relevel(droplevels(d2[[group_var]]), ref=ref_level)
  
  n_ref <- sum(d2[[group_var]]==ref_level)
  n_g   <- sum(d2[[group_var]]==g_level)
  if(n_ref < 5 || n_g < 5) return(list(p_perm=NA_real_, t_obs=NA_real_, n=nrow(d2)))
  
  y <- as.numeric(d2[[y_col]])
  
  if(length(covars) > 0){
    X0 <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse=" + "))), data=d2)
  } else {
    X0 <- matrix(1, nrow=nrow(d2), ncol=1)
    colnames(X0) <- "(Intercept)"
  }
  
  X1 <- model.matrix(as.formula(paste0("~ ",
                                       if(length(covars)>0) paste(covars, collapse=" + ") else "1",
                                       " + ", group_var)), data=d2)
  
  cn1 <- colnames(X1)
  idx_g <- which(grepl(paste0("^", group_var), cn1) & grepl(g_level, cn1))
  if(length(idx_g) != 1) return(list(p_perm=NA_real_, t_obs=NA_real_, n=nrow(d2)))
  
  fit0 <- lm.fit(X0, y)
  yhat <- as.numeric(X0 %*% fit0$coefficients)
  res  <- as.numeric(y - yhat)
  
  fit1 <- lm.fit(X1, y)
  r1   <- y - as.numeric(X1 %*% fit1$coefficients)
  df1  <- nrow(X1) - ncol(X1)
  s2   <- sum(r1^2) / df1
  XtX1_inv <- solve(crossprod(X1))
  se_g <- sqrt(s2 * XtX1_inv[idx_g, idx_g])
  t_obs <- fit1$coefficients[idx_g] / se_g
  
  if(!is.null(strata_var) && strata_var %in% colnames(d2)){
    strata <- as.character(d2[[strata_var]])
  } else {
    strata <- rep("all", nrow(d2))
  }
  strata_levels <- unique(strata)
  strata_idx <- lapply(strata_levels, function(st) which(strata == st))
  
  perm_t <- rep(NA_real_, B)
  
  for(b in 1:B){
    res_p <- res
    for(ii in strata_idx){
      if(length(ii) > 1) res_p[ii] <- res[sample(ii, length(ii), replace=FALSE)]
    }
    y_perm <- yhat + res_p
    
    fitb <- lm.fit(X1, y_perm)
    rb   <- y_perm - as.numeric(X1 %*% fitb$coefficients)
    s2b  <- sum(rb^2) / df1
    se_gb <- sqrt(s2b * XtX1_inv[idx_g, idx_g])
    perm_t[b] <- fitb$coefficients[idx_g] / se_gb
  }
  
  perm_t2 <- perm_t[is.finite(perm_t)]
  if(length(perm_t2) < max(200, floor(0.2*B))) return(list(p_perm=NA_real_, t_obs=t_obs, n=nrow(d2)))
  
  p_perm <- (1 + sum(abs(perm_t2) >= abs(t_obs))) / (1 + length(perm_t2))
  list(p_perm=p_perm, t_obs=t_obs, n=nrow(d2))
}

lm_contrast_fast <- function(d, y_col, covars, group_var, g_level, ref_level="HC"){
  need <- unique(c(y_col, covars, group_var))
  need <- need[need %in% colnames(d)]
  d2 <- d %>% select(all_of(need)) %>% drop_na()
  if(nrow(d2) < 20) return(NULL)
  
  d2[[group_var]] <- factor(d2[[group_var]])
  if(!(ref_level %in% levels(d2[[group_var]]))) return(NULL)
  if(!(g_level %in% levels(d2[[group_var]]))) return(NULL)
  
  d2 <- d2 %>% filter(.data[[group_var]] %in% c(ref_level, g_level))
  d2[[group_var]] <- relevel(droplevels(d2[[group_var]]), ref=ref_level)
  
  n_ref <- sum(d2[[group_var]]==ref_level)
  n_g   <- sum(d2[[group_var]]==g_level)
  if(n_ref < 5 || n_g < 5) return(NULL)
  
  if(length(covars)>0){
    X <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse=" + "), " + ", group_var)), data=d2)
  } else {
    X <- model.matrix(as.formula(paste0("~ ", group_var)), data=d2)
  }
  y <- as.numeric(d2[[y_col]])
  
  fit <- lm.fit(X, y)
  r   <- y - as.numeric(X %*% fit$coefficients)
  df  <- nrow(X) - ncol(X)
  s2  <- sum(r^2) / df
  XtX_inv <- solve(crossprod(X))
  se  <- sqrt(diag(s2 * XtX_inv))
  
  cn <- colnames(X)
  idx <- which(grepl(paste0("^", group_var), cn) & grepl(g_level, cn))
  if(length(idx)!=1) return(NULL)
  
  tval <- fit$coefficients[idx] / se[idx]
  pval <- 2 * pt(-abs(tval), df=df)
  
  tibble(
    beta_lm = fit$coefficients[idx],
    se_lm   = se[idx],
    t_lm    = tval,
    p_lm    = pval,
    N       = nrow(d2),
    N_ref   = n_ref,
    N_g     = n_g
  )
}

residualize_vec <- function(v, d, covars){
  d2 <- d
  d2$.v <- as.numeric(v)
  if(length(covars)==0){
    return(d2$.v - mean(d2$.v, na.rm=TRUE))
  }
  X <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse=" + "))), data=d2)
  fit <- lm.fit(X, d2$.v)
  as.numeric(d2$.v - as.numeric(X %*% fit$coefficients))
}

# Generic one-contrast runner
run_one_contrast <- function(d, roi, group_var, g_level, ref_level="HC",
                             covars=covars_main, B=n_perm, seed=999,
                             strata_var=col_sex){
  need <- unique(c(roi, covars, group_var))
  d_tmp <- d %>% select(all_of(need)) %>% drop_na()
  
  if(nrow(d_tmp) < 20) return(NULL)
  if(sum(d_tmp[[group_var]] == g_level) < min_n_group) return(NULL)
  if(sum(d_tmp[[group_var]] == ref_level) < 5) return(NULL)
  
  out_lm <- lm_contrast_fast(
    d_tmp,
    y_col=roi,
    covars=covars,
    group_var=group_var,
    g_level=g_level,
    ref_level=ref_level
  )
  if(is.null(out_lm)) return(NULL)
  
  y_res <- residualize_vec(d_tmp[[roi]], d_tmp, covars)
  g <- hedges_g(
    y_res[d_tmp[[group_var]] == g_level],
    y_res[d_tmp[[group_var]] == ref_level]
  )
  
  perm <- freedman_lane_fast(
    d_tmp,
    y_col=roi,
    covars=covars,
    group_var=group_var,
    g_level=g_level,
    ref_level=ref_level,
    B=B,
    seed=seed,
    strata_var=strata_var
  )
  
  tibble(
    ROI=roi,
    N=out_lm$N,
    N_HC=out_lm$N_ref,
    N_Group=out_lm$N_g,
    beta=out_lm$beta_lm,
    se=out_lm$se_lm,
    t=out_lm$t_lm,
    p_lm=out_lm$p_lm,
    p_perm=perm$p_perm,
    Hedges_g_adj=g
  )
}

# ==========================================================
# 5) Overall MTLE/PAT vs HC
#    Reused by BOTH P3 and P4
# ==========================================================
cat("\n=== Overall MTLE/PAT vs HC ===\n")

res_all_list <- vector("list", length(roi_cols))
for(i in seq_along(roi_cols)){
  roi <- roi_cols[i]
  
  out <- run_one_contrast(
    d=df_all,
    roi=roi,
    group_var="Group2",
    g_level="PAT",
    ref_level="HC",
    covars=covars_main,
    B=n_perm,
    seed=perm_seed + 10000 + i,
    strata_var=col_sex
  )
  
  if(is.null(out)) next
  res_all_list[[i]] <- out
  
  if(i %% print_every_roi == 0){
    cat("Overall progress:", i, "/", length(roi_cols), "ROIs\n")
  }
}

res_all <- bind_rows(res_all_list)
if(nrow(res_all)==0) stop("❌ Overall MTLE/PAT vs HC results are empty.")

res_all <- res_all %>%
  mutate(
    Contrast="AllPAT_vs_HC",
    q_lm=p.adjust(p_lm, method="BH"),
    q_perm=p.adjust(p_perm, method="BH")
  ) %>%
  select(Contrast, everything())

# ==========================================================
# 6) Pathology: HS vs HC and Other/nonHS vs HC
#    Used by P3
# ==========================================================
cat("\n=== Pathology: HS / Other vs HC ===\n")

path_levels <- c("G1_HS","G2_nonHS")
res_path_all <- list()

for(glev in path_levels){
  cat("Pathology subgroup:", glev, "\n")
  res_one <- vector("list", length(roi_cols))
  
  d_g <- df_all %>%
    filter(Group4 %in% c("HC", glev)) %>%
    droplevels()
  
  for(i in seq_along(roi_cols)){
    roi <- roi_cols[i]
    
    out <- run_one_contrast(
      d=d_g,
      roi=roi,
      group_var="Group4",
      g_level=glev,
      ref_level="HC",
      covars=covars_main,
      B=n_perm,
      seed=perm_seed + i + 1000*match(glev, path_levels),
      strata_var=col_sex
    )
    
    if(is.null(out)) next
    res_one[[i]] <- out %>% mutate(Contrast=paste0(glev, "_vs_HC"))
    
    if(i %% print_every_roi == 0){
      cat("Pathology progress:", glev, i, "/", length(roi_cols), "ROIs\n")
    }
  }
  
  res_g <- bind_rows(res_one)
  if(nrow(res_g) > 0){
    res_g <- res_g %>%
      mutate(
        q_lm=p.adjust(p_lm, method="BH"),
        q_perm=p.adjust(p_perm, method="BH")
      )
  }
  res_path_all[[glev]] <- res_g
}

res_path <- bind_rows(res_path_all)
if(nrow(res_path)==0) stop("❌ HS/Other vs HC results are empty.")

# ==========================================================
# 7) Laterality: Left MTLE vs HC and Right MTLE vs HC
#    Used by P4
#    NOTE: results are calculated but NOT separately saved.
# ==========================================================
cat("\n=== Laterality: Left MTLE / Right MTLE vs HC ===\n")

lat_levels <- c("LeftOnset","RightOnset")
res_lat_all <- list()

for(glev in lat_levels){
  cat("Laterality subgroup:", glev, "\n")
  res_one <- vector("list", length(roi_cols))
  
  d_g <- df_all %>%
    filter(Laterality %in% c("HC", glev)) %>%
    droplevels()
  
  for(i in seq_along(roi_cols)){
    roi <- roi_cols[i]
    
    out <- run_one_contrast(
      d=d_g,
      roi=roi,
      group_var="Laterality",
      g_level=glev,
      ref_level="HC",
      covars=covars_main,
      B=n_perm,
      seed=perm_seed + 20000 + i + 5000*match(glev, lat_levels),
      strata_var=col_sex
    )
    
    if(is.null(out)) next
    res_one[[i]] <- out %>% mutate(Contrast=paste0(glev, "_vs_HC"))
    
    if(i %% print_every_roi == 0){
      cat("Laterality progress:", glev, i, "/", length(roi_cols), "ROIs\n")
    }
  }
  
  res_g <- bind_rows(res_one)
  if(nrow(res_g) > 0){
    res_g <- res_g %>%
      mutate(
        q_lm=p.adjust(p_lm, method="BH"),
        q_perm=p.adjust(p_perm, method="BH")
      )
  }
  res_lat_all[[glev]] <- res_g
}

res_lat <- bind_rows(res_lat_all)
if(nrow(res_lat)==0) stop("❌ Left/Right MTLE vs HC results are empty.")


# ==========================================================
# 8) Correct ROI family definition
#    IMPORTANT: hippocampal subfields now use Left_/Right_ names,
#    so they MUST be identified BEFORE generic Left_/Right_ thalamus rules.
# ==========================================================

hippocampal_names <- c(
  "Hippocampal_tail",
  "subiculum_comb",
  "CA1_comb",
  "CA2/3_comb",
  "CA4_comb",
  "GC_ML_DG_comb",
  "molecular_layer_HP_comb",
  "HATA",
  "fimbria",
  "presubiculum_comb",
  "parasubiculum"
)


roi_family <- function(roi){
  
  r <- as.character(roi)
  
  # ----------------------------------------------------------
  # Subcortical aseg
  # ----------------------------------------------------------
  if(
    grepl("^Left-Subcort-", r) ||
    grepl("^Right-Subcort-", r)
  ){
    return("SubcorticalAseg")
  }
  
  # ----------------------------------------------------------
  # Global measures
  # ----------------------------------------------------------
  if(
    r %in% c(
      "lh_BrainSegVolNotVent",
      "rh_BrainSegVolNotVent",
      "BrainSegVolNotVent"
    )
  ){
    return("GlobalOther")
  }
  
  # ----------------------------------------------------------
  # Cortex
  # ----------------------------------------------------------
  if(
    grepl("^lh_lh_", r) ||
    grepl("^rh_rh_", r)
  ){
    return("Cortex")
  }
  
  # ----------------------------------------------------------
  # Hippocampal subfields
  # IMPORTANT:
  # identify these BEFORE generic Left_/Right_ thalamus rules
  # ----------------------------------------------------------
  if(
    r %in% c(
      paste0("Left_", hippocampal_names),
      paste0("Right_", hippocampal_names)
    )
  ){
    return("HippocampusSubfields")
  }
  
  # ----------------------------------------------------------
  # Thalamic nuclei
  # ----------------------------------------------------------
  if(
    grepl("^Left_", r) ||
    grepl("^Right_", r)
  ){
    return("ThalamusNuclei")
  }
  
  # ----------------------------------------------------------
  # Everything else
  # ----------------------------------------------------------
  return("GlobalOther")
}


# ==========================================================
# Family order
# ==========================================================

family_levels <- c(
  "SubcorticalAseg",
  "HippocampusSubfields",
  "ThalamusNuclei",
  "Cortex",
  "GlobalOther"
)


# ==========================================================
# Family-wise FDR
# ==========================================================

add_family_fdr <- function(df){
  
  df %>%
    
    mutate(
      
      Family = vapply(
        ROI,
        roi_family,
        character(1)
      ),
      
      Family = factor(
        Family,
        levels = family_levels
      )
      
    ) %>%
    
    group_by(
      Contrast,
      Family
    ) %>%
    
    mutate(
      
      q_perm_family = p.adjust(
        p_perm,
        method = "BH"
      )
      
    ) %>%
    
    ungroup() %>%
    
    mutate(
      
      star = case_when(
        
        is.finite(q_perm_family) &
          q_perm_family < 0.001 ~ "***",
        
        is.finite(q_perm_family) &
          q_perm_family < 0.01 ~ "**",
        
        is.finite(q_perm_family) &
          q_perm_family < 0.05 ~ "*",
        
        is.finite(q_perm_family) &
          q_perm_family < 0.10 ~ "·",
        
        TRUE ~ ""
        
      )
      
    )
}


# ==========================================================
# Heatmap function
#
# Figure style:
#   - Times New Roman
#   - ALL text bold
#   - ALL text black
#   - significance symbols bold + black
#   - significance symbols vertically adjusted to visual center
# ==========================================================

make_family_heatmap <- function(
    df,
    contrast_levels,
    contrast_labels,
    title_text,
    filename,
    alpha_show = 0.10,
    topN = 70
){
  
  # ----------------------------------------------------------
  # Contrast order
  # ----------------------------------------------------------
  
  df <- df %>%
    
    mutate(
      
      Contrast = factor(
        as.character(Contrast),
        levels = contrast_levels
      )
      
    )
  
  
  # ----------------------------------------------------------
  # Rank ROIs according to minimum family-wise q value
  # ----------------------------------------------------------
  
  roi_rank <- df %>%
    
    group_by(
      ROI,
      Family
    ) %>%
    
    summarise(
      
      min_qfam = suppressWarnings(
        min(
          q_perm_family,
          na.rm = TRUE
        )
      ),
      
      .groups = "drop"
      
    ) %>%
    
    filter(
      is.finite(min_qfam)
    ) %>%
    
    arrange(
      min_qfam
    )
  
  
  # ----------------------------------------------------------
  # Select ROIs for heatmap
  # ----------------------------------------------------------
  
  top_roi <- roi_rank %>%
    
    filter(
      min_qfam < alpha_show
    ) %>%
    
    slice_head(
      n = topN
    ) %>%
    
    pull(
      ROI
    )
  
  
  # ----------------------------------------------------------
  # If no ROI reaches q_family < alpha_show:
  # show the 30 lowest-q ROIs instead.
  # Stars still reflect the real q values.
  # ----------------------------------------------------------
  
  if(length(top_roi) == 0){
    
    message(
      "⚠️ No ROI with q_perm_family < ",
      alpha_show,
      "; showing the 30 lowest-q ROIs instead."
    )
    
    top_roi <- roi_rank %>%
      
      slice_head(
        n = min(
          30,
          nrow(roi_rank)
        )
      ) %>%
      
      pull(
        ROI
      )
  }
  
  
  # ----------------------------------------------------------
  # Heatmap dataframe
  # ----------------------------------------------------------
  
  hm <- df %>%
    
    filter(
      ROI %in% top_roi
    ) %>%
    
    mutate(
      
      ROI = factor(
        ROI,
        levels = rev(top_roi)
      )
      
    )
  
  
  # ----------------------------------------------------------
  # Color limit
  # ----------------------------------------------------------
  
  lim <- suppressWarnings(
    
    quantile(
      abs(hm$Hedges_g_adj),
      0.95,
      na.rm = TRUE
    )
    
  )
  
  
  if(
    !is.finite(lim) ||
    lim < 0.5
  ){
    lim <- 0.5
  }
  
  
  # ----------------------------------------------------------
  # Draw heatmap
  # ----------------------------------------------------------
  
  p <- ggplot(
    hm,
    aes(
      x = Contrast,
      y = ROI,
      fill = Hedges_g_adj
    )
  ) +
    
    # --------------------------------------------------------
  # Heatmap tiles
  # --------------------------------------------------------
  
  geom_tile(
    color = "white",
    linewidth = 0.25
  ) +
    
    # --------------------------------------------------------
  # Significance symbols
  #
  # hjust = 0.5 -> horizontal center
  # vjust = 0.72 -> visually move * slightly downward
  #                so it appears centered in each tile
  # --------------------------------------------------------
  
  geom_text(
    aes(
      label = star
    ),
    family = "Times New Roman",
    fontface = "bold",
    colour = "black",
    size = 5,
    hjust = 0.5,
    vjust = 0.72
  ) +
    
    # --------------------------------------------------------
  # ROI family panels
  # --------------------------------------------------------
  
  facet_grid(
    Family ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
    
    # --------------------------------------------------------
  # X-axis labels
  # --------------------------------------------------------
  
  scale_x_discrete(
    labels = contrast_labels
  ) +
    
    # --------------------------------------------------------
  # Effect-size color scale
  # --------------------------------------------------------
  
  scale_fill_gradient2(
    
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    
    limits = c(
      -lim,
      lim
    ),
    
    oob = scales::squish,
    
    na.value = "grey90"
  ) +
    
    # --------------------------------------------------------
  # Base theme
  # --------------------------------------------------------
  
  theme_minimal(
    base_family = "Times New Roman",
    base_size = 12
  ) +
    
    # --------------------------------------------------------
  # Font style
  # ALL text:
  #   Times New Roman
  #   bold
  #   black
  # --------------------------------------------------------
  
  theme(
    
    # Global text setting
    text = element_text(
      family = "Times New Roman",
      face = "bold",
      colour = "black"
    ),
    
    
    # ------------------------------------------------------
    # X-axis text
    # e.g. All PAT / HS / Other
    # or MTLE / Left MTLE / Right MTLE
    # ------------------------------------------------------
    
    axis.text.x = element_text(
      family = "Times New Roman",
      size = 12,
      face = "bold",
      colour = "black",
      angle = 25,
      hjust = 1,
      vjust = 1
    ),
    
    
    # ------------------------------------------------------
    # Y-axis ROI names
    # ------------------------------------------------------
    
    axis.text.y = element_text(
      family = "Times New Roman",
      size = 11,
      face = "bold",
      colour = "black"
    ),
    
    
    # ------------------------------------------------------
    # Axis titles
    # ------------------------------------------------------
    
    axis.title.x = element_text(
      family = "Times New Roman",
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    axis.title.y = element_text(
      family = "Times New Roman",
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    
    # ------------------------------------------------------
    # Family labels
    # e.g.
    # SubcorticalAseg
    # HippocampusSubfields
    # ThalamusNuclei
    # Cortex
    # ------------------------------------------------------
    
    strip.text.y = element_text(
      family = "Times New Roman",
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    
    # ------------------------------------------------------
    # Main title
    # ------------------------------------------------------
    
    plot.title = element_text(
      family = "Times New Roman",
      size = 14,
      face = "bold",
      colour = "black",
      hjust = 0
    ),
    
    
    # ------------------------------------------------------
    # Subtitle
    # ------------------------------------------------------
    
    plot.subtitle = element_text(
      family = "Times New Roman",
      size = 11,
      face = "bold",
      colour = "black",
      hjust = 0
    ),
    
    
    # ------------------------------------------------------
    # Legend title
    # ------------------------------------------------------
    
    legend.title = element_text(
      family = "Times New Roman",
      size = 11,
      face = "bold",
      colour = "black"
    ),
    
    
    # ------------------------------------------------------
    # Legend values
    # ------------------------------------------------------
    
    legend.text = element_text(
      family = "Times New Roman",
      size = 10,
      face = "bold",
      colour = "black"
    ),
    
    
    # ------------------------------------------------------
    # Optional:
    # keep panel/grid appearance clean
    # ------------------------------------------------------
    
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 15,
      b = 10,
      l = 10
    )
    
  ) +
    
    # --------------------------------------------------------
  # Labels
  # --------------------------------------------------------
  
  labs(
    
    title = title_text,
    
    subtitle = paste0(
      "Fill = Hedges' g (covariate-adjusted). ",
      "Stars = family-wise BH-FDR: ",
      "*<0.05, **<0.01, ***<0.001, ·<0.10"
    ),
    
    x = "",
    
    y = "ROI",
    
    fill = "Hedges' g"
  )
  
  
  # ----------------------------------------------------------
  # Save figure
  # ----------------------------------------------------------
  
  ggplot2::ggsave(
    
    filename = file.path(
      output_dir,
      filename
    ),
    
    plot = p,
    
    width = 7,
    
    height = max(
      10,
      length(top_roi) * 0.22
    ),
    
    dpi = 600
  )
}


# ==========================================================
# Summary function
# ==========================================================

make_summary <- function(df){
  
  df %>%
    
    mutate(
      
      sig05 =
        is.finite(q_perm_family) &
        q_perm_family < 0.05
      
    ) %>%
    
    group_by(
      Contrast,
      Family
    ) %>%
    
    summarise(
      
      n_test = sum(
        is.finite(p_perm)
      ),
      
      n_sig05 = sum(
        sig05,
        na.rm = TRUE
      ),
      
      .groups = "drop"
      
    ) %>%
    
    arrange(
      Contrast,
      Family
    )
}


# ==========================================================
# 9) P3: All PAT / HS / Other vs HC
# ==========================================================

cat(
  "\n================ P3: All PAT / HS / Other vs HC ================\n"
)


P3 <- bind_rows(
  
  # ----------------------------------------------------------
  # All patients vs HC
  # ----------------------------------------------------------
  
  res_all %>%
    
    transmute(
      
      Contrast = "AllPAT_vs_HC",
      
      ROI,
      
      p_perm,
      
      Hedges_g_adj
      
    ),
  
  
  # ----------------------------------------------------------
  # HS / Other vs HC
  # ----------------------------------------------------------
  
  res_path %>%
    
    transmute(
      
      Contrast = case_when(
        
        grepl(
          "^G1_HS",
          Contrast
        ) ~ "HS_vs_HC",
        
        grepl(
          "^G2_nonHS",
          Contrast
        ) ~ "Other_vs_HC",
        
        TRUE ~ NA_character_
        
      ),
      
      ROI,
      
      p_perm,
      
      Hedges_g_adj
      
    ) %>%
    
    filter(
      !is.na(Contrast)
    )
  
) %>%
  
  mutate(
    
    Contrast = factor(
      
      Contrast,
      
      levels = c(
        "AllPAT_vs_HC",
        "HS_vs_HC",
        "Other_vs_HC"
      )
      
    )
    
  ) %>%
  
  add_family_fdr()


# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

cat(
  "[P3 contrasts]\n"
)

print(
  table(
    P3$Contrast,
    useNA = "ifany"
  )
)


cat(
  "[P3 ROI families]\n"
)

print(
  table(
    P3$Family,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# Save full P3 results
# ------------------------------------------------------------

utils::write.csv(
  
  P3,
  
  file = file.path(
    output_dir,
    "P3_familyFDR_AllPAT_HS_Other_vs_HC.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "GBK"
)


# ------------------------------------------------------------
# P3 heatmap
# ------------------------------------------------------------

make_family_heatmap(
  
  df = P3,
  
  contrast_levels = c(
    "AllPAT_vs_HC",
    "HS_vs_HC",
    "Other_vs_HC"
  ),
  
  contrast_labels = c(
    "AllPAT_vs_HC" = "All PAT",
    "HS_vs_HC" = "HS",
    "Other_vs_HC" = "Other"
  ),
  
  title_text =
    "All PAT / HS / Other vs HC (Family-wise FDR)",
  
  filename =
    "FigP3_heatmap_AllPAT_HS_Other_familyFDR_HedgesG_star.png"
)


# ------------------------------------------------------------
# P3 summary
# ------------------------------------------------------------

P3_summary <- make_summary(
  P3
)


utils::write.csv(
  
  P3_summary,
  
  file = file.path(
    output_dir,
    "P3_familyFDR_summary_counts.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "GBK"
)


cat(
  "✅ P3 saved.\n"
)


# ==========================================================
# 10) P4: MTLE / Left MTLE / Right MTLE vs HC
# ==========================================================

cat(
  "\n================ P4: MTLE / Left MTLE / Right MTLE vs HC ================\n"
)


P4 <- bind_rows(
  
  # ----------------------------------------------------------
  # All MTLE vs HC
  # ----------------------------------------------------------
  
  res_all %>%
    
    transmute(
      
      Contrast = "MTLE_vs_HC",
      
      ROI,
      
      p_perm,
      
      Hedges_g_adj
      
    ),
  
  
  # ----------------------------------------------------------
  # Left / Right MTLE vs HC
  # ----------------------------------------------------------
  
  res_lat %>%
    
    transmute(
      
      Contrast = case_when(
        
        grepl(
          "^LeftOnset",
          Contrast
        ) ~ "Left_MTLE_vs_HC",
        
        grepl(
          "^RightOnset",
          Contrast
        ) ~ "Right_MTLE_vs_HC",
        
        TRUE ~ NA_character_
        
      ),
      
      ROI,
      
      p_perm,
      
      Hedges_g_adj
      
    ) %>%
    
    filter(
      !is.na(Contrast)
    )
  
) %>%
  
  mutate(
    
    Contrast = factor(
      
      Contrast,
      
      levels = c(
        "MTLE_vs_HC",
        "Left_MTLE_vs_HC",
        "Right_MTLE_vs_HC"
      )
      
    )
    
  ) %>%
  
  add_family_fdr()


# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------

cat(
  "[P4 contrasts]\n"
)

print(
  table(
    P4$Contrast,
    useNA = "ifany"
  )
)


cat(
  "[P4 ROI families]\n"
)

print(
  table(
    P4$Family,
    useNA = "ifany"
  )
)


# ------------------------------------------------------------
# Save full P4 results
# ------------------------------------------------------------

utils::write.csv(
  
  P4,
  
  file = file.path(
    output_dir,
    "P4_familyFDR_MTLE_LeftMTLE_RightMTLE_vs_HC.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "GBK"
)


# ------------------------------------------------------------
# P4 heatmap
# ------------------------------------------------------------

make_family_heatmap(
  
  df = P4,
  
  contrast_levels = c(
    "MTLE_vs_HC",
    "Left_MTLE_vs_HC",
    "Right_MTLE_vs_HC"
  ),
  
  contrast_labels = c(
    "MTLE_vs_HC" = "MTLE",
    "Left_MTLE_vs_HC" = "Left MTLE",
    "Right_MTLE_vs_HC" = "Right MTLE"
  ),
  
  title_text =
    "MTLE / Left MTLE / Right MTLE vs HC (Family-wise FDR)",
  
  filename =
    "FigP4_heatmap_MTLE_LeftMTLE_RightMTLE_familyFDR_HedgesG_star.png"
)


# ------------------------------------------------------------
# P4 summary
# ------------------------------------------------------------

P4_summary <- make_summary(
  P4
)


utils::write.csv(
  
  P4_summary,
  
  file = file.path(
    output_dir,
    "P4_familyFDR_summary_counts.csv"
  ),
  
  row.names = FALSE,
  
  fileEncoding = "GBK"
)


cat(
  "✅ P4 saved.\n"
)


# ==========================================================
# 11) STRICT FINAL OUTPUT CHECK
#     Final folder must contain EXACTLY six files.
# ==========================================================

.final_files <- list.files(
  output_dir,
  all.files = FALSE,
  no.. = TRUE
)


.extra_files <- setdiff(
  .final_files,
  .allowed_outputs
)


# ----------------------------------------------------------
# Delete any non-target outputs
# ----------------------------------------------------------

if(length(.extra_files) > 0){
  
  unlink(
    
    file.path(
      output_dir,
      .extra_files
    ),
    
    recursive = TRUE,
    
    force = TRUE
  )
}


# ----------------------------------------------------------
# Re-check final folder
# ----------------------------------------------------------

.final_files <- list.files(
  output_dir,
  all.files = FALSE,
  no.. = TRUE
)


.missing_files <- setdiff(
  .allowed_outputs,
  .final_files
)


# ----------------------------------------------------------
# Stop if any requested file is missing
# ----------------------------------------------------------

if(length(.missing_files) > 0){
  
  stop(
    
    "❌ Analysis finished, but these required files were not generated:\n- ",
    
    paste(
      .missing_files,
      collapse = "\n- "
    )
  )
}


# ----------------------------------------------------------
# Final message
# ----------------------------------------------------------

cat(
  
  "\n✅ STRICT CHECK PASSED. The output folder contains exactly ",
  
  length(
    .allowed_outputs
  ),
  
  " files:\n - ",
  
  paste(
    sort(.allowed_outputs),
    collapse = "\n - "
  ),
  
  "\nFolder:\n",
  
  output_dir,
  
  "\n",
  
  sep = ""
)
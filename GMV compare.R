#这部分代码是两组的分析，然后脑区只有皮层+皮层下核团16个+丘脑亚区+海马亚区

suppressPackageStartupMessages({
  library(tidyverse)
})

# ==========================================================
# 0) Required input objects (ONLY these are assumed to exist)
# ==========================================================
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found (patients)")
if(!exists("df_hc_clean"))  stop("❌ df_hc_clean not found (healthy controls)")

df_pat <- df_raw_clean
df_hc0 <- df_hc_clean

# ==========================================================
# 1) Settings (EDIT HERE)
# ==========================================================
output_dir <- "C:/Users/86150/Documents/HIPP/02011/HC_Atrophy_AB_noEdu/new_2group"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

set.seed(2026)

# --- Column names ---
col_id   <- "ID"
col_age  <- "Age"
col_sex  <- "sex"
col_etiv <- "eTIV"
col_side <- "side"                  # may be missing in HC
col_path_group <- "Pathology_Group" # optional in patients

# --- Covariates (NO education) ---
covars_candidate <- c(col_etiv, col_age, col_sex, col_side)

# --- Permutation ---
n_perm <- 5000
perm_seed <- 999

# --- Minimum subgroup size ---
min_n_group <- 8

# --- Progress print ---
print_every_roi <- 10

# --- 分成2个亚组：G1_HS vs G2_Other（统一命名）---
make_pathology_class <- function(x){
  out <- dplyr::case_when(
    x == 1 ~ "G1_HS",
    TRUE   ~ "G2_Other"
  )
  factor(out, levels = c("G1_HS","G2_Other"))
}

# ==========================================================
# 2) ROI list (YOUR PROVIDED LIST) + NEW Subcortical16 + mandatory exclusions
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
  
  
  # ---------------- Hippocampus subfields N= 26----------------
  "lh_Hippocampal_tail","lh_subiculum_comb","lh_CA1_comb","lh_CA3_comb","lh_CA4_comb",
  "lh_GC_ML_DG_comb","lh_molecular_layer_HP_comb","lh_HATA","lh_fimbria","lh_presubiculum_comb",
  "lh_parasubiculum","lh_Whole_hippocampal_head","lh_Whole_hippocampal_body",
  
  "rh_Hippocampal_tail","rh_subiculum_comb","rh_CA1_comb","rh_CA3_comb","rh_CA4_comb",
  "rh_GC_ML_DG_comb","rh_molecular_layer_HP_comb","rh_HATA","rh_fimbria","rh_presubiculum_comb",
  "rh_parasubiculum","rh_Whole_hippocampal_head","rh_Whole_hippocampal_body",
  
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
# 3) Prepare merged dataset + robust covariate selection
# ==========================================================
# Ensure ID exists
if(!(col_id %in% colnames(df_pat))) df_pat[[col_id]] <- paste0("PAT_", seq_len(nrow(df_pat)))
if(!(col_id %in% colnames(df_hc0))) df_hc0[[col_id]] <- paste0("HC_", seq_len(nrow(df_hc0)))

# Factorize where possible
for(cc in c(col_sex, col_side)){
  if(cc %in% colnames(df_pat)) df_pat[[cc]] <- as.factor(df_pat[[cc]])
  if(cc %in% colnames(df_hc0)) df_hc0[[cc]] <- as.factor(df_hc0[[cc]])
}

# Build subgroup label for patients (2-subgroup)
if(col_path_group %in% colnames(df_pat)){
  df_pat$Group4 <- as.character(make_pathology_class(df_pat[[col_path_group]]))
} else {
  df_pat$Group4 <- "G2_Other"
}
df_hc0$Group4 <- "HC"

df_all <- bind_rows(df_hc0, df_pat) %>%
  mutate(
    Group4 = factor(Group4, levels = c("HC","G1_HS","G2_Other")),
    Group2 = factor(ifelse(Group4=="HC","HC","PAT"), levels=c("HC","PAT"))
  )

cat("\n[Sanity check] Group4 counts (including NA):\n")
print(table(df_all$Group4, useNA = "ifany"))

# Covariates MUST exist in BOTH original datasets; avoids "HC side all NA" trap
covars <- covars_candidate[
  covars_candidate %in% colnames(df_pat) &
    covars_candidate %in% colnames(df_hc0)
]
cat("Covariates used (present in BOTH PAT and HC):", paste(covars, collapse=", "), "\n")
if(!all(c(col_age, col_sex, col_etiv) %in% covars)){
  stop("❌ Need Age/sex/eTIV in BOTH df_raw_clean and df_hc_clean. Please check column names.")
}

# ROI columns available in BOTH groups
roi_cols <- intersect(roi_list_user, intersect(colnames(df_pat), colnames(df_hc0)))
if(length(roi_cols) < 10) stop("❌ ROI matched too few. Check colnames consistency.")
cat("ROI used:", length(roi_cols), "\n")

# Save data availability table
roi_avail_tbl <- tibble(
  ROI = roi_cols,
  N_HC  = sapply(roi_cols, function(r) sum(is.finite(df_hc0[[r]]))),
  N_PAT = sapply(roi_cols, function(r) sum(is.finite(df_pat[[r]])))
) %>% arrange(desc(N_PAT+N_HC))
write.csv(roi_avail_tbl, file.path(output_dir, "00_ROI_Data_Availability.csv"),
          row.names=FALSE, fileEncoding="GBK")

# ==========================================================
# 4) Helper functions (FAST + robust)
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

# ==========================================================
# 5) A) Overall PAT vs HC (Group2)
# ==========================================================
cat("\n=== A) Overall PAT vs HC ===\n")

res_A_list <- vector("list", length(roi_cols))

for(i in seq_along(roi_cols)){
  roi <- roi_cols[i]
  need <- unique(c(roi, covars, "Group2"))
  d_tmp <- df_all %>% select(all_of(need)) %>% drop_na()
  if(nrow(d_tmp) < 20) next
  
  out_lm <- lm_contrast_fast(d_tmp, y_col=roi, covars=covars, group_var="Group2", g_level="PAT", ref_level="HC")
  if(is.null(out_lm)) next
  
  y_res <- residualize_vec(d_tmp[[roi]], d_tmp, covars)
  g <- hedges_g(y_res[d_tmp$Group2=="PAT"], y_res[d_tmp$Group2=="HC"])
  
  perm <- freedman_lane_fast(d_tmp, y_col=roi, covars=covars,
                             group_var="Group2", g_level="PAT", ref_level="HC",
                             B=n_perm, seed=perm_seed + i, strata_var=col_sex)
  
  res_A_list[[i]] <- tibble(
    Contrast="PAT_vs_HC",
    ROI=roi,
    N=out_lm$N, N_HC=out_lm$N_ref, N_PAT=out_lm$N_g,
    beta=out_lm$beta_lm, se=out_lm$se_lm, t=out_lm$t_lm, p_lm=out_lm$p_lm,
    p_perm=perm$p_perm,
    Hedges_g_adj=g
  )
  
  if(i %% print_every_roi == 0){
    cat("A progress:", i, "/", length(roi_cols), "ROIs done.\n")
  }
}

res_A <- bind_rows(res_A_list)
if(nrow(res_A)==0) stop("❌ res_A is empty. Check ROI/covariate availability and missingness.")
if(!("p_lm" %in% colnames(res_A))) stop("❌ p_lm missing in res_A (should never happen).")

res_A <- res_A %>%
  mutate(
    q_lm   = p.adjust(p_lm, method="BH"),
    q_perm = p.adjust(p_perm, method="BH"),
    Direction = case_when(beta < 0 ~ "Atrophy(PAT<HC)", beta > 0 ~ "Enlarge(PAT>HC)", TRUE ~ "Zero"),
    Sig_perm_Q05 = ifelse(is.finite(q_perm) & q_perm < 0.05, 1, 0)
  ) %>%
  arrange(q_perm, q_lm, p_perm, p_lm)

write.csv(res_A, file.path(output_dir, "A_overall_PAT_vs_HC_results.csv"),
          row.names=FALSE, fileEncoding="GBK")

# ==========================================================
# 6) B) Subgroups vs HC (Group4)
# ==========================================================
cat("\n=== B) Subgroups vs HC ===\n")

sub_levels <- c("G1_HS","G2_Other")
res_B_all <- list()

for(glev in sub_levels){
  cat("Subgroup:", glev, "\n")
  res_B_list <- vector("list", length(roi_cols))
  
  for(i in seq_along(roi_cols)){
    roi <- roi_cols[i]
    
    need <- unique(c(roi, covars, "Group4"))
    d_tmp <- df_all %>%
      filter(Group4 %in% c("HC", glev)) %>%
      select(all_of(need)) %>% drop_na()
    
    if(nrow(d_tmp) < 20) next
    if(sum(d_tmp$Group4==glev) < min_n_group) next
    if(sum(d_tmp$Group4=="HC") < 5) next
    
    out_lm <- lm_contrast_fast(d_tmp, y_col=roi, covars=covars, group_var="Group4", g_level=glev, ref_level="HC")
    if(is.null(out_lm)) next
    
    y_res <- residualize_vec(d_tmp[[roi]], d_tmp, covars)
    g <- hedges_g(y_res[d_tmp$Group4==glev], y_res[d_tmp$Group4=="HC"])
    
    perm <- freedman_lane_fast(d_tmp, y_col=roi, covars=covars,
                               group_var="Group4", g_level=glev, ref_level="HC",
                               B=n_perm, seed=perm_seed + i + 1000*match(glev, sub_levels),
                               strata_var=col_sex)
    
    res_B_list[[i]] <- tibble(
      Contrast=paste0(glev, "_vs_HC"),
      ROI=roi,
      N=out_lm$N, N_HC=out_lm$N_ref, N_Group=out_lm$N_g,
      beta=out_lm$beta_lm, se=out_lm$se_lm, t=out_lm$t_lm, p_lm=out_lm$p_lm,
      p_perm=perm$p_perm,
      Hedges_g_adj=g
    )
    
    if(i %% print_every_roi == 0){
      cat("B progress:", glev, ":", i, "/", length(roi_cols), "ROIs done.\n")
    }
  }
  
  res_B_g <- bind_rows(res_B_list)
  if(nrow(res_B_g) > 0){
    res_B_g <- res_B_g %>%
      mutate(
        q_lm   = p.adjust(p_lm, method="BH"),
        q_perm = p.adjust(p_perm, method="BH"),
        Direction = case_when(beta < 0 ~ "Atrophy(Group<HC)", beta > 0 ~ "Enlarge(Group>HC)", TRUE ~ "Zero"),
        Sig_perm_Q05 = ifelse(is.finite(q_perm) & q_perm < 0.05, 1, 0)
      ) %>%
      arrange(q_perm, q_lm, p_perm, p_lm)
  }
  res_B_all[[glev]] <- res_B_g
}

res_B <- bind_rows(res_B_all)
write.csv(res_B, file.path(output_dir, "B_subgroups_vs_HC_results.csv"),
          row.names=FALSE, fileEncoding="GBK")

# ==========================================================
# 7) “萎缩在哪里”清单（beta<0）
# ==========================================================
atrophy_A <- res_A %>% filter(beta < 0) %>%
  mutate(Tier = case_when(is.finite(q_perm) & q_perm < 0.05 ~ "Qperm<0.05",
                          is.finite(q_perm) & q_perm < 0.10 ~ "Qperm<0.10",
                          TRUE ~ "ns")) %>%
  arrange(q_perm)

write.csv(atrophy_A, file.path(output_dir, "A_overall_ATROPHY_list.csv"),
          row.names=FALSE, fileEncoding="GBK")

atrophy_B <- res_B %>% filter(beta < 0) %>%
  mutate(Tier = case_when(is.finite(q_perm) & q_perm < 0.05 ~ "Qperm<0.05",
                          is.finite(q_perm) & q_perm < 0.10 ~ "Qperm<0.10",
                          TRUE ~ "ns")) %>%
  arrange(Contrast, q_perm)

write.csv(atrophy_B, file.path(output_dir, "B_subgroups_ATROPHY_list.csv"),
          row.names=FALSE, fileEncoding="GBK")

# ==========================================================
# 8) Visualization (paper-friendly)
# ==========================================================
# --- Fig1: Overall volcano (beta vs -log10(q_perm)) ---
p1 <- res_A %>%
  mutate(logQ = -log10(pmax(q_perm, 1e-300))) %>%
  ggplot(aes(x=beta, y=logQ)) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(aes(color = (beta < 0 & q_perm < 0.05)), alpha=0.85) +
  theme_minimal() +
  labs(
    title="A) Overall: PAT vs HC (covariate-controlled; no Edu)",
    subtitle="x = beta (PAT-HC). beta<0 => PAT smaller than HC (atrophy). y = -log10(q_perm_BH).",
    x="Beta (PAT - HC)", y="-log10(q_perm_BH)"
  )

ggsave(file.path(output_dir, "Fig1_A_overall_volcano_beta_qperm.png"),
       p1, width=8, height=5.5, dpi=300)

# --- Fig2: Overall Top ROIs by q_perm ---
topN <- 30
topA <- res_A %>%
  filter(is.finite(q_perm)) %>%
  arrange(q_perm) %>%
  slice_head(n=topN) %>%
  mutate(ROI = factor(ROI, levels=rev(ROI)))

p2 <- ggplot(topA, aes(x=beta, y=ROI)) +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(size=2) +
  theme_minimal() +
  labs(
    title=paste0("A) Top ", topN, " ROIs by q_perm (PAT vs HC)"),
    subtitle="Negative beta suggests atrophy after covariate control",
    x="Beta (PAT - HC)", y=""
  )

ggsave(file.path(output_dir, "Fig2_A_overall_topROIs_beta.png"),
       p2, width=10, height=max(6, topN*0.28), dpi=300)

# --- Fig3: Subgroup heatmap of EFFECT SIZE (Hedges' g) ---
if(nrow(res_B) > 0){
  top_roi_any <- res_B %>%
    group_by(ROI) %>%
    summarise(min_q = suppressWarnings(min(q_perm, na.rm=TRUE)), .groups="drop") %>%
    arrange(min_q) %>%
    slice_head(n=35) %>%
    pull(ROI)
  
  hm_sub <- res_B %>%
    filter(ROI %in% top_roi_any) %>%
    mutate(ROI = factor(ROI, levels=rev(top_roi_any)))
  
  lim <- suppressWarnings(quantile(abs(hm_sub$Hedges_g_adj), 0.95, na.rm=TRUE))
  if(!is.finite(lim) || lim < 0.5) lim <- 0.5
  
  p3 <- ggplot(hm_sub, aes(x=Contrast, y=ROI, fill=Hedges_g_adj)) +
    geom_tile(color="white") +
    scale_fill_gradient2(
      low="#2166AC", mid="white", high="#B2182B",
      limits=c(-lim, lim),
      oob=scales::squish
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(
      title="B) Subgroups vs HC: Effect-size heatmap (Hedges' g, covariate-adjusted)",
      subtitle=paste0("Fill = Hedges' g on covariate-residualized ROI; color clipped at ±", round(lim,2), " (95th |g|)"),
      x="", y="ROI"
    )
  
  ggsave(file.path(output_dir, "Fig3_B_subgroups_g_heatmap.png"),
         p3, width=6, height=9, dpi=600)
}

cat("\n✅ DONE. Output dir:\n", output_dir, "\n")
cat("A overall results:", nrow(res_A), "ROIs\n")
cat("B subgroup results:", nrow(res_B), "ROI-contrasts\n")
cat("Covariates used:", paste(covars, collapse=", "), "\n")


# ===================================================================================
# Scheme 1: Primary = FDR(q_perm), Secondary = exploratory
# ===================================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

output_dir <- "C:/Users/86150/Documents/HIPP/02011/HC_Atrophy_AB_noEdu/new_2group"

have_resA <- exists("res_A")
have_resB <- exists("res_B")

if(!have_resA){
  fA <- file.path(output_dir, "A_overall_PAT_vs_HC_results.csv")
  if(!file.exists(fA)) stop("❌ Cannot find A_overall_PAT_vs_HC_results.csv in output_dir.")
  res_A <- read.csv(fA, fileEncoding="GBK")
}
if(!have_resB){
  fB <- file.path(output_dir, "B_subgroups_vs_HC_results.csv")
  if(!file.exists(fB)) stop("❌ Cannot find B_subgroups_vs_HC_results.csv in output_dir.")
  res_B <- read.csv(fB, fileEncoding="GBK")
}

alpha_primary <- 0.05
alpha_explore <- 0.05
g_cutoff      <- 0.50
require_atrophy_only <- TRUE

need_cols_A <- c("Contrast","ROI","beta","p_perm","q_perm","Hedges_g_adj")
missA <- setdiff(need_cols_A, colnames(res_A))
if(length(missA)>0) stop("❌ res_A missing columns: ", paste(missA, collapse=", "))

need_cols_B <- c("Contrast","ROI","beta","p_perm","q_perm","Hedges_g_adj")
missB <- setdiff(need_cols_B, colnames(res_B))
if(length(missB)>0) stop("❌ res_B missing columns: ", paste(missB, collapse=", "))

num_cols <- c("beta","p_perm","q_perm","Hedges_g_adj")
for(cc in intersect(num_cols, colnames(res_A))) res_A[[cc]] <- as.numeric(res_A[[cc]])
for(cc in intersect(num_cols, colnames(res_B))) res_B[[cc]] <- as.numeric(res_B[[cc]])

tag_tier <- function(df){
  df %>%
    mutate(
      is_atrophy = is.finite(beta) & beta < 0,
      primary_sig = is.finite(q_perm) & q_perm < alpha_primary,
      explore_sig = is.finite(p_perm) & p_perm < alpha_explore & is.finite(Hedges_g_adj) & abs(Hedges_g_adj) >= g_cutoff,
      Tier = case_when(
        primary_sig ~ "Primary (q_perm<0.05)",
        explore_sig ~ "Exploratory (p_perm<0.05 & |g|>=cutoff)",
        TRUE ~ "NS"
      )
    ) %>%
    { if(require_atrophy_only) filter(., is_atrophy) else . }
}

A_tier <- tag_tier(res_A)
B_tier <- tag_tier(res_B)

A_primary <- A_tier %>% filter(Tier == "Primary (q_perm<0.05)") %>% arrange(q_perm, p_perm)
B_primary <- B_tier %>% filter(Tier == "Primary (q_perm<0.05)") %>% arrange(Contrast, q_perm, p_perm)

A_explore <- A_tier %>% filter(Tier == "Exploratory (p_perm<0.05 & |g|>=cutoff)") %>%
  arrange(p_perm, desc(abs(Hedges_g_adj)))
B_explore <- B_tier %>% filter(Tier == "Exploratory (p_perm<0.05 & |g|>=cutoff)") %>%
  arrange(Contrast, p_perm, desc(abs(Hedges_g_adj)))

A_summary <- A_tier %>%
  mutate(Tier = factor(Tier, levels=c("Primary (q_perm<0.05)","Exploratory (p_perm<0.05 & |g|>=cutoff)","NS"))) %>%
  count(Tier) %>% mutate(Contrast="PAT_vs_HC") %>% select(Contrast, Tier, n)

B_summary <- B_tier %>%
  mutate(Tier = factor(Tier, levels=c("Primary (q_perm<0.05)","Exploratory (p_perm<0.05 & |g|>=cutoff)","NS"))) %>%
  count(Contrast, Tier) %>% arrange(Contrast, Tier)

write.csv(A_primary, file.path(output_dir, "S1_A_primary_qperm_lt_005.csv"), row.names=FALSE, fileEncoding="GBK")
write.csv(A_explore, file.path(output_dir, "S1_A_exploratory_candidates.csv"), row.names=FALSE, fileEncoding="GBK")
write.csv(A_summary, file.path(output_dir, "S1_A_summary_counts.csv"), row.names=FALSE, fileEncoding="GBK")

write.csv(B_primary, file.path(output_dir, "S1_B_primary_qperm_lt_005.csv"), row.names=FALSE, fileEncoding="GBK")
write.csv(B_explore, file.path(output_dir, "S1_B_exploratory_candidates.csv"), row.names=FALSE, fileEncoding="GBK")
write.csv(B_summary, file.path(output_dir, "S1_B_summary_counts.csv"), row.names=FALSE, fileEncoding="GBK")

pA_volcano <- res_A %>%
  mutate(
    beta = as.numeric(beta),
    q_perm = as.numeric(q_perm),
    p_perm = as.numeric(p_perm),
    Hedges_g_adj = as.numeric(Hedges_g_adj),
    is_atrophy = is.finite(beta) & beta < 0,
    primary_sig = is.finite(q_perm) & q_perm < alpha_primary,
    explore_sig = is.finite(p_perm) & p_perm < alpha_explore & is.finite(Hedges_g_adj) & abs(Hedges_g_adj) >= g_cutoff,
    Tier = case_when(
      primary_sig & (!require_atrophy_only | is_atrophy) ~ "Primary (q_perm<0.05)",
      explore_sig & (!require_atrophy_only | is_atrophy) ~ "Exploratory",
      TRUE ~ "NS"
    ),
    y = -log10(pmax(q_perm, 1e-300))
  ) %>%
  ggplot(aes(x=beta, y=y)) +
  geom_hline(yintercept=-log10(alpha_primary), linetype="dashed") +
  geom_vline(xintercept=0, linetype="dashed") +
  geom_point(aes(color=Tier), alpha=0.85) +
  theme_minimal() +
  labs(
    title="Scheme 1 (Overall): Primary=FDR(q_perm), Secondary=Exploratory",
    subtitle=paste0("Exploratory rule: p_perm<", alpha_explore, " & |g|>=", g_cutoff,
                    if(require_atrophy_only) " (atrophy only: beta<0)" else ""),
    x="Beta (PAT - HC)", y="-log10(q_perm_BH)"
  )

ggsave(file.path(output_dir, "S1_FigA_overall_volcano_tiers.png"),
       pA_volcano, width=8, height=5.5, dpi=600)

B_focus <- B_tier %>%
  filter(Tier != "NS") %>%
  mutate(Tier = factor(Tier, levels=c("Primary (q_perm<0.05)","Exploratory (p_perm<0.05 & |g|>=cutoff)")))

if(nrow(B_focus) > 0){
  top_roi <- B_focus %>% group_by(ROI) %>%
    summarise(min_p = min(p_perm, na.rm=TRUE), .groups="drop") %>%
    arrange(min_p) %>% slice_head(n=40) %>% pull(ROI)
  
  pB_heat <- B_focus %>%
    filter(ROI %in% top_roi) %>%
    mutate(ROI = factor(ROI, levels=rev(top_roi))) %>%
    ggplot(aes(x=Contrast, y=ROI, fill=Hedges_g_adj)) +
    geom_tile(color="white") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(
      title="Scheme 1 (Subgroups): Focus ROIs (Primary + Exploratory only)",
      subtitle="Hedges_g_adj<0 = smaller than HC (atrophy) after covariate control",
      x="", y="ROI"
    )
  
  ggsave(file.path(output_dir, "S1_FigB_subgroups_heatmap_focus.png"),
         pB_heat, width=7, height=9, dpi=600)
} else {
  message("⚠️ No subgroup ROI survived Primary or Exploratory criteria; heatmap not produced.")
}

cat("\n=== Scheme 1 outputs saved to ===\n", output_dir, "\n\n")
cat("A primary n =", nrow(A_primary), "\n")
cat("A exploratory n =", nrow(A_explore), "\n")
cat("B primary n =", nrow(B_primary), "\n")
cat("B exploratory n =", nrow(B_explore), "\n")
cat("Done.\n")


# ==========================================================
# FINAL APPENDIX (PASTE AT THE VERY END AND RUN)
# Heatmap for 3 contrasts:
#   (1) All PAT vs HC
#   (2) Left-onset vs HC   [side=1]
#   (3) Right-onset vs HC  [side=0]
# side coding corrected: Left=1, Right=0
# ==========================================================

suppressPackageStartupMessages({ library(tidyverse) })

cat("\n\n================ FINAL APPENDIX: All/Left/Right vs HC + Heatmap ================\n")

req <- c(
  "df_raw_clean","df_hc_clean","roi_cols",
  "col_id","col_age","col_sex","col_etiv","col_side",
  "n_perm","perm_seed","min_n_group","print_every_roi",
  "output_dir",
  "lm_contrast_fast","freedman_lane_fast","residualize_vec","hedges_g"
)
miss <- req[!vapply(req, exists, logical(1), inherits=TRUE)]
if(length(miss) > 0){
  stop("❌ Missing required objects/functions before running appendix:\n- ",
       paste(miss, collapse="\n- "))
}

df_pat <- df_raw_clean
df_hc0 <- df_hc_clean
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if(!(col_id %in% colnames(df_pat))) df_pat[[col_id]] <- paste0("PAT_", seq_len(nrow(df_pat)))
if(!(col_id %in% colnames(df_hc0))) df_hc0[[col_id]] <- paste0("HC_", seq_len(nrow(df_hc0)))

df_pat$Group2 <- "PAT"
df_hc0$Group2 <- "HC"

to_int_safe <- function(x){
  if(is.factor(x)) x <- as.character(x)
  suppressWarnings(as.integer(x))
}

if(!(col_side %in% colnames(df_pat))){
  stop("❌ Patients df has no 'side' column; cannot run laterality module.")
}
if(!(col_side %in% colnames(df_hc0))){
  df_hc0[[col_side]] <- NA
}

df_all_lat <- bind_rows(df_hc0, df_pat) %>%
  mutate(
    Group2 = factor(Group2, levels=c("HC","PAT")),
    side_num = to_int_safe(.data[[col_side]]),
    Laterality = case_when(
      Group2 == "HC" ~ "HC",
      is.finite(side_num) & side_num == 1 ~ "LeftOnset",
      is.finite(side_num) & side_num == 0 ~ "RightOnset",
      TRUE ~ NA_character_
    ),
    Laterality = factor(Laterality, levels=c("HC","LeftOnset","RightOnset"))
  )

cat("\n[Sanity check] Laterality counts (including NA):\n")
print(table(df_all_lat$Laterality, useNA="ifany"))

if(sum(df_all_lat$Laterality=="HC", na.rm=TRUE) < 5){
  stop("❌ Too few HC after merge; check df_hc_clean.")
}
if(sum(df_all_lat$Laterality %in% c("LeftOnset","RightOnset"), na.rm=TRUE) < min_n_group){
  stop("❌ Too few patients with valid side coding (Left=1/Right=0). Check patient 'side' values.")
}

covars_lat <- c(col_etiv, col_age, col_sex)
covars_lat <- covars_lat[covars_lat %in% colnames(df_all_lat)]
if(!all(c(col_etiv, col_age, col_sex) %in% covars_lat)){
  stop("❌ Need eTIV/Age/sex present for appendix models.")
}
cat("Covariates used (All/Left/Right models):", paste(covars_lat, collapse=", "), "\n")

run_one_contrast <- function(d, roi, group_var, g_level, ref_level, covars, B, seed, strata_var){
  need <- unique(c(roi, covars, group_var))
  d_tmp <- d %>% select(all_of(need)) %>% drop_na()
  if(nrow(d_tmp) < 20) return(NULL)
  if(sum(d_tmp[[group_var]]==g_level) < min_n_group) return(NULL)
  if(sum(d_tmp[[group_var]]==ref_level) < 5) return(NULL)
  
  out_lm <- lm_contrast_fast(d_tmp, y_col=roi, covars=covars, group_var=group_var, g_level=g_level, ref_level=ref_level)
  if(is.null(out_lm)) return(NULL)
  
  y_res <- residualize_vec(d_tmp[[roi]], d_tmp, covars)
  g <- hedges_g(y_res[d_tmp[[group_var]]==g_level], y_res[d_tmp[[group_var]]==ref_level])
  
  perm <- freedman_lane_fast(d_tmp, y_col=roi, covars=covars,
                             group_var=group_var, g_level=g_level, ref_level=ref_level,
                             B=B, seed=seed, strata_var=strata_var)
  
  tibble(
    ROI=roi,
    N=out_lm$N,
    N_ref=out_lm$N_ref,
    N_g=out_lm$N_g,
    beta=out_lm$beta_lm,
    se=out_lm$se_lm,
    t=out_lm$t_lm,
    p_lm=out_lm$p_lm,
    p_perm=perm$p_perm,
    Hedges_g_adj=g
  )
}

cat("\n=== A0) All PAT vs HC (covars: eTIV+Age+sex) ===\n")

res_all_list <- vector("list", length(roi_cols))
for(i in seq_along(roi_cols)){
  roi <- roi_cols[i]
  d <- df_all_lat %>% mutate(Group2 = droplevels(Group2))
  
  out <- run_one_contrast(
    d=d, roi=roi, group_var="Group2", g_level="PAT", ref_level="HC",
    covars=covars_lat, B=n_perm, seed=perm_seed + 10000 + i, strata_var=col_sex
  )
  if(is.null(out)) next
  
  res_all_list[[i]] <- out %>% mutate(Contrast="AllPAT_vs_HC")
  
  if(i %% print_every_roi == 0) cat("All progress:", i, "/", length(roi_cols), "ROIs done.\n")
}
res_all <- bind_rows(res_all_list) %>%
  mutate(
    q_lm   = p.adjust(p_lm, method="BH"),
    q_perm = p.adjust(p_perm, method="BH"),
    Direction = case_when(beta < 0 ~ "Smaller(PAT<HC)", beta > 0 ~ "Larger(PAT>HC)", TRUE ~ "Zero"),
    Sig_perm_Q05 = ifelse(is.finite(q_perm) & q_perm < 0.05, 1, 0),
    N_HC = N_ref, N_Group = N_g
  ) %>%
  select(Contrast, ROI, N, N_HC, N_Group, beta, se, t, p_lm, q_lm, p_perm, q_perm, Hedges_g_adj, Direction, Sig_perm_Q05) %>%
  arrange(q_perm, q_lm)

write.csv(res_all, file.path(output_dir, "A0_allPAT_vs_HC_noSide_results.csv"),
          row.names=FALSE, fileEncoding="GBK")

cat("A0 rows:", nrow(res_all), "\n")

cat("\n=== C) Laterality: Left vs HC; Right vs HC (side Left=1 Right=0) ===\n")
lat_levels <- c("LeftOnset","RightOnset")
res_lat_all <- list()

for(glev in lat_levels){
  cat("Laterality subgroup:", glev, "\n")
  res_one <- vector("list", length(roi_cols))
  
  for(i in seq_along(roi_cols)){
    roi <- roi_cols[i]
    d <- df_all_lat %>% filter(Laterality %in% c("HC", glev)) %>% droplevels()
    
    out <- run_one_contrast(
      d=d, roi=roi, group_var="Laterality", g_level=glev, ref_level="HC",
      covars=covars_lat, B=n_perm, seed=perm_seed + 20000 + i + 5000*match(glev, lat_levels),
      strata_var=col_sex
    )
    if(is.null(out)) next
    
    res_one[[i]] <- out %>% mutate(Contrast=paste0(glev, "_vs_HC"))
    
    if(i %% print_every_roi == 0) cat("Laterality progress:", glev, ":", i, "/", length(roi_cols), "ROIs done.\n")
  }
  
  res_g <- bind_rows(res_one)
  if(nrow(res_g) > 0){
    res_g <- res_g %>%
      mutate(
        q_lm   = p.adjust(p_lm, method="BH"),
        q_perm = p.adjust(p_perm, method="BH"),
        Direction = case_when(beta < 0 ~ "Smaller(Group<HC)", beta > 0 ~ "Larger(Group>HC)", TRUE ~ "Zero"),
        Sig_perm_Q05 = ifelse(is.finite(q_perm) & q_perm < 0.05, 1, 0),
        N_HC = N_ref, N_Group = N_g
      ) %>%
      select(Contrast, ROI, N, N_HC, N_Group, beta, se, t, p_lm, q_lm, p_perm, q_perm, Hedges_g_adj, Direction, Sig_perm_Q05) %>%
      arrange(q_perm, q_lm)
  }
  res_lat_all[[glev]] <- res_g
}

res_lat <- bind_rows(res_lat_all)

write.csv(res_lat, file.path(output_dir, "C_laterality_LeftRight_vs_HC_results.csv"),
          row.names=FALSE, fileEncoding="GBK")

cat("C rows:", nrow(res_lat), "\n")

cat("\n=== Heatmap: All / Left / Right vs HC (Hedges' g) ===\n")

res_heat <- bind_rows(
  res_all %>% select(Contrast, ROI, Hedges_g_adj, q_perm, p_perm, beta),
  res_lat %>% select(Contrast, ROI, Hedges_g_adj, q_perm, p_perm, beta)
) %>%
  mutate(
    Contrast = factor(Contrast, levels=c("AllPAT_vs_HC","LeftOnset_vs_HC","RightOnset_vs_HC"))
  )

topN <- 45
top_roi <- res_heat %>%
  group_by(ROI) %>%
  summarise(min_q = suppressWarnings(min(q_perm, na.rm=TRUE)), .groups="drop") %>%
  arrange(min_q) %>%
  slice_head(n=topN) %>%
  pull(ROI)

hm <- res_heat %>%
  filter(ROI %in% top_roi) %>%
  mutate(ROI = factor(ROI, levels=rev(top_roi)))

lim <- suppressWarnings(quantile(abs(hm$Hedges_g_adj), 0.95, na.rm=TRUE))
if(!is.finite(lim) || lim < 0.5) lim <- 0.5

p_hm_g <- ggplot(hm, aes(x=Contrast, y=ROI, fill=Hedges_g_adj)) +
  geom_tile(color="white") +
  scale_fill_gradient2(
    low="#2166AC", mid="white", high="#B2182B",
    limits=c(-lim, lim),
    oob=scales::squish,
    na.value="grey90"
  ) +
  theme_minimal() +
  theme(  axis.text.x = element_text(size=14, face="bold", angle=25, hjust=1),
          axis.text.y = element_text(size=14, face="plain"),
          strip.text.y = element_text(size=12, face="bold"),
          plot.title = element_text(size=14, face="bold"),
          plot.subtitle = element_text(size=11, face="plain")) +
  labs(
    title="Effect-size heatmap: All PAT / Left-onset / Right-onset vs HC",
    subtitle=paste0("Hedges' g on covariate-residualized ROI (eTIV+Age+sex). Color clipped at ±", round(lim,2),
                    ". Negative = smaller than HC."),
    x="", y="ROI"
  )

ggsave(file.path(output_dir, "FigC_heatmap_All_Left_Right_vs_HC_HedgesG.png"),
       p_hm_g, width=9.5, height=max(8, topN*0.23), dpi=600)

hm2 <- hm %>%
  mutate(logQ = -log10(pmax(q_perm, 1e-300)))

p_hm_q <- ggplot(hm2, aes(x=Contrast, y=ROI, fill=logQ)) +
  geom_tile(color="white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=35, hjust=1)) +
  labs(
    title="Significance heatmap: -log10(q_perm_BH)",
    subtitle="Higher = more significant (BH-FDR on permutation p).",
    x="", y="ROI"
  )

ggsave(file.path(output_dir, "FigD_heatmap_All_Left_Right_vs_HC_logQperm.png"),
       p_hm_q, width=9.5, height=max(8, topN*0.23), dpi=300)

cat("\n✅ DONE. Outputs saved in:\n", output_dir, "\n")
cat(" - A0_allPAT_vs_HC_noSide_results.csv\n")
cat(" - C_laterality_LeftRight_vs_HC_results.csv\n")
cat(" - FigC_heatmap_All_Left_Right_vs_HC_HedgesG.png\n")
cat(" - FigD_heatmap_All_Left_Right_vs_HC_logQperm.png\n")


# ==========================================================
# APPENDIX: Family-wise FDR (Laterality 3 contrasts)
# Families: Cortex / HippocampusSubfields / ThalamusNuclei / SubcorticalAseg / GlobalOther
# Primary significance: q_perm_family (BH within family & contrast)
# ==========================================================

suppressPackageStartupMessages({ library(tidyverse) })

cat("\n\n================ APPENDIX: Family-wise FDR + Heatmap (3 contrasts) ================\n")

if(!exists("output_dir")) stop("❌ output_dir not found.")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if(!exists("res_all")){
  f_all <- file.path(output_dir, "A0_allPAT_vs_HC_noSide_results.csv")
  if(!file.exists(f_all)) stop("❌ Cannot find A0_allPAT_vs_HC_noSide_results.csv in output_dir.")
  res_all <- read.csv(f_all, fileEncoding="GBK")
}
if(!exists("res_lat")){
  f_lat <- file.path(output_dir, "C_laterality_LeftRight_vs_HC_results.csv")
  if(!file.exists(f_lat)) stop("❌ Cannot find C_laterality_LeftRight_vs_HC_results.csv in output_dir.")
  res_lat <- read.csv(f_lat, fileEncoding="GBK")
}

need_cols <- c("Contrast","ROI","beta","p_perm","q_perm","Hedges_g_adj")
miss_all <- setdiff(need_cols, colnames(res_all))
miss_lat <- setdiff(need_cols, colnames(res_lat))
if(length(miss_all)>0) stop("❌ res_all missing columns: ", paste(miss_all, collapse=", "))
if(length(miss_lat)>0) stop("❌ res_lat missing columns: ", paste(miss_lat, collapse=", "))

num_cols <- c("beta","p_perm","q_perm","Hedges_g_adj")
for(cc in intersect(num_cols, colnames(res_all))) res_all[[cc]] <- as.numeric(res_all[[cc]])
for(cc in intersect(num_cols, colnames(res_lat))) res_lat[[cc]] <- as.numeric(res_lat[[cc]])

res_3 <- bind_rows(
  res_all %>% select(Contrast, ROI, beta, p_perm, q_perm, Hedges_g_adj),
  res_lat %>% select(Contrast, ROI, beta, p_perm, q_perm, Hedges_g_adj)
) %>%
  mutate(
    Contrast = as.character(Contrast),
    Contrast = case_when(
      Contrast %in% c("AllPAT_vs_HC","PAT_vs_HC","PAT_vs_HC_noSide","AllPATvsHC") ~ "AllPAT_vs_HC",
      grepl("^LeftOnset", Contrast) ~ "LeftOnset_vs_HC",
      grepl("^RightOnset", Contrast) ~ "RightOnset_vs_HC",
      TRUE ~ Contrast
    ),
    Contrast = factor(Contrast, levels=c("AllPAT_vs_HC","LeftOnset_vs_HC","RightOnset_vs_HC"))
  )

cat("[Sanity check] contrasts present:\n")
print(table(res_3$Contrast, useNA="ifany"))

# --- UPDATED family definition: add SubcorticalAseg ---
roi_family <- function(roi){
  r <- as.character(roi)
  
  if(grepl("^Left-Subcort-", r) || grepl("^Right-Subcort-", r)) return("SubcorticalAseg")
  
  if(r %in% c("lh_BrainSegVolNotVent","rh_BrainSegVolNotVent","BrainSegVolNotVent")) return("GlobalOther")
  if(grepl("^lh_lh_", r) || grepl("^rh_rh_", r)) return("Cortex")
  if(grepl("^Left_", r) || grepl("^Right_", r)) return("ThalamusNuclei")
  if(grepl("^lh_", r) || grepl("^rh_", r)) return("HippocampusSubfields")
  
  return("GlobalOther")
}

res_3 <- res_3 %>%
  mutate(
    Family = vapply(ROI, roi_family, character(1)),
    Family = factor(Family, levels=c("Cortex","HippocampusSubfields","ThalamusNuclei","SubcorticalAseg","GlobalOther"))
  )

cat("[Sanity check] ROI family counts:\n")
print(table(res_3$Family, useNA="ifany"))

res_3 <- res_3 %>%
  group_by(Contrast, Family) %>%
  mutate(q_perm_family = p.adjust(p_perm, method="BH")) %>%
  ungroup() %>%
  mutate(
    star = case_when(
      is.finite(q_perm_family) & q_perm_family < 0.001 ~ "***",
      is.finite(q_perm_family) & q_perm_family < 0.01  ~ "**",
      is.finite(q_perm_family) & q_perm_family < 0.05  ~ "*",
      is.finite(q_perm_family) & q_perm_family < 0.10  ~ "·",
      TRUE ~ ""
    )
  )

write.csv(res_3, file.path(output_dir, "E_familyFDR_All_Left_Right_vs_HC.csv"),
          row.names=FALSE, fileEncoding="GBK")
cat("Saved: E_familyFDR_All_Left_Right_vs_HC.csv\n")

for(ff in levels(res_3$Family)){
  tmp <- res_3 %>% filter(Family==ff) %>% arrange(Contrast, q_perm_family, p_perm)
  write.csv(tmp, file.path(output_dir, paste0("E_familyFDR_", ff, ".csv")),
            row.names=FALSE, fileEncoding="GBK")
}

alpha_show <- 0.10
topN <- 70

# ✅ BUGFIX: res3 -> res_3
top_roi <- res_3 %>%
  group_by(ROI, Family) %>%
  summarise(min_qfam = suppressWarnings(min(q_perm_family, na.rm=TRUE)), .groups="drop") %>%
  filter(is.finite(min_qfam) & min_qfam < alpha_show) %>%
  arrange(min_qfam) %>%
  slice_head(n=topN) %>%
  pull(ROI)

hm <- res_3 %>%
  filter(ROI %in% top_roi) %>%
  mutate(ROI = factor(ROI, levels=rev(top_roi)))

lim <- suppressWarnings(quantile(abs(hm$Hedges_g_adj), 0.95, na.rm=TRUE))
if(!is.finite(lim) || lim < 0.5) lim <- 0.5

p_hm <- ggplot(hm, aes(x=Contrast, y=ROI, fill=Hedges_g_adj)) +
  geom_tile(color="white") +
  geom_text(aes(label=star), size=5) +
  facet_grid(Family ~ ., scales="free_y", space="free_y") +
  scale_fill_gradient2(
    low="#2166AC", mid="white", high="#B2182B",
    limits=c(-lim, lim),
    oob=scales::squish,
    na.value="grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size=12, face="bold", angle=25, hjust=1),
    axis.text.y = element_text(size=12, face="plain"),
    strip.text.y = element_text(size=12, face="bold"),
    plot.title = element_text(size=14, face="bold"),
    plot.subtitle = element_text(size=11, face="plain")
  ) +
  labs(
    title="3-contrast heatmap (Family-wise FDR)",
    subtitle=paste0("Fill=Hedges' g (covariate-adjusted). Stars by q_perm_family (BH within Family & Contrast). ",
                    "Color clipped at ±", round(lim,2),
                    ".  *<0.05, **<0.01, ***<0.001, ·<0.10"),
    x="", y="ROI"
  )

ggsave(file.path(output_dir, "FigE_heatmap_familyFDR_HedgesG_star.png"),
       p_hm, width=8, height=max(9, topN*0.23), dpi=600)

cat("Saved: FigE_heatmap_familyFDR_HedgesG_star.png\n")

cat("\n✅ APPENDIX DONE. Output dir:\n", output_dir, "\n")
cat(" - E_familyFDR_All_Left_Right_vs_HC.csv\n")
cat(" - E_familyFDR_Cortex.csv / E_familyFDR_HippocampusSubfields.csv / E_familyFDR_ThalamusNuclei.csv / E_familyFDR_SubcorticalAseg.csv / E_familyFDR_GlobalOther.csv\n")
cat(" - FigE_heatmap_familyFDR_HedgesG_star.png\n")


# ==========================================================
# APPENDIX: Pathology 3-contrast heatmap
#   (1) AllPAT_vs_HC
#   (2) HS_vs_HC
#   (3) Other_vs_HC
# Family-wise BH-FDR (within each Family & Contrast) on p_perm
# Families include SubcorticalAseg
# ==========================================================

suppressPackageStartupMessages({ library(tidyverse) })

cat("\n\n================ APPENDIX: Pathology 3-contrast Family-FDR + Heatmap ================\n")

if(!exists("output_dir")) stop("❌ output_dir not found.")
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if(exists("res_all")){
  df_all_con <- res_all
  df_all_con$Contrast <- "AllPAT_vs_HC"
} else if(exists("res_A")){
  df_all_con <- res_A
  df_all_con$Contrast <- "AllPAT_vs_HC"
} else {
  f_new <- file.path(output_dir, "A0_allPAT_vs_HC_noSide_results.csv")
  f_old <- file.path(output_dir, "A_overall_PAT_vs_HC_results.csv")
  if(file.exists(f_new)){
    df_all_con <- read.csv(f_new, fileEncoding="GBK")
    df_all_con$Contrast <- "AllPAT_vs_HC"
  } else if(file.exists(f_old)){
    df_all_con <- read.csv(f_old, fileEncoding="GBK")
    df_all_con$Contrast <- "AllPAT_vs_HC"
  } else {
    stop("❌ Cannot find overall contrast results. Need res_all/res_A or CSV in output_dir.")
  }
}

if(!exists("res_B")){
  fB <- file.path(output_dir, "B_subgroups_vs_HC_results.csv")
  if(!file.exists(fB)) stop("❌ Cannot find B_subgroups_vs_HC_results.csv in output_dir, and res_B not in workspace.")
  res_B <- read.csv(fB, fileEncoding="GBK")
}

need_cols <- c("ROI","p_perm","Hedges_g_adj")
for(nm in need_cols){
  if(!(nm %in% colnames(df_all_con))) stop("❌ overall results missing column: ", nm)
  if(!(nm %in% colnames(res_B)))      stop("❌ res_B missing column: ", nm)
}

df_all_con$p_perm <- as.numeric(df_all_con$p_perm)
df_all_con$Hedges_g_adj <- as.numeric(df_all_con$Hedges_g_adj)

res_B$p_perm <- as.numeric(res_B$p_perm)
res_B$Hedges_g_adj <- as.numeric(res_B$Hedges_g_adj)

df_path <- res_B %>%
  mutate(
    Contrast = as.character(Contrast),
    Contrast = case_when(
      grepl("^G1_HS", Contrast) ~ "HS_vs_HC",
      grepl("^G2_Other", Contrast) ~ "Other_vs_HC",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Contrast)) %>%
  mutate(Contrast = factor(Contrast, levels=c("HS_vs_HC","Other_vs_HC")))

if(nrow(df_path)==0) stop("❌ No HS/Other contrasts found in res_B. Check Contrast naming in res_B.")

res3 <- bind_rows(
  df_all_con %>%
    transmute(Contrast = factor("AllPAT_vs_HC", levels=c("AllPAT_vs_HC","HS_vs_HC","Other_vs_HC")),
              ROI, p_perm, Hedges_g_adj),
  df_path %>%
    transmute(Contrast = factor(as.character(Contrast), levels=c("AllPAT_vs_HC","HS_vs_HC","Other_vs_HC")),
              ROI, p_perm, Hedges_g_adj)
)

cat("[Sanity check] 3 contrasts present:\n")
print(table(res3$Contrast, useNA="ifany"))

# --- UPDATED family definition: add SubcorticalAseg ---
roi_family <- function(roi){
  r <- as.character(roi)
  
  if(grepl("^Left-Subcort-", r) || grepl("^Right-Subcort-", r)) return("SubcorticalAseg")
  if(r %in% c("lh_BrainSegVolNotVent","rh_BrainSegVolNotVent","BrainSegVolNotVent")) return("GlobalOther")
  if(grepl("^lh_lh_", r) || grepl("^rh_rh_", r)) return("Cortex")
  if(grepl("^Left_", r) || grepl("^Right_", r)) return("ThalamusNuclei")
  if(grepl("^lh_", r) || grepl("^rh_", r)) return("HippocampusSubfields")
  
  return("GlobalOther")
}

res3 <- res3 %>%
  mutate(
    Family = vapply(ROI, roi_family, character(1)),
    Family = factor(Family, levels=c("Cortex","HippocampusSubfields","ThalamusNuclei","SubcorticalAseg","GlobalOther"))
  )

cat("[Sanity check] family counts:\n")
print(table(res3$Family, useNA="ifany"))

res3 <- res3 %>%
  group_by(Contrast, Family) %>%
  mutate(q_perm_family = p.adjust(p_perm, method="BH")) %>%
  ungroup() %>%
  mutate(
    star = case_when(
      is.finite(q_perm_family) & q_perm_family < 0.001 ~ "***",
      is.finite(q_perm_family) & q_perm_family < 0.01  ~ "**",
      is.finite(q_perm_family) & q_perm_family < 0.05  ~ "*",
      is.finite(q_perm_family) & q_perm_family < 0.10  ~ "·",
      TRUE ~ ""
    )
  )

write.csv(res3, file.path(output_dir, "P3_familyFDR_AllPAT_HS_Other_vs_HC.csv"),
          row.names=FALSE, fileEncoding="GBK")
cat("Saved: P3_familyFDR_AllPAT_HS_Other_vs_HC.csv\n")

alpha_show <- 0.10
topN <- 70

top_roi <- res3 %>%
  group_by(ROI, Family) %>%
  summarise(min_qfam = suppressWarnings(min(q_perm_family, na.rm=TRUE)), .groups="drop") %>%
  filter(is.finite(min_qfam) & min_qfam < alpha_show) %>%
  arrange(min_qfam) %>%
  slice_head(n=topN) %>%
  pull(ROI)

hm <- res3 %>%
  filter(ROI %in% top_roi) %>%
  mutate(ROI = factor(ROI, levels=rev(top_roi)))

lim <- suppressWarnings(quantile(abs(hm$Hedges_g_adj), 0.95, na.rm=TRUE))
if(!is.finite(lim) || lim < 0.5) lim <- 0.5

p_hm <- ggplot(hm, aes(x=Contrast, y=ROI, fill=Hedges_g_adj)) +
  geom_tile(color="white") +
  geom_text(aes(label=star), size=5) +
  facet_grid(Family ~ ., scales="free_y", space="free_y") +
  scale_fill_gradient2(
    low="#2166AC", mid="white", high="#B2182B",
    limits=c(-lim, lim),
    oob=scales::squish,
    na.value="grey90"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size=12, face="bold", angle=25, hjust=1),
    axis.text.y = element_text(size=12, face="plain"),
    strip.text.y = element_text(size=12, face="bold"),
    plot.title = element_text(size=14, face="bold"),
    plot.subtitle = element_text(size=11, face="plain")
  ) +
  labs(
    title="All PAT / HS / Other vs HC（Family-wise FDR）",
    subtitle=paste0("填充=Hedges' g（协变量校正后效应量）. 星号=q_perm_family（每个 Family×Contrast 内 BH-FDR 于 p_perm）. ",
                    "颜色剪裁 ±", round(lim,2),
                    "；*<0.05，**<0.01，***<0.001，·<0.10"),
    x="", y="ROI"
  )

ggsave(file.path(output_dir, "FigP3_heatmap_AllPAT_HS_Other_familyFDR_HedgesG_star.png"),
       p_hm, width=8, height=max(9, topN*0.22), dpi=600)
cat("Saved: FigP3_heatmap_AllPAT_HS_Other_familyFDR_HedgesG_star.png\n")

sum_tbl <- res3 %>%
  mutate(sig05 = is.finite(q_perm_family) & q_perm_family < 0.05) %>%
  group_by(Contrast, Family) %>%
  summarise(
    n_test = sum(is.finite(p_perm)),
    n_sig05 = sum(sig05, na.rm=TRUE),
    .groups="drop"
  ) %>%
  arrange(Contrast, Family)

write.csv(sum_tbl, file.path(output_dir, "P3_familyFDR_summary_counts.csv"),
          row.names=FALSE, fileEncoding="GBK")
cat("Saved: P3_familyFDR_summary_counts.csv\n")

cat("\n✅ DONE. Output dir:\n", output_dir, "\n")
cat(" - P3_familyFDR_AllPAT_HS_Other_vs_HC.csv\n")
cat(" - FigP3_heatmap_AllPAT_HS_Other_familyFDR_HedgesG_star.png\n")
cat(" - P3_familyFDR_summary_counts.csv\n")

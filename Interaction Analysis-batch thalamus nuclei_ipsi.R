# ==============================================================================
# FINAL SCRIPT (PAT-only residualized atrophy; NO PCA; batch IPSI/CONTRA thalamus nuclei × IPSI/CONTRA hip subfields)
# IPSI/CONTRA recode FIRST (RAW), then residualize in PAT (eTIV + Age) -> atrophy = -z(resid)
# Freedman–Lane permutation for interaction; small-scope FDR within (Outcome × Thal) over Hippo set size (6 or 16)
# Only IPSI tracks; LR track removed
# Plot ALL significant results (perm-FDR < 0.05): simple slopes + Cook's + LOO
# Excel-friendly outputs (UTF-8 BOM) via write_excel_csv()
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(readr)
  library(fs)
})

# =========================
# 0) SETTINGS
# =========================
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found (patients).")
df_pat0 <- df_raw_clean

# 建议用短路径，避免“链结无效”
OUTPUT_DIR <- "C:/Users/86150/Documents/ROI_Cognition_Interaction/batch_signal_IPThal/5000"
dir_create(OUTPUT_DIR)

SEED   <- 2026
N_PERM <- 5000
MIN_N  <- 20
set.seed(SEED)

TARGET_OUTCOMES <- c(
  "工作记忆指数WMI","总记忆商","听觉记忆指数","视觉记忆指数","即刻记忆指数","延迟记忆指数",
  "VA-immediate","VA-recall","VTrials I-V Total","CAVLT-slope"
  #"FA-immediate","FA-recall","FTrials I-V Total","AFLT-slope"
)

# Main covariates (Age NOT included; already adjusted in atrophy construction)
COV_IP_MAIN <- c("Edu_year", "Duration_Dis")

# Optional sensitivity (+Age in regression model)做敏感性分析的开关
DO_AGE_SENS <- TRUE
COV_IP_AGESENS <- c("Age", "Edu_year", "Duration_Dis")

ALPHA_FDR <- 0.05

# =========================
# 1) OUTPUT helper (Excel-friendly UTF-8 BOM)
# =========================
save_xcsv <- function(df, path){
  dir_create(path_dir(path))
  readr::write_excel_csv(df, path)
  invisible(path)
}

# =========================
# 2) HELPERS
# =========================
to_numeric_safely <- function(x) {
  if (is.numeric(x)) return(x)
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    x2 <- gsub("[,，]", "", x)
    x2 <- trimws(x2)
    suppressWarnings(as.numeric(x2))
  } else suppressWarnings(as.numeric(x))
}

# 你已明确：side 0/1, 0 = R, 1 = L
normalize_side_lr <- function(x){
  if(is.factor(x)) x <- as.character(x)
  
  if(is.character(x)){
    xl <- tolower(trimws(x))
    out <- rep(NA_character_, length(xl))
    out[xl %in% c("l","left","lt","lh","lhs","左","左侧")] <- "L"
    out[xl %in% c("r","right","rt","rh","rhs","右","右侧")] <- "R"
    out[xl %in% c("0")] <- "R"
    out[xl %in% c("1")] <- "L"
    return(out)
  }
  
  if(is.numeric(x) || is.integer(x)){
    ux <- sort(unique(x[is.finite(x)]))
    out <- rep(NA_character_, length(x))
    if(all(ux %in% c(0,1))){
      out[x == 0] <- "R"
      out[x == 1] <- "L"
      return(out)
    }
  }
  rep(NA_character_, length(x))
}

# PAT-only residualized atrophy with joint adjustment for eTIV + Age
# atrophy = - z(residual)
pat_resid_atrophy <- function(y_pat, etiv_pat, age_pat){
  d <- data.frame(y=y_pat, eTIV=etiv_pat, Age=age_pat)
  ok <- complete.cases(d)
  out <- rep(NA_real_, length(y_pat))
  if(sum(ok) < 10) return(out)
  
  fit <- lm(y ~ eTIV + Age, data=d[ok,])
  r_pat <- rep(NA_real_, length(y_pat))
  r_pat[ok] <- d$y[ok] - predict(fit, newdata=d[ok,])
  
  sdv <- sd(r_pat[ok], na.rm=TRUE)
  if(!is.finite(sdv) || sdv == 0) return(out)
  
  z_pat <- rep(NA_real_, length(y_pat))
  z_pat[ok] <- as.numeric(scale(r_pat[ok]))
  -z_pat
}

# 通用 ipsi/contra builder：要求 left_vars 和 right_vars 是成对的、前缀分别可指定
make_ipsi_contra <- function(df, side_col, left_vars, right_vars,
                             left_prefix, right_prefix, out_prefix){
  left_suf  <- sub(paste0("^", left_prefix),  "", left_vars)
  right_suf <- sub(paste0("^", right_prefix), "", right_vars)
  if(!setequal(left_suf, right_suf)) stop("❌ Left/Right suffix mismatch for ipsi/contra recoding.")
  
  suf_all <- sort(left_suf)
  for(suf in suf_all){
    L <- paste0(left_prefix, suf)
    R <- paste0(right_prefix, suf)
    ipsi   <- paste0("ipsi_", out_prefix, "_", suf)
    contra <- paste0("contra_", out_prefix, "_", suf)
    df[[ipsi]]   <- ifelse(df[[side_col]]=="L", df[[L]], df[[R]])
    df[[contra]] <- ifelse(df[[side_col]]=="L", df[[R]], df[[L]])
  }
  df
}

# Freedman–Lane permutation test for interaction term
freedman_lane_interaction <- function(dat, y, x, m, covars, n_perm=5000, seed=2026){
  set.seed(seed)
  
  dat <- dat %>%
    mutate(
      y_z = as.numeric(scale(.data[[y]])),
      x_z = as.numeric(scale(.data[[x]])),
      m_z = as.numeric(scale(.data[[m]]))
    )
  
  cov_str <- if(length(covars)>0) paste(covars, collapse=" + ") else "1"
  f_red  <- as.formula(paste0("y_z ~ x_z + m_z + ", cov_str))
  f_full <- as.formula(paste0("y_z ~ x_z + m_z + x_z:m_z + ", cov_str))
  
  fit_red  <- lm(f_red, data = dat)
  fit_full <- lm(f_full, data = dat)
  
  coefs_full <- summary(fit_full)$coefficients
  if(!("x_z:m_z" %in% rownames(coefs_full))) stop("❌ interaction term not found in full model.")
  
  t_obs   <- unname(coefs_full["x_z:m_z", "t value"])
  p_param <- unname(coefs_full["x_z:m_z", "Pr(>|t|)"])
  beta    <- unname(coefs_full["x_z:m_z", "Estimate"])
  se      <- unname(coefs_full["x_z:m_z", "Std. Error"])
  
  y_hat <- fitted(fit_red)
  e     <- resid(fit_red)
  
  t_perm <- numeric(n_perm)
  for(b in seq_len(n_perm)){
    y_perm <- y_hat + sample(e, replace = FALSE)
    dat$y_perm <- y_perm
    fit_p <- lm(update(f_full, y_perm ~ .), data = dat)
    t_perm[b] <- unname(summary(fit_p)$coefficients["x_z:m_z", "t value"])
  }
  
  p_perm <- (1 + sum(abs(t_perm) >= abs(t_obs), na.rm=TRUE)) / (n_perm + 1)
  list(beta_int=beta, se_int=se, t_obs=t_obs, p_param=p_param, p_perm=p_perm, n=nrow(dat))
}

# =========================
# 3) CORE: run one block (Outcome × Thal × Hippo), small-scope FDR within Outcome × Thal
# =========================
run_block <- function(df, outcomes, thal_vars, hip_vars, covars, block_name){
  out_list <- list()
  k <- 1L
  
  for(out in outcomes){
    if(!(out %in% names(df))) stop("❌ Missing outcome in df: ", out)
    
    for(t in thal_vars){
      for(h in hip_vars){
        needed <- c(out, t, h, covars)
        miss <- setdiff(needed, names(df))
        if(length(miss)>0) stop("❌ Missing columns in analysis: ", paste(miss, collapse=", "))
        
        dat <- df %>% select(all_of(needed)) %>% filter(complete.cases(.))
        
        if(nrow(dat) < MIN_N){
          out_list[[k]] <- tibble(
            Block=block_name, Outcome=out, Thal=t, Hippo=h, n=nrow(dat),
            beta_int=NA_real_, se_int=NA_real_, t_obs=NA_real_,
            p_param=NA_real_, p_perm=NA_real_,
            Covars=paste(covars, collapse=" + ")
          )
          k <- k + 1L
          next
        }
        
        ans <- freedman_lane_interaction(dat, y=out, x=t, m=h, covars=covars,
                                         n_perm=N_PERM, seed=SEED)
        
        out_list[[k]] <- tibble(
          Block=block_name, Outcome=out, Thal=t, Hippo=h, n=ans$n,
          beta_int=ans$beta_int, se_int=ans$se_int, t_obs=ans$t_obs,
          p_param=ans$p_param, p_perm=ans$p_perm,
          Covars=paste(covars, collapse=" + ")
        )
        
        cat("✅", block_name, "|", out, "|", t, "×", h,
            "| N=", ans$n, "| p_perm=", sprintf("%.4f", ans$p_perm), "\n")
        
        k <- k + 1L
      }
    }
  }
  
  res <- bind_rows(out_list) %>%
    group_by(Block, Outcome, Thal) %>%   # ✅关键：固定 Outcome + 固定 Thal，只跨 Hippo 做FDR
    mutate(
      p_perm_fdr  = p.adjust(p_perm,  method="fdr"),
      p_param_fdr = p.adjust(p_param, method="fdr"),
      n_in_fdr_set = sum(is.finite(p_perm))
    ) %>% ungroup()
  
  res
}

# =========================
# 4) Diagnostics + Plots
# =========================
influence_and_loo <- function(fit_full, term="x:h"){
  cooks <- cooks.distance(fit_full)
  hat   <- hatvalues(fit_full)
  
  mf <- model.frame(fit_full)
  n <- nrow(mf)
  beta_loo <- rep(NA_real_, n)
  for(i in seq_len(n)){
    fit_i <- lm(formula(fit_full), data = mf[-i, , drop=FALSE])
    beta_loo[i] <- coef(fit_i)[term]
  }
  list(cooks=cooks, hat=hat, beta_loo=beta_loo, beta_full=coef(fit_full)[term])
}

plot_simple_slopes <- function(df, outcome, thal_var, hip_var, covars, out_png){
  dat <- df %>%
    select(all_of(c(outcome, thal_var, hip_var, covars))) %>%
    filter(complete.cases(.)) %>%
    mutate(
      y = as.numeric(scale(.data[[outcome]])),
      x = as.numeric(scale(.data[[thal_var]])),
      h = as.numeric(scale(.data[[hip_var]]))
    )
  
  q <- quantile(dat$h, probs=c(1/3, 2/3), na.rm=TRUE)
  dat <- dat %>% mutate(h_grp = case_when(
    h <= q[1] ~ "Low hippo atrophy",
    h <= q[2] ~ "Mid hippo atrophy",
    TRUE      ~ "High hippo atrophy"
  ))
  
  make_typical <- function(v){
    if(is.numeric(v)) mean(v, na.rm=TRUE) else names(sort(table(v), decreasing=TRUE))[1]
  }
  
  cov_str <- if(length(covars)>0) paste(covars, collapse=" + ") else "1"
  f <- as.formula(paste0("y ~ x + h + x:h + ", cov_str))
  fit <- lm(f, data=dat)
  
  xseq <- seq(min(dat$x), max(dat$x), length.out=80)
  h_levels <- c("Low hippo atrophy","Mid hippo atrophy","High hippo atrophy")
  h_val <- tapply(dat$h, dat$h_grp, median, na.rm=TRUE)[h_levels]
  
  grid <- expand.grid(x=xseq, h_grp=h_levels) %>% mutate(h = as.numeric(h_val[h_grp]))
  for(cv in covars) grid[[cv]] <- make_typical(dat[[cv]])
  
  pred <- predict(fit, newdata=grid, se.fit=TRUE)
  grid$yhat <- pred$fit
  grid$se   <- pred$se.fit
  grid$lo   <- grid$yhat - 1.96*grid$se
  grid$hi   <- grid$yhat + 1.96*grid$se
  
  p <- ggplot(grid, aes(x=x, y=yhat, color=h_grp, fill=h_grp)) +
    geom_line(linewidth=1.1) +
    geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.15, color=NA) +
    theme_minimal(base_size = 12) +
    labs(
      x=paste0(thal_var, " (z; higher=worse atrophy)"),
      y=paste0(outcome, " (z)"),
      color="Hippo atrophy terciles",
      fill="Hippo atrophy terciles",
      title=paste0("Interaction: ", thal_var, " × ", hip_var),
      subtitle="Plot uses parametric prediction; inference uses Freedman–Lane permutation"
    )
  
  ggsave(out_png, p, width=10, height=7, dpi=400, bg="white")
  p
}

# =========================
# 5) COERCE CORE VARIABLES
# =========================
need_global <- c("side","Age","Edu_year","Duration_Dis","eTIV")
miss_pat <- setdiff(need_global, names(df_pat0))
if(length(miss_pat)>0) stop("❌ Missing required columns in PAT: ", paste(miss_pat, collapse=", "))

df_pat <- df_pat0 %>%
  mutate(
    Age = to_numeric_safely(Age),
    Edu_year = to_numeric_safely(Edu_year),
    Duration_Dis = to_numeric_safely(Duration_Dis),
    eTIV = to_numeric_safely(eTIV),
    side_lr = normalize_side_lr(.data[["side"]])
  )

if(any(is.na(df_pat$side_lr))){
  print(table(df_pat$side, useNA="ifany"))
  stop("❌ side cannot be normalized to L/R. Your side coding should be 0=R,1=L.")
}

for(y in TARGET_OUTCOMES){
  if(!(y %in% names(df_pat))) stop("❌ Missing outcome in PAT: ", y)
  df_pat[[y]] <- to_numeric_safely(df_pat[[y]])
}

# =========================
# 6) ROI LISTS
# =========================
thal_rois_lr <- c(
  "Left_AV","Left_VA","Left_VAmc","Left_VLa","Left_VLp","Left_VPL","Left_VM","Left_CL","Left_CeM","Left_CM","Left_Pf",
  "Left_MDm","Left_MDl","Left_LD","Left_LP","Left_Pu_Total","Left_MGN","Left_LGN","Left_L_Sg",
  "Left_MV_Re","Left_Pc","Left_Pt",
  "Right_AV","Right_VA","Right_VAmc","Right_VLa","Right_VLp","Right_VPL","Right_VM","Right_CL","Right_CeM","Right_CM","Right_Pf",
  "Right_MDm","Right_MDl","Right_LD","Right_LP","Right_Pu_Total","Right_MGN","Right_LGN","Right_L_Sg",
  "Right_MV_Re","Right_Pc","Right_Pt"
)

hip_primary_lr <- c(
  "lh_Whole_hippocampal_head","lh_Whole_hippocampal_body","lh_Hippocampal_tail",
  "rh_Whole_hippocampal_head","rh_Whole_hippocampal_body","rh_Hippocampal_tail"
)

hip_explor_lr <- c(
  "lh_subiculum_comb","lh_CA1_comb","lh_CA3_comb","lh_CA4_comb","lh_GC_ML_DG_comb",
  "lh_molecular_layer_HP_comb","lh_presubiculum_comb","lh_parasubiculum",
  "rh_subiculum_comb","rh_CA1_comb","rh_CA3_comb","rh_CA4_comb","rh_GC_ML_DG_comb",
  "rh_molecular_layer_HP_comb","rh_presubiculum_comb","rh_parasubiculum"
)

# Checks
miss_thal <- setdiff(thal_rois_lr, names(df_pat))
miss_hp1  <- setdiff(hip_primary_lr, names(df_pat))
miss_hp2  <- setdiff(hip_explor_lr, names(df_pat))
if(length(miss_thal)>0) stop("❌ Missing thalamus ROIs in PAT: ", paste(miss_thal, collapse=", "))
if(length(miss_hp1)>0)  stop("❌ Missing PRIMARY hip ROIs in PAT: ", paste(miss_hp1, collapse=", "))
if(length(miss_hp2)>0)  stop("❌ Missing EXPLOR hip ROIs in PAT: ", paste(miss_hp2, collapse=", "))

# =========================
# 7) IPSI/CONTRA RECODE FIRST (RAW volumes): HIPPO + THALAMUS
# =========================
df_pat_a <- df_pat

# ---- HIPPO ipsi/contra ----
build_ipsi_contra_raw_for_list <- function(df, lr_raw, out_prefix, left_prefix, right_prefix){
  left_vars  <- lr_raw[grepl(paste0("^", left_prefix), lr_raw)]
  right_vars <- lr_raw[grepl(paste0("^", right_prefix), lr_raw)]
  df <- make_ipsi_contra(df, "side_lr", left_vars, right_vars,
                         left_prefix=left_prefix, right_prefix=right_prefix, out_prefix=out_prefix)
  suf <- sub(paste0("^", left_prefix), "", left_vars)
  ipsi_contra <- c(paste0("ipsi_", out_prefix, "_", suf),
                   paste0("contra_", out_prefix, "_", suf))
  list(df=df, ipsi_contra=ipsi_contra)
}

tmpH1 <- build_ipsi_contra_raw_for_list(df_pat_a, hip_primary_lr, out_prefix="hipA", left_prefix="lh_", right_prefix="rh_")
df_pat_a <- tmpH1$df
HIP_IP_PRIMARY_RAW <- tmpH1$ipsi_contra

tmpH2 <- build_ipsi_contra_raw_for_list(df_pat_a, hip_explor_lr, out_prefix="hipB", left_prefix="lh_", right_prefix="rh_")
df_pat_a <- tmpH2$df
HIP_IP_EXPLOR_RAW <- tmpH2$ipsi_contra

# ---- THAL ipsi/contra ----
# 将 Left_* 与 Right_* 配对为 ipsi_thal_XXX / contra_thal_XXX
left_thal  <- thal_rois_lr[grepl("^Left_", thal_rois_lr)]
right_thal <- thal_rois_lr[grepl("^Right_", thal_rois_lr)]
# 这里用前缀 Left_/Right_ 来做 suffix 匹配
df_pat_a <- make_ipsi_contra(df_pat_a, "side_lr", left_thal, right_thal,
                             left_prefix="Left_", right_prefix="Right_", out_prefix="thal")
thal_suf <- sub("^Left_", "", left_thal)
THAL_IP_RAW <- c(paste0("ipsi_thal_", thal_suf), paste0("contra_thal_", thal_suf))

# =========================
# 8) BUILD PAT-ONLY residualized atrophy variables (eTIV + Age)
# =========================
# ensure numeric for raw LR ROIs used in recode
for(r in c(thal_rois_lr, hip_primary_lr, hip_explor_lr, HIP_IP_PRIMARY_RAW, HIP_IP_EXPLOR_RAW, THAL_IP_RAW)){
  df_pat_a[[r]] <- to_numeric_safely(df_pat_a[[r]])
}

# (A) THAL ipsi/contra atrophy
for(r in THAL_IP_RAW){
  df_pat_a[[paste0(r, "_atrophy_eTIVAge_z")]] <- pat_resid_atrophy(
    y_pat    = df_pat_a[[r]],
    etiv_pat = df_pat_a$eTIV,
    age_pat  = df_pat_a$Age
  )
}
THAL_IP_ATROPHY <- paste0(THAL_IP_RAW, "_atrophy_eTIVAge_z")

# (B) HIP ipsi/contra atrophy (primary + explor)
for(r in c(HIP_IP_PRIMARY_RAW, HIP_IP_EXPLOR_RAW)){
  df_pat_a[[paste0(r, "_atrophy_eTIVAge_z")]] <- pat_resid_atrophy(
    y_pat    = df_pat_a[[r]],
    etiv_pat = df_pat_a$eTIV,
    age_pat  = df_pat_a$Age
  )
}
HIP_IP_PRIMARY <- paste0(HIP_IP_PRIMARY_RAW, "_atrophy_eTIVAge_z") # 6
HIP_IP_EXPLOR  <- paste0(HIP_IP_EXPLOR_RAW,  "_atrophy_eTIVAge_z") # 16

# =========================
# 9) RUN MAIN (IPSI only; NO LR)
# =========================
 res_ip_primary_main <- run_block(df_pat_a, TARGET_OUTCOMES, THAL_IP_ATROPHY, HIP_IP_PRIMARY, COV_IP_MAIN, "IPSI_PRIMARY_MAIN")
save_xcsv(res_ip_primary_main, file.path(OUTPUT_DIR, "Interaction_IPSI_PRIMARY_MAIN.csv"))

res_ip_explor_main  <- run_block(df_pat_a, TARGET_OUTCOMES, THAL_IP_ATROPHY, HIP_IP_EXPLOR,  COV_IP_MAIN, "IPSI_EXPLOR_MAIN")
save_xcsv(res_ip_explor_main,  file.path(OUTPUT_DIR, "Interaction_IPSI_EXPLOR_MAIN.csv"))


res_ip_primary_main <-Interaction_IPSI_PRIMARY_MAIN
sig_ip_primary_main <- res_ip_primary_main %>%
  filter(is.finite(p_perm_fdr) & p_perm_fdr < ALPHA_FDR) %>% arrange(p_perm_fdr, p_perm)
sig_ip_explor_main  <- res_ip_explor_main  %>%
  filter(is.finite(p_perm_fdr) & p_perm_fdr < ALPHA_FDR) %>% arrange(p_perm_fdr, p_perm)

save_xcsv(sig_ip_primary_main, file.path(OUTPUT_DIR, paste0("SIG_IPSI_PRIMARY_MAIN_permFDRlt", ALPHA_FDR, ".csv")))
save_xcsv(sig_ip_explor_main,  file.path(OUTPUT_DIR, paste0("SIG_IPSI_EXPLOR_MAIN_permFDRlt",  ALPHA_FDR, ".csv")))

cat("\n==================== SUMMARY (MAIN) ====================\n")
cat("IPSI_PRIMARY_MAIN sig:", nrow(sig_ip_primary_main), "\n")
cat("IPSI_EXPLOR_MAIN  sig:", nrow(sig_ip_explor_main),  "\n")
cat("Output dir:", normalizePath(OUTPUT_DIR), "\n")

# =========================
# 10) OPTIONAL AGE-IN-MODEL SENSITIVITY (IPSI only)
# =========================
if(DO_AGE_SENS){
  res_ip_primary_age <- run_block(df_pat_a, TARGET_OUTCOMES, THAL_IP_ATROPHY, HIP_IP_PRIMARY, COV_IP_AGESENS, "IPSI_PRIMARY_AGESENS")
  res_ip_explor_age  <- run_block(df_pat_a, TARGET_OUTCOMES, THAL_IP_ATROPHY, HIP_IP_EXPLOR,  COV_IP_AGESENS, "IPSI_EXPLOR_AGESENS")
  
  save_xcsv(res_ip_primary_age, file.path(OUTPUT_DIR, "Interaction_IPSI_PRIMARY_AGESENS.csv"))
  save_xcsv(res_ip_explor_age,  file.path(OUTPUT_DIR, "Interaction_IPSI_EXPLOR_AGESENS.csv"))
}

# =========================
# 11) TABLE: All interactions (MAIN)
# =========================
tab_main <- bind_rows(
  res_ip_primary_main %>% mutate(Tier="Primary"),
  res_ip_explor_main  %>% mutate(Tier="Exploratory")
) %>%
  select(Tier, Block, Outcome, Thal, Hippo, n, n_in_fdr_set,
         beta_int, se_int, t_obs, p_perm, p_perm_fdr, p_param, p_param_fdr, Covars) %>%
  arrange(Tier, Outcome, Thal, p_perm_fdr)

save_xcsv(tab_main, file.path(OUTPUT_DIR, "TABLE_AllInteractions_MAIN.csv"))

# =========================
# 12) PLOTS for ALL significant results (perm-FDR < 0.05)
#     For EACH significant pair: simple slopes + Cook's + LOO
# =========================
plot_one_sig <- function(df_all, out, thal, hip, covars, tag_prefix){
  # 1) simple slopes
  fn1 <- file.path(OUTPUT_DIR, paste0("FIG_SimpleSlopes_", tag_prefix, "_", out, "_", thal, "_X_", hip, ".png"))
  p1 <- plot_simple_slopes(df_all, out, thal, hip, covars, fn1)
  print(p1)
  
  # 2) diagnostics (parametric)
  dat <- df_all %>%
    select(all_of(c(out, thal, hip, covars))) %>%
    filter(complete.cases(.)) %>%
    mutate(
      y = as.numeric(scale(.data[[out]])),
      x = as.numeric(scale(.data[[thal]])),
      h = as.numeric(scale(.data[[hip]]))
    )
  
  f_full <- as.formula(paste0("y ~ x + h + x:h + ", paste(covars, collapse=" + ")))
  fit_full <- lm(f_full, data=dat)
  
  inf <- influence_and_loo(fit_full, term="x:h")
  diag_df <- tibble(
    idx = seq_along(inf$cooks),
    cooks = inf$cooks,
    hat = inf$hat,
    beta_loo = inf$beta_loo,
    beta_full = inf$beta_full
  )
  
  fn_diag <- file.path(OUTPUT_DIR, paste0("DIAG_", tag_prefix, "_", out, "_", thal, "_X_", hip, ".csv"))
  save_xcsv(diag_df, fn_diag)
  
  p_cooks <- ggplot(diag_df, aes(x=idx, y=cooks)) +
    geom_point() + theme_minimal(base_size=12) +
    labs(title=paste0("Cook's distance: ", out, " | ", thal, " × ", hip),
         x="Observation index", y="Cook's distance")
  ggsave(file.path(OUTPUT_DIR, paste0("FIG_Cooks_", tag_prefix, ".png")),
         p_cooks, width=6, height=5, dpi=600, bg="white")
  
  p_loo <- ggplot(diag_df, aes(x=idx, y=beta_loo)) +
    geom_line() + geom_hline(yintercept=inf$beta_full, linetype=2) +
    theme_minimal(base_size=12) +
    labs(title=paste0("LOO beta (interaction): ", out, " | ", thal, " × ", hip),
         subtitle="Dashed line is full-sample beta",
         x="Left-out observation index", y="Interaction beta (LOO)")
  ggsave(file.path(OUTPUT_DIR, paste0("FIG_LOO_", tag_prefix, ".png")),
         p_loo, width=6, height=5, dpi=600, bg="white")
}

# Primary significant
if(nrow(sig_ip_primary_main) > 0){
  for(i in seq_len(nrow(sig_ip_primary_main))){
    out  <- sig_ip_primary_main$Outcome[i]
    thal <- sig_ip_primary_main$Thal[i]
    hip  <- sig_ip_primary_main$Hippo[i]
    tag  <- paste0("PRIMARY_sig", i)
    plot_one_sig(df_pat_a, out, thal, hip, COV_IP_MAIN, tag)
  }
}

# Exploratory significant
if(nrow(sig_ip_explor_main) > 0){
  for(i in seq_len(nrow(sig_ip_explor_main))){
    out  <- sig_ip_explor_main$Outcome[i]
    thal <- sig_ip_explor_main$Thal[i]
    hip  <- sig_ip_explor_main$Hippo[i]
    tag  <- paste0("EXPLOR_sig", i)
    plot_one_sig(df_pat_a, out, thal, hip, COV_IP_MAIN, tag)
  }
}

cat("\n✅ ALL DONE. IPSI/CONTRA thalamus+hippocampus, PAT-only residualized atrophy (eTIV+Age),\n")
cat("   Freedman–Lane permutation + small-scope FDR within (Outcome×Thal), plots for ALL significant results.\n")
cat("📁 Outputs saved to: ", normalizePath(OUTPUT_DIR), "\n", sep="")






















# ==============================================================================
# ADD-ON ONLY (append to the END of your script)另一种交互图
# Purpose:
#   1) Create a NEW subfolder under OUTPUT_DIR
#   2) Draw "SWAP simple slopes" for ALL significant results (perm-FDR < 0.05)
#   3) NO re-drawing Cook's / LOO (same model; not needed)
# Run:
#   - Assumes these already exist in your environment from the main script:
#       df_pat_a, OUTPUT_DIR, COV_IP_MAIN, sig_ip_primary_main, sig_ip_explor_main
#   - Just run this block.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# -------------------------
# 0) New subfolder
# -------------------------
PLOT_DIR_SWAP <- file.path(OUTPUT_DIR, "PLOTS_SwapSimpleSlopes")
dir.create(PLOT_DIR_SWAP, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 1) Safe filename helper (avoid invalid path / Chinese issues)
# -------------------------
safe_filename <- function(x){
  x <- as.character(x)
  x <- gsub("[/\\\\:*?\"<>|]", "_", x)   # Windows-illegal
  x <- gsub("\\s+", "_", x)
  x <- gsub("__+", "_", x)
  x <- substr(x, 1, 180)
  x
}

# -------------------------
# 2) Name cleaner for titles (shorten; remove suffix/prefix noise)
# -------------------------
clean_roi_name <- function(x){
  x <- as.character(x)
  x <- gsub("_atrophy_eTIVAge_z$", "", x)
  x <- gsub("^ipsi_thal_", "ipsi_", x)
  x <- gsub("^contra_thal_", "contra_", x)
  x <- gsub("^ipsi_hipA_", "ipsi_", x)
  x <- gsub("^contra_hipA_", "contra_", x)
  x <- gsub("^ipsi_hipB_", "ipsi_", x)
  x <- gsub("^contra_hipB_", "contra_", x)
  x <- gsub("\\s+", " ", x)
  x
}

# -------------------------
# 3) Swap simple slopes function
#    x-axis = hip atrophy (h)
#    color/group = thal atrophy terciles (x_grp)
#    NOTE: NO subtitle; title = "Interaction: <Thal> × <Hip>"
# -------------------------
plot_simple_slopes_swap <- function(df, outcome, thal_var, hip_var, covars, out_png,
                                    base_size = 14,
                                    axis_title_size = 18,
                                    axis_text_size = 14,
                                    title_size = 18,
                                    legend_inside = TRUE,
                                    legend_pos = c(0.03, 0.03)){
  
  dat <- df %>%
    dplyr::select(dplyr::all_of(c(outcome, thal_var, hip_var, covars))) %>%
    dplyr::filter(stats::complete.cases(.)) %>%
    dplyr::mutate(
      y = as.numeric(scale(.data[[outcome]])),
      x = as.numeric(scale(.data[[thal_var]])), # thal
      h = as.numeric(scale(.data[[hip_var]]))   # hip
    )
  
  # terciles of thal atrophy (for color groups)
  qx <- stats::quantile(dat$x, probs = c(1/3, 2/3), na.rm = TRUE)
  dat <- dat %>% dplyr::mutate(x_grp = dplyr::case_when(
    x <= qx[1] ~ "Low thal atrophy",
    x <= qx[2] ~ "Mid thal atrophy",
    TRUE       ~ "High thal atrophy"
  ))
  
  # fit same interaction model (parametric for visualization only)
  make_typical <- function(v){
    if(is.numeric(v)) mean(v, na.rm = TRUE) else names(sort(table(v), decreasing = TRUE))[1]
  }
  cov_str <- if(length(covars) > 0) paste(covars, collapse = " + ") else "1"
  f <- stats::as.formula(paste0("y ~ x + h + x:h + ", cov_str))
  fit <- stats::lm(f, data = dat)
  
  # prediction grid: vary hip (h) on x-axis; fix thal (x) at group medians
  hseq <- seq(min(dat$h), max(dat$h), length.out = 80)
  x_levels <- c("Low thal atrophy","Mid thal atrophy","High thal atrophy")
  x_val <- tapply(dat$x, dat$x_grp, stats::median, na.rm = TRUE)[x_levels]
  
  grid <- expand.grid(h = hseq, x_grp = x_levels) %>%
    dplyr::mutate(x = as.numeric(x_val[x_grp]))
  
  for(cv in covars) grid[[cv]] <- make_typical(dat[[cv]])
  
  pred <- stats::predict(fit, newdata = grid, se.fit = TRUE)
  grid$yhat <- pred$fit
  grid$se   <- pred$se.fit
  grid$lo   <- grid$yhat - 1.96 * grid$se
  grid$hi   <- grid$yhat + 1.96 * grid$se
  
  # --- Title: show who × who (clean names) ---
  thal_title <- clean_roi_name(thal_var)
  hip_title  <- clean_roi_name(hip_var)
  title_text <- paste0("Interaction: ", thal_title, " × ", hip_title)
  
  p <- ggplot2::ggplot(grid, ggplot2::aes(x = h, y = yhat, color = x_grp, fill = x_grp)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::labs(
      x = "Hippocampal atrophy (z; higher=worse)",
      y = paste0(outcome, " (z)"),
      color = "Thal terciles",
      fill  = "Thal terciles",
      title = title_text
      # ✅ no subtitle
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold"),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text  = ggplot2::element_text(size = axis_text_size),
      legend.title = ggplot2::element_text(size = axis_text_size + 1),
      legend.text  = ggplot2::element_text(size = axis_text_size)
    )
  
  if(isTRUE(legend_inside)){
    p <- p + ggplot2::theme(
      legend.position = legend_pos,
      legend.justification = c(0, 0),
      legend.background = ggplot2::element_rect(fill = "white", color = "grey80")
    )
  }
  
  ggplot2::ggsave(out_png, p, width = 8, height = 7.5, dpi = 600, bg = "white")
  invisible(p)
}

# -------------------------
# 4) Batch: draw swap plots for all significant results
# -------------------------
plot_swap_for_sig_table <- function(sig_tbl, tag_prefix){
  if(is.null(sig_tbl) || nrow(sig_tbl) == 0) return(invisible(NULL))
  
  for(i in seq_len(nrow(sig_tbl))){
    out  <- sig_tbl$Outcome[i]
    thal <- sig_tbl$Thal[i]
    hip  <- sig_tbl$Hippo[i]
    
    fn <- file.path(
      PLOT_DIR_SWAP,
      paste0(
        "FIG_SwapSlopes_", tag_prefix, i, "_",
        safe_filename(out), "_",
        safe_filename(thal), "_X_",
        safe_filename(hip), ".png"
      )
    )
    
    plot_simple_slopes_swap(
      df = df_pat_a,
      outcome = out,
      thal_var = thal,
      hip_var  = hip,
      covars   = COV_IP_MAIN,
      out_png  = fn,
      # 你想调就改这几行
      base_size = 14,
      axis_title_size = 14,
      axis_text_size = 18,
      title_size = 18,
      legend_inside = TRUE,
      legend_pos = c(0.03, 0.03)
    )
    
    cat("✅ SWAP plot saved:", fn, "\n")
  }
}

# -------------------------
# 5) Run (Primary + Exploratory significant tables)
# -------------------------
plot_swap_for_sig_table(sig_ip_primary_main, "PRIMARY_sig")
plot_swap_for_sig_table(sig_ip_explor_main,  "EXPLOR_sig")

cat("\n✅ DONE. Swap simple-slopes saved to:\n", normalizePath(PLOT_DIR_SWAP), "\n", sep="")

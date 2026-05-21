############################################################
## Additive models of hippocampal and thalamic volume measures

############################################################


## =========================
## 0. 基础检查与参数
## =========================

if (!exists("df_raw_clean")) stop("环境中未找到 df_raw_clean")

OUT_DIR <- "C:/Users/86150/Documents/HIPP_Thal/additive_model"
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

N_PERM <- 5000
SEED   <- 2026
set.seed(SEED)

LEFT_CODES  <- c("1", "L", "l", "Left", "left", "左")
RIGHT_CODES <- c("0", "2", "R", "r", "Right", "right", "右")

## 主 outcome
OUTCOMES_PRIMARY <- c(
  # "WAIS-总智商",
  # "言语理解指数VCI",
  # "知觉推理指数PRI",
  # "加工速度指数PSI",
  "工作记忆指数WMI",
  "总记忆商",
  "听觉记忆指数",
  "视觉记忆指数",
  "即刻记忆指数",
  "延迟记忆指数"
)

## 主分析协变量
## 与 HC-referenced burden 版本不同，这里把 eTIV 直接放进认知模型中校正，注意，这里不再进行年龄的校正，因为韦氏量表里面已经进行了
MAIN_COVARS <- c("Edu_year", "Duration_Dis_num", "sex_f", "eTIV")

df_pat <- df_raw_clean


## =========================
## 1. 工具函数
## =========================

bt <- function(x) paste0("`", x, "`")

make_rhs <- function(vars) {
  if (length(vars) == 0) return("1")
  paste(bt(vars), collapse = " + ")
}

calc_r2_from_fitted <- function(y, fitted_vals) {
  ok <- is.finite(y) & is.finite(fitted_vals)
  y <- y[ok]
  fitted_vals <- fitted_vals[ok]
  
  if (length(y) < 3) return(NA_real_)
  
  tss <- sum((y - mean(y))^2)
  rss <- sum((y - fitted_vals)^2)
  
  if (tss == 0) return(NA_real_)
  
  1 - rss / tss
}

row_sum_min_n <- function(df, vars, min_prop = 1.0) {
  vars_exist <- intersect(vars, names(df))
  
  if (length(vars_exist) == 0) {
    return(rep(NA_real_, nrow(df)))
  }
  
  x <- as.matrix(df[, vars_exist, drop = FALSE])
  x <- apply(x, 2, function(col) suppressWarnings(as.numeric(col)))
  
  if (is.null(dim(x))) {
    x <- matrix(x, ncol = 1)
  }
  
  nonmiss_n <- rowSums(!is.na(x))
  need_n <- ceiling(length(vars_exist) * min_prop)
  
  out <- rowSums(x, na.rm = TRUE)
  out[nonmiss_n < need_n] <- NA_real_
  
  out
}

recode_side_lr <- function(x, left_codes = LEFT_CODES, right_codes = RIGHT_CODES) {
  x_chr <- trimws(as.character(x))
  
  out <- ifelse(
    x_chr %in% as.character(left_codes), "L",
    ifelse(x_chr %in% as.character(right_codes), "R", NA)
  )
  
  out
}

drop_constant_vars <- function(df, vars) {
  vars <- intersect(vars, names(df))
  
  keep <- vars[sapply(vars, function(v) {
    x <- df[[v]]
    x_nonmiss <- x[!is.na(x)]
    
    if (length(x_nonmiss) == 0) return(FALSE)
    
    if (is.factor(x) || is.character(x)) {
      length(unique(x_nonmiss)) > 1
    } else {
      xx <- suppressWarnings(as.numeric(x_nonmiss))
      if (all(is.na(xx))) return(FALSE)
      sd(xx, na.rm = TRUE) > 0
    }
  })]
  
  keep
}

freedman_lane_delta_r2 <- function(df, y_var, reduced_vars, full_vars,
                                   n_perm = 5000, seed = 1234, min_n = 15) {
  
  need_vars <- unique(c(y_var, reduced_vars, full_vars))
  need_vars <- intersect(need_vars, names(df))
  d <- df[, need_vars, drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  
  if (nrow(d) < min_n) {
    return(list(
      n = nrow(d),
      obs_delta_r2 = NA_real_,
      p_perm = NA_real_
    ))
  }
  
  reduced_vars <- drop_constant_vars(d, reduced_vars)
  full_vars    <- drop_constant_vars(d, full_vars)
  full_vars    <- unique(c(reduced_vars, setdiff(full_vars, reduced_vars)))
  
  y <- d[[y_var]]
  
  Xr <- model.matrix(as.formula(paste0("~ ", make_rhs(reduced_vars))), data = d)
  Xf <- model.matrix(as.formula(paste0("~ ", make_rhs(full_vars))), data = d)
  
  fit_r <- lm.fit(x = Xr, y = y)
  fit_f <- lm.fit(x = Xf, y = y)
  
  r2_r <- calc_r2_from_fitted(y, fit_r$fitted.values)
  r2_f <- calc_r2_from_fitted(y, fit_f$fitted.values)
  
  obs_delta <- r2_f - r2_r
  
  if (!is.finite(obs_delta)) {
    return(list(
      n = nrow(d),
      obs_delta_r2 = NA_real_,
      p_perm = NA_real_
    ))
  }
  
  set.seed(seed)
  perm_stats <- numeric(n_perm)
  
  for (b in seq_len(n_perm)) {
    y_perm <- fit_r$fitted.values + sample(fit_r$residuals, replace = FALSE)
    
    fit_r_b <- lm.fit(x = Xr, y = y_perm)
    fit_f_b <- lm.fit(x = Xf, y = y_perm)
    
    perm_stats[b] <- calc_r2_from_fitted(y_perm, fit_f_b$fitted.values) -
      calc_r2_from_fitted(y_perm, fit_r_b$fitted.values)
  }
  
  p_perm <- (1 + sum(perm_stats >= obs_delta - 1e-12, na.rm = TRUE)) / (n_perm + 1)
  
  list(
    n = nrow(d),
    obs_delta_r2 = obs_delta,
    p_perm = p_perm
  )
}

get_coef_safe <- function(fit, varname) {
  cf <- coef(fit)
  
  if (varname %in% names(cf)) {
    return(unname(cf[varname]))
  }
  
  return(NA_real_)
}

write_csv_cn <- function(df, filename, out_dir = OUT_DIR) {
  write.csv(
    df,
    file = file.path(out_dir, filename),
    row.names = FALSE,
    fileEncoding = "GBK"
  )
}


## =========================
## 2. 基础预处理
## =========================

df_pat$sex_f <- factor(as.character(df_pat$sex))

df_pat$Duration_Dis_num <- suppressWarnings(as.numeric(df_pat$Duration_Dis))
df_pat$Age <- suppressWarnings(as.numeric(df_pat$Age))
df_pat$Edu_year <- suppressWarnings(as.numeric(df_pat$Edu_year))
df_pat$eTIV <- suppressWarnings(as.numeric(df_pat$eTIV))

df_pat$side_lr <- recode_side_lr(df_pat$side)

if (any(is.na(df_pat$side_lr))) {
  warning("部分患者 side 无法识别为 L/R，这些样本在 ipsi/contra 映射中会产生 NA。")
}

OUTCOMES_PRIMARY <- intersect(OUTCOMES_PRIMARY, names(df_pat))

if (length(OUTCOMES_PRIMARY) == 0) {
  stop("OUTCOMES_PRIMARY 中没有变量存在于 df_pat 中，请检查列名。")
}


## =========================
## 3. 构建 raw volume / summed volume indices
## =========================

## ---- 海马 raw whole hippocampus volume ----
hippo_left_raw_var  <- "Left-Subcort-Hippocampus"
hippo_right_raw_var <- "Right-Subcort-Hippocampus"

if (!(hippo_left_raw_var %in% names(df_pat))) {
  stop("df_pat 缺少变量：", hippo_left_raw_var)
}
if (!(hippo_right_raw_var %in% names(df_pat))) {
  stop("df_pat 缺少变量：", hippo_right_raw_var)
}

df_pat$Hippo_Left_Volume  <- as.numeric(df_pat[[hippo_left_raw_var]])
df_pat$Hippo_Right_Volume <- as.numeric(df_pat[[hippo_right_raw_var]])


## ---- 丘脑核团 raw summed volume ----
thal_left_vars <- c(
  "Left_AV", "Left_VA", "Left_VAmc", "Left_VLa", "Left_VLp", "Left_VPL", "Left_VM",
  "Left_CL", "Left_CeM", "Left_CM", "Left_Pf",
  "Left_MDm", "Left_MDl", "Left_LD", "Left_LP",
  "Left_MGN", "Left_LGN", "Left_L_Sg",
  "Left_MV_Re", "Left_Pc", "Left_Pt", "Left_Pu_Total"
)

thal_right_vars <- c(
  "Right_AV", "Right_VA", "Right_VAmc", "Right_VLa", "Right_VLp", "Right_VPL", "Right_VM",
  "Right_CL", "Right_CeM", "Right_CM", "Right_Pf",
  "Right_MDm", "Right_MDl", "Right_LD", "Right_LP",
  "Right_MGN", "Right_LGN", "Right_L_Sg",
  "Right_MV_Re", "Right_Pc", "Right_Pt", "Right_Pu_Total"
)

miss_thal_pat <- setdiff(c(thal_left_vars, thal_right_vars), names(df_pat))

if (length(miss_thal_pat) > 0) {
  stop("df_pat 缺少以下丘脑核团列：\n", paste(miss_thal_pat, collapse = "\n"))
}

thal_left_exist  <- intersect(thal_left_vars, names(df_pat))
thal_right_exist <- intersect(thal_right_vars, names(df_pat))

## min_prop = 1.0 表示要求所有丘脑核团成员齐全，否则该侧 summed volume 设为 NA
df_pat$Thal_Left_SummedVolume  <- row_sum_min_n(df_pat, thal_left_exist,  min_prop = 1.0)
df_pat$Thal_Right_SummedVolume <- row_sum_min_n(df_pat, thal_right_exist, min_prop = 1.0)


## =========================
## 4. Lesion-aligned volume mapping
## =========================

df_pat$Hippo_Ipsi_Volume <- ifelse(
  df_pat$side_lr == "L",
  df_pat$Hippo_Left_Volume,
  ifelse(df_pat$side_lr == "R", df_pat$Hippo_Right_Volume, NA_real_)
)

df_pat$Hippo_Contra_Volume <- ifelse(
  df_pat$side_lr == "L",
  df_pat$Hippo_Right_Volume,
  ifelse(df_pat$side_lr == "R", df_pat$Hippo_Left_Volume, NA_real_)
)

df_pat$Thal_Ipsi_SummedVolume <- ifelse(
  df_pat$side_lr == "L",
  df_pat$Thal_Left_SummedVolume,
  ifelse(df_pat$side_lr == "R", df_pat$Thal_Right_SummedVolume, NA_real_)
)

df_pat$Thal_Contra_SummedVolume <- ifelse(
  df_pat$side_lr == "L",
  df_pat$Thal_Right_SummedVolume,
  ifelse(df_pat$side_lr == "R", df_pat$Thal_Left_SummedVolume, NA_real_)
)

df_pat$Thal_Global_SummedVolume <- rowMeans(
  cbind(df_pat$Thal_Left_SummedVolume, df_pat$Thal_Right_SummedVolume),
  na.rm = TRUE
)

df_pat$Thal_Global_SummedVolume[!is.finite(df_pat$Thal_Global_SummedVolume)] <- NA_real_

df_pat$Structural_Model_Type <- "raw_volume_plus_covariate_adjustment"
df_pat$Thal_Composite_Method <- "sum_of_thalamic_nuclei"


## =========================
## 5. 通用运算核心函数
## =========================

check_n_fun <- function(df_pat, y, main_covars, hippo_var, thal_var) {
  vars <- c(y, main_covars, hippo_var, thal_var)
  vars <- intersect(vars, names(df_pat))
  d <- df_pat[, vars, drop = FALSE]
  
  c(
    Outcome_nonmiss = sum(!is.na(df_pat[[y]])),
    Hippo_nonmiss   = sum(!is.na(df_pat[[hippo_var]])),
    Thal_nonmiss    = sum(!is.na(df_pat[[thal_var]])),
    Complete_cases  = sum(complete.cases(d))
  )
}

check_collinearity_fun <- function(df_pat, y, main_covars, hippo_var, thal_var) {
  vars <- c(y, main_covars, hippo_var, thal_var)
  vars <- intersect(vars, names(df_pat))
  d <- df_pat[, vars, drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  
  if (nrow(d) < 3) {
    return(c(
      N = nrow(d),
      cor_Hippo_Thal = NA,
      sd_Hippo = NA,
      sd_Thal = NA
    ))
  }
  
  c(
    N = nrow(d),
    cor_Hippo_Thal = cor(d[[hippo_var]], d[[thal_var]]),
    sd_Hippo = sd(d[[hippo_var]]),
    sd_Thal  = sd(d[[thal_var]])
  )
}

run_step3_one_set <- function(df_pat, outcomes, main_covars, hippo_var, thal_var,
                              set_name, n_perm = 5000, seed = 20260407, min_n = 15) {
  
  outcomes <- intersect(outcomes, names(df_pat))
  out_list <- list()
  
  for (y_var in outcomes) {
    message("正在运行：", set_name, " | outcome = ", y_var)
    
    vars_need <- unique(c(y_var, main_covars, hippo_var, thal_var))
    dsub <- df_pat[, intersect(vars_need, names(df_pat)), drop = FALSE]
    dsub <- dsub[complete.cases(dsub), , drop = FALSE]
    
    use_covars <- drop_constant_vars(dsub, main_covars)
    
    if (nrow(dsub) < min_n) {
      out_list[[y_var]] <- data.frame(
        Set = set_name,
        Outcome = y_var,
        HippoVar = hippo_var,
        ThalVar = thal_var,
        N = nrow(dsub),
        R2_M0 = NA,
        R2_M1 = NA,
        R2_M2 = NA,
        R2_M3 = NA,
        DeltaR2_M1_vs_M0 = NA,
        DeltaR2_M2_vs_M0 = NA,
        DeltaR2_M3_vs_M1 = NA,
        DeltaR2_M3_vs_M2 = NA,
        PartialR2_Hippo_alone = NA,
        PartialR2_Thal_alone = NA,
        PartialR2_Hippo_in_M3 = NA,
        PartialR2_Thal_in_M3 = NA,
        Beta_Hippo_M3 = NA,
        Beta_Thal_M3 = NA,
        M3_Hippo_Estimable = NA,
        M3_Thal_Estimable = NA,
        Pperm_M1_vs_M0 = NA,
        Pperm_M2_vs_M0 = NA,
        Pperm_M3_vs_M1 = NA,
        Pperm_M3_vs_M2 = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    
    vars_m0 <- use_covars
    vars_m1 <- unique(c(use_covars, hippo_var))
    vars_m2 <- unique(c(use_covars, thal_var))
    vars_m3 <- unique(c(use_covars, hippo_var, thal_var))
    
    fit_m0 <- lm(
      as.formula(paste0(bt(y_var), " ~ ", make_rhs(vars_m0))),
      data = dsub
    )
    
    fit_m1 <- lm(
      as.formula(paste0(bt(y_var), " ~ ", make_rhs(vars_m1))),
      data = dsub
    )
    
    fit_m2 <- lm(
      as.formula(paste0(bt(y_var), " ~ ", make_rhs(vars_m2))),
      data = dsub
    )
    
    fit_m3 <- lm(
      as.formula(paste0(bt(y_var), " ~ ", make_rhs(vars_m3))),
      data = dsub
    )
    
    r2_m0 <- summary(fit_m0)$r.squared
    r2_m1 <- summary(fit_m1)$r.squared
    r2_m2 <- summary(fit_m2)$r.squared
    r2_m3 <- summary(fit_m3)$r.squared
    
    delta_m1_m0 <- r2_m1 - r2_m0
    delta_m2_m0 <- r2_m2 - r2_m0
    delta_m3_m1 <- r2_m3 - r2_m1
    delta_m3_m2 <- r2_m3 - r2_m2
    
    beta_hippo_m3 <- get_coef_safe(fit_m3, hippo_var)
    beta_thal_m3  <- get_coef_safe(fit_m3, thal_var)
    
    coef_names_m3 <- names(coef(fit_m3))
    m3_hippo_estimable <- hippo_var %in% coef_names_m3 && !is.na(beta_hippo_m3)
    m3_thal_estimable  <- thal_var %in% coef_names_m3 && !is.na(beta_thal_m3)
    
    perm_m1 <- freedman_lane_delta_r2(
      dsub, y_var, vars_m0, vars_m1,
      n_perm = n_perm,
      seed = seed + 1,
      min_n = min_n
    )
    
    perm_m2 <- freedman_lane_delta_r2(
      dsub, y_var, vars_m0, vars_m2,
      n_perm = n_perm,
      seed = seed + 2,
      min_n = min_n
    )
    
    perm_m3_vs_m1 <- freedman_lane_delta_r2(
      dsub, y_var, vars_m1, vars_m3,
      n_perm = n_perm,
      seed = seed + 3,
      min_n = min_n
    )
    
    perm_m3_vs_m2 <- freedman_lane_delta_r2(
      dsub, y_var, vars_m2, vars_m3,
      n_perm = n_perm,
      seed = seed + 4,
      min_n = min_n
    )
    
    out_list[[y_var]] <- data.frame(
      Set = set_name,
      Outcome = y_var,
      HippoVar = hippo_var,
      ThalVar = thal_var,
      N = nrow(dsub),
      R2_M0 = r2_m0,
      R2_M1 = r2_m1,
      R2_M2 = r2_m2,
      R2_M3 = r2_m3,
      DeltaR2_M1_vs_M0 = delta_m1_m0,
      DeltaR2_M2_vs_M0 = delta_m2_m0,
      DeltaR2_M3_vs_M1 = delta_m3_m1,
      DeltaR2_M3_vs_M2 = delta_m3_m2,
      PartialR2_Hippo_alone = delta_m1_m0,
      PartialR2_Thal_alone = delta_m2_m0,
      PartialR2_Hippo_in_M3 = delta_m3_m2,
      PartialR2_Thal_in_M3 = delta_m3_m1,
      Beta_Hippo_M3 = beta_hippo_m3,
      Beta_Thal_M3 = beta_thal_m3,
      M3_Hippo_Estimable = m3_hippo_estimable,
      M3_Thal_Estimable = m3_thal_estimable,
      Pperm_M1_vs_M0 = perm_m1$p_perm,
      Pperm_M2_vs_M0 = perm_m2$p_perm,
      Pperm_M3_vs_M1 = perm_m3_vs_m1$p_perm,
      Pperm_M3_vs_M2 = perm_m3_vs_m2$p_perm,
      stringsAsFactors = FALSE
    )
  }
  
  out_df <- do.call(rbind, out_list)
  
  out_df$Q_M1_vs_M0 <- p.adjust(out_df$Pperm_M1_vs_M0, method = "BH")
  out_df$Q_M2_vs_M0 <- p.adjust(out_df$Pperm_M2_vs_M0, method = "BH")
  out_df$Q_M3_vs_M1 <- p.adjust(out_df$Pperm_M3_vs_M1, method = "BH")
  out_df$Q_M3_vs_M2 <- p.adjust(out_df$Pperm_M3_vs_M2, method = "BH")
  
  out_df
}


## =========================
## 6. 运行4组模型
## =========================

model_sets <- list(
  list(
    set_name = "Primary_Ipsi_Global",
    hippo_var = "Hippo_Ipsi_Volume",
    thal_var = "Thal_Global_SummedVolume"
  ),
  list(
    set_name = "Primary_Contra_Global",
    hippo_var = "Hippo_Contra_Volume",
    thal_var = "Thal_Global_SummedVolume"
  ),
  list(
    set_name = "Sensitivity_Ipsi_Ipsi",
    hippo_var = "Hippo_Ipsi_Volume",
    thal_var = "Thal_Ipsi_SummedVolume"
  ),
  list(
    set_name = "Sensitivity_Contra_Ipsi",
    hippo_var = "Hippo_Contra_Volume",
    thal_var = "Thal_Ipsi_SummedVolume"
  )
)

all_sets_results <- list()

for (ms in model_sets) {
  cat("正在计算模型:", ms$set_name, "...\n")
  
  all_sets_results[[ms$set_name]] <- run_step3_one_set(
    df_pat = df_pat,
    outcomes = OUTCOMES_PRIMARY,
    main_covars = MAIN_COVARS,
    hippo_var = ms$hippo_var,
    thal_var = ms$thal_var,
    set_name = ms$set_name,
    n_perm = N_PERM,
    seed = SEED
  )
}

step3_all_sets_results <- do.call(rbind, all_sets_results)

write_csv_cn(step3_all_sets_results, "volume_primary_step3_all_sets_results.csv")
write_csv_cn(step3_all_sets_results, "volume_primary_step3_all_models_results.csv")


## =========================
## 7. N诊断与共线性诊断
## =========================

all_n_diag <- list()
all_col_diag <- list()

for (ms in model_sets) {
  tmp_n <- t(sapply(OUTCOMES_PRIMARY, function(y) {
    check_n_fun(df_pat, y, MAIN_COVARS, ms$hippo_var, ms$thal_var)
  }))
  
  tmp_n <- data.frame(
    Set = ms$set_name,
    Outcome = rownames(tmp_n),
    tmp_n,
    row.names = NULL
  )
  
  all_n_diag[[ms$set_name]] <- tmp_n
  
  tmp_col <- t(sapply(OUTCOMES_PRIMARY, function(y) {
    check_collinearity_fun(df_pat, y, MAIN_COVARS, ms$hippo_var, ms$thal_var)
  }))
  
  tmp_col <- data.frame(
    Set = ms$set_name,
    Outcome = rownames(tmp_col),
    tmp_col,
    row.names = NULL
  )
  
  all_col_diag[[ms$set_name]] <- tmp_col
}

step3_n_diagnostics_all_sets <- do.call(rbind, all_n_diag)
step3_collinearity_check_all_sets <- do.call(rbind, all_col_diag)

write_csv_cn(step3_n_diagnostics_all_sets, "volume_primary_step3_n_diagnostics_all_sets.csv")
write_csv_cn(step3_collinearity_check_all_sets, "volume_primary_step3_collinearity_check_all_sets.csv")


## =========================
## 8. 主文汇总表和宽表
## =========================

library(dplyr)
library(tidyr)
library(ggplot2)
library(flextable)
library(officer)

KEY_OUTCOMES <- intersect(c(
  "工作记忆指数WMI",
  "总记忆商",
  "听觉记忆指数",
  "视觉记忆指数",
  "即刻记忆指数",
  "延迟记忆指数"
), OUTCOMES_PRIMARY)

get_set_table <- function(df, set_name, suffix) {
  sub <- df[df$Set == set_name, c(
    "Outcome", "N",
    "R2_M0", "R2_M1", "R2_M2", "R2_M3",
    "DeltaR2_M1_vs_M0", "Q_M1_vs_M0",
    "DeltaR2_M2_vs_M0", "Q_M2_vs_M0",
    "DeltaR2_M3_vs_M1", "Q_M3_vs_M1",
    "DeltaR2_M3_vs_M2", "Q_M3_vs_M2",
    "Beta_Hippo_M3", "Beta_Thal_M3"
  ), drop = FALSE]
  
  names(sub)[names(sub) != "Outcome"] <- paste0(
    names(sub)[names(sub) != "Outcome"],
    "_",
    suffix
  )
  
  sub
}

summary_wide <- get_set_table(
  step3_all_sets_results,
  "Primary_Ipsi_Global",
  "Primary1"
)

summary_wide <- merge(
  summary_wide,
  get_set_table(step3_all_sets_results, "Primary_Contra_Global", "Primary2"),
  by = "Outcome",
  all = TRUE
)

summary_wide <- merge(
  summary_wide,
  get_set_table(step3_all_sets_results, "Sensitivity_Ipsi_Ipsi", "Sens1"),
  by = "Outcome",
  all = TRUE
)

summary_wide <- merge(
  summary_wide,
  get_set_table(step3_all_sets_results, "Sensitivity_Contra_Ipsi", "Sens2"),
  by = "Outcome",
  all = TRUE
)

write_csv_cn(summary_wide, "volume_primary_step3_summary_wide_all_sets.csv")

paper_table <- summary_wide %>%
  dplyr::filter(Outcome %in% KEY_OUTCOMES)

write_csv_cn(paper_table, "volume_primary_step3_summary_for_manuscript.csv")


## =========================
## 9. 派生指标导出
## =========================

derived_vars_to_save <- c(
  "side_lr",
  "Duration_Dis_num",
  "Structural_Model_Type",
  "Thal_Composite_Method",
  "Hippo_Left_Volume",
  "Hippo_Right_Volume",
  "Thal_Left_SummedVolume",
  "Thal_Right_SummedVolume",
  "Hippo_Ipsi_Volume",
  "Hippo_Contra_Volume",
  "Thal_Ipsi_SummedVolume",
  "Thal_Contra_SummedVolume",
  "Thal_Global_SummedVolume",
  "Age",
  "Edu_year",
  "sex_f",
  "eTIV"
)

derived_vars_to_save <- intersect(derived_vars_to_save, names(df_pat))

write_csv_cn(
  df_pat[, derived_vars_to_save, drop = FALSE],
  "volume_primary_step3_derived_indices.csv"
)


## =========================
## 10. Word三线表
## =========================

outcome_map <- c(
  # "WAIS-总智商"       = "WAIS FSIQ",
  # "言语理解指数VCI"   = "WAIS VCI",
  # "知觉推理指数PRI"   = "WAIS PRI",
  # "加工速度指数PSI"   = "WAIS PSI",
  "工作记忆指数WMI" = "WAIS WMI",
  "总记忆商"       = "WMS FSMQ",
  "听觉记忆指数"   = "WMS AMI",
  "视觉记忆指数"   = "WMS VMI",
  "即刻记忆指数"   = "WMS IMI",
  "延迟记忆指数"   = "WMS DMI"
)

set_map_table <- c(
  "Primary_Ipsi_Global"     = "Primary model set: ipsilateral hippocampal volume; global thalamic summed volume",
  "Primary_Contra_Global"   = "Primary model set: contralateral hippocampal volume; global thalamic summed volume",
  "Sensitivity_Ipsi_Ipsi"   = "Sensitivity model set: ipsilateral hippocampal volume; ipsilateral thalamic summed volume",
  "Sensitivity_Contra_Ipsi" = "Sensitivity model set: contralateral hippocampal volume; ipsilateral thalamic summed volume"
)

set_map_fig <- c(
  "Primary_Ipsi_Global"     = "Primary model set\nHippo: ipsilateral volume; Thal: global summed volume",
  "Primary_Contra_Global"   = "Primary model set\nHippo: contralateral volume; Thal: global summed volume",
  "Sensitivity_Ipsi_Ipsi"   = "Sensitivity model set\nHippo: ipsilateral volume; Thal: ipsilateral summed volume",
  "Sensitivity_Contra_Ipsi" = "Sensitivity model set\nHippo: contralateral volume; Thal: ipsilateral summed volume"
)

table_for_word <- step3_all_sets_results %>%
  mutate(
    Outcome = recode(Outcome, !!!outcome_map),
    Set = recode(Set, !!!set_map_table)
  ) %>%
  arrange(match(Set, unname(set_map_table)), Outcome) %>%
  select(
    Set,
    Outcome,
    N,
    DeltaR2_M1_vs_M0,
    Q_M1_vs_M0,
    DeltaR2_M2_vs_M0,
    Q_M2_vs_M0,
    DeltaR2_M3_vs_M1,
    Q_M3_vs_M1,
    DeltaR2_M3_vs_M2,
    Q_M3_vs_M2
  )

names(table_for_word) <- c(
  "Model Configuration",
  "Cognitive Domain",
  "N",
  "ΔR² (Hippo alone)",
  "Q-value (Hippo alone)",
  "ΔR² (Thal alone)",
  "Q-value (Thal alone)",
  "ΔR² (Thal beyond Hippo)",
  "Q-value (Thal beyond Hippo)",
  "ΔR² (Hippo beyond Thal)",
  "Q-value (Hippo beyond Thal)"
)

format_num <- function(x) {
  x_num <- suppressWarnings(as.numeric(x))
  
  ifelse(
    is.na(x_num),
    "",
    ifelse(x_num < 0.001, "< 0.001", sprintf("%.3f", x_num))
  )
}

table_for_word[, -c(1, 2, 3)] <- lapply(table_for_word[, -c(1, 2, 3)], format_num)

ft <- flextable(table_for_word) %>%
  theme_vanilla() %>%
  merge_v(j = "Model Configuration") %>%
  valign(j = "Model Configuration", valign = "top") %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  align(j = c(1, 2), align = "left", part = "body") %>%
  set_caption(caption = "Volume-based additive models across four configurations") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all")

save_as_docx(
  ft,
  path = file.path(OUT_DIR, "VolumePrimary_Comprehensive_Regression_Table.docx")
)


## =========================
## 11. 综合图：普通条形图，柔和红蓝色系
## =========================

plot_data_en <- step3_all_sets_results %>%
  mutate(
    Outcome = recode(Outcome, !!!outcome_map),
    Set = factor(recode(Set, !!!set_map_fig), levels = unname(set_map_fig))
  ) %>%
  select(
    Set,
    Outcome,
    Hippo_alone = DeltaR2_M1_vs_M0,
    Hippo_alone_P = Pperm_M1_vs_M0,
    Hippo_alone_Q = Q_M1_vs_M0,
    Thal_alone = DeltaR2_M2_vs_M0,
    Thal_alone_P = Pperm_M2_vs_M0,
    Thal_alone_Q = Q_M2_vs_M0,
    Thal_beyond_Hippo = DeltaR2_M3_vs_M1,
    Thal_beyond_Hippo_P = Pperm_M3_vs_M1,
    Thal_beyond_Hippo_Q = Q_M3_vs_M1,
    Hippo_beyond_Thal = DeltaR2_M3_vs_M2,
    Hippo_beyond_Thal_P = Pperm_M3_vs_M2,
    Hippo_beyond_Thal_Q = Q_M3_vs_M2
  ) %>%
  pivot_longer(
    cols = c(
      Hippo_alone,
      Hippo_beyond_Thal,
      Thal_alone,
      Thal_beyond_Hippo
    ),
    names_to = "Comparison",
    values_to = "DeltaR2"
  ) %>%
  mutate(
    P_value = case_when(
      Comparison == "Hippo_alone" ~ Hippo_alone_P,
      Comparison == "Hippo_beyond_Thal" ~ Hippo_beyond_Thal_P,
      Comparison == "Thal_alone" ~ Thal_alone_P,
      Comparison == "Thal_beyond_Hippo" ~ Thal_beyond_Hippo_P,
      TRUE ~ NA_real_
    ),
    Q_value = case_when(
      Comparison == "Hippo_alone" ~ Hippo_alone_Q,
      Comparison == "Hippo_beyond_Thal" ~ Hippo_beyond_Thal_Q,
      Comparison == "Thal_alone" ~ Thal_alone_Q,
      Comparison == "Thal_beyond_Hippo" ~ Thal_beyond_Hippo_Q,
      TRUE ~ NA_real_
    ),
    Comparison = recode(
      Comparison,
      "Hippo_alone" = "Hippo alone (M1−M0)",
      "Hippo_beyond_Thal" = "Hippo beyond Thal (M3−M2)",
      "Thal_alone" = "Thal alone (M2−M0)",
      "Thal_beyond_Hippo" = "Thal beyond Hippo (M3−M1)"
    ),
    Comparison = factor(
      Comparison,
      levels = c(
        "Hippo alone (M1−M0)",
        "Hippo beyond Thal (M3−M2)",
        "Thal alone (M2−M0)",
        "Thal beyond Hippo (M3−M1)"
      )
    ),
    Significance = case_when(
      is.na(Q_value)   ~ "",
      Q_value < 0.001  ~ "***",
      Q_value < 0.01   ~ "**",
      Q_value < 0.05   ~ "*",
      Q_value < 0.10   ~ "○",
      TRUE             ~ ""
    ),
    Label_text = sprintf("P=%.3f, Q=%.3f%s", P_value, Q_value, Significance)
  ) %>%
  mutate(
    Label_text = gsub("P=0.000", "P<0.001", Label_text),
    Label_text = gsub("Q=0.000", "Q<0.001", Label_text)
  )

comparison_colors <- c(
  "Hippo alone (M1−M0)" = "#D7E3F0",
  "Hippo beyond Thal (M3−M2)" = "#7FA6C9",
  "Thal alone (M2−M0)" = "#F0D9DD",
  "Thal beyond Hippo (M3−M1)" = "#C98A95"
)

max_delta <- max(plot_data_en$DeltaR2, na.rm = TRUE)

p_en_pq <- ggplot(
  plot_data_en,
  aes(x = reorder(Outcome, DeltaR2), y = DeltaR2, fill = Comparison)
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.4
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "black",
    linewidth = 0.28,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = Label_text, y = DeltaR2 + max_delta * 0.018),
    position = position_dodge(width = 0.82),
    size = 3.7,
    color = "black",
    family = "serif",
    fontface = "bold",
    hjust = 0
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.48))
  ) +
  scale_fill_manual(values = comparison_colors) +
  facet_wrap(~ Set, ncol = 2, scales = "free_x") +
  labs(
    title = "Volume-based Analysis: Incremental Variance Explained Across Model Sets",
    x = "Cognitive and Memory Domains",
    y = expression("Incremental Explained Variance (" * Delta * R^2 * ")")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11.5, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(size = 11.5, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  file.path(OUT_DIR, "VolumePrimary_Comprehensive_Incremental_R2_Panel.png"),
  plot = p_en_pq,
  width = 14,
  height = 14,
  dpi = 600
)

ggsave(
  file.path(OUT_DIR, "VolumePrimary_Comprehensive_Incremental_R2_Panel.pdf"),
  plot = p_en_pq,
  width = 14,
  height = 14
)

cat(">>> [完成] Volume-based 综合 2×2 面板图已导出。\n")


## =========================
## 12. 单独生成 Primary 结果图
## =========================

primary_sets <- c("Primary_Ipsi_Global", "Primary_Contra_Global")

plot_data_primary <- step3_all_sets_results %>%
  filter(Set %in% primary_sets) %>%
  mutate(
    Outcome = recode(Outcome, !!!outcome_map),
    Set = factor(recode(Set, !!!set_map_fig), levels = unname(set_map_fig[primary_sets]))
  ) %>%
  select(
    Set,
    Outcome,
    Hippo_alone = DeltaR2_M1_vs_M0,
    Hippo_alone_P = Pperm_M1_vs_M0,
    Hippo_alone_Q = Q_M1_vs_M0,
    Thal_alone = DeltaR2_M2_vs_M0,
    Thal_alone_P = Pperm_M2_vs_M0,
    Thal_alone_Q = Q_M2_vs_M0,
    Thal_beyond_Hippo = DeltaR2_M3_vs_M1,
    Thal_beyond_Hippo_P = Pperm_M3_vs_M1,
    Thal_beyond_Hippo_Q = Q_M3_vs_M1,
    Hippo_beyond_Thal = DeltaR2_M3_vs_M2,
    Hippo_beyond_Thal_P = Pperm_M3_vs_M2,
    Hippo_beyond_Thal_Q = Q_M3_vs_M2
  ) %>%
  pivot_longer(
    cols = c(
      Hippo_alone,
      Hippo_beyond_Thal,
      Thal_alone,
      Thal_beyond_Hippo
    ),
    names_to = "Comparison",
    values_to = "DeltaR2"
  ) %>%
  mutate(
    P_value = case_when(
      Comparison == "Hippo_alone" ~ Hippo_alone_P,
      Comparison == "Hippo_beyond_Thal" ~ Hippo_beyond_Thal_P,
      Comparison == "Thal_alone" ~ Thal_alone_P,
      Comparison == "Thal_beyond_Hippo" ~ Thal_beyond_Hippo_P,
      TRUE ~ NA_real_
    ),
    Q_value = case_when(
      Comparison == "Hippo_alone" ~ Hippo_alone_Q,
      Comparison == "Hippo_beyond_Thal" ~ Hippo_beyond_Thal_Q,
      Comparison == "Thal_alone" ~ Thal_alone_Q,
      Comparison == "Thal_beyond_Hippo" ~ Thal_beyond_Hippo_Q,
      TRUE ~ NA_real_
    ),
    Comparison = recode(
      Comparison,
      "Hippo_alone" = "Hippo alone (M1−M0)",
      "Hippo_beyond_Thal" = "Hippo beyond Thal (M3−M2)",
      "Thal_alone" = "Thal alone (M2−M0)",
      "Thal_beyond_Hippo" = "Thal beyond Hippo (M3−M1)"
    ),
    Comparison = factor(
      Comparison,
      levels = c(
        "Hippo alone (M1−M0)",
        "Hippo beyond Thal (M3−M2)",
        "Thal alone (M2−M0)",
        "Thal beyond Hippo (M3−M1)"
      )
    ),
    Significance = case_when(
      is.na(Q_value)   ~ "",
      Q_value < 0.001  ~ "***",
      Q_value < 0.01   ~ "**",
      Q_value < 0.05   ~ "*",
      Q_value < 0.10   ~ "○",
      TRUE             ~ ""
    ),
    Label_text = sprintf("P=%.3f, Q=%.3f%s", P_value, Q_value, Significance)
  ) %>%
  mutate(
    Label_text = gsub("P=0.000", "P<0.001", Label_text),
    Label_text = gsub("Q=0.000", "Q<0.001", Label_text)
  )

max_delta_primary <- max(plot_data_primary$DeltaR2, na.rm = TRUE)

p_primary <- ggplot(
  plot_data_primary,
  aes(x = reorder(Outcome, DeltaR2), y = DeltaR2, fill = Comparison)
) +
  geom_hline(
    yintercept = 0.05,
    linetype = "dashed",
    color = "gray60",
    linewidth = 0.4
  ) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.82),
    width = 0.72,
    color = "black",
    linewidth = 0.28,
    alpha = 0.95
  ) +
  geom_text(
    aes(label = Label_text, y = DeltaR2 + max_delta_primary * 0.018),
    position = position_dodge(width = 0.82),
    size = 3.7,
    color = "black",
    family = "serif",
    fontface = "bold",
    hjust = 0
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.48))
  ) +
  scale_fill_manual(values = comparison_colors) +
  facet_wrap(~ Set, ncol = 2, scales = "free_x") +
  labs(
    title = "Volume-based Primary Model Sets: Incremental Variance Explained",
    x = "Cognitive and Memory Domains",
    y = expression("Incremental Explained Variance (" * Delta * R^2 * ")")
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "serif"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(size = 11.5, color = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    strip.background = element_rect(fill = "grey92", color = "black"),
    strip.text = element_text(size = 11.5, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

ggsave(
  file.path(OUT_DIR, "VolumePrimary_Primary_Incremental_R2_Panel.png"),
  plot = p_primary,
  width = 12,
  height = 10,
  dpi = 600
)

ggsave(
  file.path(OUT_DIR, "VolumePrimary_Primary_Incremental_R2_Panel.pdf"),
  plot = p_primary,
  width = 12,
  height = 10
)

cat(">>> 已导出：VolumePrimary_Primary_Incremental_R2_Panel.png\n")


## =========================
## 13. 控制台输出
## =========================

cat("\n========================================\n")
cat("Volume-based Step 3 主分析已运行完成。\n")
cat("输出目录：", OUT_DIR, "\n")
cat("1. 完整总表：volume_primary_step3_all_sets_results.csv\n")
cat("2. 兼容总表：volume_primary_step3_all_models_results.csv\n")
cat("3. 宽表：volume_primary_step3_summary_wide_all_sets.csv\n")
cat("4. 主文表：volume_primary_step3_summary_for_manuscript.csv\n")
cat("5. N诊断：volume_primary_step3_n_diagnostics_all_sets.csv\n")
cat("6. 共线性诊断：volume_primary_step3_collinearity_check_all_sets.csv\n")
cat("7. 派生指标：volume_primary_step3_derived_indices.csv\n")
cat("8. Word表：VolumePrimary_Comprehensive_Regression_Table.docx\n")
cat("9. 综合图：VolumePrimary_Comprehensive_Incremental_R2_Panel.png\n")
cat("10. Primary图：VolumePrimary_Primary_Incremental_R2_Panel.png\n")
cat("========================================\n\n")
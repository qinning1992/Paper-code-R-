# ==============================================================================
#  TABLE 1+ figure 1 A-D
# HC + All patients + HS + Other
# - HC and All patients are descriptive only
# - Statistical comparison is HS vs Other only
# - Continuous variables: Welch t-test (fixed, not data-driven switching)
# - Categorical variables: Fisher's exact test; if RxC Fisher fails, fallback to
#   chi-square with Monte Carlo simulation
# - Continuous effect size: Hedges g with 95% CI
# - Categorical effect size: Cramér's V
# - Variable-specific n is embedded in each cell
# - Abbreviations in table body; full names provided in note
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(effectsize)
  library(officer)
  library(flextable)
})

# -------------------------
# 0) Checks + setup
# -------------------------
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found.")
if(!exists("df_hc_clean"))  stop("❌ df_hc_clean not found.")

df_pat <- df_raw_clean
df_hc  <- df_hc_clean

OUTDIR <- "C:/Users/86150/Documents/HIPP/TABLES_Demographics"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# 1) Helper functions
# -------------------------
to_num <- function(x){
  if(is.numeric(x)) return(x)
  if(is.factor(x)) x <- as.character(x)
  if(is.character(x)){
    x <- str_replace_all(x, "[,，]", "")
    x <- str_trim(x)
    suppressWarnings(as.numeric(x))
  } else {
    suppressWarnings(as.numeric(x))
  }
}

pick_first_exist <- function(cands, df){
  hit <- cands[cands %in% names(df)]
  if(length(hit) == 0) return(NA_character_)
  hit[1]
}

as_three_line <- function(ft){
  ft <- autofit(ft)
  ft <- theme_vanilla(ft)
  ft <- border_remove(ft)
  ft <- hline_top(ft, border = fp_border(width = 1.2))
  ft <- hline(ft, i = 1, part = "header", border = fp_border(width = 1.0))
  ft <- hline_bottom(ft, border = fp_border(width = 1.2))
  ft
}

fmt_p <- function(p){
  ifelse(
    is.na(p), "",
    ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
  )
}

fmt_mean_sd_n <- function(x, digits = 2){
  x <- x[is.finite(x)]
  n <- length(x)
  if(n == 0) return("")
  paste0(
    sprintf(paste0("%.", digits, "f ± %.", digits, "f"), mean(x), sd(x)),
    " (n=", n, ")"
  )
}

fmt_nN_pct <- function(n, N, digits = 1){
  if(is.na(n) || is.na(N) || N == 0) return("")
  sprintf("%d/%d (%.1f%%)", n, N, 100 * n / N)
}

fmt_cat_dist <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0) return("")
  tb <- table(x, useNA = "no")
  N  <- sum(tb)
  paste0(
    names(tb), ": ", as.integer(tb), "/", N,
    " (", sprintf("%.1f", 100 * as.integer(tb) / N), "%)",
    collapse = "; "
  )
}

cramers_v <- function(tab){
  suppressWarnings({
    r <- nrow(tab); c <- ncol(tab)
    n <- sum(tab)
    if(n == 0 || r < 2 || c < 2) return(NA_real_)
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
    sqrt(as.numeric(chi$statistic) / (n * (min(r - 1, c - 1))))
  })
}

safe_cat_test <- function(tab){
  # returns list(test, p)
  if(any(dim(tab) < 2) || sum(tab) == 0){
    return(list(test = "NA", p = NA_real_))
  }
  
  out <- tryCatch(
    fisher.test(tab),
    error = function(e) NULL
  )
  
  if(!is.null(out)){
    return(list(test = "Fisher's exact", p = out$p.value))
  }
  
  # fallback for larger / computationally difficult tables
  out2 <- tryCatch(
    chisq.test(tab, simulate.p.value = TRUE, B = 20000),
    error = function(e) NULL
  )
  
  if(!is.null(out2)){
    return(list(test = "Chi-square (Monte Carlo)", p = out2$p.value))
  }
  
  list(test = "NA", p = NA_real_)
}

safe_hedges_g <- function(x, g){
  df <- data.frame(x = x, g = g) %>%
    filter(is.finite(x), !is.na(g))
  
  if(nrow(df) < 4) {
    return(list(g = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  
  out <- tryCatch(
    effectsize::hedges_g(x ~ g, data = df, pooled_sd = FALSE, ci = 0.95),
    error = function(e) NULL
  )
  
  if(is.null(out)){
    return(list(g = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  
  g_val <- if("Hedges_g" %in% names(out)) out$Hedges_g[1] else NA_real_
  ci_low <- if("CI_low" %in% names(out)) out$CI_low[1] else NA_real_
  ci_high <- if("CI_high" %in% names(out)) out$CI_high[1] else NA_real_
  
  list(g = g_val, ci_low = ci_low, ci_high = ci_high)
}

fmt_g_ci <- function(g, lo, hi){
  if(!is.finite(g)) return("")
  if(is.finite(lo) && is.finite(hi)){
    sprintf("Hedges g=%.2f [%.2f, %.2f]", g, lo, hi)
  } else {
    sprintf("Hedges g=%.2f", g)
  }
}

# -------------------------
# 2) User-editable coding section
# IMPORTANT: verify sex coding here
# -------------------------
SEX_MALE_CODES   <- c("1", "m", "male", "男")
SEX_FEMALE_CODES <- c("2", "f", "female", "女", "0")

is_male <- function(x){
  if(is.factor(x)) x <- as.character(x)
  xl <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(xl))
  out[xl %in% tolower(SEX_MALE_CODES)]   <- 1L
  out[xl %in% tolower(SEX_FEMALE_CODES)] <- 0L
  out
}

# -------------------------
# 3) Variable labels in table body
# Use abbreviations only, as requested
# -------------------------
label_map <- c(
  "sex" = "Sex",
  "Age" = "Age",
  "Edu_year" = "Education",
  "Duration_Dis" = "Disease duration",
  "HADS-D" = "HADS-D",
  "HADS-A" = "HADS-A",
  "生活质量" = "QoL",
  "BNT-命名" = "BNT",
  "WAIS-总智商" = "WAIS FSIQ",
  "言语理解指数VCI" = "WAIS VCI",
  "知觉推理指数PRI" = "WAIS PRI",
  "工作记忆指数WMI" = "WAIS WMI",
  "加工速度指数PSI" = "WAIS PSI",
  "总记忆商" = "WMS FSMQ",
  "听觉记忆指数" = "WMS AMI",
  "视觉记忆指数" = "WMS VMI",
  "即刻记忆指数" = "WMS IMI",
  "延迟记忆指数" = "WMS DMI",
  "VA-immediate" = "VA Immediate",
  "VA-recall" = "VA Delay",
  "VTrials I-V Total" = "VTrials I–V Total",
  "CAVLT-slope" = "CAVLT Slope",
  "FA-immediate" = "FA Immediate",
  "FA-recall" = "FA Delay",
  "FTrials I-V Total" = "FTrials I–V Total",
  "AFLT-slope" = "AFLT Slope",
  "IID侧别（单=1；双=2）" = "IID laterality",
  "病灶侧IID分布（1=T ant, SP1/F7/T3；2=T front, FP1/F3；3=T5）" = "IID distribution on lesion side"
)

var_label <- function(v){
  if(v %in% names(label_map)) return(unname(label_map[[v]]))
  v
}

# -------------------------
# 4) Variables to include
# -------------------------
cog_vars <- c(
  "WAIS-总智商",
  "言语理解指数VCI",
  "知觉推理指数PRI",
  "工作记忆指数WMI",
  "加工速度指数PSI",
  "总记忆商",
  "听觉记忆指数",
  "视觉记忆指数",
  "即刻记忆指数",
  "延迟记忆指数",
  "BNT-命名",
  "VA-immediate",
  "VA-recall",
  "VTrials I-V Total",
  "CAVLT-slope",
  "FA-immediate",
  "FA-recall",
  "FTrials I-V Total",
  "AFLT-slope"
)

pat_vars <- c(
  "IID侧别（单=1；双=2）",
  "病灶侧IID分布（1=T ant, SP1/F7/T3；2=T front, FP1/F3；3=T5）",
  "Age",
  "Duration_Dis",
  "Edu_year",
  "sex",
  "HADS-D",
  "HADS-A",
  "生活质量",
  cog_vars
)
pat_vars <- unique(pat_vars)

# -------------------------
# 5) Build HS vs Other groups
# -------------------------
if(!("Pathology_Group" %in% names(df_pat))) stop("❌ df_raw_clean missing: Pathology_Group")

df_pat2 <- df_pat %>%
  mutate(Pathology_Group = to_num(.data[["Pathology_Group"]])) %>%
  mutate(PathGroup2 = case_when(
    Pathology_Group == 1 ~ "HS",
    !is.na(Pathology_Group) & Pathology_Group != 1 ~ "Other",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(PathGroup2)) %>%
  mutate(PathGroup2 = factor(PathGroup2, levels = c("HS", "Other")))

# keep variables existing in patient data
vars2_exist <- pat_vars[pat_vars %in% names(df_pat2)]
vars2_missing <- setdiff(pat_vars, vars2_exist)
if(length(vars2_missing) > 0){
  message("⚠️ These requested variables are missing and will be skipped:\n  - ",
          paste(vars2_missing, collapse = "\n  - "))
}

# -------------------------
# 6) Explicit variable-type rules
# No data-driven switching
# -------------------------
col_sex <- pick_first_exist(c("sex", "Sex"), df_pat2)
col_iid_lat <- pick_first_exist(
  c("IID侧别（单=1；双=2）", "IID laterality (1=unilateral; 2=bilateral)"),
  df_pat2
)
col_iid_dist <- pick_first_exist(
  c("病灶侧IID分布（1=T ant, SP1/F7/T3；2=T front, FP1/F3；3=T5）",
    "IID distribution on lesion side (1=T ant; 2=T front; 3=T5)"),
  df_pat2
)

categorical_vars <- c(col_sex, col_iid_lat, col_iid_dist)
categorical_vars <- categorical_vars[!is.na(categorical_vars)]
continuous_vars <- setdiff(vars2_exist, categorical_vars)

# recoding helpers
recode_iid_lat <- function(x){
  x <- to_num(x)
  out <- rep(NA_character_, length(x))
  out[x == 1] <- "Unilateral"
  out[x == 2] <- "Bilateral"
  factor(out, levels = c("Unilateral", "Bilateral"))
}

recode_iid_dist <- function(x){
  x <- to_num(x)
  out <- rep(NA_character_, length(x))
  out[x == 1] <- "T ant (SP1/F7/T3)"
  out[x == 2] <- "T front (FP1/F3)"
  out[x == 3] <- "T5"
  factor(out, levels = c("T ant (SP1/F7/T3)", "T front (FP1/F3)", "T5"))
}

# -------------------------
# 7) Build one row
# -------------------------
make_row_pat_hc <- function(v){
  
  hc_has_var <- v %in% names(df_hc)
  
  # ---- sex special case: show Male n/N (%)
  if(!is.na(col_sex) && v == col_sex){
    male_pat <- is_male(df_pat2[[col_sex]])
    g <- df_pat2$PathGroup2
    
    ok_all <- !is.na(male_pat) & !is.na(g)
    N_all <- sum(ok_all)
    n_all <- sum(male_pat[ok_all] == 1)
    
    ok_hs <- ok_all & g == "HS"
    ok_ot <- ok_all & g == "Other"
    N_hs <- sum(ok_hs); n_hs <- sum(male_pat[ok_hs] == 1)
    N_ot <- sum(ok_ot); n_ot <- sum(male_pat[ok_ot] == 1)
    
    hc_txt <- ""
    N_hc <- NA_integer_
    if(hc_has_var){
      male_hc <- is_male(df_hc[[col_sex]])
      ok_hc <- !is.na(male_hc)
      N_hc <- sum(ok_hc)
      n_hc <- sum(male_hc[ok_hc] == 1)
      hc_txt <- fmt_nN_pct(n_hc, N_hc)
    }
    
    tab2 <- table(
      g[ok_all],
      factor(male_pat[ok_all], levels = c(1, 0), labels = c("Male", "Not male"))
    )
    
    tst <- safe_cat_test(tab2)
    V <- cramers_v(tab2)
    
    return(tibble::tibble(
      Variable = "Sex (Male)",
      HC = hc_txt,
      `All patients` = fmt_nN_pct(n_all, N_all),
      HS = fmt_nN_pct(n_hs, N_hs),
      Other = fmt_nN_pct(n_ot, N_ot),
      Test = tst$test,
      `Effect size` = if(is.finite(V)) sprintf("Cramér's V=%.2f", V) else "",
      `P value` = fmt_p(tst$p),
      N_HC = N_hc,
      N_All = N_all,
      N_HS = N_hs,
      N_Other = N_ot,
      VarType = "categorical"
    ))
  }
  
  # ---- IID categorical variables
  if(v %in% categorical_vars){
    x_pat <- df_pat2[[v]]
    if(!is.na(col_iid_lat) && v == col_iid_lat)  x_pat <- recode_iid_lat(x_pat)
    if(!is.na(col_iid_dist) && v == col_iid_dist) x_pat <- recode_iid_dist(x_pat)
    
    g <- df_pat2$PathGroup2
    ok <- !is.na(x_pat) & !is.na(g)
    x_ok <- droplevels(factor(x_pat[ok]))
    g_ok <- droplevels(g[ok])
    
    hc_txt <- ""
    N_hc <- NA_integer_
    if(hc_has_var){
      x_hc <- df_hc[[v]]
      if(!is.na(col_iid_lat) && v == col_iid_lat)  x_hc <- recode_iid_lat(x_hc)
      if(!is.na(col_iid_dist) && v == col_iid_dist) x_hc <- recode_iid_dist(x_hc)
      x_hc <- x_hc[!is.na(x_hc)]
      N_hc <- length(x_hc)
      hc_txt <- fmt_cat_dist(x_hc)
    }
    
    tab <- table(g_ok, x_ok)
    tst <- safe_cat_test(tab)
    V <- cramers_v(tab)
    
    n_all <- length(x_ok)
    n_hs <- sum(g_ok == "HS")
    n_ot <- sum(g_ok == "Other")
    
    return(tibble::tibble(
      Variable = var_label(v),
      HC = hc_txt,
      `All patients` = fmt_cat_dist(x_ok),
      HS = fmt_cat_dist(x_ok[g_ok == "HS"]),
      Other = fmt_cat_dist(x_ok[g_ok == "Other"]),
      Test = tst$test,
      `Effect size` = if(is.finite(V)) sprintf("Cramér's V=%.2f", V) else "",
      `P value` = fmt_p(tst$p),
      N_HC = N_hc,
      N_All = n_all,
      N_HS = n_hs,
      N_Other = n_ot,
      VarType = "categorical"
    ))
  }
  
  # ---- continuous variables
  if(v %in% continuous_vars){
    x_pat <- to_num(df_pat2[[v]])
    g <- df_pat2$PathGroup2
    ok <- is.finite(x_pat) & !is.na(g)
    
    x_all <- x_pat[ok]
    x_hs <- x_pat[ok & g == "HS"]
    x_ot <- x_pat[ok & g == "Other"]
    
    hc_txt <- ""
    N_hc <- NA_integer_
    if(hc_has_var){
      x_hc <- to_num(df_hc[[v]])
      x_hc <- x_hc[is.finite(x_hc)]
      N_hc <- length(x_hc)
      hc_txt <- fmt_mean_sd_n(x_hc, digits = 2)
    }
    
    # fixed Welch t-test
    p <- NA_real_
    test_name <- "Welch t-test"
    if(length(x_hs) >= 2 && length(x_ot) >= 2){
      p <- tryCatch(
        t.test(x_hs, x_ot, var.equal = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }
    
    eg <- safe_hedges_g(x = x_pat[ok], g = g[ok])
    
    return(tibble::tibble(
      Variable = var_label(v),
      HC = hc_txt,
      `All patients` = fmt_mean_sd_n(x_all, digits = 2),
      HS = fmt_mean_sd_n(x_hs, digits = 2),
      Other = fmt_mean_sd_n(x_ot, digits = 2),
      Test = test_name,
      `Effect size` = fmt_g_ci(eg$g, eg$ci_low, eg$ci_high),
      `P value` = fmt_p(p),
      N_HC = N_hc,
      N_All = length(x_all),
      N_HS = length(x_hs),
      N_Other = length(x_ot),
      VarType = "continuous"
    ))
  }
  
  # fallback
  tibble::tibble(
    Variable = var_label(v),
    HC = "",
    `All patients` = "",
    HS = "",
    Other = "",
    Test = "",
    `Effect size` = "",
    `P value` = "",
    N_HC = NA_integer_,
    N_All = NA_integer_,
    N_HS = NA_integer_,
    N_Other = NA_integer_,
    VarType = "unknown"
  )
}

# build rows
t2_rows_full <- purrr::map_dfr(vars2_exist, make_row_pat_hc)

# table body for docx/csv
t2_rows <- t2_rows_full %>%
  select(Variable, HC, `All patients`, HS, Other, Test, `Effect size`, `P value`)

# detailed stats file
t2_stats <- t2_rows_full %>%
  arrange(VarType, Variable)

# -------------------------
# 8) Export
# -------------------------
n_hc  <- nrow(df_hc)
n_all <- nrow(df_pat2)
n_hs  <- sum(df_pat2$PathGroup2 == "HS")
n_ot  <- sum(df_pat2$PathGroup2 == "Other")

ft2 <- flextable(t2_rows) %>%
  set_header_labels(
    Variable = "Variable",
    HC = paste0("HC (n=", n_hc, ")"),
    `All patients` = paste0("All patients (n=", n_all, ")"),
    HS = paste0("HS (n=", n_hs, ")"),
    Other = paste0("Other (n=", n_ot, ")"),
    Test = "Test",
    `Effect size` = "Effect size",
    `P value` = "P value"
  ) %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  as_three_line()

abbr_note <- paste(
  "Abbreviations: WAIS = Wechsler Adult Intelligence Scale;",
  "FSIQ = Full Scale IQ;",
  "VCI = Verbal Comprehension Index;",
  "PRI = Perceptual Reasoning Index;",
  "WMI = Working Memory Index;",
  "PSI = Processing Speed Index;",
  "WMS = Wechsler Memory Scale;",
  "FSMQ = Full-Scale Memory Quotient;",
  "AMI = Auditory Memory Index;",
  "VMI = Visual Memory Index;",
  "IMI = Immediate Memory Index;",
  "DMI = Delayed Memory Index;",
  "BNT = Boston Naming Test;",
  "QoL = quality of life;",
  "VA = verbal association;",
  "VTrials I–V Total = total verbal learning across trials I–V;",
  "CAVLT = Chinese Auditory Verbal Learning Test;",
  "FA = figure association;",
  "FTrials I–V Total = total figural learning across trials I–V;",
  "AFLT = Aggie Figures Learning Test."
)

note_text <- paste(
  "Notes: The HC and All patients columns are descriptive only.",
  "Pathology_Group was coded as HS = 1 and Other = all non-1.",
  "Statistical comparisons were performed between HS and Other only.",
  "Continuous variables were compared using Welch's t-test, with Hedges g and 95% confidence intervals as effect sizes.",
  "Categorical variables were compared using Fisher's exact test; if exact computation failed for larger contingency tables, chi-square with Monte Carlo simulation was used.",
  "Cramér's V was reported as the effect size for categorical variables.",
  "Variable-specific effective sample sizes are embedded in each cell rather than inferred only from group header n.",
  "Observed post hoc power was not reported, because it is generally not recommended; variable-specific n and effect-size confidence intervals provide a more informative representation of statistical information.",
  abbr_note
)

doc2 <- read_docx() %>%
  body_add_par(
    "Table 2. Clinical and cognitive characteristics of HC, all patients, HS, and Other groups",
    style = "heading 1"
  ) %>%
  body_add_flextable(ft2) %>%
  body_add_par(note_text, style = "Normal")

print(doc2, target = file.path(OUTDIR, "Table2_HC_AllPatients_HS_Other_revised.docx"))

readr::write_excel_csv(
  t2_rows,
  file.path(OUTDIR, "Table2_HC_AllPatients_HS_Other_revised.csv")
)

readr::write_excel_csv(
  t2_stats,
  file.path(OUTDIR, "Table2_HC_AllPatients_HS_Other_stats_revised.csv")
)

cat("\n✅ Revised Table 2 exported to the same folder:\n",
    " - ", normalizePath(file.path(OUTDIR, "Table2_HC_AllPatients_HS_Other_revised.docx")), "\n",
    " - ", normalizePath(file.path(OUTDIR, "Table2_HC_AllPatients_HS_Other_revised.csv")), "\n",
    " - ", normalizePath(file.path(OUTDIR, "Table2_HC_AllPatients_HS_Other_stats_revised.csv")), "\n",
    sep = "")







############################################################
## figure 1 A-D  HS vs non-HS：所有认知量表组间比较柱状图
## 1）自动识别 HS 或 HS(1=with,0=non with) 列
## 2）比较 HS 与 non-HS 在所有量表上的差异
## 3）按量表类型模块化绘图
## 4）模块内共享 Y 轴，不同模块使用不同 Y 轴尺度
## 5）生成统计结果表和图片
## 6）不使用花纹，仅使用颜色区分 HS / non-HS
## 7）每个子图右上角单独显示图例
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(grid)
})

# ==========================================================
# 0) INPUT
# ==========================================================
if (!exists("df_raw_clean")) stop("❌ df_raw_clean not found in workspace.")

df_pat <- df_raw_clean

# ==========================================================
# 1) 可修改参数
# ==========================================================

# 显著性标识基于哪一个统计量：
# "p_welch"     = Welch t-test 原始 P 值
# "q_module"    = 在每个模块内 FDR 校正后的 Q 值，推荐
# "q_all"       = 所有量表统一 FDR 校正后的 Q 值，更严格
SIG_SOURCE <- "q_module"

# 是否显示 q < 0.10 的空心圆
SHOW_TREND_CIRCLE <- TRUE

# 误差线类型："SE" 或 "SD"
ERROR_TYPE <- "SE"

# 柱状图颜色
HS_COLORS <- c(
  "HS"     = "#C98A95",
  "non-HS" = "#7FA6C9"
)

# 图片尺寸
COMBINED_WIDTH  <- 7.5
COMBINED_HEIGHT <- 8
DPI <- 600

# ==========================================================
# 2) 工具函数
# ==========================================================

to_num <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  suppressWarnings(as.numeric(x))
}

fmt_p <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}

sig_label_from_p <- function(p, show_trend = TRUE) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    show_trend & p < 0.10 ~ "circle",
    TRUE ~ ""
  )
}

cohens_d_two_group <- function(x, g) {
  x1 <- x[g == "HS"]
  x0 <- x[g == "non-HS"]
  
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]
  
  n1 <- length(x1)
  n0 <- length(x0)
  
  if (n1 < 2 || n0 < 2) return(NA_real_)
  
  sd1 <- sd(x1)
  sd0 <- sd(x0)
  
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n0 - 1) * sd0^2) / (n1 + n0 - 2))
  
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  
  (mean(x1) - mean(x0)) / pooled_sd
}

safe_welch_p <- function(x, g) {
  x1 <- x[g == "HS"]
  x0 <- x[g == "non-HS"]
  
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]
  
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  if (sd(x1) == 0 && sd(x0) == 0) return(NA_real_)
  
  out <- tryCatch(
    t.test(x1, x0, var.equal = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  out
}

safe_wilcox_p <- function(x, g) {
  x1 <- x[g == "HS"]
  x0 <- x[g == "non-HS"]
  
  x1 <- x1[is.finite(x1)]
  x0 <- x0[is.finite(x0)]
  
  if (length(x1) < 2 || length(x0) < 2) return(NA_real_)
  
  out <- tryCatch(
    wilcox.test(x1, x0, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  
  out
}

# ==========================================================
# 3) HS / non-HS 分组
# ==========================================================

hs_col <- NULL
if ("HS" %in% names(df_pat)) hs_col <- "HS"
if ("HS(1=with,0=non with)" %in% names(df_pat)) hs_col <- "HS(1=with,0=non with)"
if (is.null(hs_col)) stop("❌ 未找到 HS 列或 HS(1=with,0=non with) 列")

df_pat <- df_pat %>%
  mutate(
    HS_raw = to_num(.data[[hs_col]]),
    HS_group = case_when(
      HS_raw == 1 ~ "HS",
      HS_raw == 0 ~ "non-HS",
      TRUE ~ NA_character_
    ),
    HS_group = factor(HS_group, levels = c("HS", "non-HS"))
  )

cat("✅ HS grouping:\n")
print(table(df_pat$HS_group, useNA = "ifany"))

df_pat <- df_pat %>%
  filter(!is.na(HS_group))

# ==========================================================
# 4) 量表定义与模块分组
# ==========================================================

scale_meta <- tibble::tribble(
  ~Scale, ~Is_Normed,
  "WAIS-总智商", TRUE,
  "言语理解指数VCI", TRUE,
  "知觉推理指数PRI", TRUE,
  "工作记忆指数WMI", TRUE,
  "加工速度指数PSI", TRUE,
  "总记忆商", TRUE,
  "听觉记忆指数", TRUE,
  "视觉记忆指数", TRUE,
  "即刻记忆指数", TRUE,
  "延迟记忆指数", TRUE,
  "VA-immediate", FALSE,
  "VA-recall", FALSE,
  "VTrials I-V Total", FALSE,
  "CAVLT-slope", FALSE,
  "FA-immediate", FALSE,
  "FA-recall", FALSE,
  "FTrials I-V Total", FALSE,
  "AFLT-slope", FALSE
)

scale_map <- c(
  "WAIS-总智商"        = "WAIS FSIQ",
  "言语理解指数VCI"    = "WAIS VCI",
  "知觉推理指数PRI"    = "WAIS PRI",
  "工作记忆指数WMI"    = "WAIS WMI",
  "加工速度指数PSI"    = "WAIS PSI",
  "总记忆商"          = "WMS FSMQ",
  "听觉记忆指数"      = "WMS AMI",
  "视觉记忆指数"      = "WMS VMI",
  "即刻记忆指数"      = "WMS IMI",
  "延迟记忆指数"      = "WMS DMI",
  "VA-immediate"      = "VA Immediate",
  "VA-recall"         = "VA Delay",
  "VTrials I-V Total" = "VTrials I–V Total",
  "CAVLT-slope"       = "CAVLT Slope",
  "FA-immediate"      = "FA Immediate",
  "FA-recall"         = "FA Delay",
  "FTrials I-V Total" = "FTrials I–V Total",
  "AFLT-slope"        = "AFLT Slope"
)

scale_module_map <- c(
  "WAIS-总智商"        = "WAIS index scores",
  "言语理解指数VCI"    = "WAIS index scores",
  "知觉推理指数PRI"    = "WAIS index scores",
  "工作记忆指数WMI"    = "WAIS index scores",
  "加工速度指数PSI"    = "WAIS index scores",
  
  "总记忆商"          = "WMS index scores",
  "听觉记忆指数"      = "WMS index scores",
  "视觉记忆指数"      = "WMS index scores",
  "即刻记忆指数"      = "WMS index scores",
  "延迟记忆指数"      = "WMS index scores",
  
  "VA-immediate"      = "Verbal associative learning",
  "VA-recall"         = "Verbal associative learning",
  "VTrials I-V Total" = "Verbal associative learning",
  "CAVLT-slope"       = "Verbal associative learning",
  
  "FA-immediate"      = "Figural associative learning",
  "FA-recall"         = "Figural associative learning",
  "FTrials I-V Total" = "Figural associative learning",
  "AFLT-slope"        = "Figural associative learning"
)

module_order <- c(
  "WAIS index scores",
  "WMS index scores",
  "Verbal associative learning",
  "Figural associative learning"
)

to_scale_display <- function(x) {
  out <- unname(scale_map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

scale_meta <- scale_meta %>%
  filter(Scale %in% names(df_pat)) %>%
  mutate(
    Scale_Display = to_scale_display(Scale),
    Module = unname(scale_module_map[Scale])
  ) %>%
  filter(!is.na(Module))

scale_cols <- scale_meta$Scale

if (length(scale_cols) == 0) {
  stop("❌ 未找到可分析的量表列，请检查 df_raw_clean 中的列名。")
}

cat("\n✅ Included scales:\n")
print(scale_meta %>% select(Scale, Scale_Display, Module))

# ==========================================================
# 5) 整理成长表
# ==========================================================

df_long <- df_pat %>%
  select(HS_group, all_of(scale_cols)) %>%
  mutate(across(all_of(scale_cols), to_num)) %>%
  pivot_longer(
    cols = all_of(scale_cols),
    names_to = "Scale",
    values_to = "Score"
  ) %>%
  left_join(scale_meta %>% select(Scale, Scale_Display, Module), by = "Scale") %>%
  filter(!is.na(Score), !is.na(HS_group), !is.na(Module)) %>%
  mutate(
    Module = factor(Module, levels = module_order),
    HS_group = factor(HS_group, levels = c("HS", "non-HS"))
  )

# ==========================================================
# 6) 描述统计
# ==========================================================

summary_tbl <- df_long %>%
  group_by(Module, Scale, Scale_Display, HS_group) %>%
  summarise(
    N = sum(!is.na(Score)),
    Mean = mean(Score, na.rm = TRUE),
    SD = sd(Score, na.rm = TRUE),
    SE = SD / sqrt(N),
    Median = median(Score, na.rm = TRUE),
    Q1 = quantile(Score, 0.25, na.rm = TRUE),
    Q3 = quantile(Score, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Error = ifelse(ERROR_TYPE == "SD", SD, SE),
    ymin = Mean - Error,
    ymax = Mean + Error
  )

# ==========================================================
# 7) 组间比较：HS vs non-HS
# ==========================================================

test_tbl <- df_long %>%
  group_by(Module, Scale, Scale_Display) %>%
  summarise(
    N_HS = sum(HS_group == "HS" & !is.na(Score)),
    N_nonHS = sum(HS_group == "non-HS" & !is.na(Score)),
    Mean_HS = mean(Score[HS_group == "HS"], na.rm = TRUE),
    Mean_nonHS = mean(Score[HS_group == "non-HS"], na.rm = TRUE),
    MeanDiff_HS_minus_nonHS = Mean_HS - Mean_nonHS,
    Cohen_d = cohens_d_two_group(Score, HS_group),
    p_welch = safe_welch_p(Score, HS_group),
    p_wilcox = safe_wilcox_p(Score, HS_group),
    .groups = "drop"
  ) %>%
  group_by(Module) %>%
  mutate(q_module = p.adjust(p_welch, method = "BH")) %>%
  ungroup() %>%
  mutate(
    q_all = p.adjust(p_welch, method = "BH"),
    sig_value = case_when(
      SIG_SOURCE == "p_welch"  ~ p_welch,
      SIG_SOURCE == "q_module" ~ q_module,
      SIG_SOURCE == "q_all"    ~ q_all,
      TRUE ~ q_module
    ),
    SigLabel = sig_label_from_p(sig_value, show_trend = SHOW_TREND_CIRCLE),
    p_welch_fmt = fmt_p(p_welch),
    p_wilcox_fmt = fmt_p(p_wilcox),
    q_module_fmt = fmt_p(q_module),
    q_all_fmt = fmt_p(q_all)
  )

# ==========================================================
# 8) 计算显著性标识位置
# ==========================================================

annot_tbl <- summary_tbl %>%
  group_by(Module, Scale, Scale_Display) %>%
  summarise(
    y_base = max(ymax, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(test_tbl, by = c("Module", "Scale", "Scale_Display")) %>%
  group_by(Module) %>%
  mutate(
    module_min = min(df_long$Score[df_long$Module == unique(Module)], na.rm = TRUE),
    module_max = max(df_long$Score[df_long$Module == unique(Module)], na.rm = TRUE),
    module_range = module_max - module_min,
    module_range = ifelse(
      !is.finite(module_range) | module_range == 0,
      abs(module_max),
      module_range
    ),
    module_range = ifelse(!is.finite(module_range) | module_range == 0, 1, module_range),
    y_annot = y_base + module_range * 0.08
  ) %>%
  ungroup()

# ==========================================================
# 9) 保存统计结果表
# ==========================================================

stats_export <- test_tbl %>%
  arrange(Module, Scale_Display) %>%
  select(
    Module,
    Scale,
    Scale_Display,
    N_HS,
    N_nonHS,
    Mean_HS,
    Mean_nonHS,
    MeanDiff_HS_minus_nonHS,
    Cohen_d,
    p_welch,
    p_wilcox,
    q_module,
    q_all,
    SigLabel
  )

write.csv(
  stats_export,
  file.path(output_dir, "HS_nonHS_scale_group_comparison_stats.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

write.csv(
  summary_tbl,
  file.path(output_dir, "HS_nonHS_scale_group_descriptive_stats.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

# ==========================================================
# 10) 绘图函数：每个模块单独绘图
#     无花纹版本 + 每个子图右上角图例
# ==========================================================

plot_one_module <- function(module_name) {
  
  dat_m <- summary_tbl %>%
    filter(Module == module_name) %>%
    mutate(
      Scale_Display = factor(
        Scale_Display,
        levels = scale_meta %>%
          filter(Module == module_name) %>%
          pull(Scale_Display)
      ),
      HS_group = factor(HS_group, levels = c("HS", "non-HS"))
    )
  
  ann_m <- annot_tbl %>%
    filter(Module == module_name) %>%
    mutate(
      Scale_Display = factor(
        Scale_Display,
        levels = scale_meta %>%
          filter(Module == module_name) %>%
          pull(Scale_Display)
      )
    )
  
  if (nrow(dat_m) == 0) return(NULL)
  
  y_min <- min(c(dat_m$ymin, df_long$Score[df_long$Module == module_name]), na.rm = TRUE)
  y_max <- max(c(dat_m$ymax, ann_m$y_annot), na.rm = TRUE)
  
  y_pad <- (y_max - y_min) * 0.18
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- abs(y_max) * 0.18
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 1
  
  p <- ggplot(
    dat_m,
    aes(
      x = Scale_Display,
      y = Mean,
      fill = HS_group
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.76),
      width = 0.66,
      color = "black",
      linewidth = 0.25,
      alpha = 0.95
    ) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax, group = HS_group),
      position = position_dodge(width = 0.76),
      width = 0.16,
      linewidth = 0.42,
      show.legend = FALSE
    ) +
    
    # 显著性星号
    geom_text(
      data = ann_m %>% filter(SigLabel %in% c("*", "**", "***")),
      aes(x = Scale_Display, y = y_annot, label = SigLabel),
      inherit.aes = FALSE,
      size = 4.5,
      fontface = "bold",
      color = "black"
    ) +
    
    # q < 0.10 的空心圆
    geom_point(
      data = ann_m %>% filter(SigLabel == "circle"),
      aes(x = Scale_Display, y = y_annot),
      inherit.aes = FALSE,
      shape = 1,
      size = 1.7,
      stroke = 0.8,
      color = "black"
    ) +
    
    scale_fill_manual(
      values = HS_COLORS,
      breaks = c("HS", "non-HS"),
      labels = c("HS", "non-HS"),
      drop = FALSE
    ) +
    guides(
      fill = guide_legend(
        title = NULL,
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          alpha = 0.95,
          color = "black",
          linewidth = 0.25
        )
      )
    ) +
    coord_cartesian(ylim = c(y_min - y_pad * 0.10, y_max + y_pad)) +
    labs(
      title = module_name,
      x = NULL,
      y = paste0("Score, mean ± ", ERROR_TYPE),
      fill = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      text = element_text(family = "Arial"),
      plot.title = element_text(size = 12.5, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 11, angle = 32, hjust = 1, face = "bold", color = "black"),
      axis.text.y = element_text(size = 11, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold"),
      
      # 每个子图内部右上角图例
      legend.position = c(0.98, 0.98),
      legend.justification = c(1, 1),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.background = element_rect(
        fill = scales::alpha("white", 0.78),
        color = NA
      ),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.width = unit(0.58, "cm"),
      legend.key.height = unit(0.40, "cm"),
      legend.text = element_text(size = 9.5, family = "Arial", color = "black"),
      legend.spacing.x = unit(0.15, "cm"),
      legend.margin = margin(t = 1, r = 2, b = 1, l = 2),
      
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      plot.margin = margin(t = 6, r = 8, b = 6, l = 6)
    )
  
  p
}

# ==========================================================
# 11) 生成组合图：2列 × 2行，每个子图保留右上角图例
# ==========================================================

plots <- module_order %>%
  keep(~ .x %in% unique(as.character(df_long$Module))) %>%
  map(plot_one_module)

names(plots) <- module_order[module_order %in% unique(as.character(df_long$Module))]

combined_plot <- wrap_plots(
  plots,
  ncol = 2,
  nrow = 2
)

combined_plot <- combined_plot +
  plot_annotation(
    title = "HS vs non-HS comparisons across cognitive scales",
    subtitle = paste0(
      "Bars show group mean ± ", ERROR_TYPE,
      ". Significance is based on ", SIG_SOURCE,
      " from Welch t-tests; ○<0.10, *<0.05, **<0.01, ***<0.001."
    ),
    theme = theme(
      plot.title = element_text(
        family = "Arial",
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        family = "Arial",
        size = 12.5,
        hjust = 0.5
      )
    )
  )

ggsave(
  file.path(output_dir, "HS_nonHS_all_scales_module_barplot_2x2_legendInside.png"),
  combined_plot,
  width = COMBINED_WIDTH,
  height = COMBINED_HEIGHT,
  dpi = DPI
)

ggsave(
  file.path(output_dir, "HS_nonHS_all_scales_module_barplot_2x2_legendInside.pdf"),
  combined_plot,
  width = COMBINED_WIDTH,
  height = COMBINED_HEIGHT
)

# ==========================================================
# 12) 每个模块单独保存
# ==========================================================

module_plot_dir <- file.path(output_dir, "module_plots")
if (!dir.exists(module_plot_dir)) dir.create(module_plot_dir, recursive = TRUE)

for (module_name in names(plots)) {
  ggsave(
    file.path(module_plot_dir, paste0(gsub("[^A-Za-z0-9]+", "_", module_name), "_legendInside.png")),
    plots[[module_name]],
    width = 4,
    height = 4.2,
    dpi = DPI
  )
  
  ggsave(
    file.path(module_plot_dir, paste0(gsub("[^A-Za-z0-9]+", "_", module_name), "_legendInside.pdf")),
    plots[[module_name]],
    width = 4,
    height = 4.2
  )
}

# ==========================================================
# 13) 完成提示
# ==========================================================

cat("\n✅ HS vs non-HS 量表组间比较完成。\n")
cat("输出目录：", output_dir, "\n")
cat("1. 总组合图：HS_nonHS_all_scales_module_barplot_2x2_legendInside.png / .pdf\n")
cat("2. 单模块图：", module_plot_dir, "\n")
cat("3. 组间统计表：HS_nonHS_scale_group_comparison_stats.csv\n")
cat("4. 描述统计表：HS_nonHS_scale_group_descriptive_stats.csv\n\n")
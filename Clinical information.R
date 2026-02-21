# ==============================================================================
# FIXED: Table 2 grouping + add "All patients" column
# - HS: Pathology_Group == 1
# - Other: Pathology_Group != 1  (all non-1, including 2/3/...)
# - Table 2 columns: All patients (n=71), HS, Other, Test, Effect size, P value
# - Statistical test: HS vs Other only
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(effectsize)
  library(pwr)
  library(officer)
  library(flextable)
})

# -------------------------
# 0) Checks + setup
# -------------------------
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found.")
df_pat <- df_raw_clean

OUTDIR <- "C:/Users/86150/Documents/HIPP/TABLES_Demographics"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

to_num <- function(x){
  if(is.numeric(x)) return(x)
  if(is.factor(x)) x <- as.character(x)
  if(is.character(x)){
    x <- str_replace_all(x, "[,，]", "")
    x <- str_trim(x)
    suppressWarnings(as.numeric(x))
  } else suppressWarnings(as.numeric(x))
}

# three-line style
as_three_line <- function(ft){
  ft <- autofit(ft)
  ft <- theme_vanilla(ft)
  ft <- border_remove(ft)
  ft <- hline_top(ft, border = fp_border(width = 1.2))
  ft <- hline(ft, i = 1, part = "header", border = fp_border(width = 1.0))
  ft <- hline_bottom(ft, border = fp_border(width = 1.2))
  ft
}

fmt_mean_sd <- function(x, digits=2){
  x <- x[is.finite(x)]
  if(length(x)==0) return(NA_character_)
  sprintf(paste0("%.",digits,"f ± %.",digits,"f"), mean(x), sd(x))
}
fmt_n_pct_by_level <- function(x){
  tb <- table(x, useNA="no")
  n  <- sum(tb)
  out <- paste0(names(tb), ": ", as.integer(tb), " (", sprintf("%.1f", 100*as.integer(tb)/n), "%)")
  paste(out, collapse="; ")
}
cramers_v <- function(tab){
  suppressWarnings({
    r <- nrow(tab); c <- ncol(tab)
    n <- sum(tab)
    if(n==0) return(NA_real_)
    chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
    sqrt(as.numeric(chi$statistic) / (n * (min(r-1, c-1))))
  })
}
fisher_or <- function(tab2x2){
  ft <- fisher.test(tab2x2)
  or <- unname(ft$estimate)
  ci <- unname(ft$conf.int)
  list(p = ft$p.value, or = or, lo = ci[1], hi = ci[2])
}

# continuous test: Welch default; Wilcoxon only if strong non-normality
cont_test <- function(x, g){
  df <- data.frame(x=x, g=g) %>% filter(is.finite(x) & !is.na(g))
  if(nrow(df) < 5) return(list(test="NA", p=NA_real_))
  sw1 <- if(sum(df$g==levels(df$g)[1]) >= 3) shapiro.test(df$x[df$g==levels(df$g)[1]])$p.value else 1
  sw2 <- if(sum(df$g==levels(df$g)[2]) >= 3) shapiro.test(df$x[df$g==levels(df$g)[2]])$p.value else 1
  use_wilcox <- (sw1 < 0.01) || (sw2 < 0.01)
  if(use_wilcox){
    wt <- wilcox.test(x ~ g, data=df, exact=FALSE)
    list(test="Wilcoxon rank-sum", p=wt$p.value)
  } else {
    tt <- t.test(x ~ g, data=df, var.equal=FALSE)
    list(test="Welch t-test", p=tt$p.value)
  }
}

# -------------------------
# 1) English labels mapping (EDIT KEYS if your real column names differ)
# -------------------------
label_map <- c(
  "sex" = "Sex",
  "Age" = "Age (years)",
  "Edu_year" = "Education (years)",
  "Duration_Dis" = "Disease duration (years)",
  "HADS-D" = "HADS-Depression",
  "HADS-A" = "HADS-Anxiety",
  "生活质量" = "Quality of life",
  "BNT-命名" = "BNT naming",
   "WAIS-总智商" = "WAIS Full Scale IQ",
  "言语理解指数VCI" = "Verbal Comprehension Index (VCI)",
  "知觉推理指数PRI" = "Perceptual Reasoning Index (PRI)",
  "工作记忆指数WMI" = "Working Memory Index (WMI)",
  "加工速度指数PSI" = "Processing Speed Index (PSI)",
  "总记忆商" = "General Memory Quotient",
  "听觉记忆指数" = "Auditory Memory Index",
  "视觉记忆指数" = "Visual Memory Index",
  "即刻记忆指数" = "Immediate Memory Index",
  "延迟记忆指数" = "Delayed Memory Index",
  "IID侧别（单=1；双=2）" = "IID laterality (1=unilateral; 2=bilateral)",
  "病灶侧IID分布（1=T ant, SP1/F7/T3；2=T front, FP1/F3；3=T5）" = "IID distribution on lesion side (1=T ant; 2=T front; 3=T5)"
)
var_label <- function(v){
  if(v %in% names(label_map)) return(unname(label_map[[v]]))
  v
}

# Neuropsych list requested
cog_vars <- c(
  "言语理解指数VCI","工作记忆指数WMI","总记忆商",
  "听觉记忆指数","视觉记忆指数","即刻记忆指数","延迟记忆指数",
  "VA-immediate", "VA-recall","VTrials I-V Total","CAVLT-slope",
  "FA-immediate",  "FA-recall", "FTrials I-V Total", "AFLT-slope"
)

# Patient table variables requested (plus the cognition list)
pat_vars <- c(
  "IID侧别（单=1；双=2）",
  "病灶侧IID分布（1=T ant, SP1/F7/T3；2=T front, FP1/F3；3=T5）",
  "Age","Duration_Dis","Edu_year","sex",
  "HADS-D","HADS-A","生活质量","BNT-命名","WAIS-总智商",
  "言语理解指数VCI","知觉推理指数PRI","工作记忆指数WMI","加工速度指数PSI",
  "总记忆商","听觉记忆指数","视觉记忆指数","即刻记忆指数","延迟记忆指数"
)
pat_vars <- unique(c(pat_vars, cog_vars))

# -------------------------
# 2) Build HS vs Other groups (Other = all non-1)
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
  mutate(PathGroup2 = factor(PathGroup2, levels=c("HS","Other")))

# keep only variables that exist
vars2_exist <- pat_vars[pat_vars %in% names(df_pat2)]
vars2_missing <- setdiff(pat_vars, vars2_exist)
if(length(vars2_missing)>0){
  message("⚠️ These requested variables are missing and will be skipped:\n  - ",
          paste(vars2_missing, collapse="\n  - "))
}

# =========================
# 3) Row builder (FIXED): force categorical for sex + IID variables
# =========================

# --- helper: percent string
fmt_nN_pct <- function(n, N, digits=1){
  if(is.na(n) || is.na(N) || N==0) return(NA_character_)
  sprintf("%d/%d (%.1f%%)", n, N, 100*n/N)
}

# --- safe factor
as_factor_dropna <- function(x){
  if(is.factor(x)) x <- droplevels(x) else x <- factor(x)
  x
}

# --- Identify the actual column names in your df_pat2 (either CN or EN)
pick_first_exist <- function(cands, df){
  hit <- cands[cands %in% names(df)]
  if(length(hit)==0) return(NA_character_)
  hit[1]
}

col_sex <- pick_first_exist(c("sex","Sex"), df_pat2)
col_iid_lat <- pick_first_exist(c("IID侧别（单=1；双=2）","IID laterality (1=unilateral; 2=bilateral)"), df_pat2)
col_iid_dist <- pick_first_exist(c("病灶侧IID分布（1=T ant, SP1/F7/T3；2=T front, FP1/F3；3=T5）",
                                   "IID distribution on lesion side (1=T ant; 2=T front; 3=T5)"), df_pat2)

# --- force categorical set (only add those that exist)
force_categorical <- c(col_sex, col_iid_lat, col_iid_dist)
force_categorical <- force_categorical[!is.na(force_categorical)]

# --- labeling for IID levels (optional but recommended)
recode_iid_lat <- function(x){
  x <- to_num(x)
  out <- rep(NA_character_, length(x))
  out[x==1] <- "Unilateral"
  out[x==2] <- "Bilateral"
  factor(out, levels=c("Unilateral","Bilateral"))
}
recode_iid_dist <- function(x){
  x <- to_num(x)
  out <- rep(NA_character_, length(x))
  out[x==1] <- "T ant (SP1/F7/T3)"
  out[x==2] <- "T front (FP1/F3)"
  out[x==3] <- "T5"
  factor(out, levels=c("T ant (SP1/F7/T3)","T front (FP1/F3)","T5"))
}

# --- detect sex coding and produce Male indicator
#     supports: 0/1, 1/2, "M/F", "Male/Female", "男/女"
is_male <- function(x){
  if(is.factor(x)) x <- as.character(x)
  xl <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(xl))
  out[xl %in% c("m","male","man","boy","男","1")] <- 1L
  out[xl %in% c("f","female","woman","girl","女","0","2")] <- 0L
  out
}

# --- format categorical distribution within group
fmt_cat_dist <- function(x){
  x <- x[!is.na(x)]
  if(length(x)==0) return(NA_character_)
  tb <- table(x)
  N  <- sum(tb)
  paste0(
    names(tb), ": ", as.integer(tb), " (", sprintf("%.1f", 100*as.integer(tb)/N), "%)",
    collapse="; "
  )
}

# --- the main row function (fixed)
make_row_pat2 <- function(v){
  
  # ---- special case: Sex -> show "Sex (Male)" row (binary)
  if(!is.na(col_sex) && v == col_sex){
    male <- is_male(df_pat2[[col_sex]])
    g    <- df_pat2$PathGroup2
    
    # counts (exclude NA)
    ok_all <- !is.na(male) & !is.na(g)
    N_all  <- sum(ok_all)
    n_all  <- sum(male[ok_all]==1)
    
    ok_hs <- ok_all & g=="HS"
    ok_ot <- ok_all & g=="Other"
    N_hs  <- sum(ok_hs); n_hs <- sum(male[ok_hs]==1)
    N_ot  <- sum(ok_ot); n_ot <- sum(male[ok_ot]==1)
    
    # fisher 2x2
    tab2 <- table(
      g[ok_all],
      factor(male[ok_all], levels=c(1,0), labels=c("Male","Not male"))
    )
    p <- fisher.test(tab2)$p.value
    
    # effect size: Cramér's V on 2x2
    V <- cramers_v(tab2)
    
    return(tibble::tibble(
      Variable = "Sex (Male)",
      `All patients` = fmt_nN_pct(n_all, N_all),
      HS = fmt_nN_pct(n_hs, N_hs),
      Other = fmt_nN_pct(n_ot, N_ot),
      Test = "Fisher's exact",
      `Effect size` = if(is.finite(V)) sprintf("Cramér's V=%.2f", V) else NA_character_,
      `P value` = if(is.finite(p)) sprintf("%.4f", p) else NA_character_
    ))
  }
  
  # ---- forced categorical (IID vars, or any other you add)
  if(v %in% force_categorical){
    
    x <- df_pat2[[v]]
    
    # apply recode for the two IID vars (if matched)
    if(!is.na(col_iid_lat) && v == col_iid_lat)  x <- recode_iid_lat(x)
    if(!is.na(col_iid_dist) && v == col_iid_dist) x <- recode_iid_dist(x)
    
    g <- df_pat2$PathGroup2
    ok <- !is.na(x) & !is.na(g)
    
    x <- as_factor_dropna(x[ok])
    g2 <- droplevels(g[ok])
    
    tab <- table(g2, x)
    p <- fisher.test(tab)$p.value
    V <- cramers_v(tab)
    
    return(tibble::tibble(
      Variable = var_label(v),
      `All patients` = fmt_cat_dist(x),
      HS = fmt_cat_dist(x[g2=="HS"]),
      Other = fmt_cat_dist(x[g2=="Other"]),
      Test = "Fisher's exact",
      `Effect size` = if(is.finite(V)) sprintf("Cramér's V=%.2f", V) else NA_character_,
      `P value` = if(is.finite(p)) sprintf("%.4f", p) else NA_character_
    ))
  }
  
  # ---- otherwise: numeric-like heuristic (continuous vs categorical)
  x <- df_pat2[[v]]
  x_num <- to_num(x)
  numeric_like <- sum(is.finite(x_num)) >= max(10, floor(0.6 * sum(!is.na(x))))
  
  if(numeric_like){
    df_tmp <- df_pat2 %>% transmute(g=PathGroup2, x=to_num(.data[[v]]))
    x_all <- df_tmp$x
    x_hs  <- df_tmp %>% filter(g=="HS") %>% pull(x)
    x_ot  <- df_tmp %>% filter(g=="Other") %>% pull(x)
    
    tst <- cont_test(df_tmp$x, df_tmp$g)
    eg  <- tryCatch(effectsize::hedges_g(x ~ g, data=df_tmp, pooled_sd=FALSE)$Hedges_g,
                    error=function(e) NA_real_)
    
    tibble::tibble(
      Variable = var_label(v),
      `All patients` = fmt_mean_sd(x_all, digits=2),
      HS = fmt_mean_sd(x_hs, digits=2),
      Other = fmt_mean_sd(x_ot, digits=2),
      Test = tst$test,
      `Effect size` = if(is.finite(eg)) sprintf("Hedges g=%.2f", eg) else NA_character_,
      `P value` = if(is.finite(tst$p)) sprintf("%.4f", tst$p) else NA_character_
    )
  } else {
    g <- df_pat2$PathGroup2
    x_cat <- as_factor_dropna(x)
    
    tab <- table(g, x_cat, useNA="no")
    p <- fisher.test(tab)$p.value
    V <- cramers_v(tab)
    
    tibble::tibble(
      Variable = var_label(v),
      `All patients` = fmt_n_pct_by_level(x_cat),
      HS = fmt_n_pct_by_level(x_cat[g=="HS"]),
      Other = fmt_n_pct_by_level(x_cat[g=="Other"]),
      Test = "Fisher's exact",
      `Effect size` = if(is.finite(V)) sprintf("Cramér's V=%.2f", V) else NA_character_,
      `P value` = if(is.finite(p)) sprintf("%.4f", p) else NA_character_
    )
  }
}

# rebuild Table 2 rows with the fixed function
t2_rows <- purrr::map_dfr(vars2_exist, make_row_pat2)

# -------------------------
# 4) Export to Word (three-line) + Excel-friendly CSV
# -------------------------
n_all <- nrow(df_pat2)
n_hs  <- sum(df_pat2$PathGroup2=="HS")
n_ot  <- sum(df_pat2$PathGroup2=="Other")

ft2 <- flextable(t2_rows) %>%
  set_header_labels(
    Variable="Variable",
    `All patients`=paste0("All patients (n=", n_all, ")"),
    HS=paste0("HS (n=", n_hs, ")"),
    Other=paste0("Other (n=", n_ot, ")"),
    Test="Test",
    `Effect size`="Effect size",
    `P value`="P value"
  ) %>%
  align(align="center", part="all") %>%
  align(j=1, align="left", part="all") %>%
  as_three_line()

doc2 <- read_docx() %>%
  body_add_par("Table 2. Clinical and cognitive comparison between HS and Other pathology groups (patients only)", style="heading 1") %>%
  body_add_flextable(ft2) %>%
  body_add_par(
    "Notes: Pathology_Group coded as HS=1; Other=all non-1. Statistical comparisons are HS vs Other only. Categorical variables were tested using Fisher’s exact test (effect size: Cramér’s V). Continuous variables were tested using Welch’s t-test by default; Wilcoxon rank-sum was used only when strong non-normality was detected. Effect size for continuous variables is Hedges g.",
    style="Normal"
  )

print(doc2, target = file.path(OUTDIR, "Table2_HS_vs_Other_PAT.docx"))

# Excel-friendly CSV (UTF-8 BOM)
readr::write_excel_csv(t2_rows, file.path(OUTDIR, "Table2_HS_vs_Other_PAT.csv"))

cat("\n✅ Fixed Table 2 exported:\n",
    " - ", normalizePath(file.path(OUTDIR, "Table2_HS_vs_Other_PAT.docx")), "\n",
    " - ", normalizePath(file.path(OUTDIR, "Table2_HS_vs_Other_PAT.csv")), "\n",
    sep="")
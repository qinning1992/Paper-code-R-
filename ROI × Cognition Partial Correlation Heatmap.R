suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ==========================================================
# Step 2: ROI × Cognition Partial Correlation Heatmap
# Permutation-based significance (5000 permutations)
# BH correction within each Scale column × Family: group_by(Scale, Family)
# Significance threshold: q_perm < 0.05
# ==========================================================

# ----------------------------
# 0) INPUT objects
# ----------------------------
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found in workspace.")
df_use <- df_raw_clean

# Save outputs here
output_dir <- "C:/Users/86150/Documents/HIPP/02011/ROI -Cognition Partial Correlation Heatmap/0217"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ----------------------------
# USER CONTROLS (EDIT HERE)
# ----------------------------
N_PERM <- 5000
SEED   <- 2026

# Whether to include "side" as a covariate in both ROI/Scale residualization
# (You used it before; keep TRUE for consistency with your previous runs)
INCLUDE_SIDE <- TRUE

# minimum n for correlation / permutation
MIN_N <- 15

set.seed(SEED)

# ---- Covariate pool (for existence check) ----
covar_pool <- c("sex","Duration_Dis","side","Age","Edu_year","eTIV")

# ---- Your scale metadata (as provided) ----
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
  "BNT-命名", FALSE,
  "VA-immediate", FALSE,
  "VA-recall", FALSE,
  "VTrials I-V Total", FALSE,
  "CAVLT-slope", FALSE,
  "FA-immediate", FALSE,
  "FA-recall", FALSE,
  "FTrials I-V Total", FALSE,
  "AFLT-slope", FALSE
)

# ----------------------------
# 0.1) Scale display name mapping (Chinese -> English for plotting)
# ----------------------------
scale_map <- c(
  "WAIS-总智商"       = "WAIS FSIQ",
  "言语理解指数VCI"   = "WAIS VCI",
  "知觉推理指数PRI"   = "WAIS PRI",
  "工作记忆指数WMI"   = "WAIS WMI",
  "加工速度指数PSI"   = "WAIS PSI",
  "总记忆商"         = "MQ (Total)",
  "听觉记忆指数"     = "Auditory Memory",
  "视觉记忆指数"     = "Visual Memory",
  "即刻记忆指数"     = "Immediate Memory",
  "延迟记忆指数"     = "Delayed Memory",
  "BNT-命名"         = "BNT Naming",
  "VA-immediate"     = "VA Immediate",
  "VA-recall"        = "VA Recall",
  "VTrials I-V Total"= "VTrials I-V Total",
  "CAVLT-slope"      = "CAVLT Slope",
  "FA-immediate"     = "FA Immediate",
  "FA-recall"        = "FA Recall",
  "FTrials I-V Total"= "FTrials I-V Total",
  "AFLT-slope"       = "AFLT Slope"
)

to_scale_display <- function(x){
  out <- unname(scale_map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

# ----------------------------
# 1) ROI list (as you provided)
# ----------------------------
roi_cols <- c(
  # Hippocampus subfields
  "lh_Hippocampal_tail","lh_subiculum_comb","lh_CA1_comb","lh_CA3_comb","lh_CA4_comb",
  "lh_GC_ML_DG_comb","lh_molecular_layer_HP_comb","lh_HATA","lh_fimbria","lh_presubiculum_comb",
  "lh_parasubiculum","lh_Whole_hippocampal_head","lh_Whole_hippocampal_body",
  "rh_Hippocampal_tail","rh_subiculum_comb","rh_CA1_comb","rh_CA3_comb","rh_CA4_comb",
  "rh_GC_ML_DG_comb","rh_molecular_layer_HP_comb","rh_HATA","rh_fimbria","rh_presubiculum_comb",
  "rh_parasubiculum","rh_Whole_hippocampal_head","rh_Whole_hippocampal_body",
  
  # Thalamus nuclei
  "Left_AV","Left_VA","Left_VAmc","Left_VLa","Left_VLp","Left_VPL","Left_VM","Left_CL","Left_CeM","Left_CM","Left_Pf",
  "Left_MDm","Left_MDl","Left_LD","Left_LP","Left_Pu_Total","Left_MGN","Left_LGN","Left_L_Sg",
  "Left_MV_Re","Left_Pc","Left_Pt",
  "Right_AV","Right_VA","Right_VAmc","Right_VLa","Right_VLp","Right_VPL","Right_VM","Right_CL","Right_CeM","Right_CM","Right_Pf",
  "Right_MDm","Right_MDl","Right_LD","Right_LP","Right_Pu_Total","Right_MGN","Right_LGN","Right_L_Sg",
  "Right_MV_Re","Right_Pc","Right_Pt"
)

# ----------------------------
# 2) Sanity checks
# ----------------------------
scale_cols <- scale_meta$Scale
need_cols <- unique(c(roi_cols, scale_cols, covar_pool))
missing_cols <- setdiff(need_cols, colnames(df_use))
if(length(missing_cols) > 0){
  stop("❌ Missing columns in df_use:\n- ", paste(missing_cols, collapse = "\n- "))
}

# Coerce types for covariates
df_use <- df_use %>%
  dplyr::mutate(
    sex = as.factor(sex),
    side = as.factor(side),
    Age = as.numeric(Age),
    Edu_year = as.numeric(Edu_year),
    Duration_Dis = as.numeric(Duration_Dis),
    eTIV = as.numeric(eTIV)  # 新增：eTIV转为数值型
  )

# ----------------------------
# 3) Helpers
# ----------------------------
roi_family <- function(roi){
  r <- as.character(roi)
  if(grepl("^lh_", r) || grepl("^rh_", r)) return("HippocampusSubfields")
  if(grepl("^Left_", r) || grepl("^Right_", r)) return("ThalamusNuclei")
  return("Other")
}

# residualize y ~ covars (returns residuals aligned to input vector length)
residualize_one <- function(y, d, covars){
  y <- as.numeric(y)
  if(length(covars) == 0){
    return(y - mean(y, na.rm = TRUE))
  }
  dd <- d %>%
    dplyr::select(dplyr::all_of(covars)) %>%
    dplyr::mutate(.y = y)
  
  cc <- complete.cases(dd)
  if(sum(cc) < MIN_N) return(rep(NA_real_, length(y)))
  
  X <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse = " + "))),
                    data = dd[cc, , drop=FALSE])
  fit <- lm.fit(X, dd$.y[cc])
  res <- rep(NA_real_, length(y))
  res[cc] <- dd$.y[cc] - as.numeric(X %*% fit$coefficients)
  res
}

# permutation p from residual vectors
perm_p_from_residuals <- function(x, y, n_perm = N_PERM){
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  n <- length(x)
  if(n < MIN_N) return(NA_real_)
  obs <- suppressWarnings(cor(x, y, method = "pearson"))
  if(!is.finite(obs)) return(NA_real_)
  
  perm_r <- replicate(n_perm, suppressWarnings(cor(x, sample(y), method = "pearson")))
  mean(abs(perm_r) >= abs(obs))
}

star_from_q <- function(q){
  dplyr::case_when(
    is.finite(q) & q < 0.001 ~ "***",
    is.finite(q) & q < 0.01  ~ "**",
    is.finite(q) & q < 0.05  ~ "*",
    is.finite(q) & q < 0.10  ~ "·",
    TRUE ~ ""
  )
}

# Covariates per scale (your rule; now explicitly consistent with comments)
get_covars_for_scale <- function(is_normed){
  base <- c("sex","Duration_Dis","Edu_year","eTIV")
  if(isTRUE(INCLUDE_SIDE)) base <- c(base, "side")
  if(isTRUE(is_normed)){
    # normed: do NOT include Age/Edu_year (unless you explicitly want)
    return(base)
  } else {
    return(c(base, "Age"))
  }
}

# ----------------------------
# 4) Compute ROI × Scale partial r + permutation p (long table)
# ----------------------------
cat("▶ Computing partial correlations + permutation p: ",
    length(roi_cols), " ROIs × ", length(scale_cols), " scales (",
    length(roi_cols) * length(scale_cols), " tests)\n", sep="")

res_list <- vector("list", length(roi_cols) * length(scale_cols))
k <- 1L

for(roi in roi_cols){
  fam <- roi_family(roi)
  
  for(i in seq_along(scale_cols)){
    sc <- scale_cols[i]
    is_normed <- scale_meta$Is_Normed[i]
    covars <- get_covars_for_scale(is_normed)
    
    # complete-case subset for this pair
    sub <- df_use %>%
      dplyr::select(dplyr::all_of(c(roi, sc, covars))) %>%
      tidyr::drop_na()
    n <- nrow(sub)
    
    if(n < MIN_N){
      tmp <- tibble(
        Family = fam, ROI = roi, Scale = sc,
        Scale_Display = to_scale_display(sc),
        Is_Normed = is_normed,
        CovarsUsed = paste(covars, collapse = " + "),
        n = n,
        r = NA_real_,
        perm_p = NA_real_
      )
      res_list[[k]] <- tmp
      k <- k + 1L
      next
    }
    
    # residualize (within the same complete-case subset)
    r_roi <- residualize_one(sub[[roi]], sub, covars)
    r_sc  <- residualize_one(sub[[sc]],  sub, covars)
    
    keep <- is.finite(r_roi) & is.finite(r_sc)
    n2 <- sum(keep)
    
    if(n2 < MIN_N){
      tmp <- tibble(
        Family = fam, ROI = roi, Scale = sc,
        Scale_Display = to_scale_display(sc),
        Is_Normed = is_normed,
        CovarsUsed = paste(covars, collapse = " + "),
        n = n2,
        r = NA_real_,
        perm_p = NA_real_
      )
      res_list[[k]] <- tmp
      k <- k + 1L
      next
    }
    
    # partial r
    r_obs <- suppressWarnings(cor(r_roi[keep], r_sc[keep], method = "pearson"))
    
    # permutation p (two-sided)
    perm_p <- perm_p_from_residuals(r_roi, r_sc, n_perm = N_PERM)
    
    tmp <- tibble(
      Family = fam, ROI = roi, Scale = sc,
      Scale_Display = to_scale_display(sc),
      Is_Normed = is_normed,
      CovarsUsed = paste(covars, collapse = " + "),
      n = n2,
      r = as.numeric(r_obs),
      perm_p = as.numeric(perm_p)
    )
    
    res_list[[k]] <- tmp
    k <- k + 1L
  }
}

res_long <- dplyr::bind_rows(res_list)

# ----------------------------
# 5) BH correction on permutation p within each Scale × Family
# ----------------------------
res_long <- res_long %>%
  dplyr::group_by(Scale, Family) %>%
  dplyr::mutate(q_perm = p.adjust(perm_p, method = "BH")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(star = star_from_q(q_perm))

# Save long table
write.csv(res_long,
          file.path(output_dir, "Step2_partialcorr_longtable_ROI_x_Scale_permBH.csv"),
          row.names = FALSE, fileEncoding = "GBK")

# Wide matrices (original Scale names as columns)
mat_r <- res_long %>%
  dplyr::select(ROI, Scale, r) %>%
  tidyr::pivot_wider(names_from = Scale, values_from = r)

mat_q <- res_long %>%
  dplyr::select(ROI, Scale, q_perm) %>%
  tidyr::pivot_wider(names_from = Scale, values_from = q_perm)

write.csv(mat_r, file.path(output_dir, "Step2_partialcorr_matrix_r_permBH.csv"),
          row.names = FALSE, fileEncoding = "GBK")
write.csv(mat_q, file.path(output_dir, "Step2_partialcorr_matrix_q_permBH.csv"),
          row.names = FALSE, fileEncoding = "GBK")

# Summary counts per scale & family
sum_tbl <- res_long %>%
  dplyr::mutate(sig05 = is.finite(q_perm) & q_perm < 0.05) %>%
  dplyr::group_by(Scale, Family) %>%
  dplyr::summarise(
    n_test = sum(is.finite(perm_p)),
    n_sig_q05 = sum(sig05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Scale, Family)

write.csv(sum_tbl, file.path(output_dir, "Step2_partialcorr_summary_counts_permBH.csv"),
          row.names = FALSE, fileEncoding = "GBK")

# ----------------------------
# 6) Heatmap: fill = partial r, label = stars (q_perm)
# ----------------------------
roi_order <- res_long %>%
  dplyr::distinct(ROI, Family) %>%
  dplyr::arrange(match(Family, c("HippocampusSubfields","ThalamusNuclei","Other")), ROI) %>%
  dplyr::pull(ROI)

scale_display_order <- scale_meta$Scale %>% to_scale_display()

# robust color limit
lim <- suppressWarnings(quantile(abs(res_long$r), 0.98, na.rm = TRUE))
if(!is.finite(lim) || lim < 0.2) lim <- 0.2

hm <- res_long %>%
  dplyr::mutate(
    ROI = factor(ROI, levels = rev(roi_order)),
    Scale_Display = factor(Scale_Display, levels = scale_display_order)
  )

# ----------------------------
# 6.1) Plot styling controls
# ----------------------------
BASE_SIZE          <- 12
X_TEXT_SIZE        <- 10
Y_TEXT_SIZE        <- 9
STRIP_TEXT_SIZE    <- 11
TITLE_SIZE         <- 15
SUBTITLE_SIZE      <- 10
LEGEND_TITLE_SIZE  <- 11
LEGEND_TEXT_SIZE   <- 10

X_TEXT_BOLD        <- TRUE
STRIP_TEXT_BOLD    <- TRUE
TITLE_BOLD         <- TRUE

AXIS_TEXT_COLOR    <- "black"
STRIP_TEXT_COLOR   <- "black"
TITLE_COLOR        <- "black"
STAR_COLOR         <- "black"

STAR_SIZE          <- 4
STAR_BOLD          <- TRUE

# ----------------------------
# 6.2) Draw heatmap
# ----------------------------
p_hm <- ggplot(hm, aes(x = Scale_Display, y = ROI, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = star),
    size = STAR_SIZE,
    color = STAR_COLOR,
    fontface = ifelse(STAR_BOLD, "bold", "plain"),
    hjust = 0.5,
    vjust = 0.7
  ) +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    limits = c(-lim, lim),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal(base_size = BASE_SIZE) +
  theme(
    axis.text.x = element_text(
      angle = 35, hjust = 1,
      size = X_TEXT_SIZE,
      face = ifelse(X_TEXT_BOLD, "bold", "plain"),
      color = AXIS_TEXT_COLOR
    ),
    axis.text.y = element_text(
      size = Y_TEXT_SIZE,
      color = AXIS_TEXT_COLOR
    ),
    strip.text.y = element_text(
      size = STRIP_TEXT_SIZE,
      face = ifelse(STRIP_TEXT_BOLD, "bold", "plain"),
      color = STRIP_TEXT_COLOR
    ),
    plot.title = element_text(
      size = TITLE_SIZE,
      face = ifelse(TITLE_BOLD, "bold", "plain"),
      color = TITLE_COLOR
    ),
    plot.subtitle = element_text(size = SUBTITLE_SIZE, color = AXIS_TEXT_COLOR),
    legend.title = element_text(size = LEGEND_TITLE_SIZE, face = "bold"),
    legend.text  = element_text(size = LEGEND_TEXT_SIZE)
  ) +
  labs(
    title = "ROI–Cognition partial correlations (covariate-adjusted; permutation + BH)",
    subtitle = paste0(
      "Fill = partial r (residual–residual Pearson). ",
      "Stars = q_perm (BH within Scale × Family on permutation p; Nperm=", N_PERM, "). ",
      "·<0.10, *<0.05, **<0.01, ***<0.001"
    ),
    x = "", y = "ROI", fill = "partial r"
  )

ggsave(file.path(output_dir, "Step2_partialcorr_heatmap_permBH_star.png"),
       p_hm, width = 12, height = 14, units = "in", dpi = 600)
ggsave(file.path(output_dir, "Step2_partialcorr_heatmap_permBH_star.pdf"),
       p_hm, width = 12, height = 14, units = "in")

cat("\n✅ Step 2 done (Permutation + BH within Scale×Family).\nOutputs saved in:\n", output_dir, "\n", sep="")
cat(" - Step2_partialcorr_longtable_ROI_x_Scale_permBH.csv\n")
cat(" - Step2_partialcorr_matrix_r_permBH.csv\n")
cat(" - Step2_partialcorr_matrix_q_permBH.csv\n")
cat(" - Step2_partialcorr_summary_counts_permBH.csv\n")
cat(" - Step2_partialcorr_heatmap_permBH_star.png/.pdf\n")

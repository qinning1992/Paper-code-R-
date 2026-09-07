# ==============================================================================
# Structure-cognition analyses aligned to epileptogenic laterality
# Integrated R script for manuscript data/code sharing
#
# This script implements the analyses described in Section 2.5:
#   1) fine-grained ipsilateral/contralateral hippocampal subfield and thalamic
#      nucleus-level ROI-cognition analyses;
#   2) lower-dimensional a priori summaries: hippocampal head/body/tail and
#      seven thalamic composite groups;
#   3) anatomical Left/Right ROI-cognition analyses (without epileptogenic-side
#      recoding), with seizure-onset side included as a covariate.
#
# Required input object:
#   df_raw_clean
#
# Optional standalone use:
#   Set DATA_FILE below to a de-identified CSV file containing all required
#   variables, then run this script directly.
# ==============================================================================


## Structure–cognition analyses aligned to epileptogenic laterality

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ==========================================================
# PLOT FONT
# ==========================================================
# Use Arial consistently. Register it explicitly for Windows graphics devices.
PLOT_FONT <- "Arial"
if (.Platform$OS.type == "windows") {
  suppressWarnings(
    grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
  )
}

# ==========================================================
# 0) INPUT
# ==========================================================
if(!exists("df_raw_clean")) stop("❌ df_raw_clean not found in workspace.")
df0 <- df_raw_clean

output_dir <- "C:/Users/86150/Documents/HIPP_THAL/ROI-Cognition Partial Correlation/ipsi-contra"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ==========================================================
# 1) SETTINGS
# ==========================================================
SIDE_COL <- "side"     # 病灶侧
N_PERM <- 5000
SEED   <- 2026
MIN_N  <- 15

# 主分析：ipsi/contra 对齐后一般不需要 side 协变量
SENS_INCLUDE_SIDE <- FALSE

set.seed(SEED)

# ==========================================================
# 2) SCALE META
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
  #  "加工记忆量表分",FALSE,
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
  # "加工记忆量表分"    = "WM SS",
  "BNT-命名"          = "BNT",
  "VA-immediate"      = "VA Immediate",
  "VA-recall"         = "VA Delay",
  "VTrials I-V Total" = "VTrials I–V Total",
  "CAVLT-slope"       = "CAVLT Slope",
  "FA-immediate"      = "FA Immediate",
  "FA-recall"         = "FA Delay",
  "FTrials I-V Total" = "FTrials I–V Total",
  "AFLT-slope"        = "AFLT Slope"
)

to_scale_display <- function(x){
  out <- unname(scale_map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

scale_cols <- scale_meta$Scale

# ==========================================================
# 3) ROI WHITELIST (只用你指定的海马亚区 + 丘脑核团)
# ==========================================================
# ---- Hippocampus subfields (updated Left_/Right_ column names) ----
# Current data use Left_/Right_ rather than lh_/rh_.
# CA2/3 is stored as CA2/3_comb.
hip_lh_keep <- c(
  "Left_Hippocampal_tail","Left_subiculum_comb","Left_CA1_comb","Left_CA2/3_comb","Left_CA4_comb",
  "Left_GC_ML_DG_comb","Left_molecular_layer_HP_comb","Left_HATA","Left_fimbria","Left_presubiculum_comb",
  "Left_parasubiculum"
)

hip_rh_keep <- c(
  "Right_Hippocampal_tail","Right_subiculum_comb","Right_CA1_comb","Right_CA2/3_comb","Right_CA4_comb",
  "Right_GC_ML_DG_comb","Right_molecular_layer_HP_comb","Right_HATA","Right_fimbria","Right_presubiculum_comb",
  "Right_parasubiculum"
)

# ---- Thalamus nuclei (你的清单) ----
thal_L_keep <- c(
  "Left_AV","Left_VA","Left_VAmc","Left_VLa","Left_VLp","Left_VPL","Left_VM","Left_CL","Left_CeM","Left_CM","Left_Pf",
  "Left_MDm","Left_MDl","Left_LD","Left_LP","Left_Pu_Total","Left_MGN","Left_LGN","Left_L_Sg",
  "Left_MV_Re","Left_Pc","Left_Pt"
)

thal_R_keep <- c(
  "Right_AV","Right_VA","Right_VAmc","Right_VLa","Right_VLp","Right_VPL","Right_VM","Right_CL","Right_CeM","Right_CM","Right_Pf",
  "Right_MDm","Right_MDl","Right_LD","Right_LP","Right_Pu_Total","Right_MGN","Right_LGN","Right_L_Sg",
  "Right_MV_Re","Right_Pc","Right_Pt"
)

# 检查 ROI 列是否都在数据里
need_roi_cols <- c(hip_lh_keep, hip_rh_keep, thal_L_keep, thal_R_keep)
miss_roi <- setdiff(need_roi_cols, names(df0))
if(length(miss_roi) > 0){
  stop("❌ ROI whitelist contains missing columns in df_raw_clean:\n- ",
       paste(miss_roi, collapse="\n- "))
}

# ==========================================================
# 4) SIDE NORMALIZATION -> "L"/"R"
# ==========================================================
normalize_side_lr <- function(x){
  if(is.factor(x)) x <- as.character(x)
  
  if(is.character(x)){
    xl <- tolower(trimws(x))
    out <- rep(NA_character_, length(xl))
    out[xl %in% c("1","l","left","左","左侧")] <- "L"
    out[xl %in% c("0","r","right","右","右侧")] <- "R"
    return(out)
  }
  
  if(is.numeric(x) || is.integer(x)){
    out <- rep(NA_character_, length(x))
    out[x == 1] <- "L"
    out[x == 0] <- "R"
    return(out)
  }
  
  rep(NA_character_, length(x))
}

if(!(SIDE_COL %in% names(df0))) stop("❌ side column not found: ", SIDE_COL)

df0 <- df0 %>% mutate(side_lr = normalize_side_lr(.data[[SIDE_COL]]))

if(any(is.na(df0$side_lr))){
  message("❌ side 不能完全标准化为 L/R。请检查 df_raw_clean$side 的编码。")
  print(table(df0[[SIDE_COL]], useNA="ifany"))
  stop("Fix side coding or expand normalize_side_lr().")
}

cat("✅ side normalized to L/R:\n")
print(table(df0$side_lr))

# ==========================================================
# 5) IPSI/CONTRA RECODING (白名单配对)
# ==========================================================
# 将 "Left_xxx" 和 "Right_xxx" 配对成 ipsi_hp_xxx / contra_hp_xxx
make_ipsi_contra_pairs <- function(df, side_lr_col, left_cols, right_cols,
                                   left_prefix, right_prefix, out_prefix){
  # 取 suffix（保证一一对应）
  left_suf  <- sub(paste0("^", left_prefix),  "", left_cols)
  right_suf <- sub(paste0("^", right_prefix), "", right_cols)
  
  if(!setequal(left_suf, right_suf)){
    miss_in_right <- setdiff(left_suf, right_suf)
    miss_in_left  <- setdiff(right_suf, left_suf)
    stop("❌ Left/Right suffix sets not identical.\n",
         "Missing in Right: ", paste(miss_in_right, collapse=", "), "\n",
         "Missing in Left : ", paste(miss_in_left, collapse=", "))
  }
  
  suf_all <- sort(left_suf)
  ipsi_cols <- character(0)
  contra_cols <- character(0)
  
  for(suf in suf_all){
    Lname <- paste0(left_prefix,  suf)
    Rname <- paste0(right_prefix, suf)
    
    ipsi_name   <- paste0("ipsi_", out_prefix, "_", suf)
    contra_name <- paste0("contra_", out_prefix, "_", suf)
    
    df[[ipsi_name]] <- ifelse(df[[side_lr_col]] == "L", df[[Lname]], df[[Rname]])
    df[[contra_name]] <- ifelse(df[[side_lr_col]] == "L", df[[Rname]], df[[Lname]])
    
    ipsi_cols   <- c(ipsi_cols, ipsi_name)
    contra_cols <- c(contra_cols, contra_name)
  }
  
  list(df=df, ipsi_cols=ipsi_cols, contra_cols=contra_cols)
}

# hippocampus
hp <- make_ipsi_contra_pairs(
  df0, "side_lr",
  left_cols = hip_lh_keep, right_cols = hip_rh_keep,
  left_prefix="Left_", right_prefix="Right_", out_prefix="hp"
)
df1 <- hp$df

# thalamus
th <- make_ipsi_contra_pairs(
  df1, "side_lr",
  left_cols = thal_L_keep, right_cols = thal_R_keep,
  left_prefix="Left_", right_prefix="Right_", out_prefix="thal"
)
df_use <- th$df

roi_hp_cols   <- c(hp$ipsi_cols, hp$contra_cols)       
roi_thal_cols <- c(th$ipsi_cols, th$contra_cols)      
roi_cols <- c(roi_hp_cols, roi_thal_cols)

cat("\n✅ ROI counts after recoding (expected):\n")
cat(" - Hippocampus:", length(roi_hp_cols), " (expected 22)\n")
cat(" - Thalamus   :", length(roi_thal_cols), " (expected 44)\n")
cat(" - Total      :", length(roi_cols), "\n")

# ==========================================================
# 6) COVARIATES + SANITY
# ==========================================================
covar_pool <- c("sex","Duration_Dis","Age","Edu_year","eTIV")
need_cols <- unique(c(roi_cols, scale_cols, covar_pool))
miss <- setdiff(need_cols, names(df_use))
if(length(miss) > 0){
  stop("❌ Missing required columns:\n- ", paste(miss, collapse="\n- "))
}

df_use <- df_use %>%
  mutate(
    sex = as.factor(sex),
    Age = as.numeric(Age),
    Edu_year = as.numeric(Edu_year),
    Duration_Dis = as.numeric(Duration_Dis),
    eTIV = as.numeric(eTIV)  
  )

roi_family <- function(roi){
  if(grepl("^ipsi_hp_|^contra_hp_", roi)) return("HippocampusSubfields")
  if(grepl("^ipsi_thal_|^contra_thal_", roi)) return("ThalamusNuclei")
  "Other"
}

get_covars_for_scale <- function(is_normed){
  base <- c("sex","Duration_Dis","Edu_year","eTIV")
  if(isTRUE(SENS_INCLUDE_SIDE)) base <- c(base, "side_lr")
  if(isTRUE(is_normed)) base else c(base, "Age")
}

residualize_one <- function(y, d, covars){
  y <- as.numeric(y)
  if(length(covars)==0) return(y - mean(y, na.rm=TRUE))
  dd <- d %>% select(all_of(covars)) %>% mutate(.y=y)
  cc <- complete.cases(dd)
  if(sum(cc) < MIN_N) return(rep(NA_real_, length(y)))
  X <- model.matrix(as.formula(paste0("~ ", paste(covars, collapse=" + "))),
                    data = dd[cc, , drop=FALSE])
  fit <- lm.fit(X, dd$.y[cc])
  res <- rep(NA_real_, length(y))
  res[cc] <- dd$.y[cc] - as.numeric(X %*% fit$coefficients)
  res
}

perm_p_from_residuals <- function(x, y, n_perm=N_PERM){
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  if(length(x) < MIN_N) return(NA_real_)
  obs <- suppressWarnings(cor(x, y, method="pearson"))
  if(!is.finite(obs)) return(NA_real_)
  perm_r <- replicate(n_perm, suppressWarnings(cor(x, sample(y), method="pearson")))
  mean(abs(perm_r) >= abs(obs))
}

star_from_q <- function(q){
  case_when(
    is.finite(q) & q < 0.001 ~ "***",
    is.finite(q) & q < 0.01  ~ "**",
    is.finite(q) & q < 0.05  ~ "*",
    is.finite(q) & q < 0.10  ~ "○",
    TRUE ~ ""
  )
}

# ==========================================================
# 7) STEP2: partial r + permutation p, BH within Scale×Family
# ==========================================================
cat("\n▶ Step2 IPSI/CONTRA: ",
    length(roi_cols), " ROIs × ", length(scale_cols), " scales = ",
    length(roi_cols) * length(scale_cols),
    " tests; Nperm=", N_PERM, "\n", sep="")

res_list <- vector("list", length(roi_cols) * length(scale_cols))
k <- 1L

for(roi in roi_cols){
  fam <- roi_family(roi)
  
  for(i in seq_along(scale_cols)){
    sc <- scale_cols[i]
    is_normed <- scale_meta$Is_Normed[i]
    covars <- get_covars_for_scale(is_normed)
    
    sub <- df_use %>% select(all_of(c(roi, sc, covars))) %>% drop_na()
    n <- nrow(sub)
    
    if(n < MIN_N){
      res_list[[k]] <- tibble(
        Family=fam, ROI=roi, Scale=sc, Scale_Display=to_scale_display(sc),
        Is_Normed=is_normed, CovarsUsed=paste(covars, collapse=" + "),
        n=n, r=NA_real_, perm_p=NA_real_
      )
      k <- k + 1L
      next
    }
    
    r_roi <- residualize_one(sub[[roi]], sub, covars)
    r_sc  <- residualize_one(sub[[sc]],  sub, covars)
    
    keep <- is.finite(r_roi) & is.finite(r_sc)
    n2 <- sum(keep)
    
    if(n2 < MIN_N){
      res_list[[k]] <- tibble(
        Family=fam, ROI=roi, Scale=sc, Scale_Display=to_scale_display(sc),
        Is_Normed=is_normed, CovarsUsed=paste(covars, collapse=" + "),
        n=n2, r=NA_real_, perm_p=NA_real_
      )
      k <- k + 1L
      next
    }
    
    r_obs <- suppressWarnings(cor(r_roi[keep], r_sc[keep], method="pearson"))
    perm_p <- perm_p_from_residuals(r_roi, r_sc, n_perm=N_PERM)
    
    res_list[[k]] <- tibble(
      Family=fam, ROI=roi, Scale=sc, Scale_Display=to_scale_display(sc),
      Is_Normed=is_normed, CovarsUsed=paste(covars, collapse=" + "),
      n=n2, r=as.numeric(r_obs), perm_p=as.numeric(perm_p)
    )
    k <- k + 1L
  }
}

res_long <- bind_rows(res_list)

res_long <- res_long %>%
  group_by(Scale, Family) %>%
  mutate(q_perm = p.adjust(perm_p, method="BH")) %>%
  ungroup() %>%
  mutate(star = star_from_q(q_perm))

write.csv(res_long,
          file.path(output_dir, "Step2_IPSI_longtable_permBH.csv"),
          row.names=FALSE, fileEncoding="GBK")

mat_r <- res_long %>% select(ROI, Scale, r) %>% pivot_wider(names_from=Scale, values_from=r)
mat_q <- res_long %>% select(ROI, Scale, q_perm) %>% pivot_wider(names_from=Scale, values_from=q_perm)

write.csv(mat_r, file.path(output_dir, "Step2_IPSI_matrix_r.csv"),
          row.names=FALSE, fileEncoding="GBK")
write.csv(mat_q, file.path(output_dir, "Step2_IPSI_matrix_q_permBH.csv"),
          row.names=FALSE, fileEncoding="GBK")

sum_tbl <- res_long %>%
  mutate(sig05 = is.finite(q_perm) & q_perm < 0.05) %>%
  group_by(Scale, Family) %>%
  summarise(n_test=sum(is.finite(perm_p)),
            n_sig_q05=sum(sig05, na.rm=TRUE),
            .groups="drop") %>%
  arrange(Scale, Family)

write.csv(sum_tbl, file.path(output_dir, "Step2_IPSI_summary_counts.csv"),
          row.names=FALSE, fileEncoding="GBK")

# ==========================================================
# 8) HEATMAP  face="bold",
# ==========================================================
roi_order <- res_long %>%
  distinct(ROI, Family) %>%
  arrange(match(Family, c("HippocampusSubfields","ThalamusNuclei","Other")), ROI) %>%
  pull(ROI)

scale_display_order <- scale_meta$Scale %>% to_scale_display()

lim <- suppressWarnings(quantile(abs(res_long$r), 0.98, na.rm=TRUE))
if(!is.finite(lim) || lim < 0.2) lim <- 0.2

hm <- res_long %>%
  mutate(
    ROI = factor(ROI, levels=rev(roi_order)),
    Scale_Display = factor(Scale_Display, levels=scale_display_order)
  )

p_hm <- ggplot(hm, aes(x=Scale_Display, y=ROI, fill=r)) +
  geom_tile(color="white", linewidth=0.2) +
  geom_text(aes(label=star), size=4, fontface="bold",
            hjust=0.5, vjust=0.72, color="black", family=PLOT_FONT) +
  facet_grid(Family ~ ., scales="free_y", space="free_y") +
  scale_fill_gradient2(
    low="#2166AC", mid="white", high="#B2182B",
    limits=c(-lim, lim), oob=scales::squish, na.value="grey90"
  ) +
  theme_minimal(base_size=12, base_family=PLOT_FONT) +
  theme(
    text = element_text(family=PLOT_FONT, face="bold"),
    axis.text.x = element_text(angle=35, hjust=1, size=10, face="bold", color="black"),
    axis.text.y = element_text(size=10, face="bold", color="black", family=PLOT_FONT),
    strip.text.y = element_text(size=11, face="bold"),
    plot.title = element_text(size=15, face="bold"),
    plot.subtitle = element_text(size=10, face="bold", family=PLOT_FONT)
  ) +
  labs(
    title="IPSI/CONTRA ROI–Cognition partial correlations (Permutation + BH)",
    subtitle=paste0(
      "Fill = partial r (residual–residual Pearson). ",
      "Stars = q_perm (BH within Scale × Family on permutation p; Nperm=", N_PERM, "). ",
      "○<0.10, *<0.05, **<0.01，***<0.001"),
    x="", y="ROI (ipsi/contra)", fill="partial r"
  )

ggsave(file.path(output_dir, "Step2_IPSI_heatmap_permBH.png"),
       p_hm, width=13, height=14, units="in", dpi=600)
ggsave(file.path(output_dir, "Step2_IPSI_heatmap_permBH.pdf"),
       p_hm, width=13, height=14, units="in", device=grDevices::cairo_pdf)

cat("\n✅ DONE. Output dir:\n", output_dir, "\n", sep="")



# ==========================================================
#9: 绘制并保存所有偏相关散点图 (残差 vs 残差)
# 包含拟合线、Partial r、R^2 以及 q_perm (BH调整后的p值)
# ==========================================================

cat("\n▶ 开始绘制并保存偏相关散点图...\n")


scatter_dir <- "C:/Users/86150/Documents/HIPP_THAL/ROI-Cognition Partial Correlation/ipsi-contra/scatter"
if(!dir.exists(scatter_dir)) dir.create(scatter_dir, recursive = TRUE)


plot_df_list <- res_long %>% filter(q_perm < 0.05)

total_plots <- nrow(plot_df_list)
plot_count <- 0

for(i in seq_len(total_plots)){
  row_data <- plot_df_list[i, ]
  
  if(is.na(row_data$r)) next
  
  roi <- row_data$ROI
  sc <- row_data$Scale
  sc_disp <- row_data$Scale_Display
  
  
  covars_str <- row_data$CovarsUsed
  if(covars_str == ""){
    covars <- character(0)
  } else {
    covars <- trimws(unlist(strsplit(covars_str, "\\+")))
  }
  
  sub <- df_use %>%
    dplyr::select(dplyr::all_of(c(roi, sc, covars))) %>%
    tidyr::drop_na()
  
  # 计算残差（调用你前文写好的 residualize_one 函数）
  res_roi <- residualize_one(sub[[roi]], sub, covars)
  res_sc  <- residualize_one(sub[[sc]],  sub, covars)
  
  # 组合绘图数据
  df_point <- data.frame(
    ROI_res = res_roi,
    Scale_res = res_sc
  )
  
  # 提取统计值
  r_val <- row_data$r
  r2_val <- r_val^2
  q_val <- row_data$q_perm
  star <- row_data$star
  
  # 构建副标题文本
  subtitle_txt <- sprintf(
    "Partial r = %.3f (R² = %.3f) | q_perm = %.3f %s\nCovariates: %s",
    r_val, r2_val, q_val, star, covars_str
  )
  
  # 绘制散点图
  p_scatter <- ggplot(df_point, aes(x = ROI_res, y = Scale_res)) +
    geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +        # 散点
    geom_smooth(method = "lm", formula = y ~ x, color = "#B2182B", fill = "grey70", alpha = 0.3) + # 拟合线与置信区间
    theme_minimal(base_size = 14, base_family = PLOT_FONT) +
    theme(
      text = element_text(family = PLOT_FONT, face = "bold"),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(family = PLOT_FONT, face = "bold", size = 11, color = "grey30", hjust = 0.5, margin = margin(b = 10)),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    ) +
    labs(
      title = paste(roi, "vs", sc_disp),
      subtitle = subtitle_txt,
      x = paste(roi, "(Residuals)"),
      y = paste(sc_disp, "(Residuals)")
    )
  
  # 处理文件名：替换掉 Windows 不允许的文件名字符 (例如 / \ : * ? " < > |)
  safe_roi <- gsub("[\\\\/:*?\"<>|]", "_", roi)
  safe_sc <- gsub("[\\\\/:*?\"<>|]", "_", sc)
  file_name <- paste0(safe_roi, "_vs_", safe_sc, ".png")
  
  # 保存图片
  ggsave(
    filename = file.path(scatter_dir, file_name),
    plot = p_scatter,
    width = 6, 
    height = 5, 
    dpi = 300,
    bg = "white" 
  )
  
  # 打印进度
  plot_count <- plot_count + 1
  if(plot_count %% 50 == 0){
    cat(sprintf("已完成 %d / %d 张图的绘制...\n", plot_count, total_plots))
  }
}

cat("\n✅ 所有偏相关散点图已成功绘制并保存至:\n", scatter_dir, "\n")

















suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ==========================================================
# 0) INPUT
# ==========================================================
if (!exists("df_raw_clean")) stop("❌ df_raw_clean not found in workspace.")
df0 <- df_raw_clean


output_dir <- "C:/Users/86150/Documents/HIPP_THAL/ROI-Cognition Partial Correlation/ipsi-contra/low_D"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})


# ==========================================================
# 1) SETTINGS
# ==========================================================
SIDE_COL <- "side"
N_PERM <- 5000
SEED   <- 2026
MIN_N  <- 15

# ipsi/contra 对齐后，主分析通常不再把 side 作为协变量
SENS_INCLUDE_SIDE <- FALSE


THAL_GROUP_MIN_PROP <- 1.0

set.seed(SEED)

# ==========================================================
# 2) SCALE META
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
  "加工记忆量表分", FALSE,
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
  "加工记忆量表分"    = "WM SS",
  "BNT-命名"          = "BNT Naming",
  "VA-immediate"      = "VA Immediate",
  "VA-recall"         = "VA Delay",
  "VTrials I-V Total" = "VTrials I–V Total",
  "CAVLT-slope"       = "CAVLT Slope",
  "FA-immediate"      = "FA Immediate",
  "FA-recall"         = "FA Delay",
  "FTrials I-V Total" = "FTrials I–V Total",
  "AFLT-slope"        = "AFLT Slope"
)

to_scale_display <- function(x) {
  out <- unname(scale_map[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

scale_cols <- scale_meta$Scale

# ==========================================================
# 3) ROI WHITELIST
#    海马：只保留 head/body/tail
#    丘脑：先保留核团，后续再求和为 composite groups
# ==========================================================

# ---------- 3.1 Hippocampus head/body/tail ----------
# 【核心设置】海马只保留三段，不使用高分辨率亚区
hip_segment_suffix <- c(
  "Whole_hippocampal_head",
  "Whole_hippocampal_body",
  "Hippocampal_tail"
)

hip_lh_keep <- paste0("Left_", hip_segment_suffix)
hip_rh_keep <- paste0("Right_", hip_segment_suffix)

# ---------- 3.2 Hippocampus display names ----------
# 【可修改】图和表中的海马分段名称
hip_segment_display <- c(
  "Whole_hippocampal_head" = "Hippocampal head",
  "Whole_hippocampal_body" = "Hippocampal body",
  "Hippocampal_tail"       = "Hippocampal tail"
)

hip_segment_full_name <- c(
  "Whole_hippocampal_head" = "Whole hippocampal head",
  "Whole_hippocampal_body" = "Whole hippocampal body",
  "Hippocampal_tail"       = "Hippocampal tail"
)

# ---------- 3.3 Thalamus nuclei whitelist ----------
thal_L_keep <- c(
  "Left_AV", "Left_VA", "Left_VAmc", "Left_VLa", "Left_VLp", "Left_VPL", "Left_VM",
  "Left_CL", "Left_CeM", "Left_CM", "Left_Pf",
  "Left_MDm", "Left_MDl", "Left_LD", "Left_LP", "Left_Pu_Total",
  "Left_MGN", "Left_LGN", "Left_L_Sg",
  "Left_MV_Re", "Left_Pc", "Left_Pt"
)

thal_R_keep <- c(
  "Right_AV", "Right_VA", "Right_VAmc", "Right_VLa", "Right_VLp", "Right_VPL", "Right_VM",
  "Right_CL", "Right_CeM", "Right_CM", "Right_Pf",
  "Right_MDm", "Right_MDl", "Right_LD", "Right_LP", "Right_Pu_Total",
  "Right_MGN", "Right_LGN", "Right_L_Sg",
  "Right_MV_Re", "Right_Pc", "Right_Pt"
)

need_roi_cols <- c(hip_lh_keep, hip_rh_keep, thal_L_keep, thal_R_keep)
miss_roi <- setdiff(need_roi_cols, names(df0))

if (length(miss_roi) > 0) {
  stop(
    "❌ ROI whitelist contains missing columns in df_raw_clean:\n- ",
    paste(miss_roi, collapse = "\n- ")
  )
}

# ==========================================================
# 3.5) Thalamus composite grouping
# ==========================================================
# 【可修改】丘脑 composite group 定义
thal_group_def <- list(
  AnteriorDorsalGroup = c("AV", "LD"),
  MediodorsalNuclei   = c("MDm", "MDl"),
  IntralaminarNuclei  = c("CL", "CM", "CeM", "Pf", "Pc"),
  PulvinarLPComplex   = c("Pu_Total", "L_Sg", "LP"),
  VentralNuclei       = c("VA", "VAmc", "VLa", "VLp", "VPL", "VM"),
  MidlineNuclei       = c("MV_Re", "Pt"),
  GeniculateNuclei    = c("LGN", "MGN")
)

# 【可修改】丘脑组显示名称：解剖命名为主，括号中列出组成核团
thal_group_display <- c(
  "AnteriorDorsalGroup" = "Anterior-dorsal group (AV+LD)",
  "MediodorsalNuclei"   = "Mediodorsal nuclei (MDm+MDl)",
  "IntralaminarNuclei"  = "Intralaminar nuclei (CL+CM+CeM+Pf+Pc)",
  "PulvinarLPComplex"   = "Pulvinar-LP complex (Pu_Total+L_Sg+LP)",
  "VentralNuclei"       = "Ventral nuclei (VA+VAmc+VLa+VLp+VPL+VM)",
  "MidlineNuclei"       = "Midline nuclei (MV_Re+Pt)",
  "GeniculateNuclei"    = "Geniculate nuclei (LGN+MGN)"
)

# ==========================================================
# 4) SIDE NORMALIZATION -> "L"/"R"
# ==========================================================
normalize_side_lr <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  
  if (is.character(x)) {
    xl <- tolower(trimws(x))
    out <- rep(NA_character_, length(xl))
    out[xl %in% c("1", "l", "left", "左", "左侧")] <- "L"
    out[xl %in% c("0", "r", "right", "右", "右侧")] <- "R"
    return(out)
  }
  
  if (is.numeric(x) || is.integer(x)) {
    out <- rep(NA_character_, length(x))
    out[x == 1] <- "L"
    out[x == 0] <- "R"
    return(out)
  }
  
  rep(NA_character_, length(x))
}

if (!(SIDE_COL %in% names(df0))) stop("❌ side column not found: ", SIDE_COL)

df0 <- df0 %>%
  mutate(side_lr = normalize_side_lr(.data[[SIDE_COL]]))

if (any(is.na(df0$side_lr))) {
  message("❌ side 不能完全标准化为 L/R。请检查 df_raw_clean$side 的编码。")
  print(table(df0[[SIDE_COL]], useNA = "ifany"))
  stop("Fix side coding or expand normalize_side_lr().")
}

cat("✅ side normalized to L/R:\n")
print(table(df0$side_lr))

# ==========================================================
# 5) Helper functions
# ==========================================================

row_sum_min_n <- function(df, vars, min_prop = 1.0) {
  vars_exist <- intersect(vars, names(df))
  if (length(vars_exist) == 0) return(rep(NA_real_, nrow(df)))
  
  x_df <- df[, vars_exist, drop = FALSE] %>%
    mutate(across(everything(), as.numeric))
  
  x <- as.matrix(x_df)
  nonmiss_n <- rowSums(!is.na(x))
  need_n <- ceiling(length(vars_exist) * min_prop)
  
  out <- rowSums(x, na.rm = TRUE)
  out[nonmiss_n < need_n] <- NA_real_
  out
}

make_ipsi_contra_pairs <- function(df, side_lr_col, left_cols, right_cols,
                                   left_prefix, right_prefix, out_prefix) {
  left_suf  <- sub(paste0("^", left_prefix),  "", left_cols)
  right_suf <- sub(paste0("^", right_prefix), "", right_cols)
  
  if (!setequal(left_suf, right_suf)) {
    miss_in_right <- setdiff(left_suf, right_suf)
    miss_in_left  <- setdiff(right_suf, left_suf)
    stop(
      "❌ Left/Right suffix sets not identical.\n",
      "Missing in Right: ", paste(miss_in_right, collapse = ", "), "\n",
      "Missing in Left : ", paste(miss_in_left, collapse = ", ")
    )
  }
  
  # 保持输入顺序，不再 sort，避免 head/body/tail 顺序被打乱
  suf_all <- left_suf
  
  ipsi_cols <- character(0)
  contra_cols <- character(0)
  
  for (suf in suf_all) {
    Lname <- paste0(left_prefix,  suf)
    Rname <- paste0(right_prefix, suf)
    
    ipsi_name   <- paste0("ipsi_", out_prefix, "_", suf)
    contra_name <- paste0("contra_", out_prefix, "_", suf)
    
    df[[ipsi_name]]   <- ifelse(df[[side_lr_col]] == "L", df[[Lname]], df[[Rname]])
    df[[contra_name]] <- ifelse(df[[side_lr_col]] == "L", df[[Rname]], df[[Lname]])
    
    ipsi_cols   <- c(ipsi_cols, ipsi_name)
    contra_cols <- c(contra_cols, contra_name)
  }
  
  list(df = df, ipsi_cols = ipsi_cols, contra_cols = contra_cols)
}

pretty_hp_name <- function(x) {
  out <- unname(hip_segment_display[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

pretty_thal_name <- function(x) {
  out <- unname(thal_group_display[x])
  out[is.na(out)] <- x[is.na(out)]
  out
}

pretty_roi_display <- function(roi) {
  sapply(roi, function(x) {
    if (grepl("^ipsi_hp_", x)) {
      suf <- sub("^ipsi_hp_", "", x)
      return(paste0("Ipsi hippocampus – ", pretty_hp_name(suf)))
    }
    
    if (grepl("^contra_hp_", x)) {
      suf <- sub("^contra_hp_", "", x)
      return(paste0("Contra hippocampus – ", pretty_hp_name(suf)))
    }
    
    if (grepl("^ipsi_thalgrp_", x)) {
      suf <- sub("^ipsi_thalgrp_", "", x)
      return(paste0("Ipsi thalamus – ", pretty_thal_name(suf)))
    }
    
    if (grepl("^contra_thalgrp_", x)) {
      suf <- sub("^contra_thalgrp_", "", x)
      return(paste0("Contra thalamus – ", pretty_thal_name(suf)))
    }
    
    x
  }, USE.NAMES = FALSE)
}

# ==========================================================
# 6) IPSI/CONTRA RECODING
#    6.1 海马：head/body/tail
#    6.2 丘脑：按 composite group 求和，再重编码 ipsi/contra
# ==========================================================

# ---------- 6.1 hippocampus head/body/tail ----------
hp <- make_ipsi_contra_pairs(
  df0,
  side_lr_col = "side_lr",
  left_cols = hip_lh_keep,
  right_cols = hip_rh_keep,
  left_prefix = "Left_",
  right_prefix = "Right_",
  out_prefix = "hp"
)

df1 <- hp$df

# ---------- 6.2 thalamus grouped Left / Right composites by SUM ----------
for (gname in names(thal_group_def)) {
  left_members  <- paste0("Left_",  thal_group_def[[gname]])
  right_members <- paste0("Right_", thal_group_def[[gname]])
  
  miss_l <- setdiff(left_members, names(df1))
  miss_r <- setdiff(right_members, names(df1))
  
  if (length(miss_l) > 0 || length(miss_r) > 0) {
    stop(
      "❌ Missing thalamic members for group ", gname, ":\n",
      "Left missing: ", paste(miss_l, collapse = ", "), "\n",
      "Right missing: ", paste(miss_r, collapse = ", ")
    )
  }
  
  df1[[paste0("LeftGrp_", gname)]]  <- row_sum_min_n(df1, left_members,  min_prop = THAL_GROUP_MIN_PROP)
  df1[[paste0("RightGrp_", gname)]] <- row_sum_min_n(df1, right_members, min_prop = THAL_GROUP_MIN_PROP)
}

thal_left_group_cols  <- paste0("LeftGrp_", names(thal_group_def))
thal_right_group_cols <- paste0("RightGrp_", names(thal_group_def))

th <- make_ipsi_contra_pairs(
  df1,
  side_lr_col = "side_lr",
  left_cols = thal_left_group_cols,
  right_cols = thal_right_group_cols,
  left_prefix = "LeftGrp_",
  right_prefix = "RightGrp_",
  out_prefix = "thalgrp"
)

df_use <- th$df

roi_hp_cols   <- c(hp$ipsi_cols, hp$contra_cols)
roi_thal_cols <- c(th$ipsi_cols, th$contra_cols)
roi_cols <- c(roi_hp_cols, roi_thal_cols)

cat("\n✅ ROI counts after recoding:\n")
cat(" - Hippocampus head/body/tail:", length(roi_hp_cols), " (expected 6)\n")
cat(" - Thalamus grouped SUM:", length(roi_thal_cols), " (expected 14)\n")
cat(" - Total:", length(roi_cols), " (expected 20)\n")

# ==========================================================
# 6.5) Save ROI definition tables
# ==========================================================

hip_segment_tbl <- tibble(
  Segment = names(hip_segment_display),
  Display = unname(hip_segment_display),
  Full_Name = unname(hip_segment_full_name[names(hip_segment_display)]),
  Left_Column = paste0("Left_", names(hip_segment_display)),
  Right_Column = paste0("Right_", names(hip_segment_display))
)

write.csv(
  hip_segment_tbl,
  file.path(output_dir, "Hippocampal_HeadBodyTail_Definitions.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

thal_group_tbl <- tibble(
  Group = names(thal_group_def),
  Display = unname(thal_group_display[names(thal_group_def)]),
  Members = sapply(thal_group_def, function(x) paste(x, collapse = " + ")),
  Left_Members = sapply(thal_group_def, function(x) paste0("Left_", x, collapse = "; ")),
  Right_Members = sapply(thal_group_def, function(x) paste0("Right_", x, collapse = "; "))
)

write.csv(
  thal_group_tbl,
  file.path(output_dir, "Thalamic_Composite_Group_Definitions_SUM.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

# ==========================================================
# 7) COVARIATES + SANITY
# ==========================================================
covar_pool <- c("sex", "Duration_Dis", "Age", "Edu_year", "eTIV")
need_cols <- unique(c(roi_cols, scale_cols, covar_pool))
miss <- setdiff(need_cols, names(df_use))

if (length(miss) > 0) {
  stop("❌ Missing required columns:\n- ", paste(miss, collapse = "\n- "))
}

df_use <- df_use %>%
  mutate(
    sex = as.factor(sex),
    Age = as.numeric(Age),
    Edu_year = as.numeric(Edu_year),
    Duration_Dis = as.numeric(Duration_Dis),
    eTIV = as.numeric(eTIV)
  )

roi_family <- function(roi) {
  if (grepl("^ipsi_hp_|^contra_hp_", roi)) return("Hippocampal segments")
  if (grepl("^ipsi_thalgrp_|^contra_thalgrp_", roi)) return("Thalamic composite groups")
  "Other"
}

get_covars_for_scale <- function(is_normed) {
  base <- c("sex", "Duration_Dis", "Edu_year", "eTIV")
  if (isTRUE(SENS_INCLUDE_SIDE)) base <- c(base, "side_lr")
  if (isTRUE(is_normed)) base else c(base, "Age")
}

residualize_one <- function(y, d, covars) {
  y <- as.numeric(y)
  
  if (length(covars) == 0) {
    return(y - mean(y, na.rm = TRUE))
  }
  
  dd <- d %>%
    select(all_of(covars)) %>%
    mutate(.y = y)
  
  cc <- complete.cases(dd)
  
  if (sum(cc) < MIN_N) {
    return(rep(NA_real_, length(y)))
  }
  
  X <- model.matrix(
    as.formula(paste0("~ ", paste(covars, collapse = " + "))),
    data = dd[cc, , drop = FALSE]
  )
  
  fit <- lm.fit(X, dd$.y[cc])
  
  res <- rep(NA_real_, length(y))
  res[cc] <- dd$.y[cc] - as.numeric(X %*% fit$coefficients)
  res
}

perm_p_from_residuals <- function(x, y, n_perm = N_PERM) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  
  if (length(x) < MIN_N) return(NA_real_)
  
  obs <- suppressWarnings(cor(x, y, method = "pearson"))
  if (!is.finite(obs)) return(NA_real_)
  
  perm_r <- replicate(
    n_perm,
    suppressWarnings(cor(x, sample(y), method = "pearson"))
  )
  
  mean(abs(perm_r) >= abs(obs))
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

# ==========================================================
# 8) STEP2: partial r + permutation p, BH within Scale × Family
# ==========================================================
cat(
  "\n▶ Step2 IPSI/CONTRA: ",
  length(roi_cols), " ROIs × ", length(scale_cols), " scales = ",
  length(roi_cols) * length(scale_cols),
  " tests; Nperm=", N_PERM, "\n",
  sep = ""
)

res_list <- vector("list", length(roi_cols) * length(scale_cols))
k <- 1L

for (roi in roi_cols) {
  fam <- roi_family(roi)
  
  for (i in seq_along(scale_cols)) {
    sc <- scale_cols[i]
    is_normed <- scale_meta$Is_Normed[i]
    covars <- get_covars_for_scale(is_normed)
    
    sub <- df_use %>%
      select(all_of(c(roi, sc, covars))) %>%
      drop_na()
    
    n <- nrow(sub)
    
    if (n < MIN_N) {
      res_list[[k]] <- tibble(
        Family = fam,
        ROI = roi,
        ROI_Display = pretty_roi_display(roi),
        Scale = sc,
        Scale_Display = to_scale_display(sc),
        Is_Normed = is_normed,
        CovarsUsed = paste(covars, collapse = " + "),
        n = n,
        r = NA_real_,
        perm_p = NA_real_
      )
      k <- k + 1L
      next
    }
    
    r_roi <- residualize_one(sub[[roi]], sub, covars)
    r_sc  <- residualize_one(sub[[sc]],  sub, covars)
    
    keep <- is.finite(r_roi) & is.finite(r_sc)
    n2 <- sum(keep)
    
    if (n2 < MIN_N) {
      res_list[[k]] <- tibble(
        Family = fam,
        ROI = roi,
        ROI_Display = pretty_roi_display(roi),
        Scale = sc,
        Scale_Display = to_scale_display(sc),
        Is_Normed = is_normed,
        CovarsUsed = paste(covars, collapse = " + "),
        n = n2,
        r = NA_real_,
        perm_p = NA_real_
      )
      k <- k + 1L
      next
    }
    
    r_obs <- suppressWarnings(cor(r_roi[keep], r_sc[keep], method = "pearson"))
    perm_p <- perm_p_from_residuals(r_roi, r_sc, n_perm = N_PERM)
    
    res_list[[k]] <- tibble(
      Family = fam,
      ROI = roi,
      ROI_Display = pretty_roi_display(roi),
      Scale = sc,
      Scale_Display = to_scale_display(sc),
      Is_Normed = is_normed,
      CovarsUsed = paste(covars, collapse = " + "),
      n = n2,
      r = as.numeric(r_obs),
      perm_p = as.numeric(perm_p)
    )
    k <- k + 1L
  }
}

res_long <- bind_rows(res_list)

res_long <- res_long %>%
  group_by(Scale, Family) %>%
  mutate(q_perm = p.adjust(perm_p, method = "BH")) %>%
  ungroup() %>%
  mutate(star = star_from_q(q_perm))

write.csv(
  res_long,
  file.path(output_dir, "Step2_IPSI_ThalSUM_HipHeadBodyTail_longtable_permBH.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

mat_r <- res_long %>%
  select(ROI, ROI_Display, Scale, r) %>%
  pivot_wider(names_from = Scale, values_from = r)

mat_q <- res_long %>%
  select(ROI, ROI_Display, Scale, q_perm) %>%
  pivot_wider(names_from = Scale, values_from = q_perm)

write.csv(
  mat_r,
  file.path(output_dir, "Step2_IPSI_ThalSUM_HipHeadBodyTail_matrix_r.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

write.csv(
  mat_q,
  file.path(output_dir, "Step2_IPSI_ThalSUM_HipHeadBodyTail_matrix_q_permBH.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

sum_tbl <- res_long %>%
  mutate(sig05 = is.finite(q_perm) & q_perm < 0.05) %>%
  group_by(Scale, Family) %>%
  summarise(
    n_test = sum(is.finite(perm_p)),
    n_sig_q05 = sum(sig05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Scale, Family)

write.csv(
  sum_tbl,
  file.path(output_dir, "Step2_IPSI_ThalSUM_HipHeadBodyTail_summary_counts.csv"),
  row.names = FALSE,
  fileEncoding = "GBK"
)

# ==========================================================
# 9) HEATMAP
# ==========================================================

hp_order <- c(
  paste0("ipsi_hp_", hip_segment_suffix),
  paste0("contra_hp_", hip_segment_suffix)
)

thal_order <- c(
  paste0("ipsi_thalgrp_", names(thal_group_def)),
  paste0("contra_thalgrp_", names(thal_group_def))
)

roi_order <- c(
  hp_order[hp_order %in% res_long$ROI],
  thal_order[thal_order %in% res_long$ROI]
)

scale_display_order <- scale_meta$Scale %>%
  to_scale_display()

lim <- suppressWarnings(quantile(abs(res_long$r), 0.98, na.rm = TRUE))
if (!is.finite(lim) || lim < 0.2) lim <- 0.2

roi_order_display <- pretty_roi_display(roi_order)

hm <- res_long %>%
  mutate(
    ROI_Display = pretty_roi_display(ROI),
    ROI_Display = factor(ROI_Display, levels = rev(roi_order_display)),
    Scale_Display = factor(Scale_Display, levels = scale_display_order)
  )

p_hm <- ggplot(hm, aes(x = Scale_Display, y = ROI_Display, fill = r)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = star),
    size = 4,
    fontface = "bold",
    hjust = 0.5,
    vjust = 0.72,
    color = "black",
    family = PLOT_FONT
  ) +
  facet_grid(Family ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    limits = c(-lim, lim),
    oob = scales::squish,
    na.value = "grey90"
  ) +
  theme_minimal(base_size = 12, base_family = PLOT_FONT) +
  theme(
    text = element_text(family = PLOT_FONT, face = "bold"),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 10,
      face = "bold",
      color = "black"
    ),
    axis.text.y = element_text(
      size = 9.5,
      face = "bold",
      color = "black",
      family = PLOT_FONT
    ),
    strip.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 15, face = "bold"),
    plot.subtitle = element_text(size = 10, face = "bold", family = PLOT_FONT)
  ) +
  labs(
    title = "IPSI/CONTRA ROI–Cognition partial correlations",
    subtitle = paste0(
      "Fill = partial r. ",
      "Stars = q_perm after BH correction within Scale × Family; Nperm = ", N_PERM, ". ",
      "o<0.10, *<0.05, **<0.01, ***<0.001"
    ),
    x = "",
    y = "ROI (ipsi/contra)",
    fill = "partial r"
  )

ggsave(
  file.path(output_dir, "Step2_IPSI_ThalSUM_HipHeadBodyTail_heatmap_permBH.png"),
  p_hm,
  width = 13,
  height = 8.5,
  units = "in",
  dpi = 600
)

ggsave(
  file.path(output_dir, "Step2_IPSI_ThalSUM_HipHeadBodyTail_heatmap_permBH.pdf"),
  p_hm,
  width = 13,
  height = 8.5,
  units = "in",
  device = grDevices::cairo_pdf
)

# ==========================================================
# 10) SCATTER PLOTS: ALL results
# 绘制并保存所有偏相关散点图（残差 vs 残差）
# 包含拟合线、Partial r、R^2、perm_p、q_perm
# ==========================================================
cat("\n▶ 开始绘制并保存所有偏相关散点图...\n")

scatter_dir <- file.path(output_dir, "scatter_all")
if (!dir.exists(scatter_dir)) dir.create(scatter_dir, recursive = TRUE)

# 画全部，不再筛选 q_perm < 0.05
plot_df_list <- res_long

total_plots <- nrow(plot_df_list)
plot_count <- 0

if (total_plots == 0) {
  cat("\n⚠️ res_long 为空，因此未生成散点图。\n")
} else {
  for (i in seq_len(total_plots)) {
    row_data <- plot_df_list[i, ]
    
    # 如果计算失败（样本量不足），跳过
    if (is.na(row_data$r)) next
    
    roi <- row_data$ROI
    roi_disp <- row_data$ROI_Display
    sc <- row_data$Scale
    sc_disp <- row_data$Scale_Display
    
    # 获取当前组合使用的协变量
    covars_str <- row_data$CovarsUsed
    if (covars_str == "") {
      covars <- character(0)
    } else {
      covars <- trimws(unlist(strsplit(covars_str, "\\+")))
    }
    
    # 提取完整数据集并剔除 NA
    sub <- df_use %>%
      dplyr::select(dplyr::all_of(c(roi, sc, covars))) %>%
      tidyr::drop_na()
    
    # 计算残差
    res_roi <- residualize_one(sub[[roi]], sub, covars)
    res_sc  <- residualize_one(sub[[sc]],  sub, covars)
    
    df_point <- data.frame(
      ROI_res   = res_roi,
      Scale_res = res_sc
    )
    
    # 提取统计值
    r_val <- row_data$r
    r2_val <- r_val^2
    p_val <- row_data$perm_p
    q_val <- row_data$q_perm
    star  <- row_data$star
    
    # 格式化 p / q
    p_txt <- ifelse(is.na(p_val), "NA",
                    ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)))
    q_txt <- ifelse(is.na(q_val), "NA",
                    ifelse(q_val < 0.001, "<0.001", sprintf("%.3f", q_val)))
    
    # 左上角图注：同时显示 p 和 q
    annot_txt <- sprintf(
      "Partial r = %.3f (R² = %.3f)\np_perm = %s\nq_perm = %s %s",
      r_val, r2_val, p_txt, q_txt, star
    )
    
    p_scatter <- ggplot(df_point, aes(x = ROI_res, y = Scale_res)) +
      geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        color = "#B2182B",
        fill = "grey70",
        alpha = 0.3
      ) +
      annotate(
        "text",
        x = -Inf,
        y = Inf,
        label = annot_txt,
        hjust = -0.05,
        vjust = 1.2,
        size = 4.8,
        color = "grey20",
        lineheight = 1.1,
        family = PLOT_FONT,
        fontface = "bold"
      ) +
      theme_minimal(base_size = 14, base_family = PLOT_FONT) +
      theme(
        text = element_text(family = PLOT_FONT, face = "bold"),
        plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      ) +
      labs(
        title = paste(roi_disp, "vs", sc_disp),
        x = paste(roi_disp, "(residuals)"),
        y = paste(sc_disp, "(residuals)")
      )
    
    # 安全文件名
    safe_roi <- gsub("[\\\\/:*?\"<>|]", "_", roi)
    safe_sc  <- gsub("[\\\\/:*?\"<>|]", "_", sc)
    file_name <- paste0(safe_roi, "_vs_", safe_sc, ".png")
    
    ggsave(
      filename = file.path(scatter_dir, file_name),
      plot = p_scatter,
      width = 6,
      height = 5,
      dpi = 600,
      bg = "white"
    )
    
    plot_count <- plot_count + 1
    if (plot_count %% 25 == 0) {
      cat(sprintf("已完成 %d / %d 张图的绘制...\n", plot_count, total_plots))
    }
  }
}

cat("\n✅ 所有偏相关散点图已成功绘制并保存至:\n", scatter_dir, "\n")
cat("\n✅ DONE. Output dir:\n", output_dir, "\n", sep = "")





























# ==============================================================================
# ADDITIONAL ANALYSIS: ANATOMICAL LEFT/RIGHT ROI–COGNITION ASSOCIATIONS
# ==============================================================================
# Purpose
#   Preserve the original anatomical hemisphere (Left vs Right) rather than
#   recoding ROIs relative to the epileptogenic side.
#
# Statistical principle
#   - Same partial-correlation / permutation framework as the analyses above.
#   - In the pooled MTLE sample, seizure-onset side (side_lr) is INCLUDED as a
#     covariate because anatomical Left/Right ROIs are not aligned to the
#     epileptogenic side.
#   - BH-FDR is performed within Scale × Family. Left and Right ROIs of the same
#     anatomical family therefore belong to the same multiple-comparison family.
#
# This appendix contains two parallel analyses:
#   A) Fine-grained anatomical Left/Right hippocampal subfields + thalamic nuclei
#   B) Lower-dimensional anatomical Left/Right hippocampal head/body/tail +
#      thalamic composite groups
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

if (!exists("df_raw_clean")) stop("❌ df_raw_clean not found in workspace.")

ANAT_BASE_DIR <- "C:/Users/86150/Documents/HIPP_THAL/Anatomical_LeftRight"
ANAT_FINE_DIR <- file.path(ANAT_BASE_DIR, "FineGrained")
ANAT_LOW_DIR  <- file.path(ANAT_BASE_DIR, "LowDimensional")

for (dd in c(ANAT_BASE_DIR, ANAT_FINE_DIR, ANAT_LOW_DIR)) {
  if (!dir.exists(dd)) dir.create(dd, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Shared settings / helper checks
# ------------------------------------------------------------------------------
SIDE_COL_ANAT <- "side"
N_PERM_ANAT <- if (exists("N_PERM")) N_PERM else 5000
SEED_ANAT   <- if (exists("SEED"))   SEED   else 2026
MIN_N_ANAT  <- if (exists("MIN_N"))  MIN_N  else 15
set.seed(SEED_ANAT)

if (!exists("normalize_side_lr")) {
  normalize_side_lr <- function(x) {
    if (is.factor(x)) x <- as.character(x)
    if (is.character(x)) {
      xl <- tolower(trimws(x))
      out <- rep(NA_character_, length(xl))
      out[xl %in% c("1", "l", "left", "左", "左侧")] <- "L"
      out[xl %in% c("0", "r", "right", "右", "右侧")] <- "R"
      return(out)
    }
    if (is.numeric(x) || is.integer(x)) {
      out <- rep(NA_character_, length(x))
      out[x == 1] <- "L"
      out[x == 0] <- "R"
      return(out)
    }
    rep(NA_character_, length(x))
  }
}

# Use local helper names to avoid ambiguity with earlier sections.
residualize_one_anat <- function(y, d, covars) {
  y <- as.numeric(y)
  if (length(covars) == 0) return(y - mean(y, na.rm = TRUE))
  dd <- d %>% select(all_of(covars)) %>% mutate(.y = y)
  cc <- complete.cases(dd)
  if (sum(cc) < MIN_N_ANAT) return(rep(NA_real_, length(y)))
  X <- model.matrix(
    as.formula(paste0("~ ", paste(covars, collapse = " + "))),
    data = dd[cc, , drop = FALSE]
  )
  fit <- lm.fit(X, dd$.y[cc])
  res <- rep(NA_real_, length(y))
  res[cc] <- dd$.y[cc] - as.numeric(X %*% fit$coefficients)
  res
}

perm_p_from_residuals_anat <- function(x, y, n_perm = N_PERM_ANAT) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < MIN_N_ANAT) return(NA_real_)
  obs <- suppressWarnings(cor(x, y, method = "pearson"))
  if (!is.finite(obs)) return(NA_real_)
  perm_r <- replicate(
    n_perm,
    suppressWarnings(cor(x, sample(y), method = "pearson"))
  )
  (1 + sum(abs(perm_r) >= abs(obs), na.rm = TRUE)) / (1 + n_perm)
}

star_from_q_anat <- function(q) {
  case_when(
    is.finite(q) & q < 0.001 ~ "***",
    is.finite(q) & q < 0.01  ~ "**",
    is.finite(q) & q < 0.05  ~ "*",
    is.finite(q) & q < 0.10  ~ "○",
    TRUE ~ ""
  )
}

if (!exists("to_scale_display")) {
  stop("❌ to_scale_display() not found. Run the original script sections above first.")
}
if (!exists("scale_meta")) {
  stop("❌ scale_meta not found. Run the original script sections above first.")
}

# Anatomical Left/Right pooled analysis: seizure-onset side is always controlled.
get_covars_anatomical <- function(is_normed) {
  base <- c("sex", "Duration_Dis", "Edu_year", "eTIV", "side_lr")
  if (isTRUE(is_normed)) base else c(base, "Age")
}

prepare_anatomical_df <- function(df) {
  if (!(SIDE_COL_ANAT %in% names(df))) {
    stop("❌ side column not found: ", SIDE_COL_ANAT)
  }
  
  out <- df %>%
    mutate(
      side_lr = normalize_side_lr(.data[[SIDE_COL_ANAT]]),
      side_lr = factor(side_lr, levels = c("L", "R")),
      sex = as.factor(sex),
      Age = as.numeric(Age),
      Edu_year = as.numeric(Edu_year),
      Duration_Dis = as.numeric(Duration_Dis),
      eTIV = as.numeric(eTIV)
    )
  
  if (any(is.na(out$side_lr))) {
    message("❌ side cannot be completely standardized to L/R in anatomical analysis.")
    print(table(df[[SIDE_COL_ANAT]], useNA = "ifany"))
    stop("Fix side coding or expand normalize_side_lr().")
  }
  
  out
}

run_anatomical_correlations <- function(df, roi_cols, family_fun, hemi_fun,
                                        display_fun, scale_meta, n_perm,
                                        min_n = MIN_N_ANAT) {
  results <- vector("list", length(roi_cols) * nrow(scale_meta))
  k <- 1L
  
  for (roi in roi_cols) {
    fam <- family_fun(roi)
    hemi <- hemi_fun(roi)
    
    for (ii in seq_len(nrow(scale_meta))) {
      sc <- scale_meta$Scale[ii]
      is_normed <- scale_meta$Is_Normed[ii]
      covars <- get_covars_anatomical(is_normed)
      
      required <- unique(c(roi, sc, covars))
      missing_now <- setdiff(required, names(df))
      if (length(missing_now) > 0) {
        stop("❌ Missing columns in anatomical analysis: ",
             paste(missing_now, collapse = ", "))
      }
      
      sub <- df %>% select(all_of(required)) %>% drop_na()
      n <- nrow(sub)
      
      if (n < min_n) {
        results[[k]] <- tibble(
          Family = fam, Hemisphere = hemi,
          ROI = roi, ROI_Display = display_fun(roi),
          Scale = sc, Scale_Display = to_scale_display(sc),
          Is_Normed = is_normed,
          CovarsUsed = paste(covars, collapse = " + "),
          n = n, r = NA_real_, perm_p = NA_real_
        )
        k <- k + 1L
        next
      }
      
      r_roi <- residualize_one_anat(sub[[roi]], sub, covars)
      r_sc  <- residualize_one_anat(sub[[sc]],  sub, covars)
      keep <- is.finite(r_roi) & is.finite(r_sc)
      n2 <- sum(keep)
      
      if (n2 < min_n) {
        results[[k]] <- tibble(
          Family = fam, Hemisphere = hemi,
          ROI = roi, ROI_Display = display_fun(roi),
          Scale = sc, Scale_Display = to_scale_display(sc),
          Is_Normed = is_normed,
          CovarsUsed = paste(covars, collapse = " + "),
          n = n2, r = NA_real_, perm_p = NA_real_
        )
        k <- k + 1L
        next
      }
      
      r_obs <- suppressWarnings(cor(r_roi[keep], r_sc[keep], method = "pearson"))
      p_perm <- perm_p_from_residuals_anat(
        r_roi[keep], r_sc[keep], n_perm = n_perm
      )
      
      results[[k]] <- tibble(
        Family = fam, Hemisphere = hemi,
        ROI = roi, ROI_Display = display_fun(roi),
        Scale = sc, Scale_Display = to_scale_display(sc),
        Is_Normed = is_normed,
        CovarsUsed = paste(covars, collapse = " + "),
        n = n2, r = as.numeric(r_obs), perm_p = as.numeric(p_perm)
      )
      k <- k + 1L
    }
  }
  
  bind_rows(results) %>%
    group_by(Scale, Family) %>%
    mutate(q_perm = p.adjust(perm_p, method = "BH")) %>%
    ungroup() %>%
    mutate(star = star_from_q_anat(q_perm))
}

save_anatomical_tables <- function(res, outdir, prefix) {
  write.csv(
    res,
    file.path(outdir, paste0(prefix, "_longtable_permBH.csv")),
    row.names = FALSE, fileEncoding = "GBK"
  )
  
  mat_r <- res %>%
    select(ROI, ROI_Display, Hemisphere, Scale, r) %>%
    pivot_wider(names_from = Scale, values_from = r)
  
  mat_q <- res %>%
    select(ROI, ROI_Display, Hemisphere, Scale, q_perm) %>%
    pivot_wider(names_from = Scale, values_from = q_perm)
  
  write.csv(
    mat_r,
    file.path(outdir, paste0(prefix, "_matrix_r.csv")),
    row.names = FALSE, fileEncoding = "GBK"
  )
  
  write.csv(
    mat_q,
    file.path(outdir, paste0(prefix, "_matrix_q_permBH.csv")),
    row.names = FALSE, fileEncoding = "GBK"
  )
  
  sum_tbl <- res %>%
    mutate(sig05 = is.finite(q_perm) & q_perm < 0.05) %>%
    group_by(Scale, Family) %>%
    summarise(
      n_test = sum(is.finite(perm_p)),
      n_sig_q05 = sum(sig05, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Scale, Family)
  
  write.csv(
    sum_tbl,
    file.path(outdir, paste0(prefix, "_summary_counts.csv")),
    row.names = FALSE, fileEncoding = "GBK"
  )
  
  invisible(list(mat_r = mat_r, mat_q = mat_q, summary = sum_tbl))
}

save_anatomical_heatmap <- function(res, roi_order, panel_levels,
                                    outdir, prefix, title_text,
                                    width = 13, height = 14) {
  scale_display_order <- scale_meta$Scale %>% to_scale_display()
  
  roi_display_order <- sapply(roi_order, function(x) {
    z <- res$ROI_Display[match(x, res$ROI)]
    ifelse(is.na(z), x, z)
  }, USE.NAMES = FALSE)
  
  hm <- res %>%
    mutate(
      Panel = case_when(
        Family == "HippocampusSubfields" & Hemisphere == "Left"  ~ "Left hippocampal subfields",
        Family == "HippocampusSubfields" & Hemisphere == "Right" ~ "Right hippocampal subfields",
        Family == "ThalamusNuclei" & Hemisphere == "Left"  ~ "Left thalamic nuclei",
        Family == "ThalamusNuclei" & Hemisphere == "Right" ~ "Right thalamic nuclei",
        Family == "Hippocampal segments" & Hemisphere == "Left"  ~ "Left hippocampal segments",
        Family == "Hippocampal segments" & Hemisphere == "Right" ~ "Right hippocampal segments",
        Family == "Thalamic composite groups" & Hemisphere == "Left"  ~ "Left thalamic composite groups",
        Family == "Thalamic composite groups" & Hemisphere == "Right" ~ "Right thalamic composite groups",
        TRUE ~ paste(Hemisphere, Family)
      ),
      Panel = factor(Panel, levels = panel_levels),
      ROI_Display = factor(ROI_Display, levels = rev(roi_display_order)),
      Scale_Display = factor(Scale_Display, levels = scale_display_order)
    )
  
  lim <- suppressWarnings(quantile(abs(hm$r), 0.98, na.rm = TRUE))
  if (!is.finite(lim) || lim < 0.2) lim <- 0.2
  
  p <- ggplot(hm, aes(x = Scale_Display, y = ROI_Display, fill = r)) +
    geom_tile(color = "white", linewidth = 0.2) +
    geom_text(
      aes(label = star),
      family = PLOT_FONT,
      size = 4, fontface = "bold",
      hjust = 0.5, vjust = 0.72, color = "black"
    ) +
    facet_grid(Panel ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-lim, lim), oob = scales::squish, na.value = "grey90"
    ) +
    theme_minimal(base_size = 12, base_family = PLOT_FONT) +
    theme(
      text = element_text(family = PLOT_FONT, face = "bold"),
      axis.text.x = element_text(
        angle = 35, hjust = 1, size = 10,
        face = "bold", color = "black", family = PLOT_FONT
      ),
      axis.text.y = element_text(
        size = 9.5, face = "bold", color = "black", family = PLOT_FONT
      ),
      strip.text.y = element_text(
        size = 11, face = "bold", family = PLOT_FONT
      ),
      plot.title = element_text(
        size = 15, face = "bold", family = PLOT_FONT
      ),
      plot.subtitle = element_text(
        size = 10, face = "bold", family = PLOT_FONT
      ),
      legend.title = element_text(family = PLOT_FONT, face = "bold"),
      legend.text = element_text(family = PLOT_FONT, face = "bold")
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Anatomical Left/Right ROIs; seizure-onset side included as a covariate. ",
        "Fill = partial r. Stars = q_perm after BH correction within Scale × Family; ",
        "Nperm = ", N_PERM_ANAT, ". ○<0.10, *<0.05, **<0.01, ***<0.001"
      ),
      x = "", y = "ROI (anatomical Left/Right)", fill = "partial r"
    )
  
  ggsave(
    file.path(outdir, paste0(prefix, "_heatmap_permBH.png")),
    p, width = width, height = height, units = "in", dpi = 600
  )
  ggsave(
    file.path(outdir, paste0(prefix, "_heatmap_permBH.pdf")),
    p, width = width, height = height, units = "in",
    device = grDevices::cairo_pdf
  )
  
  invisible(p)
}

save_anatomical_scatterplots <- function(res, df, outdir, only_q05 = FALSE) {
  scatter_dir <- file.path(outdir, ifelse(only_q05, "scatter_q05", "scatter_all"))
  if (!dir.exists(scatter_dir)) dir.create(scatter_dir, recursive = TRUE)
  
  plot_df <- if (only_q05) {
    res %>% filter(is.finite(q_perm) & q_perm < 0.05)
  } else {
    res
  }
  
  cat("\n▶ Anatomical Left/Right scatter plots: ", nrow(plot_df), " candidates\n", sep = "")
  if (nrow(plot_df) == 0) return(invisible(NULL))
  
  plot_count <- 0L
  
  for (ii in seq_len(nrow(plot_df))) {
    row_data <- plot_df[ii, ]
    if (!is.finite(row_data$r)) next
    
    roi <- row_data$ROI
    roi_disp <- row_data$ROI_Display
    sc <- row_data$Scale
    sc_disp <- row_data$Scale_Display
    
    covars_str <- row_data$CovarsUsed
    covars <- if (is.na(covars_str) || covars_str == "") {
      character(0)
    } else {
      trimws(unlist(strsplit(covars_str, "\\+")))
    }
    
    sub <- df %>% select(all_of(c(roi, sc, covars))) %>% drop_na()
    r_roi <- residualize_one_anat(sub[[roi]], sub, covars)
    r_sc  <- residualize_one_anat(sub[[sc]],  sub, covars)
    
    keep <- is.finite(r_roi) & is.finite(r_sc)
    if (sum(keep) < MIN_N_ANAT) next
    
    df_point <- tibble(ROI_res = r_roi[keep], Scale_res = r_sc[keep])
    
    r_val <- row_data$r
    r2_val <- r_val^2
    p_val <- row_data$perm_p
    q_val <- row_data$q_perm
    star <- row_data$star
    
    p_txt <- ifelse(is.na(p_val), "NA",
                    ifelse(p_val < 0.001, "<0.001", sprintf("%.3f", p_val)))
    q_txt <- ifelse(is.na(q_val), "NA",
                    ifelse(q_val < 0.001, "<0.001", sprintf("%.3f", q_val)))
    
    annot_txt <- sprintf(
      "Partial r = %.3f (R² = %.3f)\np_perm = %s\nq_perm = %s %s",
      r_val, r2_val, p_txt, q_txt, star
    )
    
    p_scatter <- ggplot(df_point, aes(x = ROI_res, y = Scale_res)) +
      geom_point(alpha = 0.6, size = 2, color = "#2C3E50") +
      geom_smooth(
        method = "lm", formula = y ~ x,
        color = "#B2182B", fill = "grey70", alpha = 0.3
      ) +
      annotate(
        "text", x = -Inf, y = Inf, label = annot_txt,
        hjust = -0.05, vjust = 1.2, size = 4.6,
        color = "grey20", lineheight = 1.1,
        family = PLOT_FONT, fontface = "bold"
      ) +
      theme_minimal(base_size = 14, base_family = PLOT_FONT) +
      theme(
        text = element_text(family = PLOT_FONT, face = "bold"),
        plot.title = element_text(face = "bold", size = 13.5, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      ) +
      labs(
        title = paste(roi_disp, "vs", sc_disp),
        subtitle = "Anatomical Left/Right analysis (seizure-onset side adjusted)",
        x = paste(roi_disp, "(residuals)"),
        y = paste(sc_disp, "(residuals)")
      )
    
    safe_roi <- gsub("[\\/:*?\"<>|]", "_", roi)
    safe_sc  <- gsub("[\\/:*?\"<>|]", "_", sc)
    
    ggsave(
      filename = file.path(scatter_dir, paste0(safe_roi, "_vs_", safe_sc, ".png")),
      plot = p_scatter,
      width = 6, height = 5, dpi = 600, bg = "white"
    )
    
    plot_count <- plot_count + 1L
    if (plot_count %% 25 == 0) {
      cat(sprintf("Anatomical scatter progress: %d plots saved.\n", plot_count))
    }
  }
  
  cat("✅ Anatomical scatter plots saved to: ", scatter_dir, "\n", sep = "")
  invisible(NULL)
}

# ==============================================================================
# A) FINE-GRAINED ANATOMICAL LEFT/RIGHT ANALYSIS
#    Hippocampal subfields + individual thalamic nuclei
# ==============================================================================
cat("\n\n================ ANATOMICAL LEFT/RIGHT: FINE-GRAINED ================\n")

df_anat_fine <- prepare_anatomical_df(df_raw_clean)

hip_left_fine <- c(
  "Left_Hippocampal_tail", "Left_subiculum_comb", "Left_CA1_comb",
  "Left_CA2/3_comb", "Left_CA4_comb", "Left_GC_ML_DG_comb",
  "Left_molecular_layer_HP_comb", "Left_HATA", "Left_fimbria",
  "Left_presubiculum_comb", "Left_parasubiculum"
)

hip_right_fine <- c(
  "Right_Hippocampal_tail", "Right_subiculum_comb", "Right_CA1_comb",
  "Right_CA2/3_comb", "Right_CA4_comb", "Right_GC_ML_DG_comb",
  "Right_molecular_layer_HP_comb", "Right_HATA", "Right_fimbria",
  "Right_presubiculum_comb", "Right_parasubiculum"
)

thal_left_fine <- c(
  "Left_AV", "Left_VA", "Left_VAmc", "Left_VLa", "Left_VLp", "Left_VPL",
  "Left_VM", "Left_CL", "Left_CeM", "Left_CM", "Left_Pf", "Left_MDm",
  "Left_MDl", "Left_LD", "Left_LP", "Left_Pu_Total", "Left_MGN", "Left_LGN",
  "Left_L_Sg", "Left_MV_Re", "Left_Pc", "Left_Pt"
)

thal_right_fine <- c(
  "Right_AV", "Right_VA", "Right_VAmc", "Right_VLa", "Right_VLp", "Right_VPL",
  "Right_VM", "Right_CL", "Right_CeM", "Right_CM", "Right_Pf", "Right_MDm",
  "Right_MDl", "Right_LD", "Right_LP", "Right_Pu_Total", "Right_MGN", "Right_LGN",
  "Right_L_Sg", "Right_MV_Re", "Right_Pc", "Right_Pt"
)

roi_fine_anat <- c(
  hip_left_fine, hip_right_fine,
  thal_left_fine, thal_right_fine
)

miss_fine <- setdiff(roi_fine_anat, names(df_anat_fine))
if (length(miss_fine) > 0) {
  stop("❌ Fine-grained anatomical ROI columns missing:\n- ",
       paste(miss_fine, collapse = "\n- "))
}

hip_display_map_anat <- c(
  "Hippocampal_tail" = "Hippocampal tail",
  "subiculum_comb" = "Subiculum",
  "CA1_comb" = "CA1",
  "CA2/3_comb" = "CA2/3",
  "CA4_comb" = "CA4",
  "GC_ML_DG_comb" = "GC-ML-DG",
  "molecular_layer_HP_comb" = "Molecular layer HP",
  "HATA" = "HATA",
  "fimbria" = "Fimbria",
  "presubiculum_comb" = "Presubiculum",
  "parasubiculum" = "Parasubiculum"
)

fine_family_fun <- function(roi) {
  if (roi %in% c(hip_left_fine, hip_right_fine)) return("HippocampusSubfields")
  if (roi %in% c(thal_left_fine, thal_right_fine)) return("ThalamusNuclei")
  "Other"
}

fine_hemi_fun <- function(roi) {
  if (grepl("^Left_", roi)) return("Left")
  if (grepl("^Right_", roi)) return("Right")
  NA_character_
}

fine_display_fun <- function(roi) {
  hemi <- fine_hemi_fun(roi)
  suffix <- sub("^(Left_|Right_)", "", roi)
  
  if (roi %in% c(hip_left_fine, hip_right_fine)) {
    disp <- unname(hip_display_map_anat[suffix])
    if (is.na(disp)) disp <- suffix
    return(paste(hemi, disp))
  }
  
  if (roi %in% c(thal_left_fine, thal_right_fine)) {
    return(paste(hemi, suffix))
  }
  
  roi
}

cat("Fine-grained anatomical ROI count: ", length(roi_fine_anat), "\n", sep = "")
cat("Expected: 22 hippocampal + 44 thalamic = 66 ROIs\n")

res_anat_fine <- run_anatomical_correlations(
  df = df_anat_fine,
  roi_cols = roi_fine_anat,
  family_fun = fine_family_fun,
  hemi_fun = fine_hemi_fun,
  display_fun = fine_display_fun,
  scale_meta = scale_meta,
  n_perm = N_PERM_ANAT,
  min_n = MIN_N_ANAT
)

fine_prefix <- "Step2_Anatomical_LeftRight_FineGrained"
save_anatomical_tables(res_anat_fine, ANAT_FINE_DIR, fine_prefix)

fine_roi_order <- c(
  hip_left_fine, hip_right_fine,
  thal_left_fine, thal_right_fine
)

fine_panels <- c(
  "Left hippocampal subfields",
  "Right hippocampal subfields",
  "Left thalamic nuclei",
  "Right thalamic nuclei"
)

p_anat_fine <- save_anatomical_heatmap(
  res = res_anat_fine,
  roi_order = fine_roi_order,
  panel_levels = fine_panels,
  outdir = ANAT_FINE_DIR,
  prefix = fine_prefix,
  title_text = "Anatomical Left/Right ROI–Cognition partial correlations",
  width = 13, height = 14
)

# Fine-grained scatter plots: only q_perm < 0.05 to avoid excessive output.
save_anatomical_scatterplots(
  res = res_anat_fine,
  df = df_anat_fine,
  outdir = ANAT_FINE_DIR,
  only_q05 = TRUE
)

cat("✅ Fine-grained anatomical Left/Right analysis saved to:\n", ANAT_FINE_DIR, "\n", sep = "")

# ==============================================================================
# B) LOWER-DIMENSIONAL ANATOMICAL LEFT/RIGHT ANALYSIS
#    Hippocampal head/body/tail + thalamic composite groups
# ==============================================================================
cat("\n\n================ ANATOMICAL LEFT/RIGHT: LOW-DIMENSIONAL ================\n")

df_anat_low <- prepare_anatomical_df(df_raw_clean)

hip_segment_suffix_anat <- c(
  "Whole_hippocampal_head",
  "Whole_hippocampal_body",
  "Hippocampal_tail"
)

hip_left_low <- paste0("Left_", hip_segment_suffix_anat)
hip_right_low <- paste0("Right_", hip_segment_suffix_anat)

hip_segment_display_anat <- c(
  "Whole_hippocampal_head" = "Hippocampal head",
  "Whole_hippocampal_body" = "Hippocampal body",
  "Hippocampal_tail" = "Hippocampal tail"
)

thal_group_def_anat <- list(
  AnteriorDorsalGroup = c("AV", "LD"),
  MediodorsalNuclei   = c("MDm", "MDl"),
  IntralaminarNuclei  = c("CL", "CM", "CeM", "Pf", "Pc"),
  PulvinarLPComplex   = c("Pu_Total", "L_Sg", "LP"),
  VentralNuclei       = c("VA", "VAmc", "VLa", "VLp", "VPL", "VM"),
  MidlineNuclei       = c("MV_Re", "Pt"),
  GeniculateNuclei    = c("LGN", "MGN")
)

thal_group_display_anat <- c(
  "AnteriorDorsalGroup" = "Anterior-dorsal group (AV+LD)",
  "MediodorsalNuclei"   = "Mediodorsal nuclei (MDm+MDl)",
  "IntralaminarNuclei"  = "Intralaminar nuclei (CL+CM+CeM+Pf+Pc)",
  "PulvinarLPComplex"   = "Pulvinar-LP complex (Pu_Total+L_Sg+LP)",
  "VentralNuclei"       = "Ventral nuclei (VA+VAmc+VLa+VLp+VPL+VM)",
  "MidlineNuclei"       = "Midline nuclei (MV_Re+Pt)",
  "GeniculateNuclei"    = "Geniculate nuclei (LGN+MGN)"
)

miss_hip_low <- setdiff(c(hip_left_low, hip_right_low), names(df_anat_low))
if (length(miss_hip_low) > 0) {
  stop("❌ Low-dimensional anatomical hippocampal columns missing:\n- ",
       paste(miss_hip_low, collapse = "\n- "))
}

THAL_GROUP_MIN_PROP_ANAT <- if (exists("THAL_GROUP_MIN_PROP")) THAL_GROUP_MIN_PROP else 1.0

row_sum_min_n_anat <- function(df, vars, min_prop = 1.0) {
  vars_exist <- intersect(vars, names(df))
  if (length(vars_exist) == 0) return(rep(NA_real_, nrow(df)))
  
  x_df <- df[, vars_exist, drop = FALSE] %>%
    mutate(across(everything(), as.numeric))
  
  x <- as.matrix(x_df)
  nonmiss_n <- rowSums(!is.na(x))
  need_n <- ceiling(length(vars_exist) * min_prop)
  
  out <- rowSums(x, na.rm = TRUE)
  out[nonmiss_n < need_n] <- NA_real_
  out
}

for (gname in names(thal_group_def_anat)) {
  left_members  <- paste0("Left_",  thal_group_def_anat[[gname]])
  right_members <- paste0("Right_", thal_group_def_anat[[gname]])
  
  miss_l <- setdiff(left_members, names(df_anat_low))
  miss_r <- setdiff(right_members, names(df_anat_low))
  
  if (length(miss_l) > 0 || length(miss_r) > 0) {
    stop(
      "❌ Missing anatomical thalamic members for group ", gname, ":\n",
      "Left missing: ", paste(miss_l, collapse = ", "), "\n",
      "Right missing: ", paste(miss_r, collapse = ", ")
    )
  }
  
  df_anat_low[[paste0("LeftGrp_", gname)]] <- row_sum_min_n_anat(
    df_anat_low, left_members, min_prop = THAL_GROUP_MIN_PROP_ANAT
  )
  df_anat_low[[paste0("RightGrp_", gname)]] <- row_sum_min_n_anat(
    df_anat_low, right_members, min_prop = THAL_GROUP_MIN_PROP_ANAT
  )
}

thal_left_low <- paste0("LeftGrp_", names(thal_group_def_anat))
thal_right_low <- paste0("RightGrp_", names(thal_group_def_anat))

roi_low_anat <- c(
  hip_left_low, hip_right_low,
  thal_left_low, thal_right_low
)

low_family_fun <- function(roi) {
  if (roi %in% c(hip_left_low, hip_right_low)) return("Hippocampal segments")
  if (roi %in% c(thal_left_low, thal_right_low)) return("Thalamic composite groups")
  "Other"
}

low_hemi_fun <- function(roi) {
  if (grepl("^Left", roi)) return("Left")
  if (grepl("^Right", roi)) return("Right")
  NA_character_
}

low_display_fun <- function(roi) {
  hemi <- low_hemi_fun(roi)
  
  if (roi %in% c(hip_left_low, hip_right_low)) {
    suffix <- sub("^(Left_|Right_)", "", roi)
    disp <- unname(hip_segment_display_anat[suffix])
    if (is.na(disp)) disp <- suffix
    return(paste(hemi, disp))
  }
  
  if (roi %in% c(thal_left_low, thal_right_low)) {
    suffix <- sub("^(LeftGrp_|RightGrp_)", "", roi)
    disp <- unname(thal_group_display_anat[suffix])
    if (is.na(disp)) disp <- suffix
    return(paste(hemi, "thalamus –", disp))
  }
  
  roi
}

cat("Low-dimensional anatomical ROI count: ", length(roi_low_anat), "\n", sep = "")
cat("Expected: 6 hippocampal segments + 14 thalamic composite ROIs = 20 ROIs\n")

res_anat_low <- run_anatomical_correlations(
  df = df_anat_low,
  roi_cols = roi_low_anat,
  family_fun = low_family_fun,
  hemi_fun = low_hemi_fun,
  display_fun = low_display_fun,
  scale_meta = scale_meta,
  n_perm = N_PERM_ANAT,
  min_n = MIN_N_ANAT
)

low_prefix <- "Step2_Anatomical_LeftRight_ThalSUM_HipHeadBodyTail"
save_anatomical_tables(res_anat_low, ANAT_LOW_DIR, low_prefix)

hip_def_anat <- tibble(
  Segment = names(hip_segment_display_anat),
  Display = unname(hip_segment_display_anat),
  Left_Column = paste0("Left_", names(hip_segment_display_anat)),
  Right_Column = paste0("Right_", names(hip_segment_display_anat))
)

write.csv(
  hip_def_anat,
  file.path(ANAT_LOW_DIR, "Anatomical_Hippocampal_HeadBodyTail_Definitions.csv"),
  row.names = FALSE, fileEncoding = "GBK"
)

thal_def_anat <- tibble(
  Group = names(thal_group_def_anat),
  Display = unname(thal_group_display_anat[names(thal_group_def_anat)]),
  Members = sapply(thal_group_def_anat, function(x) paste(x, collapse = " + ")),
  Left_Members = sapply(thal_group_def_anat, function(x) paste0("Left_", x, collapse = "; ")),
  Right_Members = sapply(thal_group_def_anat, function(x) paste0("Right_", x, collapse = "; "))
)

write.csv(
  thal_def_anat,
  file.path(ANAT_LOW_DIR, "Anatomical_Thalamic_Composite_Group_Definitions_SUM.csv"),
  row.names = FALSE, fileEncoding = "GBK"
)

low_roi_order <- c(
  hip_left_low, hip_right_low,
  thal_left_low, thal_right_low
)

low_panels <- c(
  "Left hippocampal segments",
  "Right hippocampal segments",
  "Left thalamic composite groups",
  "Right thalamic composite groups"
)

p_anat_low <- save_anatomical_heatmap(
  res = res_anat_low,
  roi_order = low_roi_order,
  panel_levels = low_panels,
  outdir = ANAT_LOW_DIR,
  prefix = low_prefix,
  title_text = "Anatomical Left/Right low-dimensional ROI–Cognition partial correlations",
  width = 13, height = 8.5
)

# Lower-dimensional analysis has only 20 ROIs, so save all scatter plots.
save_anatomical_scatterplots(
  res = res_anat_low,
  df = df_anat_low,
  outdir = ANAT_LOW_DIR,
  only_q05 = FALSE
)

cat("✅ Low-dimensional anatomical Left/Right analysis saved to:\n", ANAT_LOW_DIR, "\n", sep = "")

cat("\n============================================================\n")
cat("✅ ALL ANATOMICAL LEFT/RIGHT ANALYSES COMPLETED\n")
cat("Fine-grained output: ", ANAT_FINE_DIR, "\n", sep = "")
cat("Low-dimensional output: ", ANAT_LOW_DIR, "\n", sep = "")
cat("============================================================\n")


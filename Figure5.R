# =============================================================================
# Figure 5. Incremental explained variance in cognition by hippocampal coding
#
# 四种海马编码方式：
#   1. Ipsilateral to the epileptogenic hemisphere
#   2. Contralateral
#   3. Anatomical left
#   4. Anatomical right
#
# 每种海马指标均与 GLOBAL thalamic summed volume 配对。
#
# 四组嵌套模型比较：
#   M1 - M0 : Hippocampus alone
#   M3 - M2 : Hippocampus beyond thalamus
#   M2 - M0 : Thalamus alone
#   M3 - M1 : Thalamus beyond hippocampus
#
# 输出：
#   Figure_5.tif
#   Figure_5.pdf
# =============================================================================


# =============================================================================
# 0. 加载R包
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})


# =============================================================================
# 1. 字符编码设置
# =============================================================================
# 用于兼容原始CSV中的中文认知指标名称

if (!isTRUE(l10n_info()[["MBCS"]])) {
  
  for (loc in c(
    "C.UTF-8",
    "en_US.UTF-8",
    "C.utf8",
    "Chinese_China.936"
  )) {
    
    res <- suppressWarnings(
      try(Sys.setlocale("LC_CTYPE", loc), silent = TRUE)
    )
    
    if (!inherits(res, "try-error") &&
        !is.na(res) &&
        nzchar(res)) {
      break
    }
  }
}


# =============================================================================
# 2. 输入和输出路径
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)


# -----------------------------------------------------------------------------
# 输入CSV
# -----------------------------------------------------------------------------

INFILE <- if (length(args) >= 1) {
  
  args[1]
  
} else {
  
  "C:/Users/86150/Documents/HIPP_THAL/additive_model/Anatomical_LeftRight/volume_sensitivity_alignment_vs_anatomical_comparison.csv"
  
}


# -----------------------------------------------------------------------------
# 输出文件夹
# -----------------------------------------------------------------------------

OUTDIR <- "C:/Users/86150/Documents/HIPP_THAL/additive_model/Anatomical_LeftRight"


if (!dir.exists(OUTDIR)) {
  
  dir.create(
    OUTDIR,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
}


OUT_TIF <- file.path(
  OUTDIR,
  "Figure_5.tif"
)

OUT_PDF <- file.path(
  OUTDIR,
  "Figure_5.pdf"
)


# -----------------------------------------------------------------------------
# 检查输入文件是否存在
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat("Figure 5\n")
cat("============================================================\n")

cat("\nInput file:\n")
cat(INFILE, "\n")

cat("\nFile exists: ")
cat(file.exists(INFILE), "\n")


if (!file.exists(INFILE)) {
  
  stop(
    "\nInput CSV does not exist:\n",
    INFILE,
    "\n\nPlease check the folder name and CSV filename.\n"
  )
  
}


# =============================================================================
# 3. 字体
# =============================================================================
# 整张图使用 Arial
#
# 不再使用 windowsFonts()，
# PDF 使用 cairo_pdf，避免之前的 CID font 报错。
# =============================================================================

FONT <- "Arial"


if (requireNamespace("systemfonts", quietly = TRUE)) {
  
  font_info <- try(
    systemfonts::match_font(FONT),
    silent = TRUE
  )
  
  if (!inherits(font_info, "try-error")) {
    
    cat("\nRequested font:", FONT, "\n")
    
    if (!is.null(font_info$path) &&
        length(font_info$path) > 0 &&
        nzchar(font_info$path[1])) {
      
      cat("Font file found:\n")
      cat(font_info$path[1], "\n")
      
    }
  }
}


if (!capabilities("cairo")) {
  
  stop(
    "\nThis R installation does not support Cairo graphics.\n",
    "Cairo support is required for reliable Arial PDF output.\n"
  )
  
}


# =============================================================================
# 4. 读取CSV
# =============================================================================

read_csv_smart <- function(f) {
  
  if (!file.exists(f)) {
    
    stop(
      "Input file not found:\n",
      normalizePath(
        f,
        mustWork = FALSE
      )
    )
    
  }
  
  
  readers <- list(
    
    GBK = function() {
      
      read.csv(
        f,
        fileEncoding = "GBK",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
    },
    
    
    GB18030 = function() {
      
      read.csv(
        f,
        fileEncoding = "GB18030",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
    },
    
    
    UTF8_BOM = function() {
      
      read.csv(
        f,
        fileEncoding = "UTF-8-BOM",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
    },
    
    
    UTF8 = function() {
      
      read.csv(
        f,
        fileEncoding = "UTF-8",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
    },
    
    
    native = function() {
      
      read.csv(
        f,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
    }
    
  )
  
  
  for (nm in names(readers)) {
    
    x <- try(
      suppressWarnings(
        readers[[nm]]()
      ),
      silent = TRUE
    )
    
    
    if (!inherits(x, "try-error") &&
        is.data.frame(x) &&
        nrow(x) > 0) {
      
      
      # 去掉可能存在的UTF-8 BOM
      names(x) <- sub(
        "^\ufeff",
        "",
        names(x)
      )
      
      
      # 去掉列名前后的空格
      names(x) <- trimws(
        names(x)
      )
      
      
      if ("Set" %in% names(x)) {
        
        cat("\nCSV successfully read.\n")
        cat("Encoding method:", nm, "\n")
        cat("Rows:", nrow(x), "\n")
        cat("Columns:", ncol(x), "\n")
        
        return(x)
        
      }
      
    }
    
  }
  
  
  stop(
    "\nThe CSV exists, but it could not be read correctly.\n",
    "The required column 'Set' could not be identified.\n"
  )
  
}


d <- read_csv_smart(INFILE)


# =============================================================================
# 5. 检查所需变量
# =============================================================================

required_cols <- c(
  
  "Set",
  "Outcome",
  
  "DeltaR2_M1_vs_M0",
  "Pperm_M1_vs_M0",
  "Q_M1_vs_M0",
  
  "DeltaR2_M3_vs_M2",
  "Pperm_M3_vs_M2",
  "Q_M3_vs_M2",
  
  "DeltaR2_M2_vs_M0",
  "Pperm_M2_vs_M0",
  "Q_M2_vs_M0",
  
  "DeltaR2_M3_vs_M1",
  "Pperm_M3_vs_M1",
  "Q_M3_vs_M1"
  
)


missing_cols <- setdiff(
  required_cols,
  names(d)
)


if (length(missing_cols) > 0) {
  
  stop(
    "\nThe following required columns are missing:\n",
    paste(
      missing_cols,
      collapse = "\n"
    ),
    "\n"
  )
  
}


# =============================================================================
# 6. 四种海马编码方式
# =============================================================================

coding_map <- c(
  
  Primary_Ipsi_Global =
    "Ipsilateral\n(seizure-side aligned)",
  
  Primary_Contra_Global =
    "Contralateral",
  
  Sensitivity_Anatomical_Left_Global =
    "Anatomical left",
  
  Sensitivity_Anatomical_Right_Global =
    "Anatomical right"
  
)


d <- d[
  d$Set %in% names(coding_map),
  ,
  drop = FALSE
]


if (nrow(d) == 0) {
  
  stop(
    "\nNone of the expected analysis sets were found in column 'Set'.\n"
  )
  
}


d$Coding <- factor(
  
  unname(
    coding_map[
      as.character(d$Set)
    ]
  ),
  
  levels = unname(coding_map)
  
)


# =============================================================================
# 7. 认知指标名称
# =============================================================================

out_map <- c(
  
  "\u5de5\u4f5c\u8bb0\u5fc6\u6307\u6570WMI" =
    "WAIS WMI",
  
  "\u603b\u8bb0\u5fc6\u5546" =
    "WMS FSMQ",
  
  "\u542c\u89c9\u8bb0\u5fc6\u6307\u6570" =
    "WMS AMI",
  
  "\u5373\u523b\u8bb0\u5fc6\u6307\u6570" =
    "WMS IMI",
  
  "\u5ef6\u8fdf\u8bb0\u5fc6\u6307\u6570" =
    "WMS DMI",
  
  "\u89c6\u89c9\u8bb0\u5fc6\u6307\u6570" =
    "WMS VMI"
  
)


outcome_raw <- enc2utf8(
  trimws(
    as.character(d$Outcome)
  )
)


d$Outcome2 <- unname(
  out_map[
    outcome_raw
  ]
)


if (any(is.na(d$Outcome2))) {
  
  unmatched <- unique(
    outcome_raw[
      is.na(d$Outcome2)
    ]
  )
  
  cat("\nWARNING: unmatched outcome names:\n")
  print(unmatched)
  
}


if (all(is.na(d$Outcome2))) {
  
  stop(
    "\nNone of the Outcome values matched the expected memory indices.\n"
  )
  
}


# 控制Y轴上六个认知指标的排列顺序
d$Outcome2 <- factor(
  
  d$Outcome2,
  
  levels = rev(
    c(
      "WAIS WMI",
      "WMS FSMQ",
      "WMS AMI",
      "WMS IMI",
      "WMS DMI",
      "WMS VMI"
    )
  )
  
)


d <- d[
  !is.na(d$Outcome2),
  ,
  drop = FALSE
]


# =============================================================================
# 8. 四种nested contrasts转成长数据
# =============================================================================

contr <- list(
  
  c(
    "DeltaR2_M1_vs_M0",
    "Pperm_M1_vs_M0",
    "Q_M1_vs_M0",
    "Hippocampus alone (M1-M0)"
  ),
  
  c(
    "DeltaR2_M3_vs_M2",
    "Pperm_M3_vs_M2",
    "Q_M3_vs_M2",
    "Hippocampus beyond thalamus (M3-M2)"
  ),
  
  c(
    "DeltaR2_M2_vs_M0",
    "Pperm_M2_vs_M0",
    "Q_M2_vs_M0",
    "Thalamus alone (M2-M0)"
  ),
  
  c(
    "DeltaR2_M3_vs_M1",
    "Pperm_M3_vs_M1",
    "Q_M3_vs_M1",
    "Thalamus beyond hippocampus (M3-M1)"
  )
  
)


long <- do.call(
  
  rbind,
  
  lapply(
    
    contr,
    
    function(cc) {
      
      data.frame(
        
        Coding = d$Coding,
        
        Outcome2 = d$Outcome2,
        
        dR2 = suppressWarnings(
          as.numeric(
            d[[cc[1]]]
          )
        ),
        
        p = suppressWarnings(
          as.numeric(
            d[[cc[2]]]
          )
        ),
        
        q = suppressWarnings(
          as.numeric(
            d[[cc[3]]]
          )
        ),
        
        Contrast = cc[4],
        
        stringsAsFactors = FALSE
        
      )
      
    }
    
  )
  
)


long$Contrast <- factor(
  
  long$Contrast,
  
  levels = vapply(
    contr,
    `[`,
    character(1),
    4
  )
  
)


# =============================================================================
# 9. p值和q值标签
# =============================================================================

fmt_pq <- function(x) {
  
  out <- rep(
    "NA",
    length(x)
  )
  
  ok <- !is.na(x)
  
  
  out[
    ok & x < 0.001
  ] <- "<.001"
  
  
  idx <- ok & x >= 0.001
  
  
  out[idx] <- sub(
    
    "^0",
    
    "",
    
    formatC(
      x[idx],
      format = "f",
      digits = 3
    )
    
  )
  
  
  out
  
}


long$lab <- paste0(
  
  "p=",
  fmt_pq(long$p),
  
  ", q=",
  fmt_pq(long$q)
  
)


# =============================================================================
# 10. 配色
# =============================================================================

cols <- c(
  
  "Hippocampus alone (M1-M0)" =
    "#b9d0e8",
  
  "Hippocampus beyond thalamus (M3-M2)" =
    "#2f6ba8",
  
  "Thalamus alone (M2-M0)" =
    "#f6c9c6",
  
  "Thalamus beyond hippocampus (M3-M1)" =
    "#b0505a"
  
)


# =============================================================================
# 11. X轴范围
# =============================================================================
#
# 【修改X轴范围的位置】
#
# 当前固定为：
#   0 ～ 0.30
#
# 如果以后想改成0～0.4：
#   X_MIN <- 0
#   X_MAX <- 0.40
#
# =============================================================================

X_MIN <- 0
X_MAX <- 0.30


# 【修改X轴刻度间隔的位置】
#
# 当前每0.05一个刻度：
# 0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30
#
# 如果想每0.10一个刻度，将by改为0.10
#
X_BREAKS <- seq(
  X_MIN,
  X_MAX,
  by = 0.05
)


# =============================================================================
# 12. 绘图
# =============================================================================

dodge_width <- 0.82


p <- ggplot(
  
  long,
  
  aes(
    x = dR2,
    y = Outcome2,
    fill = Contrast
  )
  
) +
  
  
  # ---------------------------------------------------------------------------
# 柱形图
# ---------------------------------------------------------------------------

geom_col(
  
  position = position_dodge(
    width = dodge_width
  ),
  
  width = 0.76,
  
  colour = "black",
  
  linewidth = 0.30,
  
  na.rm = TRUE
  
) +
  
  
  # ---------------------------------------------------------------------------
# ΔR² = 0.05参考线
# ---------------------------------------------------------------------------

geom_vline(
  
  xintercept = 0.05,
  
  linetype = "dashed",
  
  colour = "grey45",
  
  linewidth = 0.60
  
) +
  
  
  # ---------------------------------------------------------------------------
# p值和q值
#
# 【修改p值/q值字体大小的位置】
#
# 当前：
#   size = 3.5
#
# 想放大可改为：
#   size = 3.8
#   size = 4.0
#
# 注意：
#   fontface = "plain"
# 表示p/q不加粗
# ---------------------------------------------------------------------------

geom_text(
  
  aes(
    label = lab
  ),
  
  position = position_dodge(
    width = dodge_width
  ),
  
  hjust = -0.06,
  
  family = FONT,
  
  fontface = "plain",
  
  colour = "black",
  
  size = 3.5,
  
  na.rm = TRUE
  
) +
  
  
  # ---------------------------------------------------------------------------
# 四个panel
# ---------------------------------------------------------------------------

facet_wrap(
  
  ~ Coding,
  
  ncol = 2
  
) +
  
  
  # ---------------------------------------------------------------------------
# 配色
# ---------------------------------------------------------------------------

scale_fill_manual(
  
  values = cols,
  
  name = NULL,
  
  drop = FALSE
  
) +
  
  
  # ---------------------------------------------------------------------------
# X轴
#
# 范围固定为0～0.30
# ---------------------------------------------------------------------------

scale_x_continuous(
  
  name = "Incremental explained variance (ΔR²)",
  
  limits = c(
    X_MIN,
    X_MAX
  ),
  
  breaks = X_BREAKS,
  
  labels = function(x) {
    
    ifelse(
      x == 0,
      "0",
      sprintf("%.2f", x)
    )
    
  },
  
  expand = expansion(
    mult = c(
      0,
      0.01
    )
  )
  
) +
  
  
  # ---------------------------------------------------------------------------
# 删除原图最下面的caption
#
# 原来的：
#
# caption = "Each bar is labelled ..."
#
# 已经完全删除。
# ---------------------------------------------------------------------------

labs(
  y = NULL
) +
  
  
  # ---------------------------------------------------------------------------
# 基础主题
#
# 【修改全图基础字体大小的位置】
#
# base_size = 14
#
# 但是下面各元素都有单独的size设置，
# 因此一般建议直接修改下面对应的位置。
# ---------------------------------------------------------------------------

theme_classic(
  
  base_size = 14,
  
  base_family = FONT
  
) +
  
  
  theme(
    
    
    # ========================================================================
    # 全局字体
    # ========================================================================
    
    text = element_text(
      
      family = FONT,
      
      face = "bold",
      
      colour = "black"
      
    ),
    
    
    # ========================================================================
    # X轴刻度数字字体大小
    #
    # 【修改这里】
    #
    # 当前 size = 12.5
    # ========================================================================
    
    axis.text.x = element_text(
      
      family = FONT,
      
      face = "bold",
      
      colour = "black",
      
      size = 12.5
      
    ),
    
    
    # ========================================================================
    # Y轴认知指标名称字体大小
    #
    # 例如：
    # WAIS WMI
    # WMS FSMQ
    # WMS AMI
    #
    # 【修改这里】
    #
    # 当前 size = 13.5
    # ========================================================================
    
    axis.text.y = element_text(
      
      family = FONT,
      
      face = "bold",
      
      colour = "black",
      
      size = 13.5,
      
      margin = margin(
        r = 7
      )
      
    ),
    
    
    # ========================================================================
    # X轴标题字体大小
    #
    # Incremental explained variance (ΔR²)
    #
    # 【修改这里】
    #
    # 当前 size = 15.5
    # ========================================================================
    
    axis.title.x = element_text(
      
      family = FONT,
      
      face = "bold",
      
      colour = "black",
      
      size = 15.5,
      
      margin = margin(
        t = 12
      )
      
    ),
    
    
    axis.title.y = element_blank(),
    
    
    # ========================================================================
    # 坐标轴
    # ========================================================================
    
    axis.line = element_line(
      
      colour = "black",
      
      linewidth = 0.55
      
    ),
    
    
    axis.ticks = element_line(
      
      colour = "black",
      
      linewidth = 0.45
      
    ),
    
    
    # ========================================================================
    # 四个panel标题字体大小
    #
    # Ipsilateral
    # Contralateral
    # Anatomical left
    # Anatomical right
    #
    # 【修改这里】
    #
    # 当前 size = 14.5
    # ========================================================================
    
    strip.text = element_text(
      
      family = FONT,
      
      face = "bold",
      
      colour = "black",
      
      size = 14.5,
      
      margin = margin(
        t = 8,
        b = 8
      )
      
    ),
    
    
    strip.background = element_rect(
      
      fill = "grey92",
      
      colour = "grey45",
      
      linewidth = 0.40
      
    ),
    
    
    # ========================================================================
    # 图例字体大小
    #
    # Hippocampus alone
    # Hippocampus beyond thalamus
    # Thalamus alone
    # Thalamus beyond hippocampus
    #
    # 【修改这里】
    #
    # 当前 size = 11.5
    # ========================================================================
    
    legend.text = element_text(
      
      family = FONT,
      
      face = "bold",
      
      colour = "black",
      
      size = 11.5
      
    ),
    
    
    legend.position = "bottom",
    
    
    legend.key.width = grid::unit(
      1.25,
      "cm"
    ),
    
    legend.key.height = grid::unit(
      0.48,
      "cm"
    ),
    
    legend.spacing.x = grid::unit(
      0.15,
      "cm"
    ),
    
    
    # ========================================================================
    # 网格线
    # ========================================================================
    
    panel.grid.major.x = element_line(
      
      colour = "grey90",
      
      linewidth = 0.35
      
    ),
    
    panel.grid.major.y = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    
    # ========================================================================
    # 四个panel之间的间隔
    #
    # 【如果希望四个panel彼此更分开，可以修改这里】
    #
    # 当前 = 1.5 lines
    # ========================================================================
    
    panel.spacing = grid::unit(
      1.5,
      "lines"
    ),
    
    
    # ========================================================================
    # 图片外围空白
    # ========================================================================
    
    plot.margin = margin(
      
      t = 18,
      
      r = 30,
      
      b = 18,
      
      l = 18
      
    )
    
  ) +
  
  
  # ---------------------------------------------------------------------------
# 图例排列为两行
# ---------------------------------------------------------------------------

guides(
  
  fill = guide_legend(
    
    nrow = 2,
    
    byrow = TRUE,
    
    override.aes = list(
      
      colour = "black",
      
      linewidth = 0.25
      
    )
    
  )
  
)


# =============================================================================
# 13. 在RStudio中显示图片
# =============================================================================

print(p)


# =============================================================================
# 14. 保存TIFF
# =============================================================================
#
# 【修改图片整体尺寸的位置】
#
# width  = 图片宽度
# height = 图片高度
#
# 当前：
#   width  = 14.5 inch
#   height = 13.5 inch
#
# 相比之前11.5 inch进一步加高，
# 因此各认知指标及p/q标签之间会更疏一些。
#
# 如果还希望更高：
#   height = 14.5
# 或：
#   height = 15
#
# =============================================================================

FIG_WIDTH  <- 14.5
FIG_HEIGHT <- 13.5


ggsave(
  
  filename = OUT_TIF,
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  units = "in",
  
  dpi = 600,
  
  device = "tiff",
  
  compression = "lzw",
  
  bg = "white"
  
)


# =============================================================================
# 15. 保存PDF
# =============================================================================

ggsave(
  
  filename = OUT_PDF,
  
  plot = p,
  
  width = FIG_WIDTH,
  
  height = FIG_HEIGHT,
  
  units = "in",
  
  device = grDevices::cairo_pdf,
  
  bg = "white"
  
)


# =============================================================================
# 16. 输出文件信息
# =============================================================================

cat("\n============================================================\n")
cat("Figure 5 successfully saved\n")
cat("============================================================\n")

cat("\nTIFF:\n")
cat(
  normalizePath(
    OUT_TIF,
    mustWork = FALSE
  ),
  "\n"
)

cat("\nPDF:\n")
cat(
  normalizePath(
    OUT_PDF,
    mustWork = FALSE
  ),
  "\n"
)

cat("\nFont:", FONT, "\n")
cat("X-axis range:", X_MIN, "to", X_MAX, "\n")
cat("Figure width:", FIG_WIDTH, "inches\n")
cat("Figure height:", FIG_HEIGHT, "inches\n")
cat("TIFF resolution: 600 dpi\n")
cat("PDF device: cairo_pdf\n\n")
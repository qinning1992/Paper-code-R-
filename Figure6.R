# =============================================================================
# Figure 6. Anatomical left/right hippocampal subfield-cognition associations
#
# Reproduces the composite main-text figure from the anatomical left/right
# fine-grained structure-cognition analysis (Section 3.4).
#
#   Panel A  Partial-correlation heatmap: 11 hippocampal subfields x left/right
#            hemisphere x 19 cognitive measures, with * for BH-FDR q < 0.05.
#   Panel B  Forest plot of the associations surviving correction (partial r
#            with Fisher-z 95% CI), grouped by cognitive measure and hemisphere.
#   Panel C  Two representative left-hemisphere verbal associations, shown as
#            residual-residual scatter plots.
#
# Inputs (relative to the repository root "code and data/"):
#   HIPP_THAL/ROI-Cognition Partial Correlation/Anatomical_LeftRight/FineGrained/
#       Step2_Anatomical_LeftRight_FineGrained_longtable_permBH.csv   (GBK)
#   df_raw_clean.xlsx
#
# The partial correlations are read from the deposited longtable produced by
# Structure_Cognition_EpileptogenicLaterality_integrated-202609.R; panel C
# residuals are recomputed from df_raw_clean.xlsx with the same covariate set.
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(readxl); library(dplyr); library(tidyr); library(cowplot)
})

## ---- paths -----------------------------------------------------------------
## Point BASE at the repository "code and data" folder.
BASE      <- "."
LONGTABLE <- file.path(BASE, "HIPP_THAL", "ROI-Cognition Partial Correlation",
                       "Anatomical_LeftRight", "FineGrained",
                       "Step2_Anatomical_LeftRight_FineGrained_longtable_permBH.csv")
DATAFILE  <- file.path(BASE, "df_raw_clean.xlsx")
OUT_TIF   <- file.path(BASE, "Figure6_Anatomical_left_right_structure_cognition.tif")
OUT_PDF   <- file.path(BASE, "Figure6_Anatomical_left_right_structure_cognition.pdf")

## ---- font: Arial, bold throughout ------------------------------------------
FONT <- "Arial"
if (.Platform$OS.type == "windows") {
  suppressWarnings(grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial")))
}
theme_fig <- function(base = 9) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      text        = element_text(family = FONT, face = "bold", colour = "black"),
      axis.text   = element_text(family = FONT, face = "bold", colour = "black"),
      axis.title  = element_text(family = FONT, face = "bold"),
      plot.title  = element_text(family = FONT, face = "bold"),
      strip.text  = element_text(family = FONT, face = "bold"),
      legend.text = element_text(family = FONT, face = "bold"),
      legend.title= element_text(family = FONT, face = "bold"),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )
}
TAG <- function(x) theme(plot.title = element_text(face = "bold", size = 13, hjust = 0))

## ---- read correlation table ------------------------------------------------
enc <- if (.Platform$OS.type == "windows") "GBK" else "GBK"
lt  <- read.csv(LONGTABLE, fileEncoding = enc, stringsAsFactors = FALSE,
                check.names = FALSE)
lt  <- lt[lt$Family == "HippocampusSubfields", ]
lt$r      <- suppressWarnings(as.numeric(lt$r))
lt$q_perm <- suppressWarnings(as.numeric(lt$q_perm))
lt$n      <- suppressWarnings(as.numeric(lt$n))

## subfield display order (head -> tail structures grouped sensibly)
sub_levels <- c("Hippocampal tail","Subiculum","Presubiculum","Parasubiculum",
                "CA1","CA2/3","CA4","GC-ML-DG","Molecular layer HP","HATA","Fimbria")
strip_hemi <- function(x, h) trimws(sub(paste0("^", h, "\\s*"), "", x))
lt$Sub <- mapply(strip_hemi, lt$ROI_Display, lt$Hemisphere)

## cognitive-measure display order and domain grouping
scale_levels <- c("WAIS VCI","BNT Naming","VA Immediate","VA Delay","VTrials I-V Total","CAVLT Slope",
                  "WAIS WMI","WAIS FSIQ","WAIS PSI",
                  "WMS FSMQ","WMS AMI","WMS IMI","WMS DMI","WMS VMI",
                  "WAIS PRI","FA Immediate","FA Delay","FTrials I-V Total","AFLT Slope")
## the deposited file uses an en-dash in some labels; normalise
lt$Scale_Display <- gsub("–", "-", lt$Scale_Display)
scale_levels     <- gsub("–", "-", scale_levels)
dom_map <- c(setNames(rep("Verbal",6), scale_levels[1:6]),
             setNames(rep("Working memory / global",3), scale_levels[7:9]),
             setNames(rep("Memory indices",5), scale_levels[10:14]),
             setNames(rep("Visuospatial",5), scale_levels[15:19]))
lt <- lt[lt$Scale_Display %in% scale_levels, ]
lt$Domain <- factor(dom_map[lt$Scale_Display],
                    levels = c("Verbal","Working memory / global","Memory indices","Visuospatial"))
lt$Scale_Display <- factor(lt$Scale_Display, levels = scale_levels)
lt$Sub  <- factor(lt$Sub,  levels = rev(sub_levels))
lt$Hemisphere <- factor(lt$Hemisphere, levels = c("Left","Right"))
lt$star <- ifelse(is.finite(lt$q_perm) & lt$q_perm < 0.05, "*", "")

## ---- Panel A: heatmap ------------------------------------------------------
pA <- ggplot(lt, aes(Scale_Display, Sub, fill = r)) +
  geom_tile(colour = "grey85", linewidth = 0.15) +
  geom_text(aes(label = star), family = FONT, fontface = "bold",
            size = 4.2, vjust = 0.78) +
  facet_grid(Hemisphere ~ Domain, scales = "free", space = "free", switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-0.45, 0.45),
                       oob = scales::squish, name = "partial r") +
  labs(x = NULL, y = NULL, title = "A") +
  theme_fig(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        panel.spacing = unit(2, "pt"),
        legend.position = "right", legend.key.height = unit(22, "pt"),
        legend.key.width = unit(9, "pt")) +
  TAG()

## ---- Panel B: forest of corrected associations -----------------------------
sig <- lt[is.finite(lt$q_perm) & lt$q_perm < 0.05, ]
## Fisher-z 95% CI for a partial correlation; k = number of covariates
k_cov <- ifelse(grepl("Age", sig$CovarsUsed), 6L, 5L)
z  <- atanh(sig$r); se <- 1 / sqrt(sig$n - k_cov - 3)
sig$lo <- tanh(z - 1.96 * se); sig$hi <- tanh(z + 1.96 * se)
sig$Measure <- factor(as.character(sig$Scale_Display),
                      levels = c("WAIS VCI","BNT Naming","WAIS WMI"))
sig <- sig[order(sig$Measure, sig$Hemisphere, sig$r), ]
sig$row <- factor(seq_len(nrow(sig)), levels = seq_len(nrow(sig)),
                  labels = paste0(sig$Hemisphere, " ", sig$Sub))
pB <- ggplot(sig, aes(r, row, colour = Hemisphere)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.7) +
  geom_point(size = 2.4) +
  facet_grid(Measure ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c(Left = "#B2182B", Right = "#2166AC"),
                      name = "Hemisphere") +
  scale_x_continuous(limits = c(-0.05, 0.75), breaks = c(0, 0.2, 0.4, 0.6)) +
  labs(x = "Partial r (95% CI)", y = NULL, title = "B") +
  theme_fig(9) +
  theme(axis.text.y = element_text(size = 7.5),
        strip.placement = "outside",
        legend.position = c(0.82, 0.12),
        panel.grid.major.x = element_line(colour = "grey93", linewidth = 0.3)) +
  TAG()

## ---- Panel C: representative scatters --------------------------------------
df <- as.data.frame(read_excel(DATAFILE))
df$side_lr <- as.numeric(df$side == 1)
resid_of <- function(y, X) {
  ok <- stats::complete.cases(cbind(y, X))
  b  <- stats::lm.fit(cbind(1, as.matrix(X[ok, , drop = FALSE])), y[ok])$coefficients
  out <- rep(NA_real_, length(y))
  out[ok] <- y[ok] - cbind(1, as.matrix(X[ok, , drop = FALSE])) %*% b
  out
}
mk_scatter <- function(roi, scale, covars, xlab, ylab, rr, qq) {
  d <- df[, c(roi, scale, covars)]
  for (cc in c(roi, scale, covars)) d[[cc]] <- suppressWarnings(as.numeric(d[[cc]]))
  d <- d[stats::complete.cases(d), ]
  xr <- resid_of(d[[roi]],   d[covars])
  yr <- resid_of(d[[scale]], d[covars])
  dd <- data.frame(x = xr, y = yr)
  ggplot(dd, aes(x, y)) +
    geom_point(size = 1.9, colour = "#B2182B", alpha = 0.75) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.8, formula = y ~ x) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
             label = sprintf("partial r = %.3f, q = %.3f, n = %d", rr, qq, nrow(dd)),
             family = FONT, fontface = "bold", size = 2.9, label.size = 0.2) +
    labs(x = xlab, y = ylab) + theme_fig(9)
}
C1 <- c("sex","Duration_Dis","Edu_year","eTIV","side_lr")
pC1 <- mk_scatter("lh_CA1_comb", "言语理解指数VCI", C1,
                  "Left CA1 volume (residual)", "WAIS VCI (residual)", 0.443, 0.048)
pC2 <- mk_scatter("lh_molecular_layer_HP_comb", "BNT-命名", c(C1,"Age"),
                  "Left molecular layer volume (residual)", "BNT Naming (residual)", 0.340, 0.028)
pC <- plot_grid(pC1, pC2, ncol = 1, labels = c("C", ""),
                label_fontfamily = FONT, label_fontface = "bold", label_size = 13,
                hjust = 0, label_x = 0.01)

## ---- assemble --------------------------------------------------------------
bottom <- plot_grid(pB, pC, ncol = 2, rel_widths = c(1.05, 1))
fig <- plot_grid(pA, bottom, ncol = 1, rel_heights = c(1.28, 1))

ggsave(OUT_TIF, fig, width = 13.5, height = 11.2, units = "in", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
ggsave(OUT_PDF, fig, width = 13.5, height = 11.2, units = "in", bg = "white")
cat("Saved:\n  ", normalizePath(OUT_TIF), "\n  ", normalizePath(OUT_PDF), "\n")

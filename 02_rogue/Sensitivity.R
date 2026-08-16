## ============================================================
## Cell-number robustness analysis for ROGUE (genotype comparison)
##
## Full, consolidated version -- combines everything run across this
## investigation (previously scattered across console commands) into
## one reproducible script:
##
##   1. Correlation checks (aggregate, within-DCM, within-control,
##      with an outlier-sensitivity check on the within-control result)
##   2. ANCOVA sensitivity checks (5-level genotype, and binary
##      control-vs-DCM), for context against the downsampling result
##   3. Downsampling-based re-analysis (the definitive test): recompute
##      ROGUE at equal cell counts per donor, repeat the primary
##      comparisons (control vs DCM; genotype independence; pairwise
##      vs control)
##
## Re-creates pathway_expr / pathway_meta / rogue.pathway from the
## saved files if they are not already in memory, so this script can
## be run standalone in a fresh session.
## ============================================================

library(ROGUE)
library(Matrix)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

set.seed(1)

data_dir    <- "/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready"
results_dir <- "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor"

## ------------------------------------------------------------
## SETUP: load / recreate objects if not already in memory
## ------------------------------------------------------------

if (!exists("pathway_expr") || !exists("pathway_meta")) {

  pathway_expr <- readMM(file.path(data_dir, "expr_lv.mtx"))
  pathway_cells <- read.csv(file.path(data_dir, "cells_lv.csv"))
  pathway_genes <- read.csv(file.path(data_dir, "genes_lv.csv"))
  pathway_meta  <- read.csv(file.path(data_dir, "meta_lv.csv"))

  rownames(pathway_expr) <- pathway_genes[[1]]
  colnames(pathway_expr) <- pathway_cells[[1]]
  pathway_expr <- as.matrix(pathway_expr)

  pathway_expr <- matr.filter(pathway_expr, min.cells = 50, min.genes = 10)
  pathway_meta <- pathway_meta[pathway_meta$X %in% colnames(pathway_expr), ]
}

if (!exists("rogue.genotype")) {
  rogue.pathway <- readRDS(file.path(results_dir, "rogue_genotype.rds"))
}

## reshape rogue.pathway (donors x genotypes, NA off-diagonal) into
## long format: one row per donor with its own genotype and score
rp <- rogue.pathway %>%
  rownames_to_column("donor") %>%
  pivot_longer(-donor, names_to = "pathway", values_to = "rogue") %>%
  filter(!is.na(rogue))

cells_per_donor <- as.data.frame(table(pathway_meta$donor_id))
names(cells_per_donor) <- c("donor", "n_cells")

rp2 <- merge(rp, cells_per_donor, by = "donor")
write.csv(rp2, file.path(results_dir, "rogue_vs_ncells.csv"), row.names = FALSE)

## ============================================================
## PART 1: CORRELATION CHECKS
## ============================================================

cat("\n============================================================\n")
cat("CHECK 1a: correlation, all donors\n")
cat("============================================================\n")
cor_all <- cor.test(rp2$rogue, rp2$n_cells, method = "spearman")
print(cor_all)

p_corr <- ggplot(rp2, aes(n_cells, rogue, colour = pathway)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40") +
  scale_x_log10() +
  labs(x = "Cells per donor (log10)", y = "ROGUE",
       title = paste0("Spearman rho = ", round(cor_all$estimate, 3),
                      ", p = ", signif(cor_all$p.value, 3))) +
  theme_classic(base_size = 12)
ggsave(file.path(results_dir, "rogue_vs_ncells.png"), p_corr, width = 6, height = 4.5, dpi = 300)

cat("\n============================================================\n")
cat("CHECK 1b: correlation, DCM donors only\n")
cat("============================================================\n")
rp2_dcm <- subset(rp2, pathway != "control")
print(cor.test(rp2_dcm$rogue, rp2_dcm$n_cells, method = "spearman"))

cat("\n============================================================\n")
cat("CHECK 1c: correlation, control donors only\n")
cat("============================================================\n")
rp2_ctrl <- subset(rp2, pathway == "control")
print(cor.test(rp2_ctrl$rogue, rp2_ctrl$n_cells, method = "spearman"))

cat("\n--- control donors sorted by cell count (check for outliers) ---\n")
print(rp2_ctrl[order(rp2_ctrl$n_cells), c("donor", "rogue", "n_cells")])

cat("\n--- CHECK 1c, sensitivity: repeat excluding the most extreme donor ---\n")
## identify the donor with lowest n_cells AND lowest rogue among controls
outlier_donor <- rp2_ctrl$donor[which.min(rp2_ctrl$n_cells)]
cat("excluding:", outlier_donor, "\n")
rp2_ctrl_noOutlier <- subset(rp2_ctrl, donor != outlier_donor)
print(cor.test(rp2_ctrl_noOutlier$rogue, rp2_ctrl_noOutlier$n_cells, method = "spearman"))

## ============================================================
## PART 2: ANCOVA SENSITIVITY CHECKS
## (context only -- Part 3's downsampling is the definitive test)
## ============================================================

rp2$log_n_cells <- log(rp2$n_cells)

cat("\n============================================================\n")
cat("ANCOVA: rogue ~ pathway (5-level) + log(n_cells)\n")
cat("============================================================\n")
fit_full <- lm(rogue ~ pathway + log_n_cells, data = rp2)
print(summary(fit_full))
fit_null <- lm(rogue ~ log_n_cells, data = rp2)
print(anova(fit_null, fit_full))

cat("\n============================================================\n")
cat("ANCOVA: rogue ~ group2 (control vs DCM, binary) + log(n_cells)\n")
cat("============================================================\n")
rp2$group2 <- ifelse(rp2$pathway == "control", "control", "DCM")
fit_bin_full <- lm(rogue ~ group2 + log_n_cells, data = rp2)
print(summary(fit_bin_full))
fit_bin_null <- lm(rogue ~ log_n_cells, data = rp2)
print(anova(fit_bin_null, fit_bin_full))

## ============================================================
## PART 3: DOWNSAMPLING-BASED ROBUSTNESS RE-ANALYSIS (definitive test)
## ============================================================

cat("\n============================================================\n")
cat("CHECK 2: downsampling\n")
cat("============================================================\n")

donor_counts <- as.integer(table(pathway_meta$donor_id))
print(summary(donor_counts))

## the raw minimum (60) is a severe outlier and would force every
## donor down to 60 cells, adding enough noise to obscure a real
## effect regardless of any confound. Instead, exclude donors below
## a defensible threshold from this sensitivity analysis only (not
## from the main ROGUE analysis), and downsample everyone else to
## that threshold. threshold = 300 excludes 4/50 donors (2 PVneg,
## 1 chromatin, 1 splicing) and leaves every genotype with >= 6
## donors -- report this exclusion explicitly in methods.
threshold <- 300
target_n  <- threshold
n_reps    <- 20

donor_ncells    <- table(pathway_meta$donor_id)
excluded_donors <- names(donor_ncells)[donor_ncells < threshold]
donors          <- setdiff(unique(as.character(pathway_meta$donor_id)), excluded_donors)

cat("donors excluded from downsampling check (< ", threshold, " cells): ",
    paste(excluded_donors, collapse = ", "), "\n", sep = "")
cat("donors retained:", length(donors), "of",
    length(unique(pathway_meta$donor_id)), "\n")

## one-time diagnostic — if any of these look wrong, that's the bug
cat("class of pathway_meta$X:", class(pathway_meta$X), "\n")
cat("class of colnames(pathway_expr):", class(colnames(pathway_expr)), "\n")
cat("n pathway_meta$X not found in colnames(pathway_expr):",
    sum(!as.character(pathway_meta$X) %in% colnames(pathway_expr)), "\n")
cat("any duplicated colnames in pathway_expr:",
    any(duplicated(colnames(pathway_expr))), "\n")

reps_list <- vector("list", n_reps)

for (i in seq_len(n_reps)) {

  keep_cells <- unlist(lapply(donors, function(d) {
    cells_d <- as.character(pathway_meta$X[as.character(pathway_meta$donor_id) == d])
    cells_d <- intersect(cells_d, colnames(pathway_expr))
    if (length(cells_d) == 0) return(character(0))
    if (length(cells_d) >= target_n) {
      sample(cells_d, target_n)
    } else {
      cells_d
    }
  }))

  keep_cells <- as.character(keep_cells)
  stopifnot(all(keep_cells %in% colnames(pathway_expr)))

  expr_sub <- pathway_expr[, keep_cells, drop = FALSE]
  meta_sub <- pathway_meta[match(keep_cells, as.character(pathway_meta$X)), ]

  rogue_rep <- rogue(expr_sub,
                     labels   = meta_sub$Genotype,
                     samples  = meta_sub$donor_id,
                     platform = "UMI",
                     filter   = FALSE,
                     span     = 0.75)

  reps_list[[i]] <- rogue_rep

  cat("rep", i, "of", n_reps, "done\n")
}

## reps_list[[i]] is a donor x genotype-style table each rep (same
## structure rogue() normally returns) — average elementwise across
## reps, ignoring NAs (a donor can be NA in a rep if it has no cells
## for a given genotype column)

sum_mat   <- Reduce(`+`, lapply(reps_list, function(x) { x[is.na(x)] <- 0; x }))
count_mat <- Reduce(`+`, lapply(reps_list, function(x) !is.na(x)))

rogue_downsampled <- sum_mat / count_mat

saveRDS(rogue_downsampled, file.path(results_dir, "rogue_genotype_downsampled.rds"))
write.csv(rogue_downsampled, file.path(results_dir, "rogue_genotype_downsampled.csv"))

## ------------------------------------------------------------
## compare downsampled result to the original (full-data) result
## ------------------------------------------------------------

rp_ds <- as.data.frame(rogue_downsampled) %>%
  rownames_to_column("donor") %>%
  pivot_longer(-donor, names_to = "pathway", values_to = "rogue") %>%
  filter(!is.na(rogue))

cat("\n--- downsampled: Kruskal-Wallis across all 5 groups (incl. control) ---\n")
print(kruskal.test(rogue ~ pathway, data = rp_ds))

cat("\n--- downsampled: control vs DCM (binary) ---\n")
rp_ds$group2 <- ifelse(rp_ds$pathway == "control", "control", "DCM")
print(wilcox.test(rogue ~ group2, data = rp_ds))

cat("\n--- downsampled: pairwise vs control, BH-corrected ---\n")
dcm_groups <- setdiff(unique(rp_ds$pathway), "control")
pvals <- sapply(dcm_groups, function(g) {
  d <- droplevels(subset(rp_ds, pathway %in% c("control", g)))
  wilcox.test(rogue ~ pathway, data = d)$p.value
})
print(data.frame(genotype = dcm_groups,
                 p = pvals,
                 p_adj = p.adjust(pvals, method = "BH")))

cat("\n--- downsampled: genotype independence, DCM only (matches original claim) ---\n")
rp_ds_dcm <- subset(rp_ds, pathway != "control")
print(kruskal.test(rogue ~ droplevels(factor(pathway)), data = rp_ds_dcm))

## ------------------------------------------------------------
## write every test result above to a single summary CSV
## ------------------------------------------------------------

cor_dcm    <- cor.test(rp2_dcm$rogue, rp2_dcm$n_cells, method = "spearman")
cor_ctrl   <- cor.test(rp2_ctrl$rogue, rp2_ctrl$n_cells, method = "spearman")
cor_ctrl_no_outlier <- cor.test(rp2_ctrl_noOutlier$rogue, rp2_ctrl_noOutlier$n_cells, method = "spearman")
kw_dcm_only <- kruskal.test(rogue ~ droplevels(factor(pathway)), data = rp_ds_dcm)
kw_all5     <- kruskal.test(rogue ~ pathway, data = rp_ds)
wilcox_bin  <- wilcox.test(rogue ~ group2, data = rp_ds)

## fixed (one-row-per-test) results
summary_fixed <- data.frame(
  test = c(
    "Check1a: correlation (rho), all donors",
    "Check1b: correlation (rho), DCM only",
    "Check1c: correlation (rho), control only",
    "Check1c: correlation (rho), control only (outlier excluded)",
    "ANCOVA: control-vs-DCM coefficient (adj. for log n_cells)",
    "Downsampled: control vs DCM (Wilcoxon W)",
    "Downsampled: genotype independence, DCM only (Kruskal-Wallis chi-sq)",
    "Downsampled: all 5 groups incl. control (Kruskal-Wallis chi-sq)"
  ),
  statistic = c(
    unname(cor_all$estimate),
    unname(cor_dcm$estimate),
    unname(cor_ctrl$estimate),
    unname(cor_ctrl_no_outlier$estimate),
    unname(coef(fit_bin_full)["group2DCM"]),
    unname(wilcox_bin$statistic),
    unname(kw_dcm_only$statistic),
    unname(kw_all5$statistic)
  ),
  p_value = c(
    cor_all$p.value,
    cor_dcm$p.value,
    cor_ctrl$p.value,
    cor_ctrl_no_outlier$p.value,
    summary(fit_bin_full)$coefficients["group2DCM", "Pr(>|t|)"],
    wilcox_bin$p.value,
    kw_dcm_only$p.value,
    kw_all5$p.value
  ),
  p_adj = NA
)

## one row per DCM genotype for the pairwise-vs-control result
summary_pairwise <- data.frame(
  test      = paste0("Downsampled: pairwise vs control (Wilcoxon) - ", dcm_groups),
  statistic = NA,
  p_value   = pvals,
  p_adj     = p.adjust(pvals, method = "BH")
)

robustness_summary <- rbind(summary_fixed, summary_pairwise)

write.csv(robustness_summary,
          file.path(results_dir, "cell_number_robustness_summary.csv"),
          row.names = FALSE)

cat("\nsaved full summary table to:",
    file.path(results_dir, "cell_number_robustness_summary.csv"), "\n")

## ============================================================
## SUMMARY -- compare these numbers directly against the originals:
##
##   Control vs DCM (Wilcoxon):
##     original    W = 99,  p = 0.0027
##     downsampled W = ?,   p = ?   (printed above)
##
##   Genotype independence, DCM only (Kruskal-Wallis):
##     original    chi-sq(3) = 2.71, p = 0.44
##     downsampled chi-sq(3) = ?,    p = ?   (printed above)
##
## If both hold up, this is your robustness evidence for methods/
## results -- report as a supplementary sensitivity-analysis table
## alongside rp2 (Part 1) and rogue_downsampled (Part 3).
## ============================================================

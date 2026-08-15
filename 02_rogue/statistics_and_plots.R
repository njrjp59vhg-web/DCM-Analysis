# ============================================================
#  ROGUE — per-donor cardiomyocyte transcriptional uniformity
#  Left ventricle, v3
#  Statistics + figures  (two-sided Wilcoxon throughout)
# ============================================================

# ------------------------------------------------------------
# 0.  Libraries
# ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(patchwork)

# ------------------------------------------------------------
# 1.  Load & reshape data
# ------------------------------------------------------------
rogue.pathway <- readRDS(
  "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_genotype.rds")

rp <- rogue.pathway |>
  as.data.frame() |>
  rownames_to_column("donor") |>
  pivot_longer(-donor, names_to = "pathway", values_to = "rogue") |>
  filter(!is.na(rogue))

lev <- c("control", "splicing", "sarcomere", "PVneg", "chromatin")
rp$pathway <- factor(rp$pathway, levels = lev)

# disease grouping (used in panels A / density)
rp$disease <- factor(
  ifelse(rp$pathway == "control", "normal", "dilated cardiomyopathy"),
  levels = c("dilated cardiomyopathy", "normal"))

cat("N donors per pathway:\n"); print(table(rp$pathway))

# ------------------------------------------------------------
# 2.  Descriptive statistics  (median / IQR = the correct
#     descriptors for a rank-based test)
# ------------------------------------------------------------
desc <- rp |>
  group_by(pathway) |>
  summarise(n = n(),
            median = median(rogue),
            IQR    = IQR(rogue),
            mean   = mean(rogue),
            sd     = sd(rogue),
            .groups = "drop")
cat("\nDescriptive statistics:\n"); print(desc)

# ------------------------------------------------------------
# 3.  Omnibus test — Kruskal-Wallis across all groups
# ------------------------------------------------------------
kw <- kruskal.test(rogue ~ pathway, data = rp)
cat("\nKruskal-Wallis omnibus:\n"); print(kw)

# ------------------------------------------------------------
# 4.  Post-hoc — each genotype vs control
#     TWO-SIDED Wilcoxon, BH-adjusted, with rank-biserial effect size
# ------------------------------------------------------------
dcm <- setdiff(lev, "control")

# helper: two-sided Wilcoxon + rank-biserial r for one genotype vs control
# r_rb sign is oriented so that NEGATIVE = genotype lower than control
wilcox_vs_control <- function(g) {
  d  <- droplevels(subset(rp, pathway %in% c("control", g)))
  w  <- wilcox.test(rogue ~ pathway, data = d, alternative = "two.sided")
  n1 <- sum(d$pathway == "control")
  n2 <- sum(d$pathway == g)
  U  <- as.numeric(w$statistic)          # W for the 'control' group
  r_rb <- (2 * U / (n1 * n2)) - 1         # >0: control higher (genotype lower)
  r_rb <- -r_rb                           # flip: negative = genotype lower
  data.frame(pathway = g, W = U, p = w$p.value, r_rb = r_rb)
}

vs  <- do.call(rbind, lapply(dcm, wilcox_vs_control))
vs$padj <- p.adjust(vs$p, method = "BH")

# assemble full results table (control row = reference, no test)
res <- data.frame(pathway = "control", W = NA, p = NA, padj = NA, r_rb = NA) |>
  rbind(vs[, c("pathway", "W", "p", "padj", "r_rb")])
res$median <- tapply(rp$rogue, rp$pathway, median)[res$pathway]
res$IQR    <- tapply(rp$rogue, rp$pathway, IQR)[res$pathway]
res$mean   <- tapply(rp$rogue, rp$pathway, mean)[res$pathway]

cat("\nPost-hoc results (two-sided Wilcoxon vs control, BH-adjusted):\n")
print(res)

write.csv(res,
          "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/tables/rogue_per_group.csv",
          row.names = FALSE)

# significance stars from adjusted p
res$star <- cut(res$padj, c(-Inf, 0.001, 0.01, 0.05, Inf),
                labels = c("***", "**", "*", "ns"))

# ------------------------------------------------------------
# 5.  Figure 1 — boxplot with BH-adjusted significance stars
# ------------------------------------------------------------
lab <- data.frame(pathway = factor(res$pathway, levels = lev),
                  star = as.character(res$star), y = 0.95)

p_box <- ggplot(rp, aes(pathway, rogue, fill = pathway)) +
  geom_boxplot(width = .55, outlier.shape = NA, alpha = .55) +
  geom_jitter(width = .12, height = 0, size = 2, alpha = .7) +
  geom_text(data = lab, aes(x = pathway, y = y, label = star),
            inherit.aes = FALSE, size = 6) +
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("Control","Splicing\n(RBM20)","Sarcomere",
                              "PVneg","Chromatin\n(LMNA)")) +
  labs(x = NULL, y = "ROGUE  (higher = more uniform)",
       title = "Per-donor cardiomyocyte transcriptional uniformity — LV, v3",
       subtitle = "Each point = one donor;  * vs control (two-sided Wilcoxon, BH-adjusted)") +
  theme_classic(base_size = 13) + theme(legend.position = "none")

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Plots/rogue_lv_v3.png",
       p_box, width = 8, height = 5.5, dpi = 300)
print(p_box)

# ------------------------------------------------------------
# 6.  Figure 2 — violin panels
#     A: DCM vs normal (overall)   B: each genotype vs control
#     Both TWO-SIDED for consistency
# ------------------------------------------------------------

## Panel A — DCM vs normal, dots coloured by genotype
pA <- ggplot(rp, aes(disease, rogue)) +
  geom_violin(aes(fill = disease), alpha = .5, colour = "grey20",
              trim = FALSE, show.legend = FALSE) +
  geom_jitter(aes(colour = pathway), width = .12, height = 0, size = 2.4, alpha = .9) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "black") +
  stat_compare_means(method = "wilcox.test", label.x = 1.3, label.y = 1.0) +
  scale_fill_manual(values = c("dilated cardiomyopathy" = "#D55E00",
                               "normal" = "#1B9E77")) +
  scale_colour_brewer(palette = "Dark2", name = "Genotype",
                      labels = c(control = "Control", splicing = "Splicing (RBM20)",
                                 sarcomere = "Sarcomere", PVneg = "PVneg",
                                 chromatin = "Chromatin (LMNA)")) +
  coord_cartesian(ylim = c(NA, 1.05)) +
  labs(x = NULL, y = "ROGUE value  (higher = more uniform)",
       title = "DCM vs normal (all donors)") +
  theme_classic(base_size = 13)

## Panel B — each genotype vs control, BH-adjusted (two-sided)
stat_df <- data.frame(
  group1 = "control",
  group2 = res$pathway[!is.na(res$padj)],
  p.adj  = signif(res$padj[!is.na(res$padj)], 2),
  y.position = seq(1.00, 1.18, length.out = sum(!is.na(res$padj))))

pB <- ggplot(rp, aes(pathway, rogue)) +
  geom_violin(aes(fill = pathway), scale = "width", alpha = .5,
              colour = "grey20", trim = FALSE, show.legend = FALSE) +
  geom_jitter(aes(colour = pathway), width = .12, height = 0, size = 1.9, alpha = .9) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "black") +
  stat_pvalue_manual(stat_df, label = "p.adj = {p.adj}", tip.length = .01, size = 3) +
  scale_fill_brewer(palette = "Set2") +
  scale_colour_brewer(palette = "Dark2", guide = "none") +
  scale_x_discrete(labels = c("Control","Splicing\n(RBM20)","Sarcomere",
                              "PVneg","Chromatin\n(LMNA)")) +
  coord_cartesian(ylim = c(NA, 1.22)) +
  labs(x = NULL, y = NULL, title = "By DCM genotype vs control") +
  theme_classic(base_size = 13)

## Combine (shared genotype legend)
fig <- pA + pB + plot_layout(widths = c(1, 2), guides = "collect") +
  plot_annotation(tag_levels = "A",
                  title = "Per-donor cardiomyocyte transcriptional uniformity (LV, v3)",
                  subtitle = "◆ = group mean;  A: two-sided Wilcoxon;  B: two-sided Wilcoxon, BH-adjusted")

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Plots/rogue_Violin.png",
       fig, width = 13, height = 6, dpi = 300)
print(fig)

# ------------------------------------------------------------
# 7.  Figure 3 — density: DCM vs normal
# ------------------------------------------------------------
pdens <- ggplot(rp, aes(rogue, fill = disease, colour = disease)) +
  geom_density(alpha = .4, linewidth = .9) +
  geom_rug(alpha = .6, linewidth = .4) +
  scale_fill_manual(values = c("dilated cardiomyopathy" = "#D55E00",
                               "normal" = "#1B9E77"), name = NULL) +
  scale_colour_manual(values = c("dilated cardiomyopathy" = "#D55E00",
                                 "normal" = "#1B9E77"), name = NULL) +
  labs(x = "ROGUE value  (higher = more uniform)", y = "Density",
       title = "Per-donor ROGUE distribution: DCM vs normal",
       subtitle = "DCM shifted lower = less transcriptional uniformity") +
  theme_classic(base_size = 13)

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Plots/rogue_density.png",
       pdens, width = 8, height = 5, dpi = 300)

# ------------------------------------------------------------
# 8.  DCM vs normal — summary + two-sided test + effect size
# ------------------------------------------------------------
summary_comb <- rp |>
  group_by(disease) |>
  summarise(n_donors = n(),
            median = median(rogue),
            IQR    = IQR(rogue),
            mean   = mean(rogue),
            sd     = sd(rogue),
            .groups = "drop")
cat("\nDCM vs normal summary:\n"); print(summary_comb)

# two-sided Wilcoxon (matches Panel A)
w_two <- wilcox.test(rogue ~ disease, data = rp, alternative = "two.sided")

# effect size (rank-biserial + Z/sqrt(N)) read from the test object
n1 <- sum(rp$disease == "dilated cardiomyopathy")
n2 <- sum(rp$disease == "normal")
U   <- as.numeric(w_two$statistic)          # from object, not hard-coded
r_rb <- 1 - 2 * U / (n1 * n2)               # rank-biserial correlation
muU  <- n1 * n2 / 2
sdU  <- sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
Z    <- (U - muU) / sdU
r_z  <- abs(Z) / sqrt(n1 + n2)              # Z/sqrt(N), as rstatix reports

#Test for any signficant differences within DCM genotypes

# DCM genotypes only
dcm_only <- droplevels(subset(rp, pathway != "control"))

# omnibus: do the DCM genotypes differ from each other at all?
kruskal.test(rogue ~ pathway, data = dcm_only)

# all pairwise two-sided Wilcoxon among DCM genotypes, BH-adjusted
pairwise.wilcox.test(dcm_only$rogue, dcm_only$pathway,
                     p.adjust.method = "BH")
ctrl_med <- median(rp$rogue[rp$pathway == "control"])  

withinp <- ggplot(dcm_only, aes(pathway, rogue, fill = pathway)) +
  geom_hline(yintercept = ctrl_med, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 0.7, y = ctrl_med + 0.006, label = "control median",
           colour = "grey40", hjust = 0, size = 3) +
  geom_boxplot(width = .55, outlier.shape = NA, alpha = .55) +
  geom_jitter(width = .12, height = 0, size = 2, alpha = .7) +
  stat_compare_means(method = "kruskal.test", label.y = 0.9) +   # prints p = 0.44
  scale_fill_brewer(palette = "Set2") +
  scale_x_discrete(labels = c("Splicing\n(RBM20)","Sarcomere\n(TTN)",
                              "PVneg","Chromatin\n(LMNA)")) +
  labs(x = NULL, y = "ROGUE  (higher = more uniform)",
       title = "Transcriptional uniformity is uniformly reduced across DCM genotypes",
       subtitle = "All genotypes sit below control; Kruskal–Wallis n.s.") +
  theme_classic(base_size = 13) + theme(legend.position = "none")

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Plots/rogue_within_dcm.png",
       withinp, width = 10, height = 5.5, dpi = 300)

library(ROGUE)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(tibble)
library(org.Hs.eg.db)
library(openxlsx)

norm_states <- readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_nor.rds")
path_states <- readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_dcm.rds")

pathstates <- rogue.boxplot(path_states) + ggtitle("A: Dilated Cardiomyopathy")
normstates <- rogue.boxplot(norm_states) + ggtitle("B: Normal")

statesfig <- (pathstates + normstates) &
  coord_cartesian(ylim = c(0.58, 1.0)) &
  theme(plot.title = element_text(face = "bold"))

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Plots/states_rogue_plot.png",statesfig, width = 15, height = 5.5, dpi = 300)

norm_med <- apply(norm_states, 2, median, na.rm = TRUE)
dcm_med  <- apply(path_states, 2, median, na.rm = TRUE)
norm_n <- colSums(!is.na(norm_states))
dcm_n  <- colSums(!is.na(path_states))

states <- intersect(names(norm_med), names(dcm_med))
drop <- tibble(
  state  = states,
  normal = norm_med[states],
  DCM    = dcm_med[states],
  n_norm = norm_n[states],
  n_dcm  = dcm_n[states]
) |>
  mutate(drop = normal - DCM) |>
  arrange(desc(drop))
print(drop)
head(drop, 3)

per_state_test <- lapply(states, function(s){
  x <- norm_states[[s]]; y <- path_states[[s]]
  x <- x[!is.na(x)];     y <- y[!is.na(y)]
  w <- wilcox.test(x, y, alternative = "two.sided")  
  data.frame(state = s, p = w$p.value)
}) |> bind_rows()
per_state_test$padj <- p.adjust(per_state_test$p, "BH")
per_state_test |> arrange(padj)


seu <- readRDS("/rds/general/user/ao225/home/CardiaFinal/Data/Post_Processed/Science/lv_Cardio.rds")

Idents(seu) <- seu$cell_states

add_sym <- function(mk) {
  ids <- sub("\\..*$", "", rownames(mk))          # drop version if present
  mk$symbol <- mapIds(org.Hs.eg.db, keys = ids,
                      keytype = "ENSEMBL", column = "SYMBOL")
  head(mk[, c("symbol","avg_log2FC","pct.1","pct.2","p_val_adj")], 25)
}

mk_5  <- FindMarkers(seu, ident.1 = "vCM5",   only.pos = TRUE,
                     logfc.threshold = 0.25, min.pct = 0.25)
mk_11 <- FindMarkers(seu, ident.1 = "vCM1.1", only.pos = TRUE,
                     logfc.threshold = 0.25, min.pct = 0.25)
mk_13 <- FindMarkers(seu, ident.1 = "vCM1.3", only.pos = TRUE,
                     logfc.threshold = 0.25, min.pct = 0.25)

add_sym(mk_5)
add_sym(mk_11)
add_sym(mk_13)

write.csv(add_sym(mk_5), "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/tables/markers_vCM5.csv")  
write.csv(add_sym(mk_11), "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/tables/markers_vCM1.1.csv")  
write.csv(add_sym(mk_13), "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/tables/markers_vCM1.3.csv")  

sheets <- list("per_state_ROGUE" = per_state_test,
               "drop"  = drop)
write.xlsx(sheets,
           "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/tables/rogue_cellstate_supp.xlsx")

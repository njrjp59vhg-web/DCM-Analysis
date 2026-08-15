library(ROGUE)
library(tidyverse)
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(patchwork)

meta <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/meta_lv.csv")
expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/expr_lv.mtx")
cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/cells_lv.csv")
genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/Science/genes_lv.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

rownames(meta) <- meta$X
colnames(expr) <- meta$X

#Seurat creation & reduction
seu <- CreateSeuratObject(counts = expr, meta.data = meta)

seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
seu <- RunHarmony(seu, group.by.vars = c("donor_id"))
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)

saveRDS(seu, "/rds/general/user/ao225/home/CardiaFinal/Data/Post_Processed/Science/lv_Cardio.rds")

#CheckPoint
seu <- readRDS("/rds/general/user/ao225/home/CardiaFinal/Data/Post_Processed/Science/lv_Cardio.rds")

#Saving UMAP data frame
umap_df <- as.data.frame(seu@reductions$umap@cell.embeddings)
umap_df$cell_states <- seu@meta.data$cell_states
umap_df$disease <- seu@meta.data$disease
umap_df$donor_id <- seu@meta.data$donor_id
umap_df$sex <- seu@meta.data$sex
umap_df$cell <- rownames(seu@reductions$umap@cell.embeddings)


write.csv(umap_df, "/rds/general/user/ao225/home/CardiaFinal/Results/Science/lv_umap.csv", 
          row.names = FALSE)

#Checkpoint
umap_df <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Results/Science/lv_umap.csv")

# UMAP by cell states
p1 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = cell_states)) +
  geom_point(size = 0.1, alpha = 0.3) +
  theme_classic() +
  labs(title = "Cell states", colour = "Cell state") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

# UMAP by disease
p2 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = disease)) +
  geom_point(size = 0.1, alpha = 0.3) +
  scale_colour_manual(values = c("normal" = "#1D9E75", 
                                 "dilated cardiomyopathy" = "#D85A30")) +
  theme_classic() +
  labs(title = "Disease condition", colour = "Condition") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1)))

p12 <- p1 + p2

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/Science/Plots/lv_umap_plots.png", p12,
       width = 14, height = 6, dpi = 300)

#On Hold
rogue.cellstate <- readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/Science/rogue_states.rds")
rogue.cellstate.long <- rogue.cellstate %>%
  as.data.frame() %>%
  rownames_to_column("donor_id") %>%
  pivot_longer(-donor_id, names_to = "cell_state", values_to = "ROGUE") %>%
  filter(!is.na(ROGUE))

# create matching keys
umap_df$key <- paste(umap_df$donor_id, umap_df$cell_states, sep = "_")
rogue.cellstate.long$key <- paste(rogue.cellstate.long$donor_id, rogue.cellstate.long$cell_state, sep = "_")

# map ROGUE score onto cells
umap_df$ROGUE_score <- rogue.cellstate.long$ROGUE[match(umap_df$key, rogue.cellstate.long$key)]

# check
cat("NAs:", sum(is.na(umap_df$ROGUE_score)), "\n")
cat("Total cells:", nrow(umap_df), "\n")

p_umap_rogue <- ggplot(umap_df, aes(x = umap_1, y = umap_2, colour = ROGUE_score)) +
  geom_point(size = 0.1, alpha = 0.3) +
  scale_colour_gradient2(
    low = "#2166AC", mid = "#F7F7B0", high = "#B2182B",   # blue–pale yellow–red
    midpoint = median(umap_df$ROGUE_score, na.rm = TRUE),
    name = "ROGUE score", na.value = "grey85") +
  theme_classic() +
  labs(title = "Transcriptional entropy by cell state (ROGUE)")

p1 + p_umap_rogue
ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/Science/lv_umap_plots_rogue.png", width = 14, height = 6, dpi = 300)

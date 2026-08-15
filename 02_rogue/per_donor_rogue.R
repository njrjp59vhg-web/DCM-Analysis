library(ROGUE)
library(Matrix)
library(tibble)
library(dplyr)

#Loading Data
pathway_expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_lv.mtx")
pathway_cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_lv.csv")
pathway_genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_lv.csv")
pathway_meta <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/meta_lv.csv")

#Merging expresion data
rownames(pathway_expr) <- pathway_genes[[1]]
colnames(pathway_expr) <- pathway_cells[[1]]

pathway_expr <- as.matrix(pathway_expr)

#Filtering expression and metadata
pathway_expr <- matr.filter(pathway_expr, min.cells = 50, min.genes = 10)
pathway_meta <- pathway_meta[pathway_meta$X %in% colnames(pathway_expr), ]

#Per-donor ROGUE score by Genotype
rogue.pathway <- rogue(
  pathway_expr,
  labels = pathway_meta$Genotype,
  samples = pathway_meta$donor_id,
  platform = "UMI",
  filter = FALSE,
  span = 0.75
)

#Per-donor ROGUE score by cell_state
rogue.states <- rogue.states <- rogue(
  pathway_expr ,
  labels = pathway_meta$cell_states,
  samples = pathway_meta$donor_id,
  platform = "UMI",
  filter = FALSE,
  span = 0.75
)

#Saving the ROGUE scores for further analysis
saveRDS(rogue.pathway, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_genotype.rds")
saveRDS(rogue.states, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_states.rds")
#-----------------------------------------------

dcm_expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_dcm.mtx")
dcm_cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_dcm.csv")
dcm_genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_dcm.csv")
dcm_meta <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/meta_dcm.csv")

rownames(dcm_expr) <- dcm_genes[[1]]
colnames(dcm_expr) <- dcm_cells[[1]]

dcm_expr <- as.matrix(dcm_expr)
dcm_expr <- matr.filter(dcm_expr, min.cells = 50, min.genes = 10)
dcm_meta <- dcm_meta[dcm_meta$X %in% colnames(dcm_expr), ]

rogue.dcm <- rogue(
  dcm_expr ,
  labels = dcm_meta$cell_states,
  samples = dcm_meta$donor_id,
  platform = "UMI",
  filter = FALSE,
  span = 0.75
)

saveRDS(rogue.dcm, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_dcm.rds")
#-----------------------------------------------
nor_expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_nor.mtx")
nor_cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_nor.csv")
nor_genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_nor.csv")
nor_meta <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/meta_nor.csv")

rownames(nor_expr) <- nor_genes[[1]]
colnames(nor_expr) <- nor_cells[[1]]

nor_expr <- as.matrix(nor_expr)
nor_expr <- matr.filter(nor_expr, min.cells = 50, min.genes = 10)
nor_meta <- nor_meta[nor_meta$X %in% colnames(nor_expr), ]

rogue.nor <- rogue(
  nor_expr ,
  labels = nor_meta$cell_states,
  samples = nor_meta$donor_id,
  platform = "UMI",
  filter = FALSE,
  span = 0.75
)

saveRDS(rogue.nor, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/Per_Donor/rogue_nor.rds")

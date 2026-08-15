library(ROGUE)
library(Matrix)
library(tibble)
library(dplyr)

#control
expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_control.mtx")
cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_control.csv")
genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_control.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_control.rds")

#sarcomere
expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_sarcomere.mtx")
cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_sarcomere.csv")
genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_sarcomere.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_sarcomere.rds")

#chromatin
expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_chromatin.mtx")
cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_chromatin.csv")
genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_chromatin.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_chromatin.rds")

#splicing
expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_splicing.mtx")
cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_splicing.csv")
genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_splicing.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_splicing.rds")

#expr_PVneg
expr <- readMM("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/expr_PVneg.mtx")
cells <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/cells_PVneg.csv")
genes <- read.csv("/rds/general/user/ao225/home/CardiaFinal/Data/Rogue_Ready/genes_PVneg.csv")

rownames(expr) <- genes[[1]]
colnames(expr) <- cells[[1]]

expr <- as.matrix(expr)
expr <- matr.filter(expr, min.cells = 50, min.genes = 10)

ent = SE_fun(expr)
saveRDS(ent,"/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_PVneg.rds")

library(dplyr)

chromatin = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_chromatin.rds")
sarcomere = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_sarcomere.rds")
splicing  = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_splicing.rds")
PVneg     = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_PVneg.rds")
control   = readRDS("/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene/ent_control.rds")

refset <- Reduce(union, list(chromatin$Gene, sarcomere$Gene,
                             splicing$Gene, PVneg$Gene, control$Gene))

ds_chromatin <- chromatin %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

ds_sarcomere <- sarcomere %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

ds_splicing <- splicing %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

ds_PVneg <- PVneg %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

ds_control<- control %>%
  filter(p.adj < 0.05) %>%
  arrange(desc(ds)) %>%
  pull(Gene)

sharedDCM <- Reduce(union, list(ds_chromatin, ds_sarcomere, ds_splicing, ds_PVneg))
DCMonly   <- setdiff(sharedDCM, ds_control)

writeLines(sharedDCM, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/sharedDCM.txt")
writeLines(DCMonly, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/DCMonly.txt")
writeLines(ds_chromatin, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/ds_chromatin.txt")
writeLines(ds_sarcomere, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/ds_sarcomere.txt")
writeLines(ds_splicing, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/ds_splicing.txt")
writeLines(ds_PVneg, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/ds_PVneg.txt")
writeLines(ds_control, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/ds_control.txt")
writeLines(refset, "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/refset.txt")

dir <- "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/per_gene"
for (g in c("chromatin","sarcomere","splicing","PVneg","control")) {
  e <- readRDS(file.path(dir, paste0("ent_", g, ".rds")))
  write.csv(e[, c("Gene","ds","p.adj")],
            file.path(dir, paste0("ds_", g, ".csv")), row.names = FALSE)
}
library(org.Hs.eg.db)
library(AnnotationDbi)

gene_dir <- "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set"

ds_chromatin <- readLines(file.path(gene_dir, "ds_chromatin.txt"))
ds_sarcomere <- readLines(file.path(gene_dir, "ds_sarcomere.txt"))
ds_splicing  <- readLines(file.path(gene_dir, "ds_splicing.txt"))
ds_PVneg     <- readLines(file.path(gene_dir, "ds_PVneg.txt"))
ds_control   <- readLines(file.path(gene_dir, "ds_control.txt"))

gene_lists <- list(chromatin = ds_chromatin, sarcomere = ds_sarcomere,
                   splicing  = ds_splicing,  PVneg     = ds_PVneg)

coreDCM <- setdiff(Reduce(intersect, gene_lists), ds_control)
cat("coreDCM n =", length(coreDCM), "\n")

core_df <- data.frame(
  Gene   = coreDCM,
  symbol = mapIds(org.Hs.eg.db, sub("\\..*$", "", coreDCM),
                  keytype = "ENSEMBL", column = "SYMBOL")
)
print(sort(core_df$symbol))
write.csv(core_df, file.path(gene_dir, "core30_annotated.csv"), row.names = FALSE)
writeLines(coreDCM, file.path(gene_dir, "coreDCM.txt"))

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggVennDiagram)

compare_ora <- function(dcm_file, ctl_file, key = "description") {
  dcm <- read.delim(dcm_file) %>% dplyr::filter(FDR < 0.05)
  ctl <- read.delim(ctl_file) %>% dplyr::filter(FDR < 0.05)
  dcm$class <- ifelse(dcm[[key]] %in% ctl[[key]], "Shared with control", "DCM-specific")
  dcm %>% dplyr::arrange(FDR) %>%
    dplyr::select(dplyr::all_of(key), enrichmentRatio, overlap, FDR, class)
}

standardise <- function(df, db) {
  names(df)[1] <- "term"          
  df$database <- db
  df
}

go   <- compare_ora("/rds/general/user/ao225/home/CardiaFinal/Data/Webgestalt/DCM/DCM_GO.txt",
                    "/rds/general/user/ao225/home/CardiaFinal/Data/Webgestalt/Control/Control_GO.txt")
kegg <- compare_ora("/rds/general/user/ao225/home/CardiaFinal/Data/Webgestalt/DCM/DCM_KEGG.txt",
                    "/rds/general/user/ao225/home/CardiaFinal/Data/Webgestalt/Control/Control_KEGG.txt")
tf   <- compare_ora("/rds/general/user/ao225/home/CardiaFinal/Data/Webgestalt/DCM/DCM_TF.txt",
                    "/rds/general/user/ao225/home/CardiaFinal/Data/Webgestalt/Control/Control_TF.txt", key = "geneSet")

all_ora <- rbind(
  standardise(go,   "GO"),
  standardise(kegg, "KEGG"),
  standardise(tf,   "TF_target")
)

col_dcm     <- "#3B6E8F"   # deep steel blue
col_shared  <- "#C6CBD1"   # neutral grey
col_dcm_lt  <- "#A7C0D2"   # light fill for Venn

mk_panel <- function(db, n_terms = 15) {
  all_ora %>%
    filter(database == db, overlap >= 5) %>%
    slice_min(FDR, n = n_terms) %>%
    mutate(term = reorder(term, -log10(FDR))) %>%
    ggplot(aes(-log10(FDR), term, fill = class)) +
    geom_col(colour = "grey30", linewidth = .3) +
    geom_text(aes(label = paste0("n=", overlap)), hjust = -0.15, size = 2.6) +
    scale_fill_manual(values = c("DCM-specific"        = col_dcm,
                                 "Shared with control" = col_shared), name = NULL)  +
    scale_x_continuous(expand = expansion(mult = c(0, .15))) +
    labs(x = expression(-log[10]~FDR), y = NULL, title = db) +
    theme_classic(base_size = 11)
}

p_go   <- mk_panel("GO")
p_kegg <- mk_panel("KEGG")
p_tf   <- mk_panel("TF_target")

write.csv(all_ora, "/rds/general/user/ao225/home/CardiaFinal/Results/Webgestalt/Table/ora_dcm_vs_control.csv", row.names = FALSE)

#----------------------------------------------------------------------------------------------------------
base <- "/rds/general/user/ao225/home/CardiaFinal/Results/ROGUE/gene_set/"
ds_control   <- readLines(paste0(base, "ds_control.txt"))
ds_chromatin <- readLines(paste0(base, "ds_chromatin.txt"))
ds_sarcomere <- readLines(paste0(base, "ds_sarcomere.txt"))
ds_splicing  <- readLines(paste0(base, "ds_splicing.txt"))
ds_PVneg     <- readLines(paste0(base, "ds_PVneg.txt"))
sharedDCM    <- readLines(paste0(base, "sharedDCM.txt"))   
DCMonly    <- readLines(paste0(base, "DCMonly.txt")) 


gene_lists <- list(chromatin = ds_chromatin, sarcomere = ds_sarcomere,
                   splicing = ds_splicing, PVneg = ds_PVneg)
all_dcm <- unique(unlist(gene_lists))
rec <- sapply(all_dcm, function(g) sum(sapply(gene_lists, function(x) g %in% x)))

pal_genotype <- c("Control" = "#66C2A5", "Splicing\n(RBM20)" = "#FC8D62",
                  "Sarcomere\n(TTN)" = "#8DA0CB", "PVneg" = "#E78AC3",
                  "Chromatin\n(LMNA)" = "#A6D854")


all_genes <- unique(c(ds_control, ds_chromatin, ds_sarcomere, ds_splicing, ds_PVneg))

pie_df <- counts %>%
  mutate(group = factor(group, levels = group)) %>%
  mutate(pct   = 100 * n / sum(n),
         ymax  = cumsum(n),
         ymin  = ymax - n,
         mid   = (ymin + ymax) / 2,
         ypos  = sum(n) - mid)          # flip for coord_polar's direction

pA <-  ggplot(pie_df, aes(x = 1, y = n, fill = group)) +
  geom_col(colour = "white", linewidth = .9, width = 1) +
  geom_text(aes(x = 1.02, y = ypos,
                label = paste0(group, "\n", n, "  (", round(pct), "%)")),
            size = 2.9, lineheight = .85) +
  scale_fill_manual(values = pal_genotype, guide = "none") +
  scale_x_continuous(limits = c(0, 1.5), expand = expansion(0)) +
  coord_polar(theta = "y") +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))



recdf <- as.data.frame(table(rec)); names(recdf) <- c("n_genotypes","n_genes")

pB <- ggplot(recdf, aes(n_genotypes, n_genes, fill = n_genotypes)) +
  geom_col(fill = col_dcm, width = .68) +
  theme_classic(base_size = 11) +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = .3))


pC <- ggVennDiagram(list(DCM = sharedDCM, Control = ds_control),
                    label = "count", label_alpha = 0, set_size = 4) +
  scale_fill_gradient(low = "#EEF2F6", high = col_dcm_lt)  +
  scale_x_continuous(expand = expansion(mult = .25)) +
  coord_cartesian(clip = "off") +
  theme(legend.position = "none", plot.margin = margin(10, 30, 10, 30))

fig <- (pA / pB / pC) | (p_go / p_kegg / p_tf)
fig <- fig + plot_layout(widths = c(1, 2.2), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 16))

ggsave("/rds/general/user/ao225/home/CardiaFinal/Results/Webgestalt/Plot/fig.png", fig, width = 22, height = 16, dpi = 300)

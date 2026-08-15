# Entropy Uncovers the Loss of Transcriptional Uniformity in Dilated Cardiomyopathy

Analysis code for an entropy-based study of cardiomyocyte transcriptional heterogeneity in
dilated cardiomyopathy (DCM), using single-nucleus RNA-seq data from Reichart et al. (2022).

## Data

Single-nucleus RNA-seq data were obtained from the CELLxGENE Discover portal
("DCM/ACM heart cell atlas: Cardiomyocytes"). No new data were generated.

Analysis was restricted to left-ventricular cardiomyocytes profiled with 10x Genomics 3' v3
chemistry. Genotypes represented by more than three donors were retained (LMNA, TTN, RBM20,
pathogenic-variant-negative) alongside non-failing controls, giving a final cohort of
91,956 cardiomyocytes from 50 donors (12 control, 38 DCM).

## External reference files

Required for pySCENIC, downloaded from https://resources.aertslab.org/cistarget/ :

- `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
- `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
- `allTFs_hg38.txt`

## Pipeline

Scripts are numbered by stage and should be run in order.

### 01_preprocessing
| Script | Purpose |
|---|---|
| `subset_and_export.ipynb` | Filters the atlas to the study cohort; exports raw UMI counts (genes x cells) plus cell/gene identifiers and metadata for ROGUE |
| `integration_umap.R` | Seurat normalisation, PCA, Harmony batch correction by donor, UMAP embedding (30 Harmony dims); exports the embedding |

### 02_rogue
| Script | Purpose |
|---|---|
| `per_donor_rogue.R` | Per-donor ROGUE scores by genotype and by cell state (UMI platform, span 0.75) |
| `per_gene_entropy.R` | Per-gene expression entropy (`SE_fun`) for control and each DCM genotype separately |
| `gene_sets.R` | Defines the reference background, per-genotype ds gene sets and the DCM-specific set |
| `core_genes.R` | Defines the genotype-independent core gene set and annotates it with gene symbols |
| `cell_states.R` | Per-cell-state ROGUE comparison, cell-state marker identification |
| `statistics_and_plots.R` | Statistical comparisons and main ROGUE figures |

### 03_enrichment
| Script | Purpose |
|---|---|
| `ora_comparison.R` | Compares WebGestalt ORA output for DCM-specific and control gene sets, classifying terms as DCM-specific or shared; builds the enrichment figure |

ORA itself was run on the WebGestalt web server using the exported gene lists, with the
`SE_fun`-tested gene set as the statistical background.

### 04_scenic
| Script | Purpose |
|---|---|
| `build_loom.ipynb` | Prepares the normalised expression matrix as `.loom` for pySCENIC |
| `jobs/step1_grn.pbs` | GRNBoost2 network inference |
| `jobs/step2_ctx.pbs` | cisTarget motif enrichment and module pruning |
| `jobs/step2b_build_regulons.py` | Builds regulons (activating modules, NES >= 3.0, direct/orthologous annotation, >= 10 targets) |
| `jobs/step3_aucell.pbs` + `jobs/aucell.py` | AUCell regulon activity scoring |
| `integrate_aucell.ipynb` | Merges AUCell scores and the Harmony UMAP into the AnnData object |
| `scenic_analysis.ipynb` | Differential regulon activity, integration with entropy (Fisher and delta-ds), per-genotype validation, figures |

pySCENIC steps were run on the Imperial College Research Computing Service.

## Paths

Scripts currently use absolute paths for the HPC environment. To run elsewhere, edit the
directory variables defined at the top of each script.

## Software

R version [X], Python version [X].
Key packages: ROGUE, Seurat, harmony, dplyr, ggplot2, patchwork, org.Hs.eg.db (R);
scanpy, anndata, loompy, pandas, numpy, scipy, pySCENIC (Python).
Full session information: `sessionInfo.txt`.

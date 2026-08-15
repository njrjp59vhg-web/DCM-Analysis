import pickle
import numpy as np
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
 
MOTIFS_FNAME = "/rds/general/user/ao225/home/CardiaFinal/Results/SCENIC/SCENIC_results/lv_regulons.motifs.csv"
REGULONS_DAT_FNAME = "/rds/general/user/ao225/home/CardiaFinal/Results/SCENIC/SCENIC_results/lv_regulons.dat"
 
# Your two actual databases (hg38 v10_clust) -- NOT the hg19 names from the
# old tutorial this code was copied from.
DB_NAMES = (
    "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings",
    "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings",
)
 
 
def contains(*elems):
    def f(context):
        return any(elem in context for elem in elems)
    return f
 
 
print("Loading motifs...")
motifs = load_motifs(MOTIFS_FNAME)
motifs.columns = motifs.columns.droplevel(0)
 
# Keep only modules built from your actual databases, and only activating modules.
mask_db = np.fromiter(map(contains(*DB_NAMES), motifs.Context), dtype=bool)
mask_activating = np.fromiter(map(contains("activating"), motifs.Context), dtype=bool)
motifs = motifs[mask_db & mask_activating]
 
print(f"Rows after database/activating filter: {len(motifs)}")
 
# Build regulons: NES >= 3.0, directly annotated or orthologous-direct, >=10 genes.
regulons = list(
    filter(
        lambda r: len(r) >= 10,
        df2regulons(
            motifs[
                (motifs["NES"] >= 3.0)
                & (
                    (motifs["Annotation"] == "gene is directly annotated")
                    | (
                        motifs["Annotation"].str.startswith("gene is orthologous to")
                        & motifs["Annotation"].str.endswith("which is directly annotated for")
                    )
                )
            ]
        ),
    )
)
 
print(f"Regulons after NES/annotation/size filter: {len(regulons)}")
 
if len(regulons) == 0:
    raise RuntimeError(
        "No regulons survived filtering -- check DB_NAMES matches your actual "
        "Context strings (print motifs.Context.iloc[0] to inspect)."
    )
 
# Rename regulons: strip suffix, keep just the TF name.
regulons = list(map(lambda r: r.rename(r.transcription_factor), regulons))
 
with open(REGULONS_DAT_FNAME, "wb") as f:
    pickle.dump(regulons, f)
 
print(f"Saved {len(regulons)} regulons to {REGULONS_DAT_FNAME}")
 
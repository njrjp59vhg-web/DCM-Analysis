import pickle
import loompy
import pandas as pd
from pyscenic.aucell import aucell

with open("/rds/general/user/ao225/home/CardiaFinal/Results/SCENIC/SCENIC_results/lv_regulons.dat", "rb") as f:
    regulons = pickle.load(f)

with loompy.connect("/rds/general/user/ao225/home/CardiaFinal/Results/SCENIC/Processed_files/lv_exp.loom") as ds:
    exp_mtx = pd.DataFrame(ds[:, :].T, index=ds.ca["CellID"], columns=ds.ra["Gene"])

auc_mtx = aucell(exp_mtx, regulons, num_workers=26)
auc_mtx.to_csv("/rds/general/user/ao225/home/CardiaFinal/Results/SCENIC/SCENIC_results/lv_auc_mtx.csv")

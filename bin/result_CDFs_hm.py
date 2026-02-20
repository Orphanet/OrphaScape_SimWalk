
import argparse
import pandas as pd
import numpy as np
import os


import path_variable as PV
from set_log import setup_logging, get_logger
from pathlib import Path
import logging



from sklearn.metrics import auc as sk_auc
from matplotlib import pyplot as plt

"""

script to calculate harmonic mean and get the CDF of each method (no group of disorder analysis here )

"""


# Logging configuration
setup_logging(level=logging.INFO,console=False,filename=f"results.log"    )
log = get_logger(Path(__file__).stem)

p = argparse.ArgumentParser()
p.add_argument("--run-name", default="default", help="Run name (used for output subfolder)")
p.add_argument("--fig-num",  default="",        help="Figure number suffix (e.g. 1 → CDF_fig1.svg)")
p.add_argument("--do-subsumed", type=int, default=0, choices=[0, 1],
               help="1 = patients_subsumed.xlsx, 0 = patients.xlsx")
p.add_argument("--product4",   default="",  help="Product4 version string (filter parquets)")
p.add_argument("--combine",    nargs="+", default=[], help="Combine methods to include")
p.add_argument("--sm-method",  nargs="+", default=[], help="SM methods to include")
p.add_argument("--vector-strs",nargs="+", default=[], help="Compact vector strings to include (e.g. 11111)")
args = p.parse_args()

out_dir = Path(PV.PATH_OUTPUT) / args.run_name
out_dir.mkdir(parents=True, exist_ok=True)

fig_suffix = f"{args.fig_num}" if args.fig_num else ""



patients = pd.read_excel(PV.get_patient_path(args.do_subsumed))
couple_patients = patients[["phenopacket","Disease"]].drop_duplicates()
# rename for the merge with sm df
couple_patients.columns = ["patients","RDs"]
 
dict_df_ra_sm = {}

list_ra= os.listdir(PV.get_dp_path(args.do_subsumed))
 

list_ra_sm = []
# Build expected filenames from config filters (if provided)
_expected = None
if args.combine and args.sm_method and args.vector_strs and args.product4:
    _expected = {
        f"{c}_{s}_{args.product4}_{v}.parquet"
        for c in args.combine
        for s in args.sm_method
        for v in args.vector_strs
    }

for ra in list_ra:
    if ("parquet" in ra) and ('~' not in ra):
        if _expected is None or ra in _expected:
            list_ra_sm.append(ra)
 
for ra in list_ra_sm:
# if ('xlsx' in ra) and ('CDF' not in ra): # RDI
    ra_rslt=  ra.rsplit('.', 1)[0] # remore the extension

    ## list of df from RA but with different vector 
    df_sm = pd.read_parquet(PV.get_dp_path(args.do_subsumed) / ra)
 
    list_sm = df_sm["patients"].drop_duplicates().tolist()
    dict_df_ra_sm[ra_rslt] = df_sm
 

    log.info(f"for ra : {ra_rslt}, nb of patient is {len(list_sm)}")


 
######################################################
#### compare patient match

list_rank_all = []
all_interactions = set()
for one_patient in list_sm:
    # log.info(one_patient)
    one_couple = couple_patients[couple_patients["patients"]==one_patient]
 
    rdi = one_couple["RDs"].values[0]


    for key,df_sm  in dict_df_ra_sm.items():
        ## extract rank from sm
        matches_sm = pd.merge(one_couple, df_sm, on=["patients", "RDs"])
        if matches_sm.empty:
            rank_sm = np.nan
        else:
            rank_sm = int(matches_sm['rank'].values[0])
        
        list_rank_all.append((one_patient,rdi,key,rank_sm))




df = pd.DataFrame(list_rank_all,columns=["patient","RD","metric",'rank'])

df_compare_rank_wide = df.pivot(index=['patient', 'RD'], 
                   columns='metric', 
                   values='rank')\
            .reset_index()



#######################################################################################
# -----------script 13_harmionic_mean_df.py only section compare rank----------------#
#######################################################################################

all_interecation = []

col_compare_rank_df = df_compare_rank_wide.columns.tolist()

for onecol in col_compare_rank_df:
    if onecol not in ['patient','RD']:
        rank_method = df_compare_rank_wide[onecol].astype(float)
        harmonic_mean = len(rank_method) / (1.0 / rank_method).sum()
        all_interecation.append((onecol,harmonic_mean))

df_hm_general = pd.DataFrame(all_interecation,columns=['method','mean_hm'])  
df_hm_general.to_excel(out_dir / f"hm_table_{fig_suffix}.xlsx")

#######################################################################################
# -----------  CDF-related to cdf and table hm in the results section---------------#
#######################################################################################
# 11 to vizualise well the top 10
rank_f = 11                              
  
methods = df_compare_rank_wide.columns.tolist()


# 3) create a “filtered” version where any rank >10 → NaN
methods = df_compare_rank_wide.columns.tolist()[2:]
df_compare_rank_filtered = df_compare_rank_wide[['patient', 'RD']].join(
    df_compare_rank_wide[methods].where(df_compare_rank_wide[methods] <= rank_f)
)
###################################################
# grab a colormap with plenty of colors
cmap   = plt.get_cmap('Set1', len(methods))
colors = cmap(np.arange(len(methods)))
 
 
###################################################
color_dict = dict(zip(methods, colors))

spacing = 0
lines_with_auc = []
plt.figure(figsize=(12,6))
for i, col in enumerate(methods):

    data_rank = df_compare_rank_filtered[col]

    n_total   = len(data_rank)
    n_missing = data_rank.isna().sum()
    missing_pct = n_missing / n_total

    x = np.sort(data_rank)
    y = np.arange(1, n_total+1) / n_total
    y_offset = y + i * spacing

    # ---- compute AUC of the empirical CDF
    data_rank_na = df_compare_rank_filtered[col].dropna()
    x_auc = np.sort(data_rank_na)
    y_auc = np.arange(1,  len(data_rank_na)+1) / len(data_rank_na)
    auc_score = sk_auc(x_auc, y_auc) if len(data_rank_na) >= 2 else float('nan')


    # draw the step curve and grab the handle
    (line_handle,) = plt.step(
        x, y,
        where='post',
        color=color_dict[col],
        label=f"{col}", #(AUC={auc_score:.3f})",
        alpha=0.8
    )
    lines_with_auc.append((line_handle, auc_score))



plt.xlabel('r = rank')
plt.ylabel('P(rank ≤ r)')

 
# now sort your handles by auc and re‐draw the legend
lines_with_auc.sort(key=lambda t: t[1],reverse=True)  # ascending AUC
handles_sorted = [h for h,_ in lines_with_auc]
labels_sorted  = [h.get_label() for h in handles_sorted]

plt.legend(handles_sorted, labels_sorted, loc='best', borderaxespad=0.5)
 
fig_cdf_nb = int(fig_suffix) + 1
plt.grid(True)
plt.tight_layout(rect=[0, 0, 0.75, 1])
plt.savefig(out_dir / f"CDF_fig_{fig_cdf_nb}.svg", dpi=300)
# plt.savefig(f"{PV.PATH_OUTPUT}/CDF.png", dpi=300)

# plt.show()


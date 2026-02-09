

import pandas as pd
import numpy as np
import os


import path_variable as PV
from set_log import setup_logging, get_logger
from pathlib import Path
import logging


 
from sklearn.metrics import auc as sk_auc
from matplotlib import pyplot as plt



setup_logging(
    level=logging.INFO,
    console=True,  # Also display in console for debug
    filename=f"{Path(__file__).stem}.log"
)
log = get_logger(Path(__file__).stem)


## script to calculate harmonic mean and get the CDF of each method (no group of disorder analysis here )

## files that will regroup all mp_sm results for each patient and disorder (from /mp_sm folder)
 
patients = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT)
list_patients = patients["phenopacket"].drop_duplicates().tolist()
couple_patients = patients[["phenopacket","Disease"]].drop_duplicates()
# rename for the merge with sm df
couple_patients.columns = ["patients","RDs"]
 
dict_df_ra_sm = {}

#"/home/maroua/Bureau/wip/my_pipeline_v2/output/mp_sm/"
list_ra= os.listdir(PV.PATH_OUTPUT_SM)
 

list_ra_sm = []
for ra in list_ra:
    if ("parquet" in ra) and ('~' not in ra): #folder_pd4 in ra:
        list_ra_sm.append(ra)
 
for ra in list_ra_sm:
# if ('xlsx' in ra) and ('CDF' not in ra): # RDI
    ra_rslt=  ra.rsplit('.', 1)[0] # remore the extension

    ## list of df from RA but with different vector 
    df_sm = pd.read_parquet(f"{PV.PATH_OUTPUT_SM}/{ra}")
    list_sm = df_sm["patients"].drop_duplicates().tolist()
    dict_df_ra_sm[ra_rslt] = df_sm


    print(f"for ra : {ra_rslt}, nb of patient is {len(list_sm)}")


 
######################################################
#### compare patient match

list_rank_all = []
all_interactions = set()
for one_patient in list_patients:
    # print(one_patient)
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
df_hm_general.to_excel(f"{PV.PATH_OUTPUT}/hm_table_1-3.xlsx")

#######################################################################################
# -----------  CDF-related to cdf and table hm in the results section---------------#
#######################################################################################
 

print(f"START  make cdf.py")

rank_f = 11                              
  
methods = df_compare_rank_wide.columns.tolist()
 



# 3) create a “filtered” version where any rank >10 → NaN
methods = df_compare_rank_wide.columns.tolist()[2:]
df_compare_rank_filtered = df_compare_rank_wide[['patient', 'RD']].join(
    df_compare_rank_wide[methods].where(df_compare_rank_wide[methods] <= rank_f)
)
###################################################
# grab a colormap with plenty of colors
cmap   = plt.get_cmap('nipy_spectral', len(methods))
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
    auc_score = sk_auc(x_auc, y_auc)


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
# plt.title('Cumulative distribution of the candidate ORPHAcode rank across GSSM')
# plt.title('Cumulative distribution of the candidate ORPHAcode rank with FunSimMaxAsym Resnik GSSM')
# plt.title('Cumulative distribution of the candidate ORPHAcode rank with the optimal GSSM\n (FunSimMaxAsym + Resnik)  and without removal   of subsum HPO terms')

 
# now sort your handles by auc and re‐draw the legend
lines_with_auc.sort(key=lambda t: t[1],reverse=True)  # ascending AUC
handles_sorted = [h for h,_ in lines_with_auc]
labels_sorted  = [h.get_label() for h in handles_sorted]

plt.legend(handles_sorted, labels_sorted, loc='best', borderaxespad=0.5)

 

plt.grid(True)
plt.tight_layout(rect=[0, 0, 0.75, 1])
plt.savefig(f"{PV.PATH_OUTPUT}/CDF.svg", dpi=300)
# plt.savefig(f"{PV.PATH_OUTPUT}/CDF.png", dpi=300)

plt.show()

print(f"END  8_cdf_auc_mrr.py")

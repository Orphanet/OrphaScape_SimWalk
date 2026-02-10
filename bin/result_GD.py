
import numpy as np
import pandas as pd

from set_log import setup_logging, get_logger
from pathlib import Path
import logging
import path_variable as PV
 
# ==============================================================================
#   FUNCTIONS
# ==============================================================================


def get_related_group(rdi,onep,rd_group_map,name_method_json_key,df_method,get_classif_method,list_rd_method,df_all_classif):
# patient,method,rd,, rank,classif_id,classif_name,type, type_id,type_name,is_rdi
    
    relationship = []
 
    for one_classif in get_classif_method:
        ## get the classif related to the pp 
        if one_classif not  in df_all_classif['root'].drop_duplicates().tolist():
            log.info(f"This classif : {one_classif} is not available in df_all_classif")
            continue

        df_cls = df_all_classif[df_all_classif["root"] == one_classif]
        classif_name = df_cls[df_cls["root"] == one_classif]['root_name'].values[0]

        for onerd in list_rd_method : #['ORPHA:610']: #['ORPHA:610'] 101685 # list_rd_method
            is_rdi = 'n'
            if onerd == rdi:
                is_rdi = "y"
            # log.info(onerd)
            try :
                rank_method = df_method[df_method['ORPHAcode'] == onerd]['rank'].values[0]
            # if get except 
            except IndexError:
                rank_method = np.nan
                log.info('no rdi rank found in sm patient')

            rd_group = rd_group_map.get(onerd)

            ## get the direct parent thus the n+1
            direct_parents = df_cls[df_cls['child_id'] == onerd]
            direct_parents_list = direct_parents['parent_id'].drop_duplicates().tolist()

            ##  Traverse hierarchy separately for each direct parent
            for parent_id in direct_parents_list: # ['ORPHA:166472']:
                visited = set()
                current_ids = [parent_id]
                n = 1
                
                while current_ids:
                    next_ids = []
                    for cid in current_ids:
                        if cid in visited:
                            continue
                        visited.add(cid)

                        direct_group = rd_group_map.get(cid)
                        # parents of the current node (i.e., grand-parents of the RD at level=2, etc.)
                        parent_of_cid = df_cls.loc[df_cls['child_id'] == cid, 'parent_id'].drop_duplicates().tolist()


                        # Case where the direct parent is the highest (before the classif)
                        if (parent_of_cid):
                            # there can be multiple; we’ll emit one row per “direct cid”
                            # and pick “undirect” as the first next parent for this row 
                            undirect_pid = parent_of_cid[0]
                            undirect_group = rd_group_map.get(undirect_pid)

                            next_ids.extend(parent_of_cid)                            
                        else:
                            undirect_pid = np.nan
                            undirect_group = np.nan

                        relationship.append({
                            'patient': onep,
                            'method': name_method_json_key,
                            'classif_id': one_classif,
                            'rd_id': onerd,
                            'is_rdi': is_rdi,
                            'rd_rank': rank_method,
                            'rd_group': rd_group,
                            'level': n,                 # 1 = direct parent, 2 = grand-parent, ...
                            'direct_parent': cid,
                            'direct_group': direct_group,
                            # 'undirect_parent': undirect_pid, # next level up if it exists
                            # 'undirect_name': undirect_name,
                            # 'undirect_type': undirect_type,
                            # 'undirect_group': undirect_group,
                        })                         

                            
                    # go up one level
                    current_ids = [nid for nid in next_ids if nid not in visited]
                    n += 1
 
    return relationship



if __name__ == "__main__":
    # ==============================================================================
    # : CONFIGURATION LOG FILE
    # ==============================================================================

    # Logging configuration
    setup_logging(level=logging.INFO,console=False,filename=f"results.log"    )  
    log = get_logger(Path(__file__).stem)


    # ==============================================================================
    #  DATA LOADING
    # ==============================================================================
    # Get path random walk with restart results 
    most_recent_file_rarw = max(PV.PATH_OUTPUT_FOLDER_RW.glob("*/*/"), key=lambda f: f.stat().st_mtime if f.is_file() else 0)
     
    # Load Semantic Similarity results
    most_recent_file_ra = max(PV.PATH_OUTPUT_SM.glob("*"), key=lambda f: f.stat().st_mtime if f.is_file() else 0)
    df_sm_all = pd.read_parquet( most_recent_file_ra)

    ## load  patient name from  df_sm_all and not from patient df output because some patient might not be there is subset was done 
    list_patient = df_sm_all['patients'].drop_duplicates().tolist()
    log.info(f"Loaded {len(list_patient)} patients")

    # -------------------------------------------------------------------------
    # Load Orphanet reference databases
    # -------------------------------------------------------------------------
    # Product 1: Contains group and type information for each ORPHAcode
    df_pd1 = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT1,index_col=0    )

    # Product 7: Contains preferential parent (PP) classification for each disease
    df7 = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT7,index_col=0    )

    # Classification hierarchy: Contains all groups and disorder relationships (excluding genetic  )
    df_classif = pd.read_excel(PV.PATH_OUTPUT_DF_CLASSIF, index_col=0   )

 
    # -------------------------------------------------------------------------
    # Load patient information with confirmed diagnoses
    # -------------------------------------------------------------------------
    patients = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT, index_col=0)


    # ==============================================================================
    #   BUILD CLASSIFICATION GROUPS FOR ANALYSIS
    # ==============================================================================
    # Create a lookup dictionary: ORPHAcode -> Group type, This avoids repeated DataFrame lookups during processing
    rd_group_map = (df_pd1.set_index("ORPHAcode")["Group"].to_dict())

    # Configuration: Number of top-ranked diseases to analyse per patient
    user_nb_top_rd = 50

    # Container for all patient results
    list_for_classif_global = []

    #  Process each patient
    for i,user_seed in enumerate(list_patient): #   list_patient
        log.info(f"--------{user_seed} ---rarw and sm part------")

        ########################################################
        ## get the RDI of the user_seed
        try:
            rdi = patients[patients['phenopacket'] == user_seed]['Disease'].drop_duplicates().values[0]
        except IndexError:
            rdi = ""
            log.info(f"  No RDI found for patient {user_seed}")

        # -------------------------------------------------------------------------
        ## Process Semantic Similarity (RA) Method Results
        # -------------------------------------------------------------------------
        df_sm_init= df_sm_all[df_sm_all['patients'] == user_seed]
        df_sm_init = df_sm_init.rename(columns={'RDs': 'ORPHAcode'})
        # Filter to top N ranked diseases (with user_nb_top_rd)
        df_sm = df_sm_init[df_sm_init['rank'].isin([*range(1,user_nb_top_rd+1)])].reset_index().sort_values(by='rank')

        # Build list of diseases to analyse (ensure RDI is included)
        list_rd_sm = df_sm['ORPHAcode'].tolist()
        if rdi not in list_rd_sm:
            list_rd_sm.append(rdi)
            log.info(f"  RA: Added RDI ({rdi}) to analysis list")
        # Identify relevant classifications for these diseases
        get_classif_sm = df_classif[(df_classif['parent_id'].isin(list_rd_sm)) | (df_classif['child_id'].isin(list_rd_sm))]['root'].drop_duplicates()

        # Extract hierarchy relationships Using df_sm_init (not df_sm) to preserve all rank information
        result_sm =  get_related_group(rdi,user_seed,rd_group_map,"RA",df_sm_init,get_classif_sm,list_rd_sm,df_classif)

        log.info(f"  RA: Top {user_nb_top_rd} RDs span {len(get_classif_sm)} classifications")


        # -------------------------------------------------------------------------
        ## Process Semantic Similarity (RA) Method Results
        # -------------------------------------------------------------------------
        # Load and prepare RARW results for this patient
        df_pg_init= pd.read_parquet(f"{most_recent_file_rarw}/{user_seed}.parquet")
        df_pg_init = df_pg_init.rename(columns={'rank_pg': 'rank'})
        # set index to ORPHAcode name 
        df_pg_init.index.name = "ORPHAcode"
        df_pg_init = df_pg_init.reset_index()

        df_pg = df_pg_init[df_pg_init['rank'].isin([*range(1,user_nb_top_rd+1)])].reset_index().sort_values(by='rank_sum_degres_pg')

        list_rd_pg = df_pg['ORPHAcode'].tolist()
        if rdi not in list_rd_pg:
            list_rd_pg.append(rdi)
            log.info(f'{user_seed} RARW :RDI: {rdi} added')

        get_classif_rarw = df_classif[(df_classif['parent_id'].isin(list_rd_pg)) | (df_classif['child_id'].isin(list_rd_pg))]['root'].drop_duplicates()

        result_pg =  get_related_group(rdi,user_seed,rd_group_map,"RARW",df_pg_init,get_classif_rarw,list_rd_pg,df_classif)

        log.info(f"  RARW: Top {user_nb_top_rd} RDs span {len(get_classif_rarw)} classifications")


        # -------------------------------------------------------------------------
        ## Combine Results from Both Methods
        # -------------------------------------------------------------------------
        concat_list = []  
        concat_list =  result_sm + result_pg
        df_classif_grp_rd = pd.DataFrame(concat_list)
        # global append for each patients 
        list_for_classif_global.append(df_classif_grp_rd)


    # Merge all patient results into a single DataFrame
    df_global_classif = pd.concat(list_for_classif_global)
    # df_global_classif.to_excel(f"{PV.PATH_OUTPUT}/global_classif.xlsx")
    log.info(f"Total records: {len(df_global_classif)}")


    # ==============================================================================
    # FROM ANALYSIS 2 IN RESULT SECTION (TABLE 5)
    # ==============================================================================
    # 
    # This analysis compares the classification groups (branch) of each candidate disease
    # with the groups of the true diagnosis (RDI), using only the Preferential
    # Parent (PP) classification from Orphanet Product 7.
    #
    # A match ("orange" condition) occurs when:
    #   - The candidate shares the same PP classification as the RDI
    #   - The candidate shares at least one "Group of disorders" ancestor with the RDI
    #
    # ==============================================================================


    # Configuration parameters
    # Maximum hierarchy level to consider
    level = 8
    # Number of top candidates to analyse per patient
    n_top = 10

    # Build lookup dictionary: ORPHAcode -> Preferential Parent classification
    orpha_to_pp = df7.dropna(subset=['ORPHAcode', 'Classif_id']) \
                        .drop_duplicates(['ORPHAcode']) \
                        .set_index('ORPHAcode')['Classif_id'].to_dict()

    results = []  
    # Process each patient-method combination
    for (patient, mtd), sub in df_global_classif.groupby(['patient', 'method']): # patient
        # Get the RDI information for this patient-method combination
        is_rdi_sub = sub[sub['is_rdi'] == 'y']
        if is_rdi_sub.empty:
            continue

        RDI_rd_id = is_rdi_sub['rd_id'].iloc[0]
        RDI_rank  = is_rdi_sub['rd_rank'].iloc[0] # rd_rank rank

        # Get the Preferential Parent classification for the RDI
        RDI_pp = orpha_to_pp.get(RDI_rd_id)
        if pd.isna(RDI_pp):
            continue

        # Filter RDI hierarchy to the specified level
        is_rdi_sub_level = is_rdi_sub[is_rdi_sub['level'] <= level]
        # Extract all "Group of disorders" ancestors of the RDI
        RDI_group_ids = (is_rdi_sub_level.loc[is_rdi_sub_level['direct_group'] == "Group of disorders", 'direct_parent']
                            .dropna().astype(str).unique().tolist())
        RDI_group_set = set(RDI_group_ids)

        # Get the top-N distinct candidate diseases by rank
        topn_rd = (sub[['rd_id', 'rd_rank']]
                    .drop_duplicates('rd_id')
                    .sort_values('rd_rank', ascending=True)
                    .head(n_top)['rd_id'].tolist())
        
        # Compare each candidate with the RDI
        for cand in topn_rd:
            cand_rows = sub[sub['rd_id'] == cand]
            # Get the Preferential Parent classification for the candidate
            cand_pp = orpha_to_pp.get(cand)
            if pd.isna(cand_pp):
                continue
            
            # Filter to only the PP classification at the specified level
            sub_cand_pp = cand_rows[cand_rows['classif_id'] == cand_pp]
            cand_rows_level = sub_cand_pp[sub_cand_pp['level'] <= level]
            # Extract all "Group of disorders" ancestors of the candidate
            cand_group_ids = (cand_rows_level.loc[cand_rows_level['direct_group'] == "Group of disorders", 'direct_parent']
                                .dropna().astype(str).unique().tolist())
            
            # Calculate overlaps
            group_overlap = RDI_group_set.intersection(cand_group_ids)
            classif_overlap = {RDI_pp}.intersection({cand_pp})

            results.append({
                'level': level,
                'patient': patient,
                'method': mtd,

                'RDI_rd_id': RDI_rd_id,
                'RDI_rank': RDI_rank,
                'RDI_classif_ids': [RDI_pp],           
                'RDI_group_ids': RDI_group_ids,        
                'candidate_rd_id': cand,
                'candidate_rank': cand_rows['rd_rank'].min(),
                'candidate_classif_ids': [cand_pp],    
                'candidate_group_ids': cand_group_ids, 

                'classif_overlap': int(len(classif_overlap)),
                'group_overlap': int(len(group_overlap)),
                'match_groups': sorted(group_overlap),
            })
    # Build summary DataFrame and identify matches
    summary_df = pd.DataFrame(results)
    # "Orange" condition: Both classification AND group overlap exist
    summary_df['is_orange'] = (summary_df['classif_overlap'] > 0) & (summary_df['group_overlap'] > 0)

 
    # Count distinct matching groups per patient-method
    rows = []
    for (patient, method), sub in summary_df.groupby(['patient', 'method']):
        # Collect all unique matching group IDs (stored as sets, no duplicates)
        groups = set()
        for gs in sub.loc[sub['is_orange'], 'match_groups']:
            groups.update(gs)
        rows.append({'patient': patient, 'method': method, 'intensity': int(len(groups)),
                    })
    long_df = pd.DataFrame(rows)


    # -------------------------------------------------------------------------
    # Create pivot table for comparison (TABLE 5)
    # -------------------------------------------------------------------------
    pivot = long_df.pivot(index='patient', columns='method', values='intensity') \
                    .fillna(0).astype(int)

    # Keep a consistent method order if present
    col_order = [m for m in ['RA','RARW'] if m in pivot.columns]
    pivot = pivot.reindex(columns=col_order)

    # Extract cases where methods disagree (TABLE 5 final result)
    if len(set(pivot['RA']).intersection(pivot['RARW'])) == len(set(pivot['RA'])):
        # condition where RA and RARW have the same values
        df_filter =pivot.copy()
    else:
        df_filter = pivot[pivot['RA'] != pivot['RARW']]


 
    df_filter.to_excel(f"{PV.PATH_OUTPUT}/nb_gp_match_rdi_table5.xlsx")

    log.info(f"\nPatients with method disagreement: {len(df_filter)}")
    log.info("\n*** TABLE 5: Group Overlap Comparison ***")
    log.info(df_filter.head(20))


    # ==============================================================================
    # FROM : ANALYSIS 1 IN THE RESULT SECTION - HARMONIC MEAN CALCULATION (TABLE 4)
    # ==============================================================================
    #
    # This analysis evaluates ranking quality by computing the harmonic mean of
    # ranks for diseases that share the same immediate  group as
    # the true diagnosis (RDI).
    #
    # ==============================================================================
    df_f = df_global_classif.drop_duplicates()
    df_f = df_f.rename(columns={'rd_rank': 'rank','type': "patient"})

    list_method = df_f['method'].unique()

    all_interaction = []
    for onep in list_patient:
        mini_df = df_f[(df_f['patient'] == onep)]
        # get the RDI
        rdi_ids = mini_df.loc[mini_df["is_rdi"] == "y", "rd_id"].unique()
        

        for onem in list_method:
            mini_df_method = mini_df[mini_df['method'] == onem]


            mini_df_method =mini_df_method[mini_df_method['rd_group'] == "Disorder"]
            mini_df_method = mini_df_method[mini_df_method['level'] == 1]
            # Get the group of disorder of the rdi 
            list_rdi_groups = mini_df_method[mini_df_method["rd_id"].isin(rdi_ids)]['direct_parent'].unique()
            # Match rd having the same group of disorder belonging to the rdi 
            rd_match_group_rdi = mini_df_method[mini_df_method['direct_parent'].isin(list_rdi_groups)]['rd_id'].unique()


            df_match_group_rdi = mini_df_method[mini_df_method['rd_id'].isin(rd_match_group_rdi)].drop_duplicates()

            # get the rank of each RDs
            df_for_hm = df_match_group_rdi[['rd_id','rank']].drop_duplicates()

            # extract ranks
            ranks_gp = df_for_hm['rank']

            h_mean_rd_same_group_rdi = len(ranks_gp) / sum(1.0 / r for r in ranks_gp)
    
            log.info(f"{onep} - {onem}\tHarmonic mean (manual calculation): {h_mean_rd_same_group_rdi}")

            all_interaction.append((onep,onem,h_mean_rd_same_group_rdi))
    

        log.info(f"{onep} - {onem}\tHarmonic mean (manual calculation): {h_mean_rd_same_group_rdi}")
 
    df_hm_general = pd.DataFrame(all_interaction,columns=['patient','method','hm'])
    # df_hm_general.to_excel('hm_each_patient.xlsx')
    df_hm_group  = df_hm_general[['patient','method','hm']]



    # Compute the average HM for each method
    method_hm_group = (
        df_hm_group
        .groupby('method', as_index=False)['hm']
        .mean()
        .rename(columns={'hm': 'mean_hm'})
    )


    method_hm_group.to_excel(f"{PV.PATH_OUTPUT}/hm_group_table4.xlsx")
    log.info("\n*** TABLE 4: Mean Harmonic Mean by Method ***")
    log.info(method_hm_group)

 
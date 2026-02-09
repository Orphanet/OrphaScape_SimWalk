"""
Aalyses the relationship between rare diseases (RDs) and their 
classification hierarchies using two methods:
    - RA (Resnik-based Semantic Similarity)
    - RARW (Random Walk approach)


"""

import numpy as np
import pandas as pd

from set_log import setup_logging, get_logger
from pathlib import Path
import logging
import path_variable as PV

# ==============================================================================
# SECTION 1: CORE FUNCTIONS
# ==============================================================================

def get_related_group(
    rdi: str,
    onep: str,
    rd_group_map: dict,
    name_method_json_key: str,
    df_method: pd.DataFrame,
    get_classif_method: pd.Series,
    list_rd_method: list,
    df_all_classif: pd.DataFrame
):
    """
    Traverse the classification hierarchy for each rare disease and extract
    parent-child relationships at each level.

    This function iterates through all relevant classifications and rare diseases,
    building a complete picture of where each RD sits within the Orphanet hierarchy.
    """
    relationship = []

    # -------------------------------------------------------------------------
    # Iterate through each relevant classification
    # -------------------------------------------------------------------------
    for one_classif in get_classif_method:
        
        # Validate that this classification exists in our reference data
        if one_classif not in df_all_classif['root'].drop_duplicates().tolist():
            log.info(f"WARNING: Classification '{one_classif}' not found in df_all_classif")
            continue

        # Extract the subset of hierarchy for this classification
        df_cls = df_all_classif[df_all_classif["root"] == one_classif]

        # ---------------------------------------------------------------------
        # Process each rare disease in the list
        # ---------------------------------------------------------------------
        for onerd in list_rd_method:
            is_rdi = 'y' if onerd == rdi else 'n'

            # Retrieve the rank of this disease from the method's results
            try:
                rank_method = df_method[df_method['ORPHAcode'] == onerd]['rank'].values[0]
            except IndexError:
                rank_method = np.nan
                log.info(f"WARNING: No rank found for {onerd} in patient's results")

            # Get the group  for this disease
            rd_group = rd_group_map.get(onerd)

            # -----------------------------------------------------------------
            # Find direct parents (level n+1 in hierarchy)
            # -----------------------------------------------------------------
            direct_parents = df_cls[df_cls['child_id'] == onerd]
            direct_parents_list = direct_parents['parent_id'].drop_duplicates().tolist()

            # -----------------------------------------------------------------
            # Traverse hierarchy upward for each direct parent branch
            # -----------------------------------------------------------------
            for parent_id in direct_parents_list:
                visited = set()
                current_ids = [parent_id]
                n = 1  # Level counter (1 = direct parent, 2 = grandparent, etc.)

                while current_ids:
                    next_ids = []

                    for cid in current_ids:
                        # Skip already visited nodes to avoid cycles
                        if cid in visited:
                            continue
                        visited.add(cid)

                        # Get group type for current ancestor
                        direct_group = rd_group_map.get(cid)

                        # Find parents of current node (next level up)
                        parent_of_cid = df_cls.loc[
                            df_cls['child_id'] == cid, 'parent_id'
                        ].drop_duplicates().tolist()

                        # Determine indirect parent info for potential future use
                        if parent_of_cid:
                            undirect_pid = parent_of_cid[0]
                            undirect_group = rd_group_map.get(undirect_pid)
                            next_ids.extend(parent_of_cid)
                        else:
                            # Reached the top of this classification branch
                            undirect_pid = np.nan
                            undirect_group = np.nan

                        # Record this relationship
                        relationship.append({
                            'patient': onep,
                            'method': name_method_json_key,
                            'classif_id': one_classif,
                            'rd_id': onerd,
                            'is_rdi': is_rdi,
                            'rd_rank': rank_method,
                            'rd_group': rd_group,
                            'level': n,
                            'direct_parent': cid,
                            'direct_group': direct_group,
                        })

                    # Move up one level, excluding already visited nodes
                    current_ids = [nid for nid in next_ids if nid not in visited]
                    n += 1
    return relationship



 

def main():
    # ==============================================================================
    # : CONFIGURATION LOG FILE
    # ==============================================================================


    setup_logging(
        level=logging.INFO,
        console=True,  # Also display in console for debug
        filename=f"{Path(__file__).stem}.log"
    )
    log = get_logger(Path(__file__).stem)


    # ==============================================================================
    #  DATA LOADING
    # ==============================================================================
    most_recent_file_rarw = max(PV.PATH_OUTPUT_FOLDER_RW.glob("*/*/"), key=lambda f: f.stat().st_mtime if f.is_file() else 0)
    # most_recent_file_rarw = "/home/maroua/Bureau/wip/my_pipeline_v2/output/rarw/0.3/3_2_2_2_1_concat_matrix/"
    
    # Load Semantic Similarity results
    most_recent_file_ra = max(PV.PATH_OUTPUT_SM.glob("*"), key=lambda f: f.stat().st_mtime if f.is_file() else 0)
    # most_recent_file_ra = "/home/maroua/Bureau/wip/my_pipeline_v2/output/mp_sm/3_2_2_2_1_rsd_resnik_n_productmai2024_all_vectors_withontologyX.xlsx"
    df_sm_all = pd.read_parquet( most_recent_file_ra)
    # df_sm_all = pd.read_excel( most_recent_file_ra)


    # -------------------------------------------------------------------------
    # Load patient information with confirmed diagnoses
    # -------------------------------------------------------------------------
    patients = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT, index_col=0)
    list_patient = patients['phenopacket'].drop_duplicates().tolist()
    log.info(f"Loaded {len(list_patient)} patients")

    # -------------------------------------------------------------------------
    # Load Orphanet reference databases
    # -------------------------------------------------------------------------

    # Product 1: Contains group and type information for each ORPHAcode
    df_pd1 = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT1,index_col=0    )

    # Product 7: Contains preferential parent (PP) classification for each disease
    df7 = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT7,index_col=0    )

    # Classification hierarchy: Contains all groups and disorder relationships (excluding genetic  )
    df_classif = pd.read_excel(PV.PATH_OUTPUT_DF_CLASSIF_NOG, index_col=0   )


    # ==============================================================================
    # SECTION 4: BUILD CLASSIFICATION GROUPS FOR ANALYSIS
    # ==============================================================================
 

    # Create a lookup dictionary: ORPHAcode -> Group type
    # This avoids repeated DataFrame lookups during processing
    rd_group_map = df_pd1.set_index("ORPHAcode")["Group"].to_dict()

    # Configuration: Number of top-ranked diseases to analyze per patient
    user_nb_top_rd = 20

    # Container for all patient results
    list_for_classif_global = []

    # -------------------------------------------------------------------------
    # Process each patient
    # -------------------------------------------------------------------------
    for i, user_seed in enumerate(list_patient):
        log.info(f"\n[{i+1}/{len(list_patient)}] Processing patient: {user_seed}")
        log.info("-" * 50)

        # ---------------------------------------------------------------------
        # Get the true diagnosis (RDI) for this patient
        # ---------------------------------------------------------------------
        try:
            rdi = patients[patients['phenopacket'] == user_seed]['Disease'].drop_duplicates().values[0]
        except IndexError:
            rdi = ""
            log.info(f"  WARNING: No RDI found for patient {user_seed}")

        # =====================================================================
        # PART A: Process Semantic Similarity (RA) Method Results
        # =====================================================================
        
        # Load and prepare SM results for this patient
        df_sm_init = df_sm_all[df_sm_all['patients'] == user_seed]
        df_sm_init = df_sm_init.rename(columns={'RDs': 'ORPHAcode'})

        # Filter to top N ranked diseases (with user_nb_top_rd)
        df_sm = df_sm_init[
            df_sm_init['rank'].isin([*range(1, user_nb_top_rd + 1)])
        ].reset_index().sort_values(by='rank')

        # Build list of diseases to analyze (ensure RDI is included)
        list_rd_sm = df_sm['ORPHAcode'].tolist()
        if rdi not in list_rd_sm:
            list_rd_sm.append(rdi)
            log.info(f"  RA: Added RDI ({rdi}) to analysis list")

        # Identify relevant classifications for these diseases
        get_classif_sm = df_classif[
            (df_classif['parent_id'].isin(list_rd_sm)) | 
            (df_classif['child_id'].isin(list_rd_sm))
        ]['root'].drop_duplicates()

        # Extract hierarchy relationships
        # Note: Using df_sm_init (not df_sm) to preserve all rank information
        result_sm = get_related_group(
            rdi, user_seed, rd_group_map, "RA",
            df_sm_init, get_classif_sm, list_rd_sm, df_classif
        )
        log.info(f"  RA: Top {user_nb_top_rd} RDs span {len(get_classif_sm)} classifications")

        # =====================================================================
        # Process Random Walk (RARW) Method Results
        # =====================================================================
        
        # Load and prepare RARW results for this patient
        df_pg_init = pd.read_excel(str(most_recent_file_rarw) + "/" + user_seed + ".xlsx")
        df_pg_init = df_pg_init.rename(columns={'rank_pg': 'rank','Unnamed: 0': 'ORPHAcode'})

        df_pg = df_pg_init[
            df_pg_init['rank'].isin([*range(1, user_nb_top_rd + 1)])
        ].reset_index().sort_values(by='rank_sum_degres_pg')

        list_rd_pg = df_pg['ORPHAcode'].tolist()
        if rdi not in list_rd_pg:
            list_rd_pg.append(rdi)
            log.info(f"  RARW: Added RDI ({rdi}) to analysis list")

        get_classif_rarw = df_classif[
            (df_classif['parent_id'].isin(list_rd_pg)) | 
            (df_classif['child_id'].isin(list_rd_pg))
        ]['root'].drop_duplicates()

        result_pg = get_related_group(
            rdi, user_seed, rd_group_map, "RARW",
            df_pg_init, get_classif_rarw, list_rd_pg, df_classif
        )
        log.info(f"  RARW: Top {user_nb_top_rd} RDs span {len(get_classif_rarw)} classifications")

        # =====================================================================
        # Combine Results from Both Methods
        # =====================================================================
        
        concat_list = result_sm + result_pg
        df_classif_grp_rd = pd.DataFrame(concat_list)
        list_for_classif_global.append(df_classif_grp_rd)

    # -------------------------------------------------------------------------
    # Merge all patient results into a single DataFrame
    # -------------------------------------------------------------------------

    df_global_classif = pd.concat(list_for_classif_global, ignore_index=True)
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
    level = 8          # Maximum hierarchy level to consider
    n_top = 10         # Number of top candidates to analyze per patient

    # Build lookup dictionary: ORPHAcode -> Preferential Parent classification
    orpha_to_pp = (
        df7.dropna(subset=['ORPHAcode', 'Classif_id'])
        .drop_duplicates(['ORPHAcode'])
        .set_index('ORPHAcode')['Classif_id']
        .to_dict()
    )

    results = []

    # -------------------------------------------------------------------------
    # Process each patient-method combination
    # -------------------------------------------------------------------------
    for (patient, mtd), sub in df_global_classif.groupby(['patient', 'method']):
        # Get the RDI information for this patient-method combination
        is_rdi_sub = sub[sub['is_rdi'] == 'y']
        if is_rdi_sub.empty:
            continue

        RDI_rd_id = is_rdi_sub['rd_id'].iloc[0]
        RDI_rank = is_rdi_sub['rd_rank'].iloc[0]

        # Get the Preferential Parent classification for the RDI
        RDI_pp = orpha_to_pp.get(RDI_rd_id)
        if pd.isna(RDI_pp):
            continue

        # Filter RDI hierarchy to the specified level
        is_rdi_sub_level = is_rdi_sub[is_rdi_sub['level'] <= level]

        # Extract all "Group of disorders" ancestors of the RDI
        RDI_group_ids = (
            is_rdi_sub_level
            .loc[is_rdi_sub_level['direct_group'] == "Group of disorders", 'direct_parent']
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )
        RDI_group_set = set(RDI_group_ids)

        # Get the top-N distinct candidate diseases by rank
        topn_rd = (
            sub[['rd_id', 'rd_rank']]
            .drop_duplicates('rd_id')
            .sort_values('rd_rank', ascending=True)
            .head(n_top)['rd_id']
            .tolist()
        )

        # -----------------------------------------------------------------
        # Compare each candidate with the RDI
        # -----------------------------------------------------------------
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
            cand_group_ids = (
                cand_rows_level
                .loc[cand_rows_level['direct_group'] == "Group of disorders", 'direct_parent']
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

            # Calculate overlaps
            group_overlap = RDI_group_set.intersection(cand_group_ids)
            classif_overlap = {RDI_pp}.intersection({cand_pp})

            # Record the comparison result
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

    # -------------------------------------------------------------------------
    # Build summary DataFrame and identify matches
    # -------------------------------------------------------------------------
    summary_df = pd.DataFrame(results)

    # "Orange" condition: Both classification AND group overlap exist
    summary_df['is_orange'] = (
        (summary_df['classif_overlap'] > 0) & 
        (summary_df['group_overlap'] > 0)
    )

    # -------------------------------------------------------------------------
    # Count distinct matching groups per patient-method
    # -------------------------------------------------------------------------
    rows = []
    for (patient, method), sub in summary_df.groupby(['patient', 'method']):
        # Collect all unique matching group IDs (stored as sets, no duplicates)
        groups = set()
        for gs in sub.loc[sub['is_orange'], 'match_groups']:
            groups.update(gs)
        rows.append({
            'patient': patient,
            'method': method,
            'intensity': int(len(groups))
        })

    long_df = pd.DataFrame(rows)

    # -------------------------------------------------------------------------
    # Create pivot table for comparison (TABLE 5)
    # -------------------------------------------------------------------------
    pivot = (
        long_df
        .pivot(index='patient', columns='method', values='intensity')
        .fillna(0)
        .astype(int)
    )

    # Ensure consistent column ordering
    col_order = [m for m in ['RA', 'RARW'] if m in pivot.columns]
    pivot = pivot.reindex(columns=col_order)


    # Extract cases where methods disagree (TABLE 5 final result)
    df_diff = pivot[pivot['RA'] != pivot['RARW']].sort_values(
        by=['RARW'], 
        ascending=False
    )

    df_diff.to_excel(f"{PV.PATH_OUTPUT}/nb_gp_match_rdi_table5.xlsx")

    log.info(f"\nPatients with method disagreement: {len(df_diff)}")
    log.info("\n*** TABLE 5: Group Overlap Comparison ***")
    log.info(df_diff.head(20))


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

            # # # # with level
            mini_df_method =mini_df_method[mini_df_method['rd_group'] == "Disorder"]
            mini_df_method = mini_df_method[mini_df_method['level'] == 1]
            list_rdi_groups = mini_df_method[mini_df_method["rd_id"].isin(rdi_ids)]['direct_parent'].unique()
            rd_match_group_rdi = mini_df_method[mini_df_method['direct_parent'].isin(list_rdi_groups)]['rd_id'].unique()


            # get the df 
            df_match_group_rdi = mini_df_method[mini_df_method['rd_id'].isin(rd_match_group_rdi)].drop_duplicates()

            # get the rank of each RDs
            df_for_hm = df_match_group_rdi[['rd_id','rank']].drop_duplicates()

            # extract ranks
            ranks_gp = df_for_hm['rank']

    
            h_mean_rd_same_group_rdi = len(ranks_gp) / sum(1.0 / r for r in ranks_gp)
    
            log.info(f"{onep} - {onem}\tHarmonic mean (manual calculation): {h_mean_rd_same_group_rdi}")



            #for oner in ranks_gp:
            all_interaction.append((onep,onem,h_mean_rd_same_group_rdi))
    

        log.info(f"{onep} - {onem}\tHarmonic mean (manual calculation): {h_mean_rd_same_group_rdi}")


    df_hm_general = pd.DataFrame(all_interaction_level,columns=['patient','method','hm'])
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


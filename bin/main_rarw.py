#!/usr/bin/env python3
 
"""
Random Walk with Restart on DD+patient graph.
Two sub-commands:
  - run: Executes PageRank for one or more patients
  - collect: Collects results and creates the RDI file
"""
from pathlib import Path
import argparse
import logging

import pandas as pd
import networkx as nx
import os

import path_variable as PV
from set_log import setup_logging, get_logger


# =============================================================================
# HELPERS
# =============================================================================

def alpha_folder(alpha: float) -> str:
    """Converts alpha to string: 0.3 -> '0.3', 0.50 -> '0.5'"""
    s = f"{float(alpha):.2f}"
    return s.rstrip("0").rstrip(".")


# =============================================================================
# RUN COMMAND - PageRank for one patient
# =============================================================================
 

def cmd_run_rarw(args, log):
    """
    Executes PageRank for one or more patients.
    Uses NetworkX as in the original version.
    """
    # Parse seeds (can be "P1 P2 P3" or "P1,P2,P3")
    if " " in args.seeds:
        seeds = args.seeds.split()
    else:
        seeds = [s.strip() for s in args.seeds.split(",") if s.strip()]
    
    alpha = args.alpha

    matrix_dir =max(PV.PATH_OUTPUT_PATIENT_ADDED.glob("*/"), key=os.path.getmtime)
    matrix_subdir_file = str(matrix_dir).split("/")[-1]
    
    log.info("Seeds: %s", seeds)
    log.info("Alpha: %s", alpha)
    log.info("Matrix subdir: %s", matrix_dir)
    
    # Directory where patient-added matrices are stored
    if not matrix_dir.exists():
        raise FileNotFoundError(f"Matrix directory not found: {matrix_dir}")
    
    # Output directory for RW results
    alpha_str = alpha_folder(alpha)
    output_dir = Path(PV.PATH_OUTPUT_FOLDER_RW) / alpha_str / matrix_subdir_file
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if patient is already done
    existing = [f.stem for f in output_dir.glob("*.parquet")]
    # loop for patient 
    for seed in seeds:
        if seed in existing:
            log.info("Patient %s already done, skipping", seed)
            continue
        
        log.info("Processing patient: %s", seed)
        

        
        joint_path = matrix_dir / f"{seed}.parquet"
        if joint_path.exists():
            df_m = pd.read_parquet(joint_path)
            log.info("Loaded matrix from Parquet: %s (shape=%s)", joint_path, df_m.shape)
 
        
        # Display patient scores to RDs (top 5)
        if seed in df_m.index:
            patient_row = df_m.loc[seed].sort_values(ascending=False)
            log.info("Top 5 scores of patient %s to RDs:", seed)
            for rd, score in patient_row.head(5).items():
                log.info("  %s -> %.4f", rd, score)
        
        # Display total sum of degrees
        log.info("Total sum of matrix: %.2f", df_m.values.sum())
        log.info("=" * 60)
         
        # Build NetworkX graph
        G_raw = nx.from_pandas_adjacency(df_m)
        G_raw.remove_edges_from(nx.selfloop_edges(G_raw))
        
        # Normalize (row-stochastic)
        A = nx.adjacency_matrix(G_raw)
        df_adj = pd.DataFrame(
            A.toarray(),
            index=df_m.index,
            columns=df_m.columns
        )
        df_adj['tot'] = df_adj.sum(axis=1)
        df_adj['tot'] = df_adj['tot'].replace(0, 1)  # Avoid division by zero
        
        df_norm = df_adj.drop(columns=['tot']).div(df_adj['tot'], axis=0)
        
        G = nx.from_pandas_adjacency(df_norm)
        
        # Sum of normalized degrees (for ranking)
        sum_degres = df_norm.sum().sort_values(ascending=False)
        
        # Personalization vector: 1 for seed, 0 elsewhere
        personalization = {n: (1.0 if n == seed else 0.0) for n in G.nodes()}
        
        # PageRank
        pr = nx.pagerank(G, personalization=personalization, alpha=alpha)
        
        # Remove seed from results
        pr.pop(seed, None)
        
        # Build result DataFrame
        df_pr = pd.DataFrame({
            'pg': pd.Series(pr),
            'sum_degres': pd.Series(pr).index.map(lambda x: sum_degres.get(x, 0))
        })
        
        df_pr['rank_pg'] = df_pr['pg'].rank(ascending=False, method='min')
        df_pr['rank_sum_degres_pg'] = df_pr['sum_degres'].rank(ascending=False, method='min')
        
        # Save
        out_path = output_dir / f"{seed}.parquet"

        df_pr.to_parquet(out_path)
        log.info("Wrote %s ", out_path)
        
 

# =============================================================================
# COLLECT COMMAND - Aggregates results into RDI file
# =============================================================================

def cmd_aggregate_rdi(args, log):
    """
    Collects all PageRank results and creates the RDI file.
    Compares with confirmed diagnoses of patients.
    """

    matrix_subdir =max(PV.PATH_OUTPUT_PATIENT_ADDED.glob("*/"), key=os.path.getmtime)
    matrix_subdir_file = str(matrix_subdir).split("/")[-1]

    alpha_str = alpha_folder(args.alpha)
    
    # Directory of RARW results
    rarw_dir = Path(PV.PATH_OUTPUT_FOLDER_RW) / alpha_str / matrix_subdir_file
    if not rarw_dir.exists():
        raise FileNotFoundError(f"RARW directory not found: {rarw_dir}")
    
    log.info("Collecting results from: %s", rarw_dir)
    
    # List of xlsx files (PageRank results)
    result_files = list(rarw_dir.glob("*.parquet"))
    
    log.info("Found %d patient result files", len(result_files))
    
    # Load patient table with their confirmed diagnoses
    # Detect patient file automatically 
    patient_dir = Path(PV.PATH_OUTPUT_DF_PATIENT).parent
    candidate_files = [
        patient_dir / "patients.xlsx",
        patient_dir / "patients_subsumed.xlsx"
    ]
    patient_file = None
    for f in candidate_files:
        if f.exists():
            patient_file = f
            break
    df_patient = pd.read_excel(patient_file, index_col=0, engine="openpyxl")


    log.info("Loaded patient table: %d rows", len(df_patient))
    
    # Expected columns
    col_patient = PV.COL_DF_PATIENT_PATIENT
    col_disease = "Disease"  
    
    if col_disease not in df_patient.columns:
        log.warning("No Disease/Orphanet column found. Will collect all rank-1 results.")
        col_disease = None
    
    all_interactions = []
    
    for pq_path in result_files:
        patient_id = pq_path.stem
        
        try:
            # Load PageRank results
            df_rarw = pd.read_parquet(pq_path)
            
            # Rename index column if necessary
            if df_rarw.index.name is None or df_rarw.index.name == "Unnamed: 0":
                df_rarw.index.name = "RD"
            
 
            ## load patient input df depending on the patientid 
            patient_data = df_patient[df_patient[col_patient] == patient_id]
            if patient_data.empty:
                log.warning("Patient %s not found in patient table", patient_id)
                # Take rank 1 by default
                best_row = df_rarw.nsmallest(1, 'rank_pg').iloc[0]
                rdi = best_row.name
                score = best_row['pg']
                rank = int(best_row['rank_pg'])
            else:
                # extract the confirmed diagnosis
                rdi = patient_data[col_disease].values[0]
                
                # TEST if rdi is NaN
                if pd.isna(rdi):
                    log.warning("Patient %s has no confirmed diagnosis", patient_id)
                    continue
                
                rdi = str(rdi)

                
                # Find rank of this RDI
                if rdi in df_rarw.index:
                    row = df_rarw.loc[rdi]
                    score = row['pg']
                    rank = int(row['rank_pg'])
                else:
                    log.warning("RDI %s not found in results for patient %s", rdi, patient_id)
                    continue
 
            
            all_interactions.append((patient_id, rdi, score, rank))
            log.info("Patient %s: RDI=%s, rank=%d", patient_id, rdi, rank)
            
        except Exception as e:
            log.error("Error processing %s: %s", pq_path, e)
            continue
    
    if not all_interactions:
        log.error("No interactions collected!")
        return
    
    # Create final DataFrame
    df_rdi = pd.DataFrame(
        all_interactions,
        columns=['patients', 'RDs', 'score', 'rank']
    )
    
    # Statistics
    log.info("=" * 60)
    log.info("RDI STATISTICS:")
    log.info("Number of patients: %d", len(df_rdi))
    log.info("Mean rank: %.2f", df_rdi['rank'].mean())
    log.info("Median rank: %.1f", df_rdi['rank'].median())
    log.info("Patients with rank=1: %d (%.1f%%)", 
             (df_rdi['rank'] == 1).sum(),
             100 * (df_rdi['rank'] == 1).mean())
    log.info("Patients with rank<=5: %d (%.1f%%)",
             (df_rdi['rank'] <= 5).sum(),
             100 * (df_rdi['rank'] <= 5).mean())
    log.info("Patients with rank<=10: %d (%.1f%%)",
             (df_rdi['rank'] <= 10).sum(),
             100 * (df_rdi['rank'] <= 10).mean())
     
    # Save
    output_path =Path(PV.PATH_OUTPUT_FOLDER_RW) / alpha_str /"RDI.xlsx"
    df_rdi.to_excel(output_path, index=False, engine="openpyxl")
    log.info("Wrote RDI file: %s", output_path)
    
 


# =============================================================================
# MAIN
# =============================================================================

def build_parser():
    p = argparse.ArgumentParser(
        description="Random Walk with Restart on DD+patient graph."
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    
    # RUN command
    run_p = sub.add_parser("run", help="Run PageRank for patient(s)")
    run_p.add_argument("--alpha","-a",default=0.3,type=float,help="Alpha value used for the run")

    run_p.add_argument("--seeds", required=True,
                       help="Patient IDs (space or comma separated): 'P1 P2' or 'P1,P2'")
    
    # COLLECT command
    col_p = sub.add_parser("collect", help="Collect results into RDI file")
    col_p.add_argument("--alpha", "-a",default=0.3,type=float,help="Alpha value used for the run")
    
    return p




def main():
    # Logging configuration
    setup_logging(level=logging.INFO,console=False,filename=f"{Path(__file__).stem}.log"    )  
    log = get_logger(Path(__file__).stem)
    
    args = build_parser().parse_args()
    
    if args.cmd == "run":
        cmd_run_rarw(args, log)
    elif args.cmd == "collect":
        cmd_aggregate_rdi(args, log)


if __name__ == "__main__":
    main()
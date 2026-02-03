#!/usr/bin/env python3
# bin/main_rarw.py 
"""
Random Walk with Restart on MM+patient graph.
Two sub-commands:
  - run: Executes PageRank for one or more patients
  - collect: Collects results and creates the RDI file
"""
from pathlib import Path
import argparse
import time
import logging

import pandas as pd
import networkx as nx
import yaml

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

def cmd_run(args, log: logging.Logger):
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
    matrix_subdir = args.matrix_subdir
    
    log.info("Seeds: %s", seeds)
    log.info("Alpha: %s", alpha)
    log.info("Matrix subdir: %s", matrix_subdir)
    
    # Directory where patient-added matrices are stored
    matrix_dir = Path(PV.PATH_OUTPUT_PATIENT_ADDED) / matrix_subdir
    if not matrix_dir.exists():
        raise FileNotFoundError(f"Matrix directory not found: {matrix_dir}")
    
    # Output directory for RW results
    alpha_str = alpha_folder(alpha)
    output_dir = Path(PV.PATH_OUTPUT_FOLDER_RW) / alpha_str / matrix_subdir
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if patient is already done
    existing = [f.stem for f in output_dir.glob("*.xlsx")]
    
    for seed in seeds:
        if seed in existing:
            log.info("Patient %s already done, skipping", seed)
            continue
        
        log.info("Processing patient: %s", seed)
        t0 = time.perf_counter()
        
        # Search for patient CSV file
        # Can be either {seed}.csv or in by_patient/{seed}_to_RDs.parquet
        parquet_path = matrix_dir / "by_patient" / f"{seed}_to_RDs.parquet"
        
        joint_path = matrix_dir / f"{seed}.parquet"
        if joint_path.exists():
            df_m = pd.read_parquet(joint_path)
            log.info("Loaded matrix from Parquet: %s (shape=%s)", joint_path, df_m.shape)

        elif parquet_path.exists():
            # Parquet mode: patient vector -> must load MM and add patient
            log.info("Loading patient vector from parquet: %s", parquet_path)
            
            # Load patient vector
            df_vec = pd.read_parquet(parquet_path)
            if "RD" in df_vec.columns:
                df_vec = df_vec.set_index("RD")
            
            # Load MM matrix from PROVENANCE
            prov_path = matrix_dir / "PROVENANCE.yaml"
            if not prov_path.exists():
                raise FileNotFoundError(f"PROVENANCE.yaml not found in {matrix_dir}")
            
            prov = yaml.safe_load(prov_path.read_text(encoding="utf-8"))
            mm_path = prov["MM"]["path"]
            
            # Load MM
            if mm_path.endswith(".parquet"):
                df_mm = pd.read_parquet(mm_path)
            else:
                df_mm = pd.read_excel(mm_path, index_col=0, engine="openpyxl")
            
            log.info("Loaded MM matrix: %s (shape=%s)", mm_path, df_mm.shape)
            
            # Build matrix with patient added
            df_m = df_mm.copy()
            
            # Align patient vector on matrix RDs
            patient_scores = df_vec["score"].reindex(df_m.index).fillna(0.0)
            
            # Add patient as new row and column
            df_m[seed] = patient_scores.values
            df_m.loc[seed] = patient_scores.reindex(df_m.columns).fillna(0.0).values
            df_m.at[seed, seed] = 0.0  # No self-loop
            
            log.info("Built patient-added matrix (shape=%s)", df_m.shape)
        else:
            log.error("No matrix found for patient %s", seed)
            continue
        
        # # ===== TEMPORARY MATRIX DISPLAY =====
        # log.info("=" * 60)
        # log.info("PATIENT-ADDED MATRIX (preview):")
        # log.info("Shape: %s", df_m.shape)
        # log.info("Index (first 10): %s", list(df_m.index[:10]))
        # log.info("Columns (first 10): %s", list(df_m.columns[:10]))
        # log.info("Patient '%s' present in index: %s", seed, seed in df_m.index)
        # log.info("Patient '%s' present in columns: %s", seed, seed in df_m.columns)
        
        # Display patient scores to RDs (top 5)
        if seed in df_m.index:
            patient_row = df_m.loc[seed].sort_values(ascending=False)
            log.info("Top 5 scores of patient %s to RDs:", seed)
            for rd, score in patient_row.head(5).items():
                log.info("  %s -> %.4f", rd, score)
        
        # Display total sum of degrees
        log.info("Total sum of matrix: %.2f", df_m.values.sum())
        log.info("=" * 60)
        # ===== END TEMPORARY DISPLAY =====
        
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
        out_path = output_dir / f"{seed}.xlsx"
        df_pr.to_excel(out_path, engine='openpyxl')
        
        log.info("Wrote %s (%.1fs)", out_path, time.perf_counter() - t0)
        
        # Display top 5 results
        top5 = df_pr.nsmallest(5, 'rank_pg')
        log.info("Top 5 RDs for %s:", seed)
        for rd, row in top5.iterrows():
            log.info("  Rank %d: %s (score=%.6f)", int(row['rank_pg']), rd, row['pg'])


# =============================================================================
# COLLECT COMMAND - Aggregates results into RDI file
# =============================================================================

def cmd_collect(args, log: logging.Logger):
    """
    Collects all PageRank results and creates the RDI file.
    Compares with confirmed diagnoses of patients.
    """
    matrix_subdir = args.matrix_subdir
    alpha_str = alpha_folder(args.alpha)
    
    # Directory of RARW results
    rarw_dir = Path(PV.PATH_OUTPUT_FOLDER_RW) / alpha_str / matrix_subdir
    if not rarw_dir.exists():
        raise FileNotFoundError(f"RARW directory not found: {rarw_dir}")
    
    log.info("Collecting results from: %s", rarw_dir)
    
    # List of xlsx files (PageRank results)
    xlsx_files = list(rarw_dir.glob("*.xlsx"))
    xlsx_files = [f for f in xlsx_files if not f.name.startswith("RDI_") and not f.name.startswith("RARW_")]
    
    log.info("Found %d patient result files", len(xlsx_files))
    
    # Load patient table with their confirmed diagnoses
    df_patient = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT, index_col=0, engine="openpyxl")
    log.info("Loaded patient table: %d rows", len(df_patient))
    
    # Expected columns
    col_patient = PV.COL_DF_PATIENT_PATIENT
    col_disease = "Disease" if "Disease" in df_patient.columns else "Orphanet"
    
    if col_disease not in df_patient.columns:
        log.warning("No Disease/Orphanet column found. Will collect all rank-1 results.")
        col_disease = None
    
    all_interactions = []
    
    for xlsx_path in xlsx_files:
        patient_id = xlsx_path.stem
        
        try:
            # Load PageRank results
            df_rarw = pd.read_excel(xlsx_path, index_col=0, engine="openpyxl")
            
            # Rename index column if necessary
            if df_rarw.index.name is None or df_rarw.index.name == "Unnamed: 0":
                df_rarw.index.name = "RD"
            
            # Find confirmed diagnosis for this patient
            if col_disease:
                patient_data = df_patient[df_patient[col_patient] == patient_id]
                if patient_data.empty:
                    log.warning("Patient %s not found in patient table", patient_id)
                    # Take rank 1 by default
                    best_row = df_rarw.nsmallest(1, 'rank_pg').iloc[0]
                    rdi = best_row.name
                    score = best_row['pg']
                    rank = int(best_row['rank_pg'])
                else:
                    rdi = patient_data[col_disease].values[0]
                    
                    # Normalize RDI identifier
                    if pd.isna(rdi):
                        log.warning("Patient %s has no confirmed diagnosis", patient_id)
                        continue
                    
                    rdi = str(rdi)
                    if not rdi.startswith("ORPHA:") and rdi.replace(".", "").isdigit():
                        rdi = f"ORPHA:{int(float(rdi))}"
                    
                    # Find rank of this RDI
                    if rdi in df_rarw.index:
                        row = df_rarw.loc[rdi]
                        score = row['pg']
                        rank = int(row['rank_pg'])
                    else:
                        log.warning("RDI %s not found in results for patient %s", rdi, patient_id)
                        continue
            else:
                # No diagnosis column: take rank 1
                best_row = df_rarw.nsmallest(1, 'rank_pg').iloc[0]
                rdi = best_row.name
                score = best_row['pg']
                rank = int(best_row['rank_pg'])
            
            all_interactions.append((patient_id, rdi, score, rank))
            log.info("Patient %s: RDI=%s, rank=%d", patient_id, rdi, rank)
            
        except Exception as e:
            log.error("Error processing %s: %s", xlsx_path, e)
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
    log.info("=" * 60)
    
    # Save
    output_path = rarw_dir / f"RDI_{matrix_subdir}.xlsx"
    df_rdi.to_excel(output_path, index=False, engine="openpyxl")
    log.info("Wrote RDI file: %s", output_path)
    
    # Also in RARW_ format for compatibility
    output_path2 = rarw_dir.parent / f"RARW_{matrix_subdir}.xlsx"
    df_rdi.to_excel(output_path2, index=False, engine="openpyxl")
    log.info("Wrote RARW file: %s", output_path2)


# =============================================================================
# MAIN
# =============================================================================

def build_parser():
    p = argparse.ArgumentParser(
        description="Random Walk with Restart on MM+patient graph."
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    
    # RUN command
    run_p = sub.add_parser("run", help="Run PageRank for patient(s)")
    run_p.add_argument("--seeds", required=True,
                       help="Patient IDs (space or comma separated): 'P1 P2' or 'P1,P2'")
    run_p.add_argument("--alpha", type=float, required=True,
                       help="Damping factor for PageRank (e.g., 0.85)")
    run_p.add_argument("--matrix_subdir", required=True,
                       help="Subdirectory in patient_added containing the matrices")
    
    # COLLECT command
    col_p = sub.add_parser("collect", help="Collect results into RDI file")
    col_p.add_argument("--matrix_subdir", required=True,
                       help="Subdirectory name")
    col_p.add_argument("--alpha", type=float, required=True,
                       help="Alpha value used for the run")
    
    return p


def main():
    setup_logging(
        level=logging.INFO,
        console=True,  # Also display in console for debug
        filename=f"{Path(__file__).stem}.log"
    )
    log = get_logger(Path(__file__).stem)
    
    args = build_parser().parse_args()
    
    if args.cmd == "run":
        cmd_run(args, log)
    elif args.cmd == "collect":
        cmd_collect(args, log)


if __name__ == "__main__":
    main()
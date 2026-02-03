#!/usr/bin/env python3

"""
Script to add patients to an MM matrix.
Joint mode: creates a complete parquet matrix with patient included (for RARW)
Split mode: creates individual parquet vectors
"""
from pathlib import Path
import argparse
import os
import time
import datetime as dt
import logging

import pandas as pd
import yaml



from classes.utils import (
    abbr_combine, abbr_method, parse_agg_filename,
    split_csv_list, load_mm_matrix, load_mp_sm
)
import path_variable as PV
from set_log import setup_logging, get_logger




def resolve_input_path(given: str, base_dir: Path) -> str:
    """
    Resolves the input path:
    - If absolute and exists: uses as is
    - Otherwise: searches in base_dir
    - Otherwise: recursive search by filename
    """
    if os.path.isabs(given) and os.path.exists(given):
        return given
    
    direct = base_dir / given
    if direct.exists():
        return str(direct)
    
    # Recursive search
    matches = list(base_dir.rglob(Path(given).name))
    matches = [p for p in matches if p.is_file() and p.name == Path(given).name]
    if matches:
        matches.sort(key=lambda p: p.stat().st_mtime, reverse=True)
        return str(matches[0])
    
    raise FileNotFoundError(f"File not found: {given} (searched in {base_dir})")


def main():
    ap = argparse.ArgumentParser(
        description="Adds patients to an MM matrix (Disease-Disease)."
    )
    ap.add_argument("--output-mode", choices=["split", "joint"], default="joint",
                    help="joint: 1 CSV/patient with complete matrix (for RARW), split: parquet vectors")
    ap.add_argument("--pp-combine", choices=["funSimMax", "funSimAvg", "BMA", "BUMS", "rsd"], 
                    default=None, help="Combination method for Patient-Patient")
    ap.add_argument("--pp-method", choices=["resnik", "lin", "jiang", "graphic", "ic", "rel"], 
                    default=None, help="Similarity method for Patient-Patient")
    ap.add_argument("--outcsv", default=None, help="Output path (joint mode)")
    ap.add_argument("--mm-file", required=True, help="MM matrix file")
    ap.add_argument("--sm-file", required=True, help="MP similarity file (patients-RDs)")
    ap.add_argument("--patients", default=None, help="List of patients 'P1,P2' or '@file.txt'")
    ap.add_argument("--out-root", default=str(PV.PATH_OUTPUT_PATIENT_ADDED))
    
    args = ap.parse_args()
    
    # Logging with console=True to see matrices
    setup_logging(level=logging.INFO, console=True, filename=f"{Path(__file__).stem}.log")
    log = get_logger(Path(__file__).stem)
    
    log.info("=" * 70)
    log.info("ADD PATIENTS TO MM MATRIX")
    log.info("=" * 70)
    log.info("Mode: %s", args.output_mode)
    
    # Path resolution
    mm_path = resolve_input_path(args.mm_file, PV.PATH_OUTPUT_MM)
    sm_path = resolve_input_path(args.sm_file, PV.PATH_OUTPUT_SM)
    
    log.info("MM file: %s", mm_path)
    log.info("SM file: %s", sm_path)
    
    if not os.path.exists(mm_path):
        raise FileNotFoundError(f"MM file not found: {mm_path}")
    if not os.path.exists(sm_path):
        raise FileNotFoundError(f"SM file not found: {sm_path}")
    
    # Metadata parsing from filenames
    mm_meta = parse_agg_filename(Path(mm_path).stem)
    sm_meta = parse_agg_filename(Path(sm_path).stem)
    
    log.info("MM metadata: %s", mm_meta)
    log.info("SM metadata: %s", sm_meta)
    
    # Patient-Patient settings (for joint mode with multiple patients)
    pp_combine = args.pp_combine or sm_meta["combine"]
    pp_method = args.pp_method or sm_meta["method"]
    
    # Tags for output directory
    mm_tag = f"{abbr_combine(mm_meta['combine'])}_{abbr_method(mm_meta['method'])}_{mm_meta['vector_compact']}"
    mp_tag = f"{abbr_combine(sm_meta['combine'])}_{abbr_method(sm_meta['method'])}_{sm_meta['vector_compact']}"
    
    run_dir = Path(args.out_root) / f"mm_{mm_tag}_mp_{mp_tag}"
    run_dir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", run_dir)
    
    # Loading data
    log.info("-" * 50)
    log.info("Loading MM matrix...")
    df_mm = load_mm_matrix(mm_path, log)
    log.info("MM shape: %s", df_mm.shape)
    
    # ===== TEMPORARY DISPLAY OF MM MATRIX =====
    log.info("=" * 60)
    log.info("MM MATRIX PREVIEW:")
    log.info("Index (first 10): %s", list(df_mm.index[:10]))
    log.info("Columns (first 10): %s", list(df_mm.columns[:10]))
    
    # Mini table display
    sample_size = min(5, len(df_mm))
    sample_mm = df_mm.iloc[:sample_size, :sample_size]
    log.info("Mini table %dx%d:", sample_size, sample_size)
    for idx in sample_mm.index:
        row_str = " | ".join([f"{v:.4f}" for v in sample_mm.loc[idx]])
        log.info("  %s: %s", str(idx)[:20].ljust(20), row_str)
    
    log.info("Total MM sum: %.2f", df_mm.values.sum())
    log.info("Min/Max: %.4f / %.4f", df_mm.values.min(), df_mm.values.max())
    log.info("=" * 60)
    # ===== END MM DISPLAY =====
    
    log.info("Loading MP scores (patients-RDs)...")
    df_sm = load_mp_sm(sm_path, log)
    log.info("SM shape: %s", df_sm.shape)
    log.info("SM unique patients: %d", df_sm["patients"].nunique())
    log.info("SM unique RDs: %d", df_sm["RDs"].nunique())
    
    # Pivot: patients × RDs
    pivot = df_sm.pivot(index="patients", columns="RDs", values="score").fillna(0.0)
    log.info("Pivot shape: %s", pivot.shape)
    
    # Patient selection
    sel = split_csv_list(args.patients) or pivot.index.tolist()
    log.info("Selected patients: %d", len(sel))
    if len(sel) <= 10:
        log.info("Patients: %s", sel)
    else:
        log.info("Patients (first 10): %s ...", sel[:10])
    
    # RD verification
    rd_labels = list(df_mm.columns)
    unknown = [c for c in pivot.columns if c not in rd_labels]
    if unknown:
        log.warning("SM contains %d RD(s) absent from MM (ignored)", len(unknown))
    
    skipped, written = [], []
    
    # ==========================================================================
    # JOINT MODE: One CSV matrix per patient (for RARW)
    # ==========================================================================
    if args.output_mode == "joint":
        log.info("=" * 50)
        log.info("JOINT MODE: Creating one CSV matrix per patient")
        log.info("=" * 50)
        
        for pat in sel:
            if pat not in pivot.index:
                log.warning("Patient '%s' not found in SM, ignored", pat)
                skipped.append(pat)
                continue
            
            t0 = time.perf_counter()
            
            # Copy of MM matrix
            mat = df_mm.copy()
            
            # Patient scores towards RDs
            scores = pivot.loc[pat]
            scores_aligned = scores.reindex(mat.index).fillna(0.0)
            
            # Add column (RDs -> patient) and row (patient -> RDs)
            mat[pat] = scores_aligned.values
            mat.loc[pat] = scores_aligned.reindex(mat.columns).fillna(0.0).values
            mat.at[pat, pat] = 0.0  # No self-loop
            
            # ===== DISPLAY OF PATIENT-ADDED MATRIX =====
            log.info("-" * 50)
            log.info("PATIENT: %s", pat)
            log.info("Patient-added matrix shape: %s", mat.shape)
            log.info("Max patient->RDs score: %.4f", scores.max())
            log.info("Mean patient->RDs score: %.4f", scores.mean())
            log.info("RDs with score > 0: %d", (scores > 0).sum())
            
            # Top 5 RDs for this patient
            top5 = scores.nlargest(5)
            log.info("Top 5 RDs for %s:", pat)
            for rd, sc in top5.items():
                log.info("  %s: %.4f", rd, sc)
            
            # Verify that patient is in the matrix
            log.info("Patient '%s' in index: %s", pat, pat in mat.index)
            log.info("Patient '%s' in columns: %s", pat, pat in mat.columns)
            # ===== END DISPLAY =====
            
            # # Save as CSV
            # out_csv = run_dir / f"{pat}.csv"
            # mat.to_csv(out_csv)
            #  written.append(str(out_csv))
            out_pq = run_dir / f"{pat}.parquet"
            mat.to_parquet(out_pq)   # keeps index/columns
            written.append(str(out_pq))

            log.info("Saved: %s (%.2fs)", out_pq.name, time.perf_counter() - t0)
    
    # ==========================================================================
    # SPLIT MODE: One parquet vector per patient
    # ==========================================================================
    else:
        log.info("=" * 50)
        log.info("SPLIT MODE: Creating one parquet vector per patient")
        log.info("=" * 50)
        
        (run_dir / "by_patient").mkdir(parents=True, exist_ok=True)
        
        for pat in sel:
            if pat not in pivot.index:
                log.warning("Patient '%s' not found in SM, ignored", pat)
                skipped.append(pat)
                continue
            
            t0 = time.perf_counter()
            
            # Patient vector -> RDs (aligned on MM order)
            scores = pivot.loc[pat].reindex(rd_labels).fillna(0.0)
            
            # ===== VECTOR DISPLAY =====
            log.info("-" * 40)
            log.info("Patient: %s", pat)
            log.info("  Max score: %.4f", scores.max())
            log.info("  Mean score: %.4f", scores.mean())
            log.info("  RDs with score > 0: %d", (scores > 0).sum())
            top5 = scores.nlargest(5)
            log.info("  Top 5 RDs:")
            for rd, sc in top5.items():
                log.info("    %s: %.4f", rd, sc)
            # ===== END DISPLAY =====
            
            # Save as parquet
            out_vec = run_dir / "by_patient" / f"{pat}_to_RDs.parquet"
            df_vec = scores.rename("score").to_frame()
            df_vec.index.name = "RD"
            df_vec.to_parquet(out_vec)
            
            log.info("  Saved: %s (%.2fs)", out_vec.name, time.perf_counter() - t0)
            written.append(str(out_vec))
    
    # ==========================================================================
    # METADATA FILES
    # ==========================================================================
    
    # INDEX.csv
    index_data = []
    for p in written:
        pat_name = Path(p).stem.replace("_to_RDs", "")
        index_data.append({"patient": pat_name, "path": str(Path(p).resolve())})
    
    pd.DataFrame(index_data).to_csv(run_dir / "INDEX.csv", index=False)
    log.info("INDEX.csv saved: %s", run_dir / "INDEX.csv")
    
    # PROVENANCE.yaml
    prov = {
        "created_at": dt.datetime.now().astimezone().isoformat(),
        "output_dir": str(run_dir.resolve()),
        "output_mode": args.output_mode,
        "patients_requested": args.patients or "ALL",
        "patients_written": [Path(p).stem.replace("_to_RDs", "") for p in written],
        "patients_skipped": skipped,
        "MM": {
            "path": mm_path,
            **mm_meta,
            "abbr": {
                "combine": abbr_combine(mm_meta["combine"]), 
                "method": abbr_method(mm_meta["method"])
            },
        },
        "MP": {
            "path": sm_path,
            **sm_meta,
            "abbr": {
                "combine": abbr_combine(sm_meta["combine"]), 
                "method": abbr_method(sm_meta["method"])
            },
        },
        "PP": {
            "combine": pp_combine, 
            "method": pp_method
        },
    }
    
    (run_dir / "PROVENANCE.yaml").write_text(
        yaml.safe_dump(prov, sort_keys=False, allow_unicode=True),
        encoding="utf-8"
    )
    log.info("PROVENANCE.yaml saved")
    
    log.info("=" * 70)
    log.info("DONE: %d patients written, %d skipped", len(written), len(skipped))
    log.info("Output directory: %s", run_dir)
    log.info("=" * 70)


if __name__ == "__main__":
    main()
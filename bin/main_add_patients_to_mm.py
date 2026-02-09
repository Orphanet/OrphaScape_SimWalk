#!/usr/bin/env python3

"""
Script to add patients to an MM matrix.
Automatically uses the latest files from PATH_OUTPUT_MM and PATH_OUTPUT_SM.
"""
from pathlib import Path
import datetime as dt
import logging

import yaml

from classes.utils import (
    abbr_combine, abbr_method, parse_agg_filename,
    load_mm_matrix, load_mp_sm
)
import path_variable as PV
from set_log import setup_logging, get_logger



def get_latest_file(directory: Path, pattern: str = "*.parquet") -> Path:
    """
    Retrieves the most recent file in a directory.
    
    Args:
        directory: Directory to scan
        pattern: File pattern (default: *.parquet)
    
    Returns:
        Path of the most recent file
    
    Raises:
        FileNotFoundError: If no file found
    """
    directory = Path(directory)
    
    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")
    
    # Find all files matching the pattern
    files = [f for f in directory.glob(pattern) if f.is_file()]
    
    if not files:
        raise FileNotFoundError(f"No {pattern} files found in {directory}")
    
    # Sort by modification date (most recent first)
    files.sort(key=lambda f: f.stat().st_mtime, reverse=True)
    
    return files[0]


# =============================================================================
# MAIN
# =============================================================================

def main():
    # Logging
    setup_logging(level=logging.INFO, console=True, filename=f"{Path(__file__).stem}.log")
    log = get_logger(Path(__file__).stem)
    
    log.info("=" * 70)
    log.info("ADD PATIENTS TO MM MATRIX")
    log.info("=" * 70)
    
    # ==========================================================================
    # AUTO-DETECTION OF LATEST FILES
    # ==========================================================================
    
    try:
        mm_path = get_latest_file(PV.PATH_OUTPUT_MM)
        log.info("MM file (latest): %s", mm_path.name)
    except FileNotFoundError as e:
        log.error("MM error: %s", e)
        return
    
    try:
        sm_path = get_latest_file(PV.PATH_OUTPUT_SM)
        log.info("SM file (latest): %s", sm_path.name)
    except FileNotFoundError as e:
        log.error("SM error: %s", e)
        return
    
    # Parse metadata from filenames
    mm_meta = parse_agg_filename(mm_path.stem)
    sm_meta = parse_agg_filename(sm_path.stem)
    
    log.info("MM metadata: %s", mm_meta)
    log.info("SM metadata: %s", sm_meta)
    
    # Patient-Patient settings (inherited from SM)
    pp_combine = sm_meta["combine"]
    pp_method = sm_meta["method"]
    
    # Build output path
    mm_tag = f"{abbr_combine(mm_meta['combine'])}_{abbr_method(mm_meta['method'])}_{mm_meta['vector_compact']}"
    mp_tag = f"{abbr_combine(sm_meta['combine'])}_{abbr_method(sm_meta['method'])}_{sm_meta['vector_compact']}"
    
    run_dir = Path(PV.PATH_OUTPUT_PATIENT_ADDED) / f"mm_{mm_tag}_mp_{mp_tag}"
    run_dir.mkdir(parents=True, exist_ok=True)
    log.info("Output directory: %s", run_dir)
    
    # ==========================================================================
    # DATA LOADING
    # ==========================================================================
    
    log.info("-" * 50)
    log.info("Loading MM matrix...")
    df_mm = load_mm_matrix(str(mm_path), log)
    log.info("MM shape: %s", df_mm.shape)
    log.info("Total MM sum: %.2f", df_mm.values.sum())
    log.info("Min/Max: %.4f / %.4f", df_mm.values.min(), df_mm.values.max())
    
    log.info("Loading MP scores (patients-RDs)...")
    df_sm = load_mp_sm(str(sm_path), log)
    log.info("SM shape: %s", df_sm.shape)
    log.info("SM unique patients: %d", df_sm["patients"].nunique())
    log.info("SM unique RDs: %d", df_sm["RDs"].nunique())
    
    # Pivot: patients × RDs
    pivot = df_sm.pivot(index="patients", columns="RDs", values="score").fillna(0.0)
    log.info("Pivot shape: %s", pivot.shape)
    
    # Patient selection
    patient_sel = pivot.index.tolist()
    log.info("Selected patients: %d", len(patient_sel))
    if len(patient_sel) <= 10:
        log.info("Patients: %s", patient_sel)
    else:
        log.info("Patients (first 10): %s ...", patient_sel[:10])
    
    # RD verification
    rd_labels = list(df_mm.columns)
    unknown = [c for c in pivot.columns if c not in rd_labels]
    if unknown:
        log.warning("SM contains %d RD(s) absent from MM (ignored)", len(unknown))
    
    skipped, written = [], []
    
    # ==========================================================================
    # MATRIX CREATION (1 per patient)
    # ==========================================================================
    
    log.info("=" * 50)
    log.info("Creating one parquet matrix per patient")
    log.info("=" * 50)
    
    for pat in patient_sel:
        if pat not in pivot.index:
            log.warning("Patient '%s' not found in SM, ignored", pat)
            skipped.append(pat)
            continue
        
        # Copy MM matrix
        mat = df_mm.copy()
        
        # Patient scores towards RDs
        scores = pivot.loc[pat]
        scores_aligned = scores.reindex(mat.index).fillna(0.0)
        
        # Add column and row
        mat[pat] = scores_aligned.values
        mat.loc[pat] = scores_aligned.reindex(mat.columns).fillna(0.0).values
        mat.at[pat, pat] = 0.0  # No self-loop
        
        log.info("-" * 50)
        log.info("PATIENT: %s", pat)
        log.info("Matrix shape: %s", mat.shape)
        log.info("Max patient->RDs score: %.4f", scores.max())
        log.info("Mean patient->RDs score: %.4f", scores.mean())
        
        out_pq = run_dir / f"{pat}.parquet"
        mat.to_parquet(out_pq)
        written.append(str(out_pq))
        log.info("Saved: %s", out_pq.name)
    
    # ==========================================================================
    # METADATA FILE
    # ==========================================================================
    
    prov = {
        "created_at": dt.datetime.now().astimezone().isoformat(),
        "output_dir": str(run_dir.resolve()),
        "patients_written": [Path(p).stem for p in written],
        "patients_skipped": skipped,
        "MM": {
            "path": str(mm_path),
            **mm_meta,
        },
        "MP": {
            "path": str(sm_path),
            **sm_meta,
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
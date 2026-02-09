# classes/concat_sm.py — Optimized
"""
Aggregation of SM/MM files by RD into a single workbook.
- Unified glob patterns
- Parquet as preferred format
- Removal of commented code
"""
import glob
import logging
from pathlib import Path

import pandas as pd

from classes.utils import safe_to_excel


class ConcatSm:
    """
    Aggregates SM/RDI files by RD into a single workbook.
    
    For MP: base_dir/combine/method/pd4/vector_str/*.parquet
    For MM: same structure
    """
    
    def __init__(
        self,
        vector_str: str,
        col1: str,
        col2: str,
        base_dir: Path | str,
        pd4: str,
        combine: str | None,
        method: str | None,
        logger: logging.Logger,
        out_path: Path | str | None = None
    ):
        self.vector_str = str(vector_str)
        self.col1 = col1
        self.col2 = col2
        self.base_dir = Path(base_dir)
        self.pd4 = pd4
        self.combine = combine
        self.method = method
        self.log = logger
        self.out_path = Path(out_path) if out_path else None
    
    # =========================================================================
    # HELPERS
    # =========================================================================
    
    def _find_vector_dir(self) -> Path:
        """Finds the directory containing files for a given vector."""
        
        pattern = (
        self.base_dir 
        / (self.combine or "*") 
        / (self.method or "*") 
        / self.pd4 
        / self.vector_str
        )
        
        if "*" in str(pattern):
            matches = sorted(pattern.parent.glob(pattern.name))
            if not matches:
                raise RuntimeError(f"No directory matches pattern: {pattern}")
            return matches[0]
        
        if not pattern.exists():
            raise RuntimeError(f"Directory does not exist: {pattern}")
        return pattern
    
    # =========================================================================
    # MP AGGREGATION (Patients × RDs)
    # =========================================================================
    
    def concat_mp(self,path_patient,) -> None:
        """
        Aggregates MP files into a single Excel + Parquet.
        Also generates RDI file (rank 1 per patient).
        """
        
        in_dir = self._find_vector_dir()

        parquet_files = sorted(glob.glob(str(in_dir / "*.parquet")))
        
        if parquet_files:
            files, mode = parquet_files, "parquet"
        else:
            raise RuntimeError(f"No SM files found in {in_dir}")
        
        dfs_sm = []
        for fn in files:
            try:
                if mode == "parquet":
                    df = pd.read_parquet(fn)
                if df is not None and not df.empty:
                    dfs_sm.append(df)
            except Exception as e:
                self.log.warning("Error reading %s: %s", fn, e)
        
        if not dfs_sm:
            raise RuntimeError(f"No valid SM files in {in_dir}")
        
        # Concat, sort and ranking
        df_sm = (
            pd.concat(dfs_sm, ignore_index=True)
            .sort_values([self.col1, "score"], ascending=[True, False])
        )
        df_sm["rank"] = df_sm.groupby(self.col1)["score"].rank(
            ascending=False, method="dense"
        )
        
        # Output path
        if self.out_path:
            out_sm_xlsx = self.out_path
        else:
            tag = f"{self.combine}_{self.method}_{self.pd4}_{self.vector_str}.xlsx"
            out_sm_xlsx = self.base_dir / tag
        
        out_sm_xlsx = Path(out_sm_xlsx)
        
        # # Export Excel
        # safe_to_excel(df_sm, out_sm_xlsx, sheet_base="SM", log=self.log)
        
        # Export Parquet
        try:
            out_parquet = out_sm_xlsx.with_suffix(".parquet")
            df_sm.to_parquet(out_parquet, index=False)
            self.log.info("Wrote SM parquet → %s", out_parquet)
        except Exception as e:
            self.log.warning("Failed to write SM parquet: %s", e)
        
        
                
        # # RDI: best RD per patient (rank == 1)
        # df_rdi = df_sm[df_sm['rank'] == 1].copy()
        # RDI get the orpha of the patient
        # load input patient to get the RDI
        df_p = pd.read_excel(path_patient,  index_col=0)
        df_p.columns = ["patients","hpo_id","RDs",'type','Group']
        df_p = df_p[[ "patients","RDs"]].drop_duplicates()
        df_rdi = df_sm.merge(df_p, on=["patients","RDs"], how="inner")
    
        rdi_xlsx = out_sm_xlsx.with_name(f"RDI_{out_sm_xlsx.name}")
        safe_to_excel(df_rdi, rdi_xlsx, sheet_base="RDI", log=self.log)
    
        self.log.info("Wrote RDI to %s", rdi_xlsx)
    
    # =========================================================================
    # MM AGGREGATION (RDs × RDs - Matrix)
    # =========================================================================
    
 
    def concat_mm(self) -> None:
        """
        Aggregates partial MM matrices into a symmetric square matrix.
        Reads Parquet, pivots, symmetrizes, and writes Parquet + Excel.
        """

        in_dir = self._find_vector_dir()
        files = sorted(glob.glob(str(in_dir)+"/*.parquet"))
        
        if not files:
            raise RuntimeError(f"No MM matrices found with pattern {in_dir}")
        
        # Load and pivot each file
        mats = []
        for f in files:
            try:
                df = pd.read_parquet(f, columns=[self.col1, self.col2, "score"])
                if df is not None and not df.empty:
                    pivoted = df.pivot(
                        index=self.col1, 
                        columns=self.col2, 
                        values="score"
                    )
                    mats.append(pivoted)
            except Exception as e:
                self.log.warning("Skipping %s: %s", f, e)
        
        if not mats:
            raise RuntimeError(f"No valid MM matrices found")
        
        # Stack all partial rows
        M = pd.concat(mats, axis=0)
        
        # Global reindex
        labels = sorted(set(M.index) | set(M.columns))
        M = M.reindex(index=labels, columns=labels)
        
        # Symmetrize (fill the other half)
        M = M.combine_first(M.T).fillna(0.0)
        
        # Output paths
        if self.out_path:
            out_xlsx = Path(self.out_path)
        else:
            tag = f"{self.combine or 'comb'}_{self.method or 'sm'}_{self.pd4}_{self.vector_str}.xlsx"
            out_xlsx = self.base_dir / tag
        
        out_xlsx.parent.mkdir(parents=True, exist_ok=True)
        
        # Write Parquet  
        out_parquet = out_xlsx.with_suffix(".parquet")
        M.to_parquet(out_parquet)
        self.log.info("Wrote MM matrix (parquet) → %s", out_parquet)
        
        # # Write Excel
        # M.to_excel(out_xlsx, engine="openpyxl")
        # self.log.info("Wrote MM matrix → %s (shape=%s)", out_xlsx, M.shape)
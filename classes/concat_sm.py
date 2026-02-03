# classes/concat_sm.py — Optimized
"""
Aggregation of SM/MM files by RD into a single workbook.
- Unified glob patterns
- Parquet as preferred format
- Removal of commented code
"""
import glob
import time
import logging
from pathlib import Path

import pandas as pd

from classes.utils import safe_to_excel


class ConcatSm:
    """
    Aggregates SM/RDI files by RD into a single workbook.
    
    For MP: base_dir/combine/sm/{y|n}/product4/vector_str/*.parquet
    For MM: same structure
    """
    
    def __init__(
        self,
        vector_str: str,
        col1: str,
        col2: str,
        base_dir: Path | str,
        product4: str,
        combine: str | None,
        sm: str | None,
        logger: logging.Logger,
        out_path: Path | str | None = None
    ):
        self.vector_str = str(vector_str)
        self.col1 = col1
        self.col2 = col2
        self.base_dir = Path(base_dir)
        self.product4 = product4
        self.combine = combine
        self.sm = sm
        self.log = logger
        self.out_path = Path(out_path) if out_path else None
    
    # =========================================================================
    # HELPERS
    # =========================================================================
    
    def _find_vector_dir(self, weight_flag: str = "n") -> Path:
        """Finds the directory containing files for a given vector."""
        pattern = (
            self.base_dir 
            / (self.combine or "*") 
            / (self.sm or "*") 
            / weight_flag 
            / self.product4 
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
    
    def _load_files(self, in_dir: Path) -> tuple[list[pd.DataFrame], str]:
        """
        Loads all files from a directory (Parquet preferred, otherwise Excel).
        
        Returns:
            Tuple (list of DataFrames, mode: 'parquet' or 'xlsx')
        """
        parquet_files = sorted(glob.glob(str(in_dir / "*.parquet")))
        xlsx_files = sorted(glob.glob(str(in_dir / "*.xlsx")))
        
        if parquet_files:
            files, mode = parquet_files, "parquet"
        elif xlsx_files:
            files, mode = xlsx_files, "xlsx"
        else:
            raise RuntimeError(f"No SM files found in {in_dir}")
        
        dfs = []
        for fn in files:
            try:
                if mode == "parquet":
                    df = pd.read_parquet(fn)
                else:
                    df = pd.read_excel(fn, engine="openpyxl")
                
                if df is not None and not df.empty:
                    dfs.append(df)
            except Exception as e:
                self.log.warning("Error reading %s: %s", fn, e)
        
        return dfs, mode
    
    # =========================================================================
    # MP AGGREGATION (Patients × RDs)
    # =========================================================================
    
    def process_similarity(self) -> None:
        """
        Aggregates MP files into a single Excel + Parquet.
        Also generates RDI file (rank 1 per patient).
        """
        t0 = time.perf_counter()
        
        in_dir = self._find_vector_dir("n")
        dfs_sm, mode = self._load_files(in_dir)
        
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
            tag = f"{self.combine}_{self.sm}_n_{self.product4}_{self.vector_str}.xlsx"
            out_sm_xlsx = self.base_dir / tag
        
        out_sm_xlsx = Path(out_sm_xlsx)
        
        # Export Excel
        safe_to_excel(df_sm, out_sm_xlsx, sheet_base="SM", log=self.log)
        
        # Export Parquet
        try:
            out_parquet = out_sm_xlsx.with_suffix(".parquet")
            df_sm.to_parquet(out_parquet, index=False)
            self.log.info("Wrote SM parquet → %s", out_parquet)
        except Exception as e:
            self.log.warning("Failed to write SM parquet: %s", e)
        
        self.log.info("Wrote SM to %s (%.1fs)", out_sm_xlsx, time.perf_counter() - t0)
        
        # RDI: best RD per patient (rank == 1)
        df_rdi = df_sm[df_sm['rank'] == 1].copy()
        rdi_xlsx = out_sm_xlsx.with_name(f"RDI_{out_sm_xlsx.name}")
        safe_to_excel(df_rdi, rdi_xlsx, sheet_base="RDI", log=self.log)
        self.log.info("Wrote RDI to %s", rdi_xlsx)
    
    # =========================================================================
    # MM AGGREGATION (RDs × RDs - Matrix)
    # =========================================================================
    
    def concat_matrix(self) -> None:
        """
        Aggregates partial MM matrices into a symmetric square matrix.
        Reads Parquet, pivots, symmetrizes, and writes Parquet + Excel.
        """
        # Pattern to find files
        pattern = str(
            self.base_dir 
            / (self.combine or "*") 
            / (self.sm or "*") 
            / "*" 
            / self.product4 
            / self.vector_str 
            / "*.parquet"
        )
        files = sorted(glob.glob(pattern))
        
        if not files:
            raise RuntimeError(f"No MM matrices found with pattern {pattern}")
        
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
            tag = f"{self.combine or 'comb'}_{self.sm or 'sm'}_n_{self.product4}_{self.vector_str}.xlsx"
            out_xlsx = self.base_dir / tag
        
        out_xlsx.parent.mkdir(parents=True, exist_ok=True)
        
        # Write Parquet (recommended)
        out_parquet = out_xlsx.with_suffix(".parquet")
        M.to_parquet(out_parquet)
        self.log.info("Wrote MM matrix (parquet) → %s", out_parquet)
        
        # Write Excel
        M.to_excel(out_xlsx, engine="openpyxl")
        self.log.info("Wrote MM matrix → %s (shape=%s)", out_xlsx, M.shape)
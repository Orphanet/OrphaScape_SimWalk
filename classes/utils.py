
"""
Centralized utility module to eliminate redundancies.
Contains: abbreviations, ORPHA normalization,
          path management, file reading/writing, etc.
"""

import re

import logging
from pathlib import Path
from typing import Any
import pandas as pd
import math
# =============================================================================
# CENTRALIZED CONSTANTS
# =============================================================================

# Mapping HPO frequencies -> scores
FREQ_MAP: dict[str, float] = {
    "Obligate": 1.0,
    "Very frequent": 0.8,
    "Frequent": 0.6,
    "Occasional": 0.4,
    "Very rare": 0.2,
    "Excluded": 0.0,
}

# Abbreviations for combine methods
COMB_ABBR: dict[str, str] = {
    "funsimmax": "fsm",
    "funsimavg": "fsa",
    "bma": "bma",
    "rsd": "rsd",
    "bums": "bums",
}

# Abbreviations for similarity methods
METH_ABBR: dict[str, str] = {
    "resnik": "r",
    "lin": "l",
    "jiang": "jc",
    "jc": "jc",
    "graphic": "g",
    "ic": "ic",
    "rel": "rel",
}
 

# =============================================================================
# ABBREVIATION FUNCTIONS
# =============================================================================

def abbr_combine(name: str) -> str:
    """Converts a combine method name to abbreviation."""
    return COMB_ABBR.get(str(name).lower(), str(name)[:3])


def abbr_method(name: str) -> str:
    """Converts a similarity method name to abbreviation."""
    return METH_ABBR.get(str(name).lower(), str(name)[:2])


# =============================================================================
# PARSING FUNCTIONS
# =============================================================================
 


def alpha_folder(alpha: float) -> str:
    """Converts alpha to folder string: 0.3 -> '0.3', 0.50 -> '0.5'"""
    s = f"{float(alpha):.2f}"
    return s.rstrip("0").rstrip(".")


def slug(x: str) -> str:
    """Cleans a string for use in a path."""
    return re.sub(r"[^A-Za-z0-9._-]+", "-", str(x))


def resolve_product4_segment(dd_p4: str | None, dp_p4: str | None) -> str:
    """Resolves the product4 segment for combined DD and DP."""
    dd_p4 = (dd_p4 or "").strip()
    dp_p4 = (dp_p4 or "").strip()
    if dd_p4 and dp_p4 and dd_p4 == dp_p4:
        return dd_p4
    if dd_p4 or dp_p4:
        return f"dd{slug(dd_p4)}__dp{slug(dp_p4)}"
    return "product4_unknown"


# =============================================================================
# NORMALIZATION FUNCTIONS
# =============================================================================

def norm_orpha(x: Any) -> str:
    """
    Normalizes an ORPHA identifier:
    - '100985' -> 'ORPHA:100985'
    - 'ORPHA_1234' -> 'ORPHA:1234'
    - '100985.0' (Excel float) -> 'ORPHA:100985'
    - 'ORPHA:100985' -> 'ORPHA:100985' (unchanged)
    """
    if x is None:
        return ""
    s = str(x).strip()
    
    # Handle Excel floats: "100985.0"
    if re.fullmatch(r"\d+(\.0+)?", s):
        s = s.split(".")[0]
    
    # ORPHA_1234 -> ORPHA:1234
    s = s.replace("ORPHA_", "ORPHA:")
    
    # "1234" -> "ORPHA:1234"
    if s.isdigit():
        s = f"ORPHA:{s}"
    
    return s


def orpha_num(x: str) -> int:
    """Extracts the number from an ORPHA identifier."""
    m = re.search(r"(\d+)$", str(x))
    return int(m.group(1)) if m else 0


def freq_to_score(text: str | None) -> float | None:
    """Converts HPO frequency text to numeric score."""
    if not text:
        return None
    for key, val in FREQ_MAP.items():
        if key in text:
            return val
    return None


# =============================================================================
# WEIGHT VECTOR PARSING FUNCTIONS
# =============================================================================

def parse_weights(vector_str: Any, expected: int = 5, default: float = 1.0) -> list[float]:
    """
    Parses a weight vector from different formats:
    - "3_2_2_2_1" (underscore-separated)
    - "32221" (compact)
    - ["3_2_2_2_1"] (YAML list)
    
    Returns: List of floats of length `expected`
    """
    # Unwrap list/tuple
    if isinstance(vector_str, (list, tuple)):
        if len(vector_str) != 1:
            raise ValueError(f"vector_str list must have length 1, got {vector_str}")
        vector_str = vector_str[0]
    
    s = str(vector_str).strip()
    
    # Underscore format
    if "_" in s:
        parts = [p for p in s.split("_") if p]
    # Compact digit format
    elif s.isdigit():
        s = s.zfill(expected)[:expected]
        parts = list(s)
    else:
        raise ValueError(
            f"vector_str='{vector_str}' invalid. "
            f"Use 'a_b_c_d_e', 'abcde', or ['a_b_c_d_e']."
        )
    
    # Convert to floats
    try:
        vec = [float(p) for p in parts]
    except ValueError as e:
        raise ValueError(f"Invalid weight value in vector_str={vector_str}: {parts}") from e
    
    # Pad/trim to expected length
    if len(vec) < expected:
        vec += [default] * (expected - len(vec))
    elif len(vec) > expected:
        vec = vec[:expected]
    
    return vec


def to_compact(v: str) -> str:
    """Converts '1_1_1_1_1' -> '11111'"""
    return v.replace("_", "")


# =============================================================================
# FILE READING/WRITING FUNCTIONS
# =============================================================================

def split_csv_list(s: str | None) -> list[str]:
    """
    Parses a CSV list or file:
    - 'A,B,C' -> ['A', 'B', 'C']
    - '@file.txt' -> lines from file
    """
    if not s:
        return []
    s = str(s).strip()
    if s.startswith("@"):
        path = s[1:]
        with open(path, "r", encoding="utf-8") as f:
            return [line.strip() for line in f if line.strip()]
    return [x.strip() for x in s.split(",") if x.strip()]


def ensure_parent(p: Path) -> None:
    """Creates the parent directory if it doesn't exist."""
    p.parent.mkdir(parents=True, exist_ok=True)


def pick_excel_sheet(xl: pd.ExcelFile, preferred: str) -> str:
    """
    Robustly chooses an Excel sheet:
    - Exact match (case-insensitive)
    - Or unique sheet if only one exists
    - Otherwise error with available names
    """
    pref_low = preferred.lower()
    for s in xl.sheet_names:
        if s.strip().lower() == pref_low:
            return s
    if len(xl.sheet_names) == 1:
        return xl.sheet_names[0]
    raise ValueError(
        f"Worksheet '{preferred}' not found. Available: {xl.sheet_names}"
    )


def load_dataframe(path: str, log: logging.Logger | None = None) -> pd.DataFrame:
    """
    Loads a DataFrame from Parquet, Excel, or CSV.
    Automatically detects the format.
    """
    path_lower = str(path).lower()
    
    try:
        if path_lower.endswith(".parquet"):
            return pd.read_parquet(path)
        elif path_lower.endswith(".csv"):
            return pd.read_csv(path)
        elif path_lower.endswith(".csv.gz"):
            return pd.read_csv(path, compression="gzip")
        elif path_lower.endswith((".xlsx", ".xls")):
            return pd.read_excel(path, engine="openpyxl")
        else:
            # Try parquet first, then excel
            try:
                return pd.read_parquet(path)
            except Exception:
                return pd.read_excel(path, engine="openpyxl")
    except Exception as e:
        if log:
            log.error("Cannot read '%s': %s", path, e)
        raise


def load_dd_matrix(path: str, log: logging.Logger) -> pd.DataFrame:
    """Loads an DD matrix (square format) from Parquet or Excel."""
    path_lower = str(path).lower()
    
    try:
        if path_lower.endswith(".parquet"):
            df = pd.read_parquet(path)
            if df is None or df.empty:
                raise ValueError(f"Parquet empty: {path}")
            return df
        
        # Fallback to Excel
        xl = pd.ExcelFile(path, engine="openpyxl")
        sheet = pick_excel_sheet(xl, "DD")
        df = xl.parse(sheet, index_col=0)
        if df is None or df.empty:
            raise ValueError(f"Sheet '{sheet}' is empty in {path}")
        return df
    except Exception as e:
        log.error("Cannot read DD matrix '%s': %s", path, e)
        raise


def load_dp_sm(path: str, log: logging.Logger) -> pd.DataFrame:
    """
    Loads an SM long table (patients, RDs, score).
    Supports Parquet, CSV, Excel (multi-sheet SM_*).
    """
    path_lower = str(path).lower()
    needed_cols = {"patients", "RDs", "score"}
    
    try:
        # CSV
        if path_lower.endswith(".csv") or path_lower.endswith(".csv.gz"):
            df = pd.read_csv(path)
            if df is None or df.empty:
                raise ValueError(f"CSV empty: {path}")
            missing = needed_cols - set(df.columns)
            if missing:
                raise ValueError(f"CSV missing columns: {sorted(missing)}")
            return df[list(needed_cols)].copy()
        
        # Parquet
        if path_lower.endswith(".parquet"):
            df = pd.read_parquet(path)
            if df is None or df.empty:
                raise ValueError(f"Parquet empty: {path}")
            missing = needed_cols - set(df.columns)
            if missing:
                raise ValueError(f"Parquet missing columns: {sorted(missing)}")
            return df[list(needed_cols)].dropna(how="any").drop_duplicates()
        
        # Excel multi-sheet
        xl = pd.ExcelFile(path, engine="openpyxl")
        
        # Look for exact 'SM'
        for s in xl.sheet_names:
            if s.strip().lower() == "sm":
                df = xl.parse(s)
                break
        else:
            # Concat SM, SM_1, SM_2... sheets
            sm_sheets = [s for s in xl.sheet_names 
                        if re.match(r"^sm(_\d+)?$", s.strip().lower())]
            if sm_sheets:
                frames = [xl.parse(s) for s in sm_sheets]
                frames = [f for f in frames if f is not None and not f.empty]
                if not frames:
                    raise ValueError(f"SM sheets empty in {path}")
                df = pd.concat(frames, ignore_index=True)
            elif len(xl.sheet_names) == 1:
                df = xl.parse(xl.sheet_names[0])
            else:
                raise ValueError(f"'SM' sheet not found. Available: {xl.sheet_names}")
        
        if df is None or df.empty:
            raise ValueError(f"SM sheet(s) empty in {path}")
        
        missing = needed_cols - set(df.columns)
        if missing:
            raise ValueError(f"SM missing columns: {sorted(missing)}")
        
        return df[list(needed_cols)].dropna(how="any").drop_duplicates()
    
    except Exception as e:
        log.error("Cannot read DP SM '%s': %s", path, e)
        raise


# =============================================================================
# EXPORT FUNCTIONS
# =============================================================================

def safe_to_excel(
    df: pd.DataFrame,
    out_xlsx: Path,
    sheet_base: str = "Sheet",
    max_rows: int = 1_000_000,
    log: logging.Logger | None = None
) -> None:
    """
    Writes a DataFrame to Excel with large file handling (multi-sheet).
    Also writes a CSV.gz version in parallel.
    """
    
    
    out_xlsx = Path(out_xlsx)
    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    
    n = len(df)
    with pd.ExcelWriter(out_xlsx, engine="openpyxl") as writer:
        if n == 0:
            pd.DataFrame(columns=df.columns).to_excel(
                writer, index=False, sheet_name=sheet_base
            )
        else:
            chunks = math.ceil(n / max_rows)
            for i in range(chunks):
                start = i * max_rows
                end = min((i + 1) * max_rows, n)
                sheet_name = f"{sheet_base}_{i+1}" if chunks > 1 else sheet_base
                df.iloc[start:end].to_excel(writer, index=False, sheet_name=sheet_name)



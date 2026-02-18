
"""
Calculation of semantic similarity measures between HPO.
- Global cache for HPO objects and similarity scores
- Removal of unnecessary S matrices (streaming max)
- Triangularization for DD (symmetric calculation)
- Type hints and docstrings
"""
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from pyhpo import Ontology, HPOSet

import path_variable as PV
from classes.utils import orpha_num
from set_log import get_logger

# Load HPO ontology only once at module level
Ontology(str(PV.PATH_INPUT_HPO), transitive=True)
log = get_logger("sim_measure_cal")
log.info("HPO Ontology version: %s", Ontology.version())


# =============================================================================
# GLOBAL CACHE
# =============================================================================
_HPO_OBJ_CACHE: dict[str, Optional[object]] = {}  # hpo_id -> Ontology object
_SIM_CACHE: dict[tuple[str, str, str], float] = {}  # (hpoA, hpoB, method) -> score


def _get_hpo_obj(hpo_id: str):
    """Gets an HPO object with caching.
        exemple : hpo_id = "HP:0001250 """
    if hpo_id in _HPO_OBJ_CACHE:
        return _HPO_OBJ_CACHE[hpo_id]
    try:
        obj = Ontology.get_hpo_object(hpo_id)
    except RuntimeError:
        obj = None
    # store it in the cache in dict format 
    _HPO_OBJ_CACHE[hpo_id] = obj
    return obj


def _sim_cached(hpo_a: str, hpo_b: str, method: str) -> float:
    """Computes similarity between two HPO with caching (symmetric).
    example : hpo_a = "HP:0001250"  hpo_b = "HP:0001251 method = resnik"""
    # Canonicalization of pair (symmetry)
    key = (hpo_a, hpo_b, method) if hpo_a <= hpo_b else (hpo_b, hpo_a, method)
    
    if key in _SIM_CACHE:
        return _SIM_CACHE[key]
    
    obj_a = _get_hpo_obj(key[0])
    obj_b = _get_hpo_obj(key[1])
    
    if obj_a is None or obj_b is None:
        _SIM_CACHE[key] = 0.0
        return 0.0
    
    try:
        v = float(obj_a.similarity_score(obj_b, 'orpha', method))
    except RuntimeError:
        v = 0.0
    
    _SIM_CACHE[key] = v
    return v


# =============================================================================
# MAIN CLASS
# =============================================================================

class Sim_measure:
    """Semantic similarity calculation between entities (RD-RD, RD-Patient, Patient-Patient)."""
    
    # Mapping frequency -> index in weight vector
    FREQ_TO_INDEX = {
        1.0: 0,   # Obligate
        0.8: 1,   # Very frequent
        0.6: 2,   # Frequent
        0.4: 3,   # Occasional
        0.2: 4,   # Very rare
    }
    
    def __init__(
        self,
        df_group1: pd.DataFrame,
        df_group2: pd.DataFrame,
        colname_1: str,
        colname_2: str,
        logger: Optional[logging.Logger] = None
    ):
        """
        Args:
            df_group1: DataFrame of group 1 (patients or RDs)
            df_group2: DataFrame of group 2 (RDs)
            colname_1: Column name identifying entities in group 1
            colname_2: Column name identifying entities in group 2
        """
        self.df_group1 = df_group1
        self.df_group2 = df_group2
        self.colname_1 = colname_1
        self.colname_2 = colname_2
        self.logger = logger or log
    
    # =========================================================================
    # FREQUENCY WEIGHTS
    # =========================================================================
    
    def _get_freq_weight(
        self,
        hpo_obj,
        mini_df: pd.DataFrame,
        vector_weight: list[float]
    ) -> float:
        """
        Computes weight based on HPO frequency.
        
        Args:
            hpo_obj: HPO object from hpo3 eg Ontology.get_hpo_object("HP:0001250")
            mini_df: Filtered DataFrame for an entity
            eg mini_df = pd.DataFrame({
                            'hpo_id': ['HP:0001250'],
                            'hpo_frequency': [0.8]
                            })
            vector_weight: Weight vector [obligate, very_freq, freq, occ, very_rare] [1.0, 1.0, 1.0, 1.0, 1.0]
        """
        try:
            freq = mini_df[mini_df['hpo_id'] == hpo_obj.id]['hpo_frequency'].values[0]
            weight_idx = self.FREQ_TO_INDEX.get(freq, None)
            add_weight = vector_weight[weight_idx] if weight_idx is not None else 1.0
        except (KeyError, IndexError):
            freq = 1.0
            add_weight = 1.0
        
        return  add_weight
    
    # =========================================================================
    # SCORE COMBINATION
    # =========================================================================
    
    @staticmethod
    def _combine_scores(
        max_row: np.ndarray,
        max_col: np.ndarray,
        len_a: int,
        len_b: int,
        combine: str
    ) -> float:
        """
        Combines max scores of rows/columns according to chosen method.
        
        max_row: best match for each HPO of B (size len_b) 
        eg : np.array([0.8, 0.6, 0.9])  # 3 HPOs len b
        max_col: best match for each HPO of A (size len_a)
        eg :  np.array([0.7, 0.85])  # 2 HPOs len a
        """
        if combine == "funSimAvg":
            avg_row = max_row.mean() if len_b else 0.0
            avg_col = max_col.mean() if len_a else 0.0
            return (avg_row + avg_col) / 2.0
        
        elif combine == "funSimMax":
            avg_row = max_row.mean() if len_b else 0.0
            avg_col = max_col.mean() if len_a else 0.0
            return max(avg_row, avg_col)
        
        elif combine == "BMA":
            denom = (len_b + len_a) or 1
            return (max_row.sum() + max_col.sum()) / denom
        
        elif combine == "rsd":
            # RSD asymmetric: sum(max_col) / len_a
            return max_col.sum() / (len_a or 1)
        
        else:
            raise ValueError(f"Unknown combine method: {combine}")
    
    # =========================================================================
    # DD CALCULATION (RD × RD) - OPTIMIZED
    # =========================================================================
    
    def run_dd_freq(
        self,
        element2: str,
        rd_id_list: list[str],
        combine: str,
        method: str,
        vector_weight: list[float]
    ) -> pd.DataFrame:
        """
        Calculates RD (element2) ↔ RDs (rd_id_list) similarity.
        eg:element2 = "ORPHA:123", rd_id_list = ["ORPHA:123", "ORPHA:456", "ORPHA:789"]
        combine="funSimMax", method="resnik", vector_weight=[1.0, 1.0, 1.0, 1.0, 1.0]

        - Cache of HPO similarities
        - Triangularization (computes only element1 >= element2)
        """
        interactions = []
        
        # HPOs of element2
        minidf_2 = self.df_group2[self.df_group2[self.colname_2] == element2]
        hpo_el_2 = minidf_2['hpo_id'].drop_duplicates().tolist()
        len_b = len(hpo_el_2)
        
        # Numeric index for triangularization 
        e2n = orpha_num(element2) # get the numeric part of ORPHAcode
        
        for element1 in rd_id_list:
            # Triangularization: compute only if element1 >= element2 because matrix is symmetric
            if orpha_num(element1) < e2n:
                continue
            
            # Diagonal
            if element2 == element1:
                score = 0.0
            else:
                minidf_1 = self.df_group1[self.df_group1[self.colname_1] == element1]
                hpo_el_1 = minidf_1['hpo_id'].drop_duplicates().tolist()
                len_a = len(hpo_el_1)
                
                # Streaming max (no S matrix)
                max_row = np.zeros(len_b, dtype=float)
                max_col = np.zeros(len_a, dtype=float)
                
                for i, h2 in enumerate(hpo_el_2):
                    obj2 = _get_hpo_obj(h2)
                    if obj2 is None:
                        continue
                    
                    w2 = self._get_freq_weight(obj2, minidf_2, vector_weight)
                    best_i = 0.0
                    
                    for j, h1 in enumerate(hpo_el_1):
                        # compute sm symmetric
                        ic = _sim_cached(h2, h1, method)
                        s = ic * w2
                        
                        if s > best_i:
                            best_i = s
                        if s > max_col[j]:
                            max_col[j] = s
                    
                    max_row[i] = best_i
                
                score = self._combine_scores(max_row, max_col, len_a, len_b, combine)
            
            interactions.append((element2, element1, float(score)))
        
        return pd.DataFrame(interactions, columns=['RDs', 'patients', 'score'])
    
    # =========================================================================
    # DP CALCULATION (Patient × RD) - OPTIMIZED
    # =========================================================================
    
    def run_sm_freq(
        self,
        element2: str,
        patient_id_list: list[str],
        combine: str,
        method: str,
        vector_weight: list[float]
    ) -> pd.DataFrame:
        """
        Calculates RD (element2) ↔ patients (patient_id_list) similarity.
        eg : element2 = "ORPHA:123", patient_id_list = ["P1", "P2", "P3"]
        combine="funSimMax", method="resnik", vector_weight=[1.0, 1.0, 1.0, 1.0, 1.0]

        Handles exclusion rule (HPO with frequency = 0), 
        if patient have hpo_frequency=0 the result is 0 between patient and RD.(avoid useless calculation)
        """
        interactions = []
        
        # HPOs of element2 (RD)
        minidf_2 = self.df_group2[self.df_group2[self.colname_2] == element2]
        hpo_el_2 = minidf_2['hpo_id'].drop_duplicates().tolist()
        hpo_excluded = set(
            minidf_2[minidf_2['hpo_frequency'] == 0]['hpo_id'].drop_duplicates()
        )
        
        for element1 in patient_id_list:
            if element2 == element1:
                score = 0.0
            else:
                minidf_1 = self.df_group1[self.df_group1[self.colname_1] == element1]
                hpo_el_1 = minidf_1['hpo_id'].drop_duplicates().tolist()
                
                # Exclusion rule
                if hpo_excluded.intersection(hpo_el_1):
                    score = 0.0
                else:
                    len_a = len(hpo_el_1)
                    len_b = len(hpo_el_2)
                    
                    # Matrix S (necessary for patients as not symmetric)
                    S = np.zeros((len_b, len_a), dtype=float)
                    
                    for i, h2 in enumerate(hpo_el_2):
                        try:
                            obj2 = Ontology.get_hpo_object(h2)
                            w2 = self._get_freq_weight(obj2, minidf_2, vector_weight)
                        except RuntimeError:
                            continue
                        
                        for j, h1 in enumerate(hpo_el_1):
                            try:
                                obj1 = Ontology.get_hpo_object(h1)
                                ic = obj2.similarity_score(obj1, 'orpha', method)
                                S[i, j] = ic * w2
                            except RuntimeError:
                                pass
                    
                    max_row = np.amax(S, axis=1) if S.size else np.array([])
                    max_col = np.amax(S, axis=0) if S.size else np.array([])
                    
                    # Special handling for BUMS
                    score = self._combine_scores(max_row, max_col, len_a, len_b, combine)
            
            interactions.append((element2, element1, float(score)))
        
        return pd.DataFrame(interactions, columns=['RDs', 'patients', 'score'])


    # =========================================================================
    # PP CALCULATION (Patient × Patient)
    # =========================================================================
    
    def run_pp(
        self,
        patients_list: list[str],
        combine: str,
        method: str
    ) -> pd.DataFrame:
        """
        Calculates Patient ↔ Patient similarity (without frequencies).
        eg :patients_list = ["P001", "P002", "P003"]
        Returns a symmetric square matrix.
        """
        # Unique while preserving order
        pats = list(dict.fromkeys(patients_list))
        n = len(pats)
        
        # Map patient -> list of HPO objects
        hpo_by_pat: dict[str, list] = {}
        for p in pats:
            tmp = self.df_group1[self.df_group1[self.colname_1] == p]
            hpo_ids = tmp['hpo_id'].drop_duplicates().tolist()
            hpo_by_pat[p] = [_get_hpo_obj(h) for h in hpo_ids]
            hpo_by_pat[p] = [o for o in hpo_by_pat[p] if o is not None]
        
        M = np.zeros((n, n), dtype=float)
        
        for i, p1 in enumerate(pats):
            h1_list = hpo_by_pat.get(p1, [])
            len_a = len(h1_list)
            
            for j in range(i, n):
                if i == j:
                    M[i, j] = 0.0
                    continue
                
                p2 = pats[j]
                h2_list = hpo_by_pat.get(p2, [])
                len_b = len(h2_list)
                
                if not len_a or not len_b:
                    score = 0.0
                else:
                    # Similarity matrix
                    S = np.zeros((len_b, len_a), dtype=float)
                    for r, hb in enumerate(h2_list):
                        for c, ha in enumerate(h1_list):
                            try:
                                S[r, c] = hb.similarity_score(ha, 'orpha', method)
                            except Exception:
                                S[r, c] = 0.0
                    
                    max_row = np.amax(S, axis=1)
                    max_col = np.amax(S, axis=0)
                    
                    # Symmetrize rsd
                    if combine == "rsd":
                        rsd_ba = max_col.sum() / (len_a or 1)
                        rsd_ab = max_row.sum() / (len_b or 1)
                        score = (rsd_ab + rsd_ba) / 2.0
                    else:
                        score = self._combine_scores(max_row, max_col, len_a, len_b, combine)
                
                M[i, j] = score
                M[j, i] = score  # Symmetric
        
        return pd.DataFrame(M, index=pats, columns=pats)


    
    # =========================================================================
    # EXPORT  
    # =========================================================================
    
    def export_sm(self, df: pd.DataFrame, path_output: str) -> None:
        """Exports to Parquet (converts .xlsx to .parquet if necessary)."""
        if not path_output.endswith(".parquet"):
            path_output = path_output.replace(".xlsx", ".parquet")
        
        Path(path_output).parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(path_output, index=False)
    
    def compute_sm(
        self,
        index: int,
        rd: str,
        rd_id_list: list[str],
        combine: str,
        method: str,
        weights: list[float],
        out_dir: str
    ) -> None:
        """Wrapper for Snakefile: computes and also  exports DD.
        eg     out_dir="./output/dd"""
        rd_file = rd.replace(":", "-")
        df_sm = self.run_dd_freq(rd, rd_id_list, combine, method,  weights)
        df_sm.rename(columns={'patients': 'OC2', 'RDs': 'OC1'}, inplace=True)
        
        sm_path = f"{out_dir}/{index}_{rd_file}.parquet"
        self.export_sm(df_sm, sm_path)
        self.logger.info("Exported SM to %s", sm_path)
    
    def compute_sm_RDI(
        self,
        index: int,
        rd: str,
        patients: list[str],
        combine: str,
        method: str,
        weights: list[float],
        out_dir: str
    ) -> None:
        """Wrapper for Snakefile: computes and also exports DP."""
        rd_file = rd.replace(":", "-")
        df_sm = self.run_sm_freq(rd, patients, combine, method,  weights)
        
        sm_path = f"{out_dir}/{index}_{rd_file}.parquet"
        self.export_sm(df_sm, sm_path)
        self.logger.info("Exported SM to %s", sm_path)
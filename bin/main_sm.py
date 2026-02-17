#!/usr/bin/env python3
# bin/main_sm.py  
"""
Main script for semantic similarity measure (MM and MP).
"""
from pathlib import Path
import argparse
import logging

import pandas as pd

import path_variable as PV
from classes.dataset import DataSet
from classes.sim_measure_cal import Sim_measure
from classes.utils import norm_orpha, parse_weights, split_csv_list, ensure_parent
from set_log import setup_logging, get_logger


def _load_pd(mini_rd: str | None) -> pd.DataFrame:
    """Loads product4 and optionally filters by mini_rd."""
    # Load product4 with all orphacodes
    df = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT4, engine="openpyxl", index_col=0)
    df.columns = [c.strip() for c in df.columns]

    
    df["ORPHAcode"] = df["ORPHAcode"].map(norm_orpha)
    # if mini_rd is provided a filtration depending on the list(in string format) is made 
    if mini_rd:
        wanted = {norm_orpha(x) for x in mini_rd.split(",")}
        wanted.discard("")
        df = df[df["ORPHAcode"].isin(wanted)]
    
    return df


def _write_yaml_from_df(df: pd.DataFrame,) -> dict:
    """Get ORPHAcode position id from the yaml file ."""
    ds = DataSet(str(PV.PATH_YAML_PRODUCT4), "")
    cfg = ds.build_yaml_rds(df, "ORPHAcode")
    return cfg


def run_dd(
    index: int, param_rd: str, combine: str, method: str, pd4: str, vector_str: str,
    mini_rd: str | None,
    log: logging.Logger = None
) -> None:
    """Executes MM calculation (RD × RD)."""

    log = log or get_logger(__name__)
    
    ## load RDs 
    df1 = _load_pd(mini_rd)

    cfg = _write_yaml_from_df(df1)

    
    if cfg.get("n", 0) == 0:
        raise RuntimeError("YAML empty: no RD to compute. Check mini_rd.")
    
    weights = parse_weights(vector_str)
    log.info("vector_str=%r -> weights=%s", vector_str, weights)
    
    # create output directories
    out_dir = Path(PV.PATH_OUTPUT_DD) / combine / method  / pd4 / str(vector_str)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    sim = Sim_measure(df1, df1, "ORPHAcode", "ORPHAcode", logger=log)
    rd_list = df1["ORPHAcode"].drop_duplicates().tolist()
    
    param_RD = cfg.get("param_RD", {})
    for k in sorted(param_RD.keys(), key=lambda x: int(x)):
        one_param = param_RD[k]
        if one_param not in set(rd_list):
            empty = pd.DataFrame(columns=["OC1", "OC2", "score"], index=[0])
            sim.export_sm(empty, f"{out_dir}/{k}_{one_param.replace(':','-')}.parquet")
        else:
            sim.compute_sm( int(k), one_param, rd_list, combine, method, weights, str(out_dir))



def run_dp(
    index: int, param_rd: str, combine: str, method: str, pd4: str, vector_str: str,
    mini_rd: str | None,  
    mini_patient: str | None = None, log: logging.Logger = None
) -> None:
    """Executes MP calculation (Patient × RD)."""
    log = log or get_logger(__name__)
    
    # Load patients
    df_p = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT, engine="openpyxl", index_col=0)
    
    # Filter patients if requested 
    keep_patients = split_csv_list(mini_patient)
    if keep_patients:
        df_p = df_p[df_p[PV.COL_DF_PATIENT_PATIENT].isin(keep_patients)].copy()
        log.info("Filtered patients: %d IDs", len(keep_patients))
    
    log.info("[MP] Patients after filter: %d uniques", 
             df_p[PV.COL_DF_PATIENT_PATIENT].nunique())
    
    # Load RDs
    df_m = _load_pd(mini_rd)
    cfg = _write_yaml_from_df(df_m)
    
    weights = parse_weights(vector_str)
    log.info("vector_str=%r -> weights=%s", vector_str, weights)
    
    out_dir = Path(PV.PATH_OUTPUT_DP) / combine / method / pd4 / str(vector_str)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    patients_ids = df_p[PV.COL_DF_PATIENT_PATIENT].drop_duplicates().tolist()
    rd_list = df_m["ORPHAcode"].drop_duplicates().tolist()
    
    # load for each RD 
    for i in range(1, cfg["n"] + 1):
        one_param = cfg["param_RD"][i]
        sim = Sim_measure(df_p, df_m, PV.COL_DF_PATIENT_PATIENT, "ORPHAcode", logger=log)
        if one_param not in rd_list:
            empty = pd.DataFrame(columns=["RDs", "patients", "score"], index=[0])
            sim.export_sm(empty, f"{out_dir}/{i}_{one_param.replace(':','-')}.parquet")
        else:
            sim.compute_sm_RDI(i, one_param, patients_ids, combine, method, weights, str(out_dir))

 


def _add_common_args(sp: argparse.ArgumentParser) -> None:
    """Adds common arguments to sub-commands."""
    ## optimal parameters add value only forr the snakefile when used 
    sp.add_argument("--index",  default=0, help="index of RD in the yaml file")
    sp.add_argument("--param_rd",  default='ORPHA:O', help="ORPHAcode of the RD")


    sp.add_argument("--combine","-c",default='funSimMax', choices=["funSimMax", "funSimAvg", "BMA", "rsd"],help="Combination method default='funSimMax'")
    sp.add_argument("--sm", "-m", default='resnik',
                    choices=["resnik", "lin", "jc", "graphic", "ic", "relevance"],
                    help="Similarity measure default='resnik'")
    sp.add_argument("--vector-str", "-v", default='1_1_1_1_1',
                    help="Weight vector string default '1_1_1_1_1'")
    sp.add_argument("--pd4","-pd4",required=True, help="Version of Orphanet product4 named provided by user ")

    sp.add_argument("--mini-patient", default=None,
                    help="List of patient IDs 'P1,P2' ")
    
    sp.add_argument("--mini-rd", default=None,
                    help="List of RDs e.g: 'ORPHA:610,ORPHA:100985'")


def main():
    # Logging configuration
    setup_logging(level=logging.INFO,console=False,filename=f"{Path(__file__).stem}.log"    )  
    log = get_logger(Path(__file__).stem)
    
    p = argparse.ArgumentParser(prog="main_sm")
    sub = p.add_subparsers(dest="mode", required=True)
    
    sp_dd = sub.add_parser("dd", help="RD × RD similarity")
    _add_common_args(sp_dd)
    
    sp_dp = sub.add_parser("dp", help="Patient × RD similarity")
    _add_common_args(sp_dp)
    
    args = p.parse_args()
    
    common_kwargs = dict(
        index=args.index, 
        param_rd=args.param_rd, 
        combine=args.combine,
        method=args.sm, 
        pd4=args.pd4,
        vector_str=args.vector_str, 
        mini_rd=args.mini_rd,
        log=log
    )

    if args.mode == "dd":
        run_dd(**common_kwargs)
    else:
        run_dp(**common_kwargs, mini_patient=getattr(args, "mini_patient", None))
 

if __name__ == "__main__":
    main()
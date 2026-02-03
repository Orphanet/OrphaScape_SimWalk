#!/usr/bin/env python3
# bin/main_sm.py  
"""
Main script for semantic similarity measure (MM and MP).
"""
from pathlib import Path
import argparse
import logging

import pandas as pd
import yaml

import path_variable as PV
from classes.dataset import DataSet
from classes.sim_measure_cal import Sim_measure
from classes.utils import norm_orpha, parse_weights, split_csv_list, ensure_parent
from set_log import setup_logging, get_logger


def _load_pd4(mini_rd_csv: str | None) -> pd.DataFrame:
    """Loads product4 and optionally filters by mini_rd_csv."""
    df = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT4, engine="openpyxl", index_col=0)
    df.columns = [c.strip() for c in df.columns]
    
    if "ORPHACode" in df.columns and "ORPHAcode" not in df.columns:
        df = df.rename(columns={"ORPHACode": "ORPHAcode"})
    
    df["ORPHAcode"] = df["ORPHAcode"].map(norm_orpha)
    
    if mini_rd_csv:
        wanted = {norm_orpha(x) for x in mini_rd_csv.split(",")}
        wanted.discard("")
        df = df[df["ORPHAcode"].isin(wanted)]
    
    return df


def _write_yaml_from_df(df: pd.DataFrame, log: logging.Logger) -> dict:
    """Generates and writes YAML config from DataFrame."""
    ds = DataSet(str(PV.PATH_YAML_PRODUCT4), "")
    cfg = ds.build_yaml_rds(df, "ORPHAcode")
    
    ensure_parent(Path(PV.PATH_YAML_PRODUCT4))
    with open(PV.PATH_YAML_PRODUCT4, "w", encoding="utf-8") as f:
        yaml.dump(cfg, f, default_flow_style=False)
    
    log.info("YAML saved → %s", PV.PATH_YAML_PRODUCT4)
    return cfg


def run_mm(
    index: int, param_rd: str, combine: str, method: str,
    is_freq: str, pd4: str, vector_str: str,
    mini_rd_csv: str | None, only_yaml: bool, run_all: bool,
    mini_patient_csv: str | None = None, log: logging.Logger = None
) -> None:
    """Executes MM calculation (RD × RD)."""
    log = log or get_logger(__name__)
    
    df1 = _load_pd4(mini_rd_csv)
    df2 = df1.copy()
    cfg = _write_yaml_from_df(df1, log)
    
    if only_yaml:
        return
    
    if cfg.get("n", 0) == 0:
        raise RuntimeError("YAML empty: no RD to compute. Check mini_rd_csv.")
    
    weights = parse_weights(vector_str)
    log.info("vector_str=%r -> weights=%s", vector_str, weights)
    
    out_dir = Path(PV.PATH_OUTPUT_MM) / combine / method / is_freq / pd4 / str(vector_str)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    sim = Sim_measure(df2, df1, "ORPHAcode", "ORPHAcode", logger=log)
    rd_list_2 = df1["ORPHAcode"].drop_duplicates().tolist()
    rd_set_2 = set(rd_list_2)
    
    def _compute(one_i: int, one_param: str):
        if one_param not in rd_set_2:
            empty = pd.DataFrame(columns=["OC1", "OC2", "score"], index=[0])
            sim.export_sm(empty, f"{out_dir}/{one_i}_{one_param.replace(':','-')}.parquet")
        else:
            sim.compute_sm(one_i, one_param, rd_list_2, combine, method, is_freq, weights, str(out_dir))
    
    if run_all:
        param_RD = cfg.get("param_RD", {})
        for k in sorted(param_RD.keys(), key=lambda x: int(x)):
            _compute(int(k), param_RD[k])
    else:
        _compute(int(index), param_rd)


def run_mp(
    index: int, param_rd: str, combine: str, method: str,
    is_freq: str, pd4: str, vector_str: str,
    mini_rd_csv: str | None, only_yaml: bool, run_all: bool,
    mini_patient_csv: str | None = None, log: logging.Logger = None
) -> None:
    """Executes MP calculation (Patient × RD)."""
    log = log or get_logger(__name__)
    
    # Load patients
    df_p = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT, engine="openpyxl", index_col=0)
    
    # Filter patients if requested
    keep_patients = split_csv_list(mini_patient_csv)
    if keep_patients:
        df_p = df_p[df_p[PV.COL_DF_PATIENT_PATIENT].isin(keep_patients)].copy()
        log.info("Filtered patients: %d IDs", len(keep_patients))
    
    log.info("[MP] Patients after filter: %d uniques", 
             df_p[PV.COL_DF_PATIENT_PATIENT].nunique())
    
    # Load RDs
    df_m = _load_pd4(mini_rd_csv)
    cfg = _write_yaml_from_df(df_m, log)
    
    if only_yaml:
        return
    
    weights = parse_weights(vector_str)
    log.info("vector_str=%r -> weights=%s", vector_str, weights)
    
    out_dir = Path(PV.PATH_OUTPUT_SM) / combine / method / is_freq / pd4 / str(vector_str)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    patients_ids = df_p[PV.COL_DF_PATIENT_PATIENT].drop_duplicates().tolist()
    
    def _compute(one_i: int, one_param: str):
        sim = Sim_measure(df_p, df_m, PV.COL_DF_PATIENT_PATIENT, "ORPHAcode", logger=log)
        if one_param not in df_m["ORPHAcode"].tolist():
            empty = pd.DataFrame(columns=["RDs", "patients", "score"], index=[0])
            sim.export_sm(empty, f"{out_dir}/{one_i}_{one_param.replace(':','-')}.parquet")
        else:
            sim.compute_sm_RDI(one_i, one_param, patients_ids, combine, method, is_freq, weights, str(out_dir))
    
    if run_all:
        for i in range(1, cfg["n"] + 1):
            _compute(i, cfg["param_RD"][i])
    else:
        _compute(int(index), param_rd)


def _add_common_args(sp: argparse.ArgumentParser) -> None:
    """Adds common arguments to sub-commands."""
    sp.add_argument("index")
    sp.add_argument("param_rd")
    sp.add_argument("combine", choices=["funSimMax", "funSimAvg", "BMA", "rsd"])
    sp.add_argument("method", choices=["resnik", "lin", "jiang"])
    sp.add_argument("is_freq", choices=["y", "n"])
    sp.add_argument("pd4")
    sp.add_argument("vector_str")
    sp.add_argument("--mini-patient-csv", default=None,
                    help="List of patient IDs 'P1,P2' or '@path.txt'")
    
    mx = sp.add_mutually_exclusive_group()
    mx.add_argument("--only-yaml", action="store_true")
    mx.add_argument("--run-all", action="store_true")
    sp.add_argument("--mini-rd-csv", default=None,
                    help="e.g: 'ORPHA:610,ORPHA:100985'")


def main():
    setup_logging(level=logging.INFO, console=False, filename=f"{Path(__file__).stem}.log")
    log = get_logger(Path(__file__).stem)
    
    p = argparse.ArgumentParser(prog="main_sm")
    sub = p.add_subparsers(dest="mode", required=True)
    
    sp_mm = sub.add_parser("mm", help="RD × RD similarity")
    _add_common_args(sp_mm)
    
    sp_mp = sub.add_parser("mp", help="Patient × RD similarity")
    _add_common_args(sp_mp)
    
    args = p.parse_args()
    
    kwargs = dict(
        index=args.index, param_rd=args.param_rd, combine=args.combine,
        method=args.method, is_freq=args.is_freq, pd4=args.pd4,
        vector_str=args.vector_str, mini_rd_csv=args.mini_rd_csv,
        only_yaml=args.only_yaml, run_all=args.run_all,
        mini_patient_csv=getattr(args, "mini_patient_csv", None),
        log=log
    )
    
    if args.mode == "mm":
        run_mm(**kwargs)
    else:
        run_mp(**kwargs)


if __name__ == "__main__":
    main()
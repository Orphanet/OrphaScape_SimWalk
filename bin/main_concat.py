#!/usr/bin/env python3
# bin/main_concat.py  
"""
Script for concatenating SM/MM files.
"""

from pathlib import Path
import argparse
import logging
import path_variable as PV
from classes.concat_sm import ConcatSm
from set_log import setup_logging, get_logger

def main():
    p = argparse.ArgumentParser(description="Concat SM/MM exports into one file per combo")
    sub = p.add_subparsers(dest="cmd", required=True)
    
    # MM: ORPHA×ORPHA matrices
    mm = sub.add_parser("concat_matrix", help="Pivot & concatenate MM score matrices")
    mm.add_argument("-v", "--vector_str", required=True)
    mm.add_argument("--col1", required=True, help="row index col (e.g. OC1)")
    mm.add_argument("--col2", required=True, help="col index col (e.g. OC2)")
    mm.add_argument("--product4", required=True)
    mm.add_argument("--combine", default=None)
    mm.add_argument("--sm", default=None)
    mm.add_argument("--outxlsx", default=None)
    
    # MP: patients ↔ RDs
    sm = sub.add_parser("process_similarity", help="Concat & rank SM/CDF outputs")
    sm.add_argument("-v", "--vector_str", required=True)
    sm.add_argument("--col1", required=True, help="sample col (e.g. patients)")
    sm.add_argument("--col2", required=True, help="label col (e.g. RDs)")
    sm.add_argument("--product4", required=True)
    sm.add_argument("--combine", default=None)
    sm.add_argument("--sm", default=None)
    sm.add_argument("--outxlsx", default=None)
    
    args = p.parse_args()
    
    setup_logging(level=logging.INFO, console=False, filename=f"{Path(__file__).stem}.log")
    log = get_logger(Path(__file__).stem)
    
    base_dir = PV.PATH_OUTPUT_MM if args.cmd == "concat_matrix" else PV.PATH_OUTPUT_SM
    
    cx = ConcatSm(
        vector_str=args.vector_str,
        col1=args.col1,
        col2=args.col2,
        base_dir=base_dir,
        product4=args.product4,
        combine=args.combine,
        sm=args.sm,
        logger=log,
        out_path=args.outxlsx,
    )
    
    if args.cmd == "concat_matrix":
        cx.concat_matrix()
    else:
        cx.process_similarity()

if __name__ == "__main__":
    main()
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
    p = argparse.ArgumentParser(description="Concat DP/DD exports into one file per combo")
    sub = p.add_subparsers(dest="cmd", required=True)
    
    # MM: ORPHA×ORPHA matrices
    mm = sub.add_parser("concat_dd", help="Pivot & concatenate DD  matrices")
    mm.add_argument("--combine", "-c",default='funSimMax', choices=["funSimMax", "funSimAvg", "BMA", "rsd"],help="Combination method default='funSimMax'")
    mm.add_argument("--sm", "-m",default='resnik',choices=["resnik", "lin", "jc","graphic","ic","relevance"], help="Similarity measure default='resnik'")
    mm.add_argument("--vector-str", "-v",default='1_1_1_1_1', help="Weight vector string, default '1_1_1_1_1'")
    mm.add_argument("--pd4",required=True, help="Version of Orphanet product4 named provided by user ")

    
    mm.add_argument("--col1",default='OC1',  help="name of the label col1 by default RDs")
    mm.add_argument("--col2",default='OC2', help="name of the label col 2 by default RDs")



    # MP: patients ↔ RDs
    sm = sub.add_parser("concat_dp", help="Concat & rank SM/CDF outputs")
    sm.add_argument("--combine", "-c",default='funSimMax', choices=["funSimMax", "funSimAvg", "BMA", "rsd"],help="Combination method default='funSimMax'")
    sm.add_argument("--sm", "-m",default='resnik',choices=["resnik", "lin", "jc","graphic","ic","relevance"], help="Similarity measure default='resnik'")
    sm.add_argument("--vector-str", "-v",default='1_1_1_1_1', help="Weight vector string, default '1_1_1_1_1'")
    sm.add_argument("--pd4","-pd4",required=True, help="Version of Orphanet product4 named provided by user ")

 
    ## parameter col option 
    sm.add_argument("--col1",default='patients',  help="name of the sample col by default patients")
    sm.add_argument("--col2",default='RDs', help="name of the label col by default RDs")


    
    args = p.parse_args()
    
    # Logging configuration
    setup_logging(level=logging.INFO,console=False,filename=f"{Path(__file__).stem}.log"    )  
    log = get_logger(Path(__file__).stem)
    
    base_dir = PV.PATH_OUTPUT_DD if args.cmd == "concat_dd" else PV.PATH_OUTPUT_DP
    
    cx = ConcatSm(
        vector_str=args.vector_str,
        col1=args.col1,
        col2=args.col2,
        base_dir=base_dir,
        pd4=args.pd4,
        combine=args.combine,
        method=args.sm,
        logger=log,
        # out_path=args.outxlsx,
    )

    
    if args.cmd == "concat_dd":
        cx.concat_dd()
    else:
        cx.concat_dp(PV.PATH_OUTPUT_DF_PATIENT)

if __name__ == "__main__":
    main()



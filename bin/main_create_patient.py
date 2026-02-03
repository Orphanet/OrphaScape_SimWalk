#!/usr/bin/env python3
# bin/main_create_patient.py 
"""
Script to create patients DataFrame from Phenopacket files.
"""

from pathlib import Path
import os
import json
import logging
import warnings
import pandas as pd
import path_variable as PV
from set_log import setup_logging, get_logger

def main():
    # Logging configuration
    warnings.filterwarnings(
        "ignore",
        message="Data Validation extension is not supported and will be removed",
        category=UserWarning,
        module="openpyxl",
    )
    logging.captureWarnings(True)
    setup_logging(level=logging.INFO, console=False, filename=f"{Path(__file__).stem}.log")
    log = get_logger(Path(__file__).stem)
    
    # Create output directory
    Path(PV.PATH_OUTPUT_PATIENT_SOLVERD).mkdir(parents=True, exist_ok=True)
    
    # Load patients
    patients_raw = os.listdir(PV.PATH_INPUT_PATIENTS_FOLDER)
    log.info("%d patients in the folder study_population", len(patients_raw))
    
    list_case_HPO = []
    for onefile in patients_raw:
        filepath = Path(PV.PATH_INPUT_PATIENTS_FOLDER) / onefile
        try:
            with open(filepath, 'r', encoding='utf8') as f:
                data = json.load(f)
                patient_id = data['id']
                phenotypes = data.get('phenotypes', [])
                for pheno in phenotypes:
                    pheno_type = pheno.get('type', {})
                    list_case_HPO.append((
                        patient_id,
                        pheno_type.get('id', ''),
                        pheno_type.get('label', '')
                    ))
        except (json.JSONDecodeError, KeyError, TypeError) as e:
            log.warning("Error reading %s: %s", onefile, e)
            continue
    
    df_raw_info = pd.DataFrame(
        list_case_HPO,
        columns=["phenopacket", "hpo_id", "hpo_label"]
    )
    log.info("%d Patients solved confirmed with ontologyX", 
        df_raw_info["phenopacket"].nunique())
    
    df_raw_info.to_excel(PV.PATH_OUTPUT_DF_PATIENT)
    log.info("Saved patients to %s", PV.PATH_OUTPUT_DF_PATIENT)

if __name__ == "__main__":
    main()
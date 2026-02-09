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
from classes.dataset import DataSet   

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
    build_patient = DataSet(PV.PATH_INPUT_PATIENTS_FOLDER  ,"") 
    ## patient solved and unsolved
    df_patient = build_patient.build_patients_df()
    df_raw_info  = df_patient[['phenopacket','hpo_id','hpo_label']]

    ########################################################
    ## Replace subtype by disorder in the patient df by disorder parent
    ########################################################
    
    orpha_map = build_patient.parse_orpha_file(PV.PATH_ORPHA_RDI_FILE)
    df_raw_info["Disease"] = df_raw_info["phenopacket"].map(orpha_map)

    df_parent_child_noclassif = pd.read_excel(PV.PATH_OUTPUT_DF_PC_v2,index_col=0)


    df1 = pd.read_excel(PV.PATH_OUTPUT_DF_PRODUCT1,index_col=0)
    df1 = df1[['ORPHAcode',"Type","Group"]]
    df1.columns = ['Disease','Type','Group']

    valid_types = df1[df1["Group"] == "Disorder"]['Type'].drop_duplicates().tolist()
    all_interactions= []
    dict_df_patient_f =  df_raw_info.to_dict('index')
    for value in dict_df_patient_f.values():
        patient = value['phenopacket']
        hpoid = value['hpo_id']
        disease = value['Disease']

        # test the disease type 
        type_disease = ""
        df_get_type = df_parent_child_noclassif[(df_parent_child_noclassif['parent_id'] ==disease ) | (df_parent_child_noclassif['child_id'] ==disease) ]
        # print(disease)
        # print(df_get_type)
        dict_df_get_type =  df_get_type.to_dict('index')
        for value in dict_df_get_type.values():
            parent_id = value['parent_id']
            parent_type = value['parent_type'].strip()
            child_id = value['child_id']
            child_type = value['child_type'].strip()
            
            foundin = ""
            # si la maladie est dans la colonne child 
            if disease in child_id:
                if child_type in valid_types:
                    type_disease = child_type
                    foundin = "child_id"
                # si non trouver alors il est dans la col parent  
                elif disease in parent_id: 
                    type_disease = parent_type
                    foundin = "child_id but parent "
                # sinon c'est un subtype meaning the parent are disorder
                elif parent_type in valid_types:
                    foundin = "child_id_subtype"
                    # get the related parent
                    new_disease = parent_id 
                    disease = parent_id 

                    print(f"before {disease}, after {new_disease}")
                            
            elif disease in parent_id:
                if parent_type in valid_types:
                    type_disease = parent_type
                    foundin = "parent_id"
                    print(f"before {disease}, after {new_disease}")
                elif child_type in valid_types:
                    foundin = "parent_id_is_category"
                    new_disease = child_id 
                    disease = child_id 
                    print(f"before {disease}, after {new_disease}")
            # else :
            #     print(disease)
            

        print(f"{disease} type is {type_disease} found in {foundin}")
        all_interactions.append((patient,hpoid,disease))

    df_patient_only_disorder = pd.DataFrame(all_interactions,columns=["phenopacket",'hpo_id','Disease'])
    df_final = df_patient_only_disorder.merge(df1,on = ['Disease'], how='inner').dropna(subset=['phenopacket'])

    df_final.to_excel(PV.PATH_OUTPUT_DF_PATIENT)
    print("Saved patients to %s", PV.PATH_OUTPUT_DF_PATIENT)



if __name__ == "__main__":
    main()
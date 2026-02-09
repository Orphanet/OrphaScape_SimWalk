from pathlib import Path
import sys

# Find /projet by going up until we see "bin" and "classes"
p = Path.cwd().resolve()
while p != p.parent and not ((p / "bin").exists() and (p / "classes").exists()):
    p = p.parent

PROJECT_ROOT = p
sys.path.insert(0, str(PROJECT_ROOT))

print("Detected project root:", PROJECT_ROOT)

"""
#### for main_sm.py
common_kwargs = dict(

index=0

param_rd="ORPHA:00"

combine="funSimMax"

method="resnik"

pd4="args.pd4"

vector_str="1_1_1_1_1"

mini_rd="ORPHA:610,ORPHA:100985"

mini_patient="P0001068"
log=log 


######## for main_concat.py


cmd = "concat_matrix_mp" #process_similarity
base_dir = PV.PATH_OUTPUT_MM if cmd == "concat_matrix_mm" else PV.PATH_OUTPUT_SM

vector_str="1_1_1_1_1"
col1="OC1" # patients
col2="OC2" # RDs
base_dir=base_dir
product4="args.pd4"
combine="funSimMax"
sm="resnik"
out_path="./output/final_MM_matrix.xlsx"

col1="patients"  
col2="RDs"  
path_patient = 
    df_p = pd.read_excel(PV.PATH_OUTPUT_DF_PATIENT, engine="openpyxl", index_col=0)

logger=log,


## pour rarw main
seeds = "P0001068"

alpha = 'O.3'

matrix_subdir =max(PV.PATH_OUTPUT_PATIENT_ADDED.glob("*/"), key=os.path.getmtime)
matrix_subdir_file = str(matrix_subdir).split("/")[-1]

"""
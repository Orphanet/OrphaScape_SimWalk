# path_variable.py 
"""
Variables de chemins centralisées pour le pipeline.
IMPORTANT: Ne contient PAS de code exécutable au chargement (pas de pd.read_excel, etc.)
"""
from pathlib import Path

# =============================================================================
# RACINE DU PROJET
# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent

# =============================================================================
# DOSSIERS PRINCIPAUX
# =============================================================================
PATH_OUTPUT = PROJECT_ROOT / "output"
PATH_INPUT = PROJECT_ROOT / "input"
PATH_INPUT_HPO = PATH_INPUT / "hpo"

# =============================================================================
# COLONNES (constantes string)
# =============================================================================
COL_DF_PATIENT_PATIENT = "phenopacket"

# =============================================================================
# DOSSIERS DE SORTIE SM/MM
# =============================================================================
PATH_OUTPUT_SM = PATH_OUTPUT / "mp_sm"
PATH_OUTPUT_MM = PATH_OUTPUT / "mm_sm"

# =============================================================================
# DOSSIERS PATIENTS
# =============================================================================
PATH_OUTPUT_PATIENT_SOLVERD = PATH_OUTPUT / "patient_solverd"
PATH_OUTPUT_PATIENT_SOLVERD_UNSOLVED = PATH_OUTPUT_PATIENT_SOLVERD / "patient_solverd_unsolved.xlsx"
PATH_OUTPUT_DF_PATIENT = PATH_OUTPUT_PATIENT_SOLVERD / "patients.xlsx"

# =============================================================================
# RANDOM WALK
# =============================================================================
PATH_OUTPUT_FOLDER_RW = PATH_OUTPUT / "rarw"

# =============================================================================
# PATIENT ADDED (matrices MM avec patients)
# =============================================================================
PATH_OUTPUT_PATIENT_ADDED = PATH_OUTPUT / "patient_added"

# =============================================================================
# ORPHANET PRODUCTS
# =============================================================================
PATH_OUTPUT_PRODUCT = PATH_OUTPUT / "pd_orphanet"
PATH_OUTPUT_PRODUCT_CLASSIF = PATH_OUTPUT_PRODUCT / "Classifications"
PATH_OUTPUT_DF_PC_CLASSIF_v2 = PATH_OUTPUT_PRODUCT / "parent_child_classif_v2.xlsx"
PATH_OUTPUT_DF_PC_v2 = PATH_OUTPUT_PRODUCT / "parent_child_noclassif_v2.xlsx"


PATH_OUTPUT_DF_CLASSIF_NOG = PATH_OUTPUT_PRODUCT / "classifs_nogenetic.xlsx"
PATH_OUTPUT_DF_CLASSIF   = PATH_OUTPUT_PRODUCT / "classifs.xlsx"

# JSON et Excel de sortie
PATH_OUTPUT_PRODUCT4_JSON = PATH_OUTPUT_PRODUCT / "all_enpd_mai_2025.json"
PATH_OUTPUT_DF_PRODUCT4 = PATH_OUTPUT_PRODUCT / "all_enpd_mai_2025.xlsx"
PATH_YAML_PRODUCT4 = PATH_OUTPUT_PRODUCT / "config_from_product4.yaml"
PATH_CLASSIFICATION_JSON = PATH_OUTPUT_PRODUCT / "classif_orpha.json"

PATH_OUTPUT_PRODUCT7_JSON = PATH_OUTPUT_PRODUCT / "en_product7.json"
PATH_OUTPUT_DF_PRODUCT7 = PATH_OUTPUT_PRODUCT / "en_product7.xlsx"

PATH_OUTPUT_PRODUCT1_JSON = PATH_OUTPUT_PRODUCT / "en_product1.json"
PATH_OUTPUT_DF_PRODUCT1 = PATH_OUTPUT_PRODUCT / "en_product1.xlsx"

PATH_OUTPUT_PRODUCT6_JSON = PATH_OUTPUT_PRODUCT / "en_product6.json"
PATH_OUTPUT_DF_PRODUCT6 = PATH_OUTPUT_PRODUCT / "en_product6.xlsx"

# =============================================================================
# INPUT PRODUCTS
# =============================================================================
PATH_INPUT_PRODUCT = PATH_INPUT / "pd_orphanet"
PATH_INPUT_PRODUCT4_XML = PATH_INPUT_PRODUCT / "en_product4.xml"
PATH_INPUT_PRODUCT1_XML = PATH_INPUT_PRODUCT / "en_product1.xml"
PATH_INPUT_PRODUCT7_XML = PATH_INPUT_PRODUCT / "en_product7.xml"
PATH_INPUT_PRODUCT6_XML = PATH_INPUT_PRODUCT / "en_product6.xml"
PATH_INPUT_PRODUCTCLASSIF_XML = PATH_INPUT_PRODUCT / "Classifications"

# =============================================================================
# INPUT PATIENTS
# =============================================================================
PATH_INPUT_PATIENTS_FOLDER = PATH_INPUT / "patient" 
PATH_ORPHA_RDI_FILE = PATH_INPUT /"patient_RDI.txt"
# PATH_INPUT_PATIENTS_FC = PATH_INPUT / "patient" / "PATIENTS_SOLVED_FC_v2.xlsx"


# =============================================================================
# CONVERSION AUTOMATIQUE EN STRING POUR COMPATIBILITÉ
# =============================================================================
def _ensure_str(p: Path) -> str:
    """Convertit un Path en string."""
    return str(p)

# Pour la rétro-compatibilité avec os.path.join et scripts existants,
# on peut accéder aux chemins comme strings via:
#   str(PATH_OUTPUT_MM) ou PATH_OUTPUT_MM.as_posix()

#!/usr/bin/env python3
# bin/main_load_data.py 
"""
Script to load and convert Orphanet data.
Generates DataFrames and YAML needed for the pipeline.
"""

from pathlib import Path
import logging

import pandas as pd
import yaml

from classes.dataset import DataSet
from classes.utils import ensure_parent
import path_variable as PV
from set_log import setup_logging, get_logger


def convert_and_build_dataframe(
    path_xml: Path,
    path_json: Path,
    nb_type: int,
    log: logging.Logger
) -> pd.DataFrame:
    """
    Converts XML -> JSON -> DataFrame.
    Uses JSON cache if available and up to date.
    """
    path_xml = Path(path_xml)
    path_json = Path(path_json)
    
    # JSON cache: rebuild if missing or XML is more recent
    need_rebuild = (
        not path_json.exists() or 
        (path_xml.exists() and path_xml.stat().st_mtime > path_json.stat().st_mtime)
    )
    
    if need_rebuild:
        log.info("[BUILD] %s => %s", path_xml.name, path_json.name)
        build_json = DataSet(str(path_xml), "")
        json_obj = build_json.from_xml_to_json()
        ensure_parent(path_json)
        build_json.save_json(str(path_json), json_obj)
    else:
        log.info("[CACHE] %s", path_json.name)
    
    # JSON -> DataFrame
    build_df = DataSet(str(path_json), "")
    
    type_handlers = {
        4: build_df.build_orpha_df,
        44: build_df.from_rsd_build_orpha_df,
        6: build_df.df_pd6,
        1: build_df.df_pd1,
        7: build_df.df_pd7,
    }
    
    handler = type_handlers.get(nb_type)
    if not handler:
        raise ValueError(f"Unknown nb_type: {nb_type}")
    
    return handler()


def build_classif(
    path_input_classif: Path,
    path_output_classif: Path,
    log: logging.Logger
) -> None:
    """Converts all classification files XML -> JSON -> XLSX."""
    pin = Path(path_input_classif)
    pout = Path(path_output_classif)
    pout.mkdir(parents=True, exist_ok=True)
    
    xml_files = sorted(p for p in pin.iterdir() if p.suffix.lower() == ".xml")
    
    for xml_path in xml_files:
        motif = xml_path.stem
        
        # XML -> JSON
        ds_json = DataSet(str(xml_path), "")
        json_obj = ds_json.from_xml_to_json()
        json_out = pout / f"{motif}.json"
        ds_json.save_json(str(json_out), json_obj)
        
        # JSON -> XLSX
        ds_cls = DataSet(str(json_out), "")
        df_pd_classif = ds_cls.df_classif()
        df_pd_classif.to_excel(pout / f"{motif}.xlsx")
        log.info("[CLASSIF] %s.xlsx", motif)


def process_classifications(log: logging.Logger) -> None:
    """Processes and merges classifications."""
    classif_dir = Path(PV.PATH_OUTPUT_PRODUCT_CLASSIF)
    out_dir = Path(PV.PATH_OUTPUT_PRODUCT)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Load and merge all classifications
    ds_cls = DataSet(str(classif_dir), "")
    df_all = ds_cls.merge_all_classif(str(classif_dir))
    df_all.to_excel(out_dir / "classifs.xlsx")
    
    # Classifications without genetics
    if "root" in df_all.columns:
        df_no_gen = df_all[df_all["root"] != "ORPHA:98053"]
        log.info("Nb classif without genetic: %d", df_no_gen["root"].nunique())
        df_no_gen.to_excel(out_dir / "classifs_nogenetic.xlsx")
        
        # Only neuro
        df_neuro = df_all[df_all["root"] == "ORPHA:98006"]
        df_neuro.to_excel(out_dir / "classifs_onlyneuro.xlsx")
    
    # Extract parent-child
    xlsx_files = [p for p in classif_dir.glob("*.xlsx") if not p.name.startswith("~$")]
    log.info("Nb classif: %d", len(xlsx_files))
    
    dfs = [pd.read_excel(p, index_col=0, engine="openpyxl") 
           for p in xlsx_files]
    dfs = [df for df in dfs if not df.empty]
    log.info("Nb non-empty classif: %d", len(dfs))
    
    df_global = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    
    # Valid types
    valid_parent_types = [
        "Disease", "Morphological anomaly", "Malformation syndrome",
        "Biological anomaly", "Clinical syndrome",
        "Particular clinical situation in a disease or syndrome",
        "Clinical subtype", "Etiological subtype", "Histopathological subtype"
    ]
    
    valid_child_subtypes = ["Etiological subtype", "Histopathological subtype", "Clinical subtype"]
    no_valid_child_type = [
        "Disease", "Morphological anomaly", "Malformation syndrome",
        "Biological anomaly", "Clinical syndrome",
        "Particular clinical situation in a disease or syndrome"
    ]
    
    # Parents with subtypes
    df_with_subtypes = df_global[
        (df_global["parent_type"].isin(valid_parent_types)) &
        (df_global["child_type"].isin(valid_child_subtypes))
    ]
    parent_ids_with_subtype = df_with_subtypes["parent_id"].unique()
    log.info("%d have a subtype", len(parent_ids_with_subtype))
    
    # Parents without subtype
    df_no_subtype = df_global[
        (~df_global["child_id"].isin(no_valid_child_type)) &
        (~df_global["parent_id"].isin(parent_ids_with_subtype))
    ]
    log.info("%d have no subtype", df_no_subtype["parent_id"].nunique())
    
    # Diseases linked to other diseases
    df_subtype_disease = df_global[
        (df_global["parent_type"].isin(valid_parent_types)) &
        (~df_global["parent_id"].isin(parent_ids_with_subtype))
    ]
    log.info("%d diseases related to another disease", 
             df_subtype_disease["parent_id"].nunique())
    
    # Concat and export
    df_parent_child = pd.concat(
        [df_with_subtypes, df_no_subtype, df_subtype_disease],
        ignore_index=True
    ).sort_values(by="parent_id")
    
    df_parent_child.to_excel(PV.PATH_OUTPUT_DF_PC_CLASSIF_v2)
    
    df_parent_child_noclassif = df_parent_child[
        ["parent_id", "parent_type", "child_id", "child_type"]
    ].drop_duplicates()
    df_parent_child_noclassif.to_excel(PV.PATH_OUTPUT_DF_PC_v2)


def main():
    setup_logging(level=logging.INFO, console=False, filename=f"{Path(__file__).stem}.log")
    log = get_logger(Path(__file__).stem)
    
    # PD4
    df4 = convert_and_build_dataframe(
        PV.PATH_INPUT_PRODUCT4_XML, PV.PATH_OUTPUT_PRODUCT4_JSON, 4, log
    )
    ensure_parent(Path(PV.PATH_OUTPUT_DF_PRODUCT4))
    df4.to_excel(PV.PATH_OUTPUT_DF_PRODUCT4)
    
    # PD6
    df6 = convert_and_build_dataframe(
        PV.PATH_INPUT_PRODUCT6_XML, PV.PATH_OUTPUT_PRODUCT6_JSON, 6, log
    )
    ensure_parent(Path(PV.PATH_OUTPUT_DF_PRODUCT6))
    df6.to_excel(PV.PATH_OUTPUT_DF_PRODUCT6)
    
    # PD1
    df1 = convert_and_build_dataframe(
        PV.PATH_INPUT_PRODUCT1_XML, PV.PATH_OUTPUT_PRODUCT1_JSON, 1, log
    )
    ensure_parent(Path(PV.PATH_OUTPUT_DF_PRODUCT1))
    df1.to_excel(PV.PATH_OUTPUT_DF_PRODUCT1)
    
    # PD7
    df7 = convert_and_build_dataframe(
        PV.PATH_INPUT_PRODUCT7_XML, PV.PATH_OUTPUT_PRODUCT7_JSON, 7, log
    )
    ensure_parent(Path(PV.PATH_OUTPUT_DF_PRODUCT7))
    df7.to_excel(PV.PATH_OUTPUT_DF_PRODUCT7)
    
    # YAML for Snakefile
    log.info("Create YAML for snakefile (mp/mm)")
    ds_yaml = DataSet(str(PV.PATH_YAML_PRODUCT4), "")
    config_yaml = ds_yaml.build_yaml_rds(df4, df4.columns[0])
    ensure_parent(Path(PV.PATH_YAML_PRODUCT4))
    with open(PV.PATH_YAML_PRODUCT4, "w", encoding="utf-8") as f:
        yaml.dump(config_yaml, f, default_flow_style=False)
    log.info("YAML saved → %s", PV.PATH_YAML_PRODUCT4)
    
    # Classifications
    build_classif(PV.PATH_INPUT_PRODUCTCLASSIF_XML, PV.PATH_OUTPUT_PRODUCT_CLASSIF, log)
    log.info("All classification files exported (xml→json→xlsx).")
    
    process_classifications(log)


if __name__ == "__main__":
    main()
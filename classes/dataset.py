# classes/dataset.py  
"""
DataSet class for manipulation of Orphanet and patient data.
- Use of utils module for common functions
- Removal of duplicate code (freq_to_score)
- Better error handling
- Type hints
"""
from pathlib import Path
import os
import json
import glob
from typing import Any

import pandas as pd
import numpy as np

from classes.datagenerate import DataGenerate
from classes.utils import freq_to_score
from set_log import get_logger

log = get_logger(__name__)


class DataSet(DataGenerate):
    """Class for construction and manipulation of Orphanet/patient datasets."""
    
    def __init__(self, input_path: str | Path, output_path: str | Path = ""):
        super().__init__(input_path, output_path)
    
    # =========================================================================
    # PATIENTS
    # =========================================================================
    
    def build_patients_df(self) -> pd.DataFrame:
        """
        Builds a DataFrame of patients from phenopacket files.
        Handles variable JSON formats and extracts HPO, diseases, genes.
        """
        pheno_with_invalid_id = set()
        list_case_HPO = []
        
        patients_raw = os.listdir(self.input_path)
        
        for onefile in patients_raw:
            filepath = self.input_path / onefile
            try:
                with open(filepath, 'r', encoding='utf8') as f:
                    data = json.load(f)
                
                patient_id = data['id']
                
                # Progress status
                progress_status = (
                    data.get('interpretations', [{}])[0].get('progressStatus', 'no_info')
                )
                
                # ERN
                ern = data.get('metaData', {}).get('externalReferences', [{}])[0].get('id', 'no_info')
                
                # Genetic information
                type_gene, gene, variant = 'no_info', 'no_info', 'no_info'
                interp = data.get('interpretations', [{}])[0]
                if 'diagnosis' in interp:
                    gene_section = interp['diagnosis'].get('genomicInterpretations', [])
                    if gene_section:
                        gs = gene_section[0]
                        type_gene = gs.get('interpretationStatus', 'no_info')
                        if 'gene' in gs:
                            gene = gs['gene'].get('symbol', 'no_info')
                        elif 'variantInterpretation' in gs:
                            vi = gs['variantInterpretation']
                            variant = vi.get('acmgPathogenicityClassification', 'no_info')
                            gene = vi.get('variationDescriptor', {}).get('geneContext', {}).get('symbol', 'no_info')
                
                # Diseases
                disease_omim, disease_orphanet = 'no_info', 'no_info'
                for disease in data.get('diseases', []):
                    disease_id = disease.get('term', {}).get('id', '')
                    if 'OMIM' in disease_id:
                        disease_omim = disease_id
                    elif 'Orphanet' in disease_id:
                        parts = disease_id.split(':')
                        disease_orphanet = f"ORPHA:{parts[1]}" if len(parts) > 1 else disease_id
                
                # HPO phenotypes
                phenotypes = data.get('phenotypicFeatures', [])
                if not phenotypes:
                    continue
                
                for pheno in phenotypes:
                    if not pheno:
                        continue
                    
                    # Handle two possible formats
                    try:
                        if 'type.label' in pheno:
                            type_data = pheno['type.label']
                            if type_data == 'Invalid id':
                                pheno_with_invalid_id.add(patient_id)
                                continue
                            if 'type.negated' in pheno:
                                continue
                            id_hpo = type_data.get('id', '')
                            label_hpo = type_data.get('label', '')
                        else:
                            type_data = pheno.get('type', {})
                            if type_data.get('label') == 'Invalid id':
                                pheno_with_invalid_id.add(patient_id)
                                continue
                            if pheno.get('negated'):
                                continue
                            id_hpo = type_data.get('id', '')
                            label_hpo = type_data.get('label', '')
                        
                        if id_hpo and patient_id not in pheno_with_invalid_id:
                            list_case_HPO.append((
                                patient_id, progress_status, disease_orphanet, disease_omim,
                                ern, gene, type_gene, variant, id_hpo, label_hpo
                            ))
                    except (KeyError, TypeError):
                        continue
                        
            except (json.JSONDecodeError, UnicodeDecodeError) as e:
                log.warning("Cannot read %s: %s", onefile, e)
                continue
        
        return pd.DataFrame(
            list_case_HPO,
            columns=['phenopacket', 'status', 'Orphanet', 'OMIM', 'ern', 
                    'gene', 'type_gene', 'variant', 'hpo_id', 'hpo_label']
        )
    
    def filter_df_keep_confirmed_only(
        self, 
        path_input_p: str | Path, 
        df_input_p: pd.DataFrame, 
        col_patient: str
    ) -> pd.DataFrame:
        """Filters to keep only confirmed patients."""
        df_confirmed = pd.read_excel(path_input_p, engine='openpyxl', sheet_name='Feuil2')
        df_confirmed = df_confirmed[df_confirmed['Result'] == 'yes']
        
        patients_c = df_confirmed['Patient ID'].drop_duplicates().tolist()
        df_confirmed = df_confirmed[['Patient ID', 'Gene', 'Disease found ORPHA']]
        df_confirmed.columns = [col_patient, 'Gene_p', 'Disease']
        
        df_filtered = df_input_p[df_input_p[col_patient].isin(patients_c)]
        df_merged = pd.merge(df_filtered, df_confirmed, on=col_patient, how='outer')
        
        return df_merged.dropna(subset=['hpo_id'])
    
    # =========================================================================
    # ORPHANET PRODUCTS
    # =========================================================================
    
    def build_orpha_df(self) -> pd.DataFrame:
        """Builds RD DataFrame with HPO from product4 JSON."""
        with open(self.input_path, 'r') as f:
            root = json.load(f)
        
        all_interactions = []
        diseases = root['JDBOR']['HPODisorderSetStatusList']['HPODisorderSetStatus']
        
        for disorder in diseases:
            the_dict = disorder['Disorder']
            disease_id = the_dict['OrphaCode']
            
            hpo_section = the_dict.get('HPODisorderAssociationList', {})
            if 'HPODisorderAssociation' not in hpo_section:
                continue
            
            hpo_list = hpo_section['HPODisorderAssociation']
            if not isinstance(hpo_list, list):
                hpo_list = [hpo_list]
            
            for hpo_entry in hpo_list:
                try:
                    hpo_id = hpo_entry['HPO']['HPOId']
                    freq_text = hpo_entry['HPOFrequency']['Name']['#text']
                    freq_score = freq_to_score(freq_text)
                    
                    if freq_score is None:
                        log.warning("Unknown frequency: %s for %s", freq_text, disease_id)
                    
                    all_interactions.append((
                        f"ORPHA:{disease_id}",
                        hpo_id,
                        freq_score
                    ))
                except (KeyError, TypeError, AttributeError) as e:
                    log.debug("Error parsing HPO for %s: %s", disease_id, e)
                    continue
        
        return pd.DataFrame(
            all_interactions,
            columns=['ORPHAcode', 'hpo_id', 'hpo_frequency']
        )
    
    def from_rsd_build_orpha_df(self) -> pd.DataFrame:
        """Builds RD DataFrame from RSD format."""
        with open(self.input_path, 'r') as f:
            data = json.load(f)
        
        rows = []
        missing = []
        
        for disorder in data['JDBOR']['DisorderList']['Disorder']:
            orpha = disorder.get('OrphaCode')
            assoc = disorder.get('HPODisorderAssociationList', {}).get('HPODisorderAssociation', [])
            
            if isinstance(assoc, dict):
                assoc = [assoc]
            
            if not assoc:
                missing.append(f"ORPHA:{orpha}")
                continue
            
            for a in assoc:
                hpo = a.get('HPO', {})
                freq_text = a.get('HPOFrequency', {}).get('Name', {}).get('#text')
                freq = freq_to_score(freq_text)
                
                rows.append({
                    'ORPHAcode': f"ORPHA:{orpha}",
                    'hpo_id': hpo.get('HPOId'),
                    'hpo_frequency': freq
                })
        
        df = pd.DataFrame(rows)
        log.info("%d RDs with HPO associations.", df['ORPHAcode'].nunique())
        
        if missing:
            log.info("Warning: %d disorders have no HPO associations.", len(missing))
        
        return df
    
    def df_pd1(self) -> pd.DataFrame:
        """Builds product1 DataFrame (nomenclature)."""
        with open(self.input_path, 'r') as f:
            root = json.load(f)
        
        interactions = set()
        for disorder in root['JDBOR']['DisorderList']['Disorder']:
            interactions.add((
                f"ORPHA:{disorder['OrphaCode']}",
                disorder['Name']['#text'],
                disorder['DisorderType']['Name']['#text'],
                disorder['DisorderGroup']['Name']['#text']
            ))
        
        return pd.DataFrame(
            interactions,
            columns=['ORPHAcode', 'Name', 'Type', 'Group']
        )
    
    def df_pd6(self) -> pd.DataFrame:
        """Builds product6 DataFrame (genes)."""
        with open(self.input_path, 'r') as f:
            root = json.load(f)
        
        interactions = set()
        for disorder in root['JDBOR']['DisorderList']['Disorder']:
            disease_id = f"ORPHA:{disorder['OrphaCode']}"
            disease_name = disorder['Name']['#text']
            
            genelist = disorder.get('DisorderGeneAssociationList', {}).get('DisorderGeneAssociation', [])
            if not isinstance(genelist, list):
                genelist = [genelist]
            
            for g in genelist:
                gene = g.get('Gene', {})
                interactions.add((
                    disease_id,
                    disease_name,
                    gene.get('@id', ''),
                    gene.get('Symbol', '')
                ))
        
        return pd.DataFrame(
            interactions,
            columns=['ORPHAcode', 'Name', 'Gene_orphaid', 'Symbol']
        )
    
    def df_pd7(self) -> pd.DataFrame:
        """Builds product7 DataFrame (hierarchy)."""
        with open(self.input_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        disorders = data['JDBOR']['DisorderList']['Disorder']
        interactions = set()
        
        for disorder in disorders:
            disease_id = disorder.get('OrphaCode')
            assoc_list = disorder.get('DisorderDisorderAssociationList', {})
            count = assoc_list.get('@count', '0')
            
            if count == '0':
                log.debug("%s has no association classif", disease_id)
                continue
            
            associations = assoc_list.get('DisorderDisorderAssociation', [])
            if not isinstance(associations, list):
                associations = [associations]
            
            for assoc in associations:
                target = assoc.get('TargetDisorder', {})
                if target:
                    interactions.add((
                        f"ORPHA:{disease_id}",
                        f"ORPHA:{target.get('OrphaCode', '')}",
                        target.get('Name', {}).get('#text', ''),
                        assoc.get('DisorderDisorderAssociationType', {}).get('Name', {}).get('#text', '')
                    ))
        
        return pd.DataFrame(
            interactions,
            columns=['ORPHAcode', 'Classif_id', 'Classif_name', 'Classif_type']
        )
    
    # =========================================================================
    # CLASSIFICATIONS
    # =========================================================================
    
    def traverse_node(
        self, 
        node: dict, 
        root_orpha_id: str, 
        root_orpha_name: str, 
        interactions: set
    ) -> None:
        """Recursively traverses a classification node."""
        current = node.get('Disorder', {})
        parent_id = f"ORPHA:{current.get('OrphaCode', '')}"
        parent_type = current.get('DisorderType', {}).get('Name', {}).get('#text', '')
        
        child_list = node.get('ClassificationNodeChildList', {})
        children = child_list.get('ClassificationNode', [])
        
        if children:
            if not isinstance(children, list):
                children = [children]
            
            for child in children:
                target = child.get('Disorder', {})
                child_id = f"ORPHA:{target.get('OrphaCode', '')}"
                child_type = target.get('DisorderType', {}).get('Name', {}).get('#text', '')
                
                interactions.add((
                    root_orpha_id, root_orpha_name,
                    parent_id, parent_type,
                    child_id, child_type
                ))
                
                # Recursion
                self.traverse_node(child, root_orpha_id, root_orpha_name, interactions)
    
    def df_classif(self) -> pd.DataFrame:
        """Builds classification DataFrame from a JSON file."""
        with open(self.input_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        
        interactions = set()
        
        classification = data['JDBOR']['ClassificationList']['Classification']
        root_node = classification['ClassificationNodeRootList']['ClassificationNode']
        root_disorder = root_node.get('Disorder', {})
        root_id = f"ORPHA:{root_disorder.get('OrphaCode', '')}"
        root_name = root_disorder.get('Name', {}).get('#text', '')
        
        child_nodes = root_node.get('ClassificationNodeChildList', {}).get('ClassificationNode', [])
        if not isinstance(child_nodes, list):
            child_nodes = [child_nodes]
        
        for node in child_nodes:
            self.traverse_node(node, root_id, root_name, interactions)
        
        return pd.DataFrame(
            list(interactions),
            columns=['root', 'root_name', 'parent_id', 'parent_type', 'child_id', 'child_type']
        )
    
    def merge_all_classif(self, path_classif: str | Path) -> pd.DataFrame:
        """Merges all classification files into one DataFrame."""
        path_classif = Path(path_classif)
        dfs = []
        
        for excel_file in path_classif.glob("*.xlsx"):
            if excel_file.name.startswith("~$"):
                continue
            log.info("Loading %s", excel_file.name)
            df = pd.read_excel(excel_file, index_col=0, engine='openpyxl')
            if not df.empty:
                dfs.append(df)
        
        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    
    # =========================================================================
    # YAML GENERATION
    # =========================================================================
    
    def build_yaml_rds(self, df_rd: pd.DataFrame, col_name: str) -> dict[str, Any]:
        """
        Extracts ORPHA codes from DataFrame and builds Snakemake config dict.
        """
        raw_codes = df_rd[col_name].dropna().astype(str)
        rds = {code.strip() for code in raw_codes if code.strip().startswith('ORPHA')}
        
        return {
            'n': len(rds),
            'param_RD': {i + 1: code for i, code in enumerate(sorted(rds))}
        }
    
    # =========================================================================
    #  OTHER
    # =========================================================================
    
    def df_omim_orpha(self, root: dict) -> pd.DataFrame:
        """Builds OMIM <-> ORPHA mapping."""
        interactions = []
        
        for disorder in root['JDBOR']['DisorderList']['Disorder']:
            disease_id = disorder['OrphaCode']
            refs = disorder.get('ExternalReferenceList', {}).get('ExternalReference', [])
            
            if not isinstance(refs, list):
                refs = [refs]
            
            for ref in refs:
                if ref.get('Source') == 'OMIM':
                    interactions.append((disease_id, ref.get('Reference', '')))
        
        log.info("Built df pd1 OMIM-ORPHA mapping")
        return pd.DataFrame(interactions, columns=['ORPHAcode', 'OMIM'])
    
    def build_df_prevalence(self) -> pd.DataFrame:
        """Builds DataFrame of prevalences."""
        with open(self.input_path, 'r', encoding='ISO-8859-1') as f:
            root = json.load(f)
        
        all_prev = [(d['orpha'], d['preval']) for d in root]
        return pd.DataFrame(all_prev, columns=['ORPHAcode', 'Estimated_prevalence'])
    
    def from_dict_to_df(self, dict_patient: dict, prefix: str) -> pd.DataFrame:
        """Converts a patient dict->HPOs into DataFrame."""
        interactions = [
            (f"{prefix}{key}", hpo)
            for key, hpos in dict_patient.items()
            for hpo in hpos
        ]
        return pd.DataFrame(interactions, columns=['Patient', 'hpo_id'])
    
    def build_noisy_patient(
        self, 
        list_random_hpo: list[str], 
        dict_omim: dict
    ) -> dict:
        """Generates noisy patients (half original HPO + half random)."""
        result = {}
        for key, value in dict_omim.items():
            half_stop = len(value) // 2
            half_value = list(value[half_stop:])
            
            while len(half_value) != len(value):
                random_hpo = list_random_hpo[np.random.randint(0, len(list_random_hpo))].strip()
                if random_hpo not in value:
                    half_value.append(random_hpo)
            
            result[key] = half_value
        return result
    
    def load_build_mm(self, json_to_extract: list[str]) -> pd.DataFrame:
        """Loads and builds MM matrix from Excel files."""
        all_files = glob.glob(str(self.input_path / "*.xlsx"))
        
        filtered_files = [
            f for f in all_files
            if any(motif in f for motif in json_to_extract)
        ]
        
        df_all = pd.concat(
            (pd.read_excel(f) for f in set(filtered_files)),
            ignore_index=True
        )
        
        df_matrix = df_all.pivot(index='j', columns='i', values='resnik')
        return df_matrix.loc[json_to_extract, json_to_extract]
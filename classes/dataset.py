
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


from pyhpo import Ontology   

def minimal_hpo_ids(hp_ids):
    """
    Remove redundant ancestors:
    if an HPO term is an ancestor of another term in the set, drop it.
    Keep only the most specific terms.
    """

    ## first loop to extract all parents of all hpo in the hpo list 
    all_parents = set()
    for hp in hp_ids:
        pyhpo_term = Ontology.get_hpo_object(hp)

        # all ancestors (parents, grandparents, ...)
        ancestors = pyhpo_term.parents
        # print(f'{hp} : \t {ancestors}')
        for anc in ancestors:
            # from hpo term type to str 
            anc_id = anc.id
            all_parents.add(anc_id)

    # extract the hpo that are parents inside the hpo list
    hp_ids_that_are_parent = all_parents.intersection(hp_ids)

    ## keep only hpo that are not in hp_ids_that_are_parent
    # minimal_set = []
    # for hp in hp_ids:
    #     if hp not in hp_ids_that_are_parent:
    #         minimal_set.append(hp)
    
    minimal_set = set(hp_ids).difference(hp_ids_that_are_parent)
    return minimal_set


class DataSet(DataGenerate):
    """Class for construction and manipulation of Orphanet/patient datasets."""
    
    def __init__(self, input_path: str | Path, output_path: str | Path = ""):
        super().__init__(input_path, output_path)

    # =========================================================================
    # PATIENTS
    # =========================================================================
    
    def build_patients_df(self,do_subsumed) -> pd.DataFrame:
        """
        Builds a DataFrame of patients from phenopacket files.
        Handles variable JSON formats and extracts HPO, diseases, genes.
        """
        # input_path = PV.PATH_INPUT_PATIENTS_FOLDER
        pheno_with_invalid_id = set()
        list_case_HPO_f = []

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
                    pass
                
                patient_rows = []  # <-- collect only this patient's rows here

                for pheno in phenotypes:
                    # don't take into account the excluded phenotype
                    if 'excluded' in pheno.keys():
                        # excluded == True : on ignore
                        pass
                    else:
                        try:
                            if 'type.label' in pheno:
                                type_data = pheno['type.label']
                                if type_data == 'Invalid id':
                                    pheno_with_invalid_id.add(patient_id)
                                    pass
                                if 'type.negated' in pheno:
                                    pass
                                id_hpo = type_data.get('id', '')
                                label_hpo = type_data.get('label', '')
                            else:
                                type_data = pheno.get('type', {})
                                if type_data.get('label') == 'Invalid id':
                                    pheno_with_invalid_id.add(patient_id)
                                    pass
                                if pheno.get('negated'):
                                    pass

                                id_hpo = type_data.get('id', '')
                                label_hpo = type_data.get('label', '')

                            if id_hpo and patient_id not in pheno_with_invalid_id:
                                patient_rows.append((
                                    patient_id, progress_status, disease_orphanet, disease_omim,
                                    ern, gene, type_gene, variant, id_hpo, label_hpo
                                ))
                        except (KeyError, TypeError):
                            pass
                if do_subsumed == 1:
                    # ---- REMOVE REDUNDANT TERMS (minimal set) ----
                    hp_ids_in_order = [row[8] for row in patient_rows]
                    keep_ids = minimal_hpo_ids(hp_ids_in_order)
                    
                    for one_row in patient_rows:
                        if one_row[8] in keep_ids:
                            list_case_HPO_f.append(one_row)
                else:
                    list_case_HPO_f.extend(patient_rows)


            except (json.JSONDecodeError, UnicodeDecodeError) as e:
                log.warning("Cannot read %s: %s", onefile, e)
                continue
        
        return pd.DataFrame(
            list_case_HPO_f,
            columns=['phenopacket', 'status', 'Orphanet', 'OMIM', 'ern', 
                    'gene', 'type_gene', 'variant', 'hpo_id', 'hpo_label']
        )
    


    # get the orphanet confirmed (RDI) for each patient
    def parse_orpha_file(self, txt_path)-> pd.DataFrame:
        """
        Parse a file of the form:
        P1,ORPHA:610
        P2,ORPHA:35689
        """
        orpha_map = {}

        with open(txt_path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                patient_id, orpha = line.split(",")
                orpha_map[patient_id] = orpha

        return orpha_map
        
    # =========================================================================
    # ORPHANET PRODUCTS
    # =========================================================================
    # ok 
    def df_pd4(self) -> pd.DataFrame:
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
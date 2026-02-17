# Phenotype-based Rare Disease Prioritisation Pipeline

## Project overview

This repository provides a phenotype-based computational pipeline for prioritising rare diseases (ORPHAcodes) from patient phenotypes. The method combines semantic similarity between clinical sign from the Human Phenotype Ontology, with network-based propagation over the Orphanet nomenclature.

The pipeline computes disease-disease and disease-patient semantic similarity using asymmetric aggregation of HPO terms, accounting for term subsumption and Orphanet frequency annotations. Candidate diseases are then ranked by integrating patients into a disease similarity network and applying Random Walk with Restart to improve clinical consistency among ranked hypotheses.

The approach was evaluated on expert-curated cases from the Solve-RD project and demonstrated improved diagnostic prioritisation and more clinically coherent candidate lists compared to baseline phenotype-matching methods. This pipeline is designed to support researchers and developers working on rare disease diagnostics hypothesis and phenotypic data integration.

See MM section for explanation of semantic similarity,how aggregation methods differ and how RARW work for this project.

---

## Target audience

This project is intended for:
- Researchers working on rare disease diagnostics, semantic similarity, or phenotype-driven analysis
- Developers building bioinformatics pipelines involving HPO, Orphanet, and network-based methods

The pipeline assumes familiarity with:
- Command-line environments
- Snakemake workflows
- YAML configuration files
- Basic concepts in semantic similarity and graph-based analysis

---

## Requirements
**Core environment**
Python 3.12.9
Snakemake 9.3

**Python dependencies**
- pandas 2.2.3 
- numpy 1.26.4
- xmltodict 0.14.2
- hpo3 1.4.0 (install with pip)
- networkx 3.4.2
- matplotlib 3.10.1  
- openpyxl 3.1.5 
- scikit-learn 1.8.0
- pyarrow
- fastparquet

**Standard library usage**
The pipeline also relies on Python standard libraries, including:
math, os, json, re, logging, warnings, pathlib, argparse, time, datetime, glob, sys, yaml.

## Pipeline overview

1. Input data   
2. Disease-Disease semantic similarity   
3. Disease-Patient semantic similarity   
4. Patient integration into the disease network 
5. Random Walk with Restart 
6. Results part 

### Prerequisites
Before running any step, ensure the following:
**Step 1 Requirements:**
-  Python 3.12.9 and Snakemake 9.3 installed
-  All Python dependencies installed (see Requirements section)
-  HPO data files present in `input/hpo/`
-  Orphanet product files present in `input/pd_orphanet/`
-  All patients are added in the `patient/` 
-  The text file is filled to get the ORPHAcode confirmed for each patient `patient_RDI.txt`
Example : 
```text
P1,ORPHA:610
P2,ORPHA:100985
P3,ORPHA:35689
```
**Step 2-4 Requirements:**
-  Step 1 completed successfully
-  Output files from Step 1 exist in `output/`
-  Complete the config file `config_sm_mm.yaml` and `config_sm_mp.yaml`
**Step 5-6 Requirements:**
-  Steps 2 and 3 completed
-  Complete the config file `config_add_rw.yaml`

---

## Project architecture 
```text
[project_name]/
├── bin/
│
├── classes/
│
├── configs/
│
├── logs/
│
├── input/
│ ├── hpo/
│ ├── patient/
│ └── pd_orphanet/
│   └── Classifications/
├── output/           # Auto-generated
│   ├── mm_sm/
│   ├── mp_sm/
│   ├── patient_solverd/
│   ├── pd_orphanet/
│   ├── patient_added/
│   └── rarw/
```

All commands must be executed from the `[project_name]/` directory.

### `bin/`
Contains all executable elements of the pipeline:
- Python scripts
- YAML configuration files

### `input/`
The input directory **must strictly follow the expected architecture**.  
Do **not rename, move, or restructure** folders or files.

#### `input/hpo/`
- Contains three HPO  files downloaded from the official HPO website
- Required by the `hpo3` Python library

#### `input/patient/study_population/`
- Patients should be added in this folder
- Contains patient files in **Phenopacket format**
- One or multiple patients
- Used for disease-patient similarity and random walk steps

#### `input/pd_orphanet/`
- Contains Orphanet product DDL files
- Downloaded from Orphadata
- Required to compute disease-related similarity matrices

#### `[project_name]`
- Snakemake workflow files (`Snakefile.*`)
- readme.md
- path_variable.py which set all path of the projet (no need to touch)
- set_log.py set the log part 

---



## Step-by-step usage

### Step 1 - Load and normalise input data

```bash
snakemake -s Snakefile.load_input --cores 8  # Replace 8 with your desired number of cores

# python commands if you don't want to use snakemake
python -m bin.main_load_data
python -m bin.main_create_patient
 

```
 


Executes Snakefile.load_input
Runs main_load_data.py and main_create_patient.py
Converts raw input files into the internal standardised format used by the pipeline

### Patients format type 
```json
{
  "id": "P1",
  "subject": {
    "id": "P1",
    "alternateIds": [
      "xxx"
    ],
    "vitalStatus": {
      ...
    },
    ...
  },

  "phenotypicFeatures": [
    {
      "type": {
        "id": "HP:xxx",
        "label": "xxxx"
      },
      "excluded": true
    },
    {
      "type": {
        "id": "HP:yyy",
        "label": "yyy"
      }
    },
    ...
    
  ],

  "interpretations": [
    {
      ...
    }
  ],
  "metaData": {
    "resources": [
      {
        ...
      },
      ...
    ],
    ...
  }
}

```
The most important part here is `phenotypicFeatures` and `id`, others sections are not mandaroty.


### Step 2 - Build the Disease-Disease  similarity matrix
Configuration file to change : `config_sm_mm.yaml`

```yaml
mode: mm  # do not change

# list format can add multiple items
combine: ["aggregation_method_name"]
sm_method: ["similarity_measure_name"]
vector_strs: ["weight_vector"] # format 1_1_1_1_1

product4: name-to-define-which-pd4-used #  Ex: pd4may2025exejan2026 avoid '_'

# optional
mini_rd_csv: "ORPHA:xxx,ORPHA:xxx..." 
# example "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949"
```
- `combine`: semantic aggregation function (multiple allowed)
- `sm_method`: pairwise semantic similarity measure (multiple allowed)
- `vector_strs`: weighting vector definitions
- `product4`: identifier for the Orphanet product version
- `mini_rd_csv` (optional): restricts computation to a subset of ORPHAcodes instead of taking all ORPHAcodes from the product4

**Execution**
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mm.yaml --cores[n]

#or with python command example 
python -m bin.main_sm mp  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985" --mini-patient "P1"

python -m bin.main_sm mm  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985" 

python -m bin.main_concat concat_mp -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name 

python -m bin.main_concat concat_mm -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name

# Go get more info related to python parameters :
# python -m bin.main_sm  --help
# python -m bin.main_sm mp --help 
# python -m bin.main_sm mm --help 
# python -m bin.main_concat --help 
# python -m bin.main_concat concat_mp --help 
# python -m bin.main_concat concat_mm --help 


```

Output files are saved in `output/mm_sm/` with the following structure:
- Individual parquet files per disease into folder [aggregation_method]/[similariy_measure]/[n]/[name_to_define_which_pd4_used]/[weight_vector] (format: `{index}_{ORPHA_code}.parquet`)
- Concatenated file: `{combine}_{method}_{product4}_{vector}.parquet` 
  (Contains all disease-disease similarity scores)


### Weight Vector Explanation

The weight vector format is: `obligate_veryfreq_freq_occ_veryrare`

Example: `2_2_2_2_1`
- Position 1: Weight for Obligate frequency terms (default: 1.0)
- Position 2: Weight for Very Frequent terms (default: 1.0)
- Position 3: Weight for Frequent terms (default: 1.0)
- Position 4: Weight for Occasional terms (default: 1.0)
- Position 5: Weight for Very Rare terms (default: 1.0)

Higher values = emphasise that frequency category


### Step 3 - Build the Disease-Patient  similarity vectors
Configuration file to change: `config_sm_mp.yaml`

```yaml
mode: mp  # do not change

# list format can add multiple items
combine: ["aggregation_method_name"]
sm_method: ["similarity_measure_name"]
vector_strs: ["weight_vector"] # format 1_1_1_1_1


product4: name_to_define_which_pd4_used

# optional (one or multiple items)
mini_rd_csv: "ORPHA:xxx,ORPHA:xxx..." 
mini_patient_csv: "P1,P2..."

```
Same configuration logic as MM
Computed for one or multiple patients, if command mini_rd_csv restrict to only a subset of patients if not take all patients.
Restricted to diseases present in the MM matrix 

**Execution**
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mp.yaml --cores 8  # Replace 8 with your desired number
```

Output files are saved in `output/mp_sm/` with the following structure:
- Individual parquet files per disease into folder [aggregation_method]/[similariy_measure]/[n]/[name_to_define_which_pd4_used]/[weight_vector] (format: `{index}_{ORPHA_code}.parquet`)
- Concatenated file: `{combine}_{method}_{product4}_{vector}.parquet` 
  (Contains all patient-disease similarity scores)

Additional file :`RDI_{combine}_{method}_{product4}_{vector}.xlsx`,  contains, for each patient, the disease with the highest disease-patient similarity score

### Step 4 and 5 - Add patients to the MM matrix and apply Random Walk with Restart
Configuration file to change: `config_add_rw.yaml`

```yaml
alphas: [value_alpha]  

```

Set the value of alpha for the random walk, multiple value for multiple walk can be done, per defaults is 0.3.

**Execution**
```bash
snakemake -s Snakefile.add_rw --configfile configs/config_add_rw.yaml --cores 8  # Replace 8 with your desired number

#or with python commands
python -m bin.main_add_patients_to_mm
python -m bin.main_rarw run -a 0.3 --seeds "P1"
python -m bin.main_rarw collect -a 0.3

#python help for parameters
#python -m bin.main_rarw --help
#python -m bin.main_rarw run --help
#python -m bin.main_rarw collect --help



```

Executes main_add_patients_to_mm.py and main_rarw.py
Integrates patient nodes into the disease similarity matrix
Applies Random Walk with Restart using NetworkX library.


Output are both in patient_added/ and rarw/. in parquet format.
patient_added/ have the matrix DD with one patient added for each patient a parquet file is created. A `PROVENANCE.yaml` file is generated which describe the config done.
rarw/ contain the results of the ranking of each patients in `rarw/[value_alpha]/[config_name]/*.parquet`

Note :Only patients available on `mp_mm` are selected.
 

### Step 6 - Results execution and provenance

The Results section of the paper is reproducible using the scripts orchestrated by the Snakefile.rslt
These analyses are computed from the final outputs of the semantic similarity (RA) and random-walk  (RARW) pipelines.

The first results script (`results_CDFs_hm.py`) aggregates all ranking outputs across patients and methods, extracts the rank of the true diagnosis (RDI), and computes global performance summaries, including the harmonic mean of ranks and empirical cumulative distribution functions (CDFs). These analyses correspond to the global ranking comparison presented in the Results section (Tables 1–3 and Figure CDF).

The second results script (`results_GDs.py`) produces the classification-based comparisons reported in the Results section (Tables 4 and 5).
**Execution**
```bash
snakemake -s Snakefile.rslt --cores 8  # Replace 8 with your desired number

#or with python commands
python -u bin/results_CDFs_hm.py 
python -u bin/results_GDs.py 

```

----


### Computational considerations
Building full disease-disease similarity matrices can be computationally expensive.
For development, testing, or exploratory analyses, it is strongly recommended to restrict computations using the `mini_rd_csv` parameter.
Example : 
```yaml
mini_rd_csv: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949"

```
The Parquet output format is recommended for large-scale analyses due to improved performance and reduced disk usage.

## Example Workflow
### Quick Start (Subset of 10 Diseases)

1. **Load data**
```bash
snakemake -s Snakefile.load_input --cores 4
```

2. **Build MM matrix for 10 diseases**

Edit `config_sm_mm.yaml`:
```yaml
mode: mm
combine: ["funSimMax"]
sm_method: ["resnik"]
vector_strs: ["1_1_1_1_1"]
product4: pd4_may_2025_exejan2026
mini_rd_csv: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949,ORPHA:610"
```
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mm.yaml --cores 4
```

3. **Build MP vector**

Edit `config_sm_mp.yaml`:
```yaml
mode: mp
combine: ["rsd"]
sm_method: ["resnik"]
vector_strs: ["2_2_2_2_1"]
product4: pd4_may_2025_exejan2026
mini_rd_csv: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949,ORPHA:610"
mini_patient_csv: "P1"
```
```bash
snakemake -s Snakefile.sim  --configfile configs/config_sm_mp.yaml --cores 4
```


4-5. **Add patient to MM matrix and run RARW**

Edit `config_add_rw.yaml`:
```yaml
mm_file: "funSimMax_resnik_n_pd4_may_2025_exejan2026_11111.parquet"
sm_file: "rsd_resnik_n_pd4_may_2025_exejan2026_11111.parquet"

alpha: 0.3

output_mode: joint
```
```bash
snakemake -s Snakefile.add_rw --configfile configs/config_add_rw.yaml --cores 4
```

6. **Results**
```bash
snakemake -s Snakefile.rslt --cores 4
```

## Troubleshooting

### Common Errors

**Error: "HPO files not found"**
- Solution: Verify `input/hpo/` contains three HPO files from official source

**Error: "No patients found"**
- Solution: Ensure Phenopacket files are in `input/patient/study_population/`

**Error: "YAML parsing error"**
- Solution: Check YAML indentation (use spaces, not tabs)

**Memory issues on large datasets**
- Solution: Use `mini_rd_csv` and/or `mini_patient_csv` to restrict to smaller subset
- Increase `--cores` and available RAM

## Citation
If you use this pipeline in an academic or scientific context, please cite the corresponding publication describing the methodology and its evaluation on Solve-RD data.

A reference to the article will be added here upon publication.

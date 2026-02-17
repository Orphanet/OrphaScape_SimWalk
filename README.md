# Phenotype-based Rare Disease Prioritisation Pipeline

> **Reproducibility notice**  
> This repository contains the code used to produce all results, tables, and figures in:  
> **"Combining phenotypic similarity and network propagation to improve performance and clinical consistency of rare disease diagnosis"** - Chahdil, M. et al.  
>
> The original patient data used in the study cannot be shared due to patient confidentiality constraints. Three **simulated patient datasets** are provided as substitutes to allow reviewers to run and verify the pipeline end-to-end.

---

## Project overview

This repository provides a phenotype-based computational pipeline for prioritising rare diseases (ORPHAcodes) from patient phenotypes. The method combines semantic similarity between clinical signs from the Human Phenotype Ontology (HPO) with network-based propagation over the Orphanet nomenclature.

The pipeline computes disease-disease and disease-patient semantic similarity using asymmetric aggregation of HPO terms, accounting for term subsumption and Orphanet frequency annotations. Candidate diseases are then ranked by integrating patients into a disease similarity network and applying Random Walk with Restart (RWR) to improve clinical consistency among ranked hypotheses.

The approach was evaluated on expert-curated cases from the Solve-RD project and demonstrated improved diagnostic prioritisation and more clinically coherent candidate lists compared to baseline phenotype-matching methods.

> See the Materials and Methods section of the paper for a detailed explanation of semantic similarity measures, aggregation methods, and how RWR is applied in this pipeline.

---

## Target audience

This project is intended for:
- Researchers working on rare disease diagnostics, semantic similarity, or phenotype-driven analysis
- Developers building bioinformatics pipelines involving HPO, Orphanet, and network-based methods
- Reviewers wishing to reproduce the results of the above publication

The pipeline assumes familiarity with:
- Command-line environments
- Snakemake workflows
- YAML configuration files
- Basic concepts in semantic similarity and graph-based analysis

---

## About the provided simulated data

Because the original Solve-RD patient data are confidential and cannot be redistributed, **three simulated patient datasets** are included in this repository to enable full pipeline execution. These datasets are synthetic: they were generated to be structurally consistent with real Phenopacket files and are suitable for reproducing the pipeline steps and verifying code correctness.

> **Note:** Results obtained with simulated data will differ from those reported in the paper, which used real expert-curated cases. The simulated data are provided solely for code verification and reproducibility review purposes.

---

## Requirements
**Core environment**
- Python 3.12.9
- Snakemake 9.3

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

## Pipeline overview

1. Input data loading and normalisation
2. Disease-Disease semantic similarity
3. Disease-Patient semantic similarity
4. Patient integration into the disease network
5. Random Walk with Restart
6. Results

### Prerequisites

Before running any step, ensure the following:

**Step 1 requirements:**
- Python 3.12.9 and Snakemake 9.3 are installed
- All Python dependencies are installed (see Requirements section)
- HPO data files are present in `input/hpo/`
- Orphanet product files are present in `input/pd_orphanet/`
- Patient files (real or simulated) are placed in `input/patient/`
- The confirmed ORPHAcode per patient are declared in `patient_RDI.txt`

Example `patient_RDI.txt`:
```text
P1,ORPHA:610
P2,ORPHA:100985
P3,ORPHA:35689
```

**Steps 2-4 requirements:**
- Step 1 completed successfully
- Output files from Step 1 present in `output/`
- Configuration files `config_sm_mm.yaml` and `config_sm_mp.yaml` filled in

**Steps 5-6 requirements:**
- Steps 2 and 3 completed
- Configuration file `config_add_rw.yaml` filled in

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
│   ├── hpo/
│   ├── patient/
│   └── pd_orphanet/
│       └── Classifications/
├── output/           # Auto-generated
│   ├── mm_sm/
│   ├── mp_sm/
│   ├── patient_solverd/
│   ├── pd_orphanet/
│   ├── patient_added/
│   └── RWR/
```

All commands must be executed from the `[project_name]/` directory.

### `bin/`
Contains all executable elements of the pipeline: Python scripts and YAML configuration files.

### `input/`
The input directory **must strictly follow the expected architecture**.  
Do **not rename, move, or restructure** folders or files.

#### `input/hpo/`
Contains three HPO files downloaded from the official HPO website. Required by the `hpo3` Python library.

#### `input/patient/study_population/`
Contains patient files in **Phenopacket format** (one file per patient). The three simulated patient files are already placed here. Used for disease-patient similarity and random walk steps.

#### `input/pd_orphanet/`
Contains Orphanet product DDL files downloaded from Orphadata. Required to compute disease-related similarity matrices.

#### `[project_name]/` (root)
- Snakemake workflow files (`Snakefile.*`)
- `README.md`
- `path_variable.py` - sets all project paths (do not edit)
- `set_log.py` - configures logging

---

## Step-by-step usage

### Step 1 - Load and normalise input data

```bash
snakemake -s Snakefile.load_input --cores n  #n The number of core Laurent: Meaning AT MOST 8 cores. Otherwise max. Can be removed

# python commands if you don't want to use snakemake
python -m bin.main_load_data
python -m bin.main_create_patient
```

Converts raw input files (HPO, Orphanet, patients) into the internal standardised format used by the pipeline.

### Patient file format

Patient files must follow the **Phenopacket** format (JSON). The most important fields are `id` and `phenotypicFeatures`. All other sections are optional.

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

---

### Step 2 - Build the Disease-Disease (DD) similarity matrix

Configuration file: `configs/config_sm_mm.yaml`

```yaml
mode: mm  # do not change

# List format - multiple values allowed
combine: ["aggregation_method_name"]
sm_method: ["similarity_measure_name"]
vector_strs: ["weight_vector"]  # format: 1_1_1_1_1

product4: name-to-define-which-pd4-used  # e.g. pd4may2025exejan2026 - avoid underscores '_'

# Optional: restrict to a subset of ORPHAcodes
mini_rd_csv: "ORPHA:xxx,ORPHA:xxx..."

# Example: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949"
```

**Parameter descriptions:**
- `combine`: semantic aggregation function (e.g. `funSimMax`, `BMA`, `rsd`)
- `sm_method`: pairwise semantic similarity measure (e.g. `resnik`)
- `vector_strs`: weighting vector (see Weight Vector section below)
- `product4`: identifier for the Orphanet product version used
- `mini_rd_csv` *(optional)*: restricts computation to a subset of ORPHAcodes instead of the full Orphanet catalogue `product4`

**Execution:**
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mm.yaml --cores n

# Alternative: run Python scripts directly
python -m bin.main_sm mp  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985" --mini-patient "P1"
python -m bin.main_sm mm  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985" 
python -m bin.main_concat concat_mp -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name 
python -m bin.main_concat concat_mm -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name

# Help
# python -m bin.main_sm  --help
# python -m bin.main_sm mp --help 
# python -m bin.main_sm mm --help 
# python -m bin.main_concat --help 
# python -m bin.main_concat concat_mp --help 
# python -m bin.main_concat concat_mm --help 
```

**Output** (saved in `output/mm_sm/`):
- Individual parquet files per disease into folder [aggregation_method]/[similariy_measure]/[n]/[name_to_define_which_pd4_used]/[weight_vector] (format: `{index}_{ORPHA_code}.parquet`)
- Concatenated file: `{combine}_{method}_{product4}_{vector}.parquet` (all disease-disease similarity scores)

---

### Weight Vector Explanation

The weight vector format is: `obligate_veryfreq_freq_occ_veryrare`

Example: `2_2_2_2_1`
- Position 1: Weight for **Obligate** frequency terms
- Position 2: Weight for **Very Frequent** terms
- Position 3: Weight for **Frequent** terms
- Position 4: Weight for **Occasional** terms
- Position 5: Weight for **Very Rare** terms

Default weight is `1.0` for all positions. Higher values increase emphasis on that frequency category.

---

### Step 3 - Build the Disease-Patient (DP) similarity vectors

Configuration file to change: `config_sm_mp.yaml`

```yaml
mode: mp  # do not change

# list format can add multiple items
combine: ["aggregation_method_name"]
sm_method: ["similarity_measure_name"]
vector_strs: ["weight_vector"] # format 1_1_1_1_1

product4: name_to_define_which_pd4_used

# Optional: restrict to a subset of diseases and/or patients
mini_rd_csv: "ORPHA:xxx,ORPHA:xxx..." 
mini_patient_csv: "P1,P2..."
```

Same configuration logic as Step 2. Computation is restricted to diseases present in the DD matrix.

**Execution:**
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mp.yaml --cores n
```

**LAURENT A CHECKER PAR MAROUA:**
```
# Alternative: run Python scripts directly
python -m bin.main_sm mp  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985" --mini-patient "P1"
python -m bin.main_concat concat_mp -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name

# Help
python -m bin.main_sm mp --help
python -m bin.main_concat concat_mp --help
```

**Output** (saved in `output/mp_sm/`):
- Individual parquet files per disease (same folder structure as DD)
- Concatenated file: `{combine}_{method}_{product4}_{vector}.parquet` (all patient-disease similarity scores)
- Additional file: `RDI_{combine}_{method}_{product4}_{vector}.xlsx` - for each patient, the disease with the highest similarity score

---

### Steps 4 & 5 - Patient integration and Random Walk with Restart

Configuration file to change: `config_add_rw.yaml`

```yaml
alphas: [0.3]  # Alpha value(s) for random walk. Multiple values run multiple walks.

```

**Execution:**
```bash
snakemake -s Snakefile.add_rw --configfile configs/config_add_rw.yaml --cores 8

# Alternative: run Python scripts directly
python -m bin.main_add_patients_to_mm
python -m bin.main_rarw run -a 0.3 --seeds "P1"
python -m bin.main_rarw collect -a 0.3

# Help
#python -m bin.main_rarw --help
#python -m bin.main_rarw run --help
#python -m bin.main_rarw collect --help
```

This step integrates patient nodes into the disease similarity matrix, then applies Random Walk with Restart using the NetworkX library.

**Output:**
- `output/patient_added/` - one parquet file per patient containing the DD matrix extended with that patient. A `PROVENANCE.yaml` file records the configuration used.
- `output/RWR/[alpha_value]/[config_name]/*.parquet` - ranked disease candidates per patient

> **Note:** Only patients present in both the DD and DP matrices are processed.

---

### Step 6 - Results

The results reported in the paper are reproducible using the scripts orchestrated by `Snakefile.rslt`, computed from the final outputs of Steps 2–5.

`results_CDFs_hm.py` aggregates all ranking outputs across patients and methods, extracts the rank of the true diagnosis (RDI), and computes global performance summaries including the harmonic mean of ranks and empirical cumulative distribution functions (CDFs). These correspond to Tables 1–3 and the CDF figure in the paper.

`results_GDs.py` produces the classification-based comparisons reported in Tables 4 and 5.

**Execution:**
```bash
snakemake -s Snakefile.rslt --cores 8

# Alternative: run Python scripts directly
python -u bin/results_CDFs_hm.py
python -u bin/results_GDs.py
```

----

## Computational considerations

Building full disease-disease similarity matrices can be computationally expensive. For development, testing, or exploratory analyses (including reproducibility review), it is strongly recommended to restrict computations using the `mini_rd_csv` parameter.

Example restriction:
```yaml
mini_rd_csv: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949"
```

The Parquet output format is recommended for large-scale analyses due to improved I/O performance and reduced disk usage.

---

## Example workflow (subset of 10 diseases)

This quick-start example uses a small subset of diseases and the provided simulated patients to verify the pipeline runs correctly end-to-end.

**1. Load data**
```bash
snakemake -s Snakefile.load_input --cores 4
```

**2. Build DD matrix (10 diseases)**

Edit `configs/config_sm_mm.yaml`:
```yaml
mode: mm
combine: ["funSimMax"]
sm_method: ["resnik"]
vector_strs: ["1_1_1_1_1"]
product4: pd4may2025exejan2026
mini_rd_csv: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949,ORPHA:610"
```
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mm.yaml --cores 4
```

**3. Build DP vectors**

Edit `configs/config_sm_mp.yaml`:
```yaml
mode: mp
combine: ["rsd"]
sm_method: ["resnik"]
vector_strs: ["2_2_2_2_1"]
product4: pd4may2025exejan2026
mini_rd_csv: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949,ORPHA:610"
mini_patient_csv: "P1"
```
```bash
snakemake -s Snakefile.sim --configfile configs/config_sm_mp.yaml --cores 4
```

**4–5. Patient integration and RWR**

Edit `configs/config_add_rw.yaml`:
```yaml
alphas: [0.3]
```
```bash
snakemake -s Snakefile.add_rw --configfile configs/config_add_rw.yaml --cores 4
```

**6. Results**
```bash
snakemake -s Snakefile.rslt --cores 4
```

---

## Troubleshooting

**"HPO files not found"**  
Verify that `input/hpo/` contains the three required HPO files from the official HPO website.

**"No patients found"**  
Ensure Phenopacket files are present in `input/patient/study_population/`.

**"YAML parsing error"**  
Check YAML indentation - use spaces, never tabs.

**Memory issues on large datasets**  
Use `mini_rd_csv` and/or `mini_patient_csv` to restrict computation to a smaller subset, and increase available RAM.

---

## Citation

If you use this pipeline in an academic or scientific context, please cite:

> Chahdil, M. et al. *Combining phenotypic similarity and network propagation to improve performance and clinical consistency of rare disease diagnosis.* (Publication details to be added upon journal release.)

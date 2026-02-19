# Phenotype-based Rare Disease Prioritisation Pipeline

> **Reproducibility notice**  
> This repository contains the code used to produce all results, tables, and figures in:  
> **"Combining phenotypic similarity and network propagation to improve performance and clinical consistency of rare disease diagnosis"** - Chahdil, M. et al.  
>
> The original patient data used in the study cannot be shared due to patient confidentiality constraints. Three **simulated patient datasets** are provided as substitutes to allow reviewers to run and verify the pipeline end-to-end.

---
## Table of contents

- [Project overview](#project-overview)
- [Target audience](#target-audience)
- [About the provided simulated data](#about-the-provided-simulated-data)
- [Requirements](#requirements)
- [Install](#install)
	- [Prerequisites](#prerequisites)
	- [Application](#application)
- [TL;DR:](#tl-dr)


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
- Python 3.12+

**Python dependencies**
- pandas 2.2-3.0+
- numpy 1.26-2.4+
- xmltodict 0.14-1.0+
- hpo3 1.4-1.5+
- networkx 3.4-3.6+
- matplotlib 3.10
- openpyxl 3.1
- scikit-learn 1.8
- snakemake 9.x
- pyarrow 23.0
- fastparquet 2025.12

## Install

### Prerequisites

You have to use a python 3.12 environment. The easiest way is to use `conda` which belong to the *Anaconda/Miniconda* distribution. We advise to use *Miniconda*: `https://www.anaconda.com/docs/getting-started/miniconda/main`.

Once miniconda is installed, from a *Terminal/Console*, type:

```bash
conda create --name simwalk -c conda-forge -c bioconda python=3.12 fastparquet pandas matplotlib networkx openpyxl pyarrow scikit-learn snakemake xmltodict
conda activate simwalk
pip install hpo3
```

### Application

From a terminal, you should `git clone` or download and unzip the application, than change your working directory to `OrphaScape_SimWalk`.


## TL;DR:

You don't have much time? Jump directly to section `Example workflow (subset of 10 diseases)` below.

## Pipeline overview

1. Input data loading and normalisation
2. Disease-Disease semantic similarity
3. Disease-Patient semantic similarity
4. Patient integration into the disease network
5. Random Walk with Restart
6. Results

Before running any step, ensure the following:

**Step 1 requirements:**
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

**Steps 5-6 requirements:**
- Steps 2 and 3 completed


## Run Configuration for all steps

The pipeline consists of 6 steps, all configured through a single YAML file located in `configs/`. Each time the pipeline runs, this one file drives all pipeline steps (`Snakefile.sim`, `Snakefile.rslt`, `Snakefile.add_rw`), so there's no need to edit multiple config files between runs.

### Config file structure
Here an exemple of how to run all pipeline steps (see run4.yaml in section 'Example workflow (subset of 10 diseases)' )
```bash
# configs/run_example.yaml

run_name: run_example            # output subfolder where figures and tables are: output/run_example/
fig_num: 1                       # figure suffix: CDF_fig{fig_num+1}.svg, hm_table{fig_num}.xlsx
do_subsumed: 0                   # 0 = raw HPO terms | 1 = subsumed HPO terms thus removing ancestors terms in the patient profil 

product4: pd4may2025exefev2026   # This refer to the orphanet product version identifier named depending on the user.
mini_rd: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:610,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949,ORPHA:35689"
   # optional: restricts computation to a subset of ORPHAcodes instead of the full Orphanet catalogue 

# optional: restricts computation to a subset of patients instead of the full set
mini_patient: "P1,P2" 

dd:                              # Disease-Disease similarity - omit section to skip
  combine: ["funSimMax"] #  aggregation function 
  sm_method: ["resnik"] # pairwise similarity measure
  vector_strs: ["1_1_1_1_1"]  # weighting vector

dp:                              # Disease-Patient similarity -- omit section to skip
  combine: ["rsd"] #  aggregation function 
  sm_method: ["resnik"] # pairwise similarity measure
  vector_strs: ["1_1_1_1_1", "2_2_2_2_1"]  # weighting vector
  # Restriction to certain patients (optional dp only )
  


alpha: 0.3                       #  Damping factor for PageRank

steps: ["cdf"]                   # which result steps to run:
                                 #   "cdf" → CDF plot + harmonic mean table
                                 #   "gd"  → clinical consistency -based classification (group of disorder) analysis (Tables 4 & 5)
```

This table map each config parameter to the Snakefile(s) and steps that use it:

| Parameter | Step 1 - Load data<br>`Snakefile.load_input` | Steps 2-3 - Similarity matrices<br>`Snakefile.sim` | Steps 4-5 - Patient integration + RWR<br>`Snakefile.add_rw` | Step 6 - Results<br>`Snakefile.rslt` |
|-----------|:---:|:---:|:---:|:---:|
| `do_subsumed`  | ✓ | ✓ | ✓ | ✓ |
| `run_name`     | x | x | x | ✓ |
| `fig_num`      | x | x | x | ✓ |
| `product4`     | x | ✓ | x | ✓ (filter) |
| `mini_rd`      | ✓ | ✓ (optional) | x | x |
| `mini_patient` | x | ✓ (Step 3 only optional) | x | ✓ |
| `dd` (section) | x | ✓ (Step 2 only) | x | x |
| `dp` (section) | x | ✓ (Step 3 only) | x | x |
| `vector_strs`        | x | ✓ | x | x |
| `sm_method`        | x | ✓ | x | x |
| `combine`        | x | ✓ | x | x |
| `alpha`        | x | x | ✓ | x |
| `steps`        | x | x | x | ✓ |



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
│   └── orphanet_data/
│       └── Classifications/
├── output/           # Auto-generated
│   ├── dd_sm/
│   ├── dp_sm/
│   ├── patient_solverd/
│   ├── pd_orphanet/
│   ├── patient_added/
│   ├── RWR/
│   └── run*/ # contain results depending on the run (depending on the configuration file )
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
Converts raw input files (HPO, Orphanet, patients) into the internal standardised format used by the pipeline.

```bash
# Requires --configfile to set do_subsumed (0 = raw HPO terms, 1 = subsumed HPO terms)
# Use any run config that contains do_subsumed, or the dedicated config:
snakemake -s Snakefile.load_input --configfile configs/run_example.yaml --cores all # The number of core  Meaning AT MOST 8 cores. Otherwise max. Can be removed

 
# Alternative: run Python scripts directly
python -m bin.main_load_data
python -m bin.main_create_patient --do-subsumed 0   # → output/patient_solved/patients.xlsx
python -m bin.main_create_patient --do-subsumed 1   # → output/patient_solved/patients_subsumed.xlsx

```


Parameter mandaroty in the **config** for this step: `do_subsumed`    


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

### Step 2 - Build the Disease-Disease (DD) similarity matrix (modif Maroua)
Computes   similarity measure between all ORPHAcodes based on their HPO phenotype annotations, then concatenates the per-ORPHAcode results into a single parquet file per  combination.

```bash
snakemake -s Snakefile.sim --configfile configs/run_example.yaml --cores all

# Alternative: run Python scripts directly
python -m bin.main_sm dd  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985"
python -m bin.main_concat concat_dd -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name

# Help
# python -m bin.main_sm  --help
# python -m bin.main_sm dp --help 
# python -m bin.main_sm dd --help 
# python -m bin.main_concat --help 
# python -m bin.main_concat concat_dp --help 
# python -m bin.main_concat concat_dd --help 
```
Parameter mandaroty in the **config** for this step: `do_subsumed`,`product4`,`dd`,`mini_rd`(optional),`vector_strs`,`sm_method`,`combine`.

**Output** (saved in `output/dd_sm/`):
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

### Step 3 - Build the Disease-Patient (DP) similarity vectors (modif Maroua)
Computes   similarity measure between each patient's HPO phenotype profile and all ORPHAcodes , then concatenates the per-ORPHAcode results into a single parquet file per  combination.


```bash
snakemake -s Snakefile.sim --configfile configs/run_example.yaml --cores all


# Alternative: run Python scripts directly
python -m bin.main_sm dp  -c BMA -m resnik -v "1_1_1_1_1" -pd4 pd4name --mini-rd "ORPHA:610,ORPHA:100985" --mini-patient "P1"
python -m bin.main_concat concat_dp -c BMA -m resnik -v "1_1_1_1_1" --pd4 pd4name

# Help
python -m bin.main_sm dp --help
python -m bin.main_concat concat_dp --help
```

Parameter mandaroty in the **config** for this step: `do_subsumed`,`product4`,`dp`,`mini_rd`(optional),`mini_patient` (optional),`vector_strs`,`sm_method`,`combine`.

**Output** (saved in `output/dp_sm/`):
- Individual parquet files per disease (same folder structure as DD)
- Concatenated file: `{sub0/1}_{combine}_{method}_{product4}_{vector}.parquet` (all patient-disease similarity scores)
- Additional file: `RDI_{combine}_{method}_{product4}_{vector}.xlsx` - for each patient, the disease with the highest similarity score

---

### Steps 4 & 5 - Patient integration and Random Walk with Restart (modif Maroua)
This step integrates patient nodes into the disease similarity matrix, then applies Random Walk with Restart using the NetworkX library.

**Execution:**
```bash
snakemake -s Snakefile.add_rw --configfile configs/run_example.yaml --cores all

# Alternative: run Python scripts directly
python -m bin.main_add_patients_to_dd
python -m bin.main_rarw run -a 0.3 --seeds "P1"
python -m bin.main_rarw collect -a 0.3

# Help
#python -m bin.main_rarw --help
#python -m bin.main_rarw run --help
#python -m bin.main_rarw collect --help
```
Parameter mandaroty in the **config** for this step: `do_subsumed`,`alpha`.


**Output:**
- `output/patient_added/` - one parquet file per patient containing the DD matrix extended with that patient. A `PROVENANCE.yaml` file records the configuration used.
- `output/RWR/[alpha_value]/[config_name]/*.parquet` - ranked disease candidates per patient

> **Note:** Only patients present in both the DD and DP matrices are processed.

---

### Step 6 - Results (modif Maroua)

The results reported in the paper are reproducible using the scripts orchestrated by `Snakefile.rslt`, computed from the final outputs of Steps 2-5.

```bash
snakemake -s Snakefile.rslt --configfile configs/run_example.yaml --cores all  

# Alternative: run Python scripts directly
python -u bin/results_CDFs_hm.py
python -u bin/results_GDs.py
```
Parameter mandaroty in the **config** for this step: `do_subsumed`,`alpha`.

**Output:**
`results_CDFs_hm.py` aggregates all ranking outputs across patients and methods, extracts the rank of the true diagnosis (RDI), and computes global performance summaries including the harmonic mean of ranks and empirical cumulative distribution functions (CDFs). These correspond to Tables 1-3 and the CDF figure in the paper.

`results_GDs.py` produces the classification-based comparisons reported in Tables 4 and 5.


----


## Computational considerations

Building full disease-disease similarity matrices can be computationally expensive. For development, testing, or exploratory analyses (including reproducibility review), it is strongly recommended to restrict computations using the `mini_rd` parameter.

Example restriction:
```yaml
mini_rd: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949"
```

The Parquet output format is used because it is recommended for large-scale analyses due to improved I/O performance and reduced disk usage.

---

## Example workflow (subset of 10 diseases) (modif Maroua)

This quick-start example uses a small subset of diseases and the provided simulated patients to verify the pipeline runs correctly end-to-end.

Four ready-to-use config files are provided in `configs/`, each corresponding to a section of the paper. They can be run independently.

---

**run1** - *"Resnik + FunSimMaxAsym outperforms other groupwise semantic similarity measures"*
```bash
snakemake -s Snakefile.load_input --configfile configs/run1.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run1.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run1.yaml --cores all
```

**Output:** `output/run1/`  Contains examples of the CDF  and the harmonic mean ranks Table related the to first results.

---

**run2** - *"Effect of subsumed HPO terms removal on ranking"*

This result requires two separate runs differing only in `do_subsumed` (0 = raw HPO terms, 1 = subsumed). 
```bash
snakemake -s Snakefile.load_input --configfile configs/run2_raw.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run2_raw.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run2_raw.yaml --cores all

snakemake -s Snakefile.load_input --configfile configs/run2_sub.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run2_sub.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run2_sub.yaml --cores all
```

**Output:** `output/run2_raw/` and `output/run2_sub/`   Contain both  examples of the CDF  and the harmonic mean ranks Table related the to seconds results.In the paper the two outputs were merged; here they are produced independently.

---

**run3** - *"Effect of frequency-based weighting on ranking performance"*
```bash
snakemake -s Snakefile.load_input --configfile configs/run3.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run3.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run3.yaml --cores all
```
**Output:** `output/run3/`  Contains examples of the CDF  and the harmonic mean ranks Table related the to third results.

---

**run4** - *"Effect of network propagation on clinical consistency of top-ranked candidates"*

Full pipeline including patient integration in the network and RWR.
```bash
snakemake -s Snakefile.load_input --configfile configs/run4.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run4.yaml --cores all
snakemake -s Snakefile.add_rw     --configfile configs/run4.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run4.yaml --cores all
```

**Output:** `output/run4/`  Contains examples of tables (table 4 and 5) related the to fourth results.

---

## Troubleshooting

**"HPO files not found"**  
Verify that `input/hpo/` contains the three required HPO files from the official HPO website.

**"No patients found"**  
Ensure Phenopacket files are present in `input/patient/`.

**"YAML parsing error"**  
Check YAML indentation - use spaces, never tabs.

**Memory issues on large datasets**  
Use `mini_rd` and/or `mini_patient` to restrict computation to a smaller subset, and increase available RAM.

---

## How to cite

If you use this pipeline in an academic or scientific context, please cite:

> Chahdil, M. et al. *Combining phenotypic similarity and network propagation to improve performance and clinical consistency of rare disease diagnosis.* (Publication details to be added upon journal release.)

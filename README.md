# Phenotype-based Rare Disease Prioritisation Pipeline

> **Reproducibility notice**  
> This repository contains the code used to produce all results, tables, and figures in:  
> **"Combining phenotypic similarity and network propagation to improve performance and clinical consistency of rare disease diagnosis"** - Chahdil, M. et al.  
>
> SimWalk was benchmarked on three patient cohorts: SHEPHERD (simulated), Phenopacket Store (published real cases) and Solve-RD (expert-curated real cases). The patient datasets derived from SHEPHERD and Phenopacket Store are provided in the [`dataset/`](dataset/) folder of this repository. The Solve-RD data cannot be shared due to patient confidentiality constraints. In addition, three small simulated patients (P1, P2, P3) are provided to run and verify the pipeline end-to-end in a few minutes.

---

<a name="top"></a>
## Table of contents

- [Project overview](#project-overview)
- [Target audience](#target-audience)
- [Datasets used in the paper](#datasets-used-in-the-paper)
- [About the provided simulated data](#about-the-provided-simulated-data)
- [Requirements](#requirements)
- [Install](#install)
	- [Prerequisites](#prerequisites)
	- [Application](#application)
- [TL;DR](#tldr)
- [How to: Example workflow (subset of 10 diseases)](#tldr) 
- [Pipeline overview](#pipeline-overview)
- [Run Configuration for all steps](#run-configuration-for-all-steps)
	- [Config file structure](#config-file-structure)
- [Project architecture](#project-architecture)
- [Step-by-step usage](#step-by-step-usage)
	- [Step 1 - Load and normalise input data](#step-1---load-and-normalise-input-data)
	- [Patient file format](#patient-file-format)
	- [Step 2 - Build the Disease-Disease (DD) similarity matrix](#step-2---build-the-disease-disease-dd-similarity-matrix)
	- [Weight Vector Explanation](#weight-vector-explanation)
	- [Step 3 - Build the Disease-Patient (DP) similarity vectors](#step-3---build-the-disease-patient-dp-similarity-vectors)
	- [Steps 4 & 5 - Patient integration and Random Walk with Restart](#steps-4--5---patient-integration-and-random-walk-with-restart)
	- [Step 6 - Results](#step-6---results)
- [Computational considerations](#computational-considerations)
- [Troubleshooting](#troubleshooting)
- [How to cite](#how-to-cite)

## Project overview

This repository provides SimWalk, a phenotype-based computational pipeline for prioritising rare diseases (ORPHAcodes) from patient phenotypes. The method combines semantic similarity between clinical signs from the Human Phenotype Ontology (HPO) with network-based propagation over the Orphanet nomenclature.

The pipeline computes disease-disease and disease-patient semantic similarity using asymmetric aggregation of HPO terms, accounting for term subsumption and Orphanet frequency annotations. Candidate diseases are then ranked by integrating patients into a disease similarity network and applying Random Walk with Restart (RWR) to improve clinical consistency among ranked hypotheses.

The approach was benchmarked on 6,183 simulated patients derived from the SHEPHERD dataset, and evaluated on two real-patient cohorts: 4,167 patients from Phenopacket Store and 117 expert-curated cases from the Solve-RD project. In a head-to-head comparison with four phenotype-driven diagnostic tools (Phenomizer, Phrank, RunSolveRD and PubCaseFinder), SimWalk achieved the highest performance on SHEPHERD and Phenopacket Store, and RWR further improved the clinical consistency of the candidate lists.


> See the Materials and Methods section of the paper for a detailed explanation of semantic similarity measures, aggregation methods, and how RWR is applied in this pipeline.

[⬆ Back to top](#top)


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

[⬆ Back to top](#top)

## Datasets used in the paper

All analyses use the HPO release of 2026-09-01 and the Orphadata "Phenotypes Associated with Rare Disorders" dataset of June 2026 (4,355 ORPHAcodes annotated with 8,758 distinct HPO terms). All diagnoses are harmonised to the Disorder level of the Orphanet nomenclature. Only patients whose diagnosis is annotated in this Orphadata dataset are analysed.



### SHEPHERD simulated cohort

Source: the simulated patients used to train and validate SHEPHERD (Alsentzer et al., *npj Digital Medicine* 2025), available on Harvard Dataverse: `https://doi.org/10.7910/DVN/TZTPFL`.

- 42,680 simulated patients have a diagnosis mapped to 2,134 ORPHAcodes. From the 6402 patients, 219 were removed because their diagnosis is not annotated in the Orphadata phenotype dataset. 6,183 patients (on average 3 per ORPHAcode) were sampled, preserving the ORPHAcode distribution of the original dataset.
  
- Three sparser versions were generated by randomly removing 10% (SHEPHERD90), 30% (SHEPHERD70) and 50% (SHEPHERD50) of the HPO terms of each patient. SHEPHERD50 is the validation set used for most analyses of the paper.

### Phenopacket Store cohort

Source: GA4GH Phenopacket Store (Danis et al., *HGG Advances* 2024), release 0.1.27 (June 2026), 10,377 curated cases:`https://github.com/monarch-initiative/phenopacket-store/releases/download/0.1.27/all_phenopackets.zip`

- the OMIM diagnosis is mapped to ORPHAcodes using Orphanet data, Exact Match (E) only; cases without an exact mapping are removed. 4,167 of them have a diagnosis annotated in the Orphadata phenotype dataset and are analysed.

### Solve-RD cohort

Ssolved Solve-RD patients had a confirmed diagnosis mapped to an ORPHAcode; an Orphanet medical expert selected 117 patients whose phenotype was consistent with their diagnosis. The phenopackets are deposited at the EGA (EGAD00001009767, EGAD00001009768, EGAD00001009769 and EGAD00001009770, study EGAS00001003851) but cannot be redistributed.

[⬆ Back to top](#top)





	
## About the provided simulated data

Because the original Solve-RD patient data are confidential and cannot be redistributed, **three simulated patient datasets** are included in this repository to enable full pipeline execution. These datasets are synthetic: they were generated to be structurally consistent with real Phenopacket files and are suitable for reproducing the pipeline steps and verifying code correctness.

To reproduce the results of the paper on SHEPHERD and Phenopacket Store, copy the Phenopacket files of one cohort from `dataset/` into `input/patient/` and its diagnoses into `input/patient_RDI.txt` (see [Pipeline overview](#pipeline-overview)).

>  **Note:** Results obtained with the three example patients will differ from those reported in the paper. They are provided solely for code verification. The results on Solve-RD cannot be reproduced because these data are confidential.
	
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

[⬆ Back to top](#top)
	
## Install

### Prerequisites

You have to use a python 3.12 environment. The easiest way is to use `conda` which belong to the *Anaconda/Miniconda* distribution. We advise to use *Miniconda*: `https://www.anaconda.com/docs/getting-started/miniconda/main`.

Once miniconda is installed, from a *Terminal/Console*, type:

```bash
conda create --name simwalk -c conda-forge -c bioconda python=3.12 fastparquet pandas matplotlib networkx openpyxl pyarrow scikit-learn snakemake xmltodict
conda activate simwalk
pip install hpo3
```


### Build OrphaScape_SimWalk environment

In a terminal (or Windows Powershell), clone repo or download and unzip it, then move your path to the working directory "OrphaScape_SimWalk".
 
### /!\ WARNING FOR WINDOWS USERS /!\ 

If you use the github clone facility, before proceeding it, be sure to desactivate the autocrlf option:

```bash
git config --global core.autocrlf false
git clone https://github.com/Orphanet/OrphaScape_SimWalk.git

``` 

You can also use the download and unzip project solution. 


## TL;DR:

You don't have much time? Jump directly to section `How to: Example workflow (subset of 10 diseases)` below.

	
	
## How to: Example workflow (subset of 10 diseases)

This quick-start example uses a small subset of diseases and the provided simulated patients to verify the pipeline runs correctly end-to-end.

Four ready-to-use config files are provided in `configs/`, each corresponding to a section of the paper. They can be run independently.

Be sure that you are already running a virtual environnement (e.g. conda activate my_env)

---

### Clean all outputs  

**snakemake to clean all output and start from scratch**
```bash
snakemake -s Snakefile.clean clean --cores all
```
---

### Start the example workflow  

#### **run1** - *"Jiang-Conrath + FunSimAvg outperforms other groupwise semantic similarity measures"*
```bash
snakemake -s Snakefile.load_input --configfile configs/run1.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run1.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run1.yaml --cores all
```

Output: `output/run1/`  Contains examples of the CDF  and the harmonic mean ranks Table related to the first results (Table 2 and Figure 8 of the paper).
 

#### **run2** - *"Effect of subsumed HPO terms removal on ranking"*

This result requires two separate runs differing only in `do_subsumed` (0 = raw HPO terms, 1 = subsumed). 
```bash
snakemake -s Snakefile.load_input --configfile configs/run2_raw.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run2_raw.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run2_raw.yaml --cores all

snakemake -s Snakefile.load_input --configfile configs/run2_sub.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run2_sub.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run2_sub.yaml --cores all
```

Output: `output/run2_raw/` and `output/run2_sub/`   Contain both  examples of the CDF  and the harmonic mean ranks table related to the second results (Figure 9 of the paper: SimWalk WR = `run2_raw`, SimWalk WOR = `run2_sub`). 
In the paper the two outputs were merged; here they are produced independently.
 

#### **run3** - *"Frequency-based weighting as an optional tuning parameter""*
```bash
snakemake -s Snakefile.load_input --configfile configs/run3.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run3.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run3.yaml --cores all
```
Output: `output/run3/`  Contains examples of the CDF and the harmonic mean ranks table comparing the weighting vectors `1_1_1_1_1`, `2_2_2_2_1` and `1_1_1_1_0` (Discussion section of the paper).

#### **run4** - *"Effect of network propagation on clinical consistency of top-ranked candidates"*

Full pipeline including patient integration in the network and RWR.
```bash
snakemake -s Snakefile.load_input --configfile configs/run4.yaml --cores all
snakemake -s Snakefile.sim        --configfile configs/run4.yaml --cores all
snakemake -s Snakefile.add_rw     --configfile configs/run4.yaml --cores all
snakemake -s Snakefile.rslt       --configfile configs/run4.yaml --cores all
```

Output: `output/run4/`  Contains examples of tables related to the fourth results (Tables 3 and 4 of the paper on SHEPHERD50; Tables 6 and 7 for Solve-RD and Phenopacket Store).

[⬆ Back to top](#top)


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
- Patient files (example patients, or one cohort from `dataset/`) are placed in `input/patient/`
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

[⬆ Back to top](#top)


## Run Configuration for all steps

The pipeline consists of 6 steps, all configured through a single YAML file located in `configs/`. Each time the pipeline runs, this one file drives all pipeline steps (`Snakefile.sim`, `Snakefile.rslt`, `Snakefile.add_rw`), so there's no need to edit multiple config files between runs.

### Config file structure
Here an example of how to run all pipeline steps with the configuration of the paper (see run4.yaml in section 'Example workflow (subset of 10 diseases)' )
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
                                 #   "gd"  → clinical consistency -based classification (group of disorder) analysis (Tables 3 & 4)
```

Available values: `combine` = `BMA`, `funSimAvg`, `funSimMax`, `rsd` (FunSimMaxAsym); `sm_method` = `resnik`, `lin`, `jc` (Jiang-Conrath), `graphic` (GraphIC), `ic` (SimIC), `relevance`.

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


[⬆ Back to top](#top)



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

### `dataset/`
Contains the patient cohorts derived for the paper from SHEPHERD (SHEPHERD100, SHEPHERD90, SHEPHERD70, SHEPHERD50) and from Phenopacket Store, as Phenopacket JSON files, with one table per cohort giving the confirmed ORPHAcode of each patient. See [Datasets used in the paper](#datasets-used-in-the-paper).

### `input/`
The input directory **must strictly follow the expected architecture**.  
Do **not rename, move, or restructure** folders or files.

#### `input/hpo/`
Contains three HPO files downloaded from the official HPO website (release 2026-09-01 in the paper). Required by the `hpo3` Python library.

#### `input/patient/`
Contains patient files in **Phenopacket format** (one file per patient). The three example patients are already placed here. To run a cohort of the paper, replace them with the files of that cohort from `dataset/` (one cohort at a time). Used for disease-patient similarity and random walk steps.

#### `input/pd_orphanet/`
Contains Orphanet product DDL files downloaded from Orphadata (June 2026 release in the paper). Required to compute disease-related similarity matrices.

#### `[project_name]/` (root)
- Snakemake workflow files (`Snakefile.*`)
- `README.md`
- `path_variable.py` - sets all project paths (do not edit)
- `set_log.py` - configures logging

[⬆ Back to top](#top)


## Step-by-step usage


### Clean all outputs

To delete all generated files and start from scratch, run the following command from the project root:

```bash
snakemake -s Snakefile.clean clean --cores all
```

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


Mandatory parameter in the **config** for this step: `do_subsumed`    


### Patient file format

Patient files must follow the **Phenopacket** format (JSON). The most important fields are `id` and `phenotypicFeatures`. All other sections are optional.
The files of `dataset/` already follow this format and can be used as is.

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

[⬆ Back to top](#top)
	
### Step 2 - Build the Disease-Disease (DD) similarity matrix
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
**Mandatory parameters** in the **config** for this step: `do_subsumed`,`product4`,`dd`,`mini_rd`(optional),`vector_strs`,`sm_method`,`combine`

**Output** (saved in `output/dd_sm/`):
- Individual parquet files per disease into folder [aggregation_method]/[similariy_measure]/[n]/[name_to_define_which_pd4_used]/[weight_vector] (format: `{index}_{ORPHA_code}.parquet`)
- Concatenated file: `{combine}_{method}_{product4}_{vector}.parquet` (all disease-disease similarity scores)

The DD matrix contains no patient: it is built once and reused for every patient cohort.

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

[justif]

[⬆ Back to top](#top)

### Step 3 - Build the Disease-Patient (DP) similarity vectors
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

Mandatory parameters in the **config** for this step: `do_subsumed`,`product4`,`dp`,`mini_rd`(optional),`mini_patient` (optional),`vector_strs`,`sm_method`,`combine`.

**Output** (saved in `output/dp_sm/`):
- Individual parquet files per disease (same folder structure as DD)
- Concatenated file: `{sub0/1}_{combine}_{method}_{product4}_{vector}.parquet` (all patient-disease similarity scores)
- Additional file: `RDI_{combine}_{method}_{product4}_{vector}.xlsx` - for each patient, the disease with the highest similarity score

The best configuration depends on the cohort: Jiang-Conrath × FunSimAvg is the best on SHEPHERD50, whereas Resnik × FunSimMaxAsym (`rsd`) with vector `2_2_2_2_1` performs best on Solve-RD (Supporting Information SF1).

[⬆ Back to top](#top)



### Steps 4 & 5 - Patient integration and Random Walk with Restart
This step integrates patient nodes into the disease similarity matrix, then applies Random Walk with Restart using the NetworkX library.

**Execution:**
```bash
snakemake -s Snakefile.add_rw --configfile configs/run_example.yaml --cores all

# Alternative: run Python scripts directly
python -m bin.main_add_patients_to_dd
python -m bin.main_rarw run -a 0.2 --seeds "P1"
python -m bin.main_rarw collect -a 0.2

# Help
#python -m bin.main_rarw --help
#python -m bin.main_rarw run --help
#python -m bin.main_rarw collect --help
```
Mandatory parameters in the **config** for this step: `do_subsumed`,`alpha`.

In the paper, α between 0.05 and 0.3 improves both the rank of the confirmed diagnosis (optimum α ≈ 0.2, Table 5) and clinical consistency (optimum α = 0.23, Table 3).

**Output:**
- `output/patient_added/` - one parquet file per patient containing the DD matrix extended with that patient. A `PROVENANCE.yaml` file records the configuration used.
- `output/RWR/[alpha_value]/[config_name]/*.parquet` - ranked disease candidates per patient

> **Note:** Only patients present in both the DD and DP matrices are processed.

[⬆ Back to top](#top)

	
### Step 6 - Results

The results reported in the paper are reproducible using the scripts orchestrated by `Snakefile.rslt`, computed from the final outputs of Steps 2-5.

```bash
snakemake -s Snakefile.rslt --configfile configs/run_example.yaml --cores all  

# Alternative: run Python scripts directly
python -u bin/results_CDFs_hm.py
python -u bin/results_GDs.py
```
Mandatory parameters in the **config** for this step: `do_subsumed`,`alpha`.

**Output:**
`results_CDFs_hm.py` aggregates all ranking outputs across patients and methods, extracts the rank of the true diagnosis (RDI), and computes global performance summaries including the harmonic mean of ranks and empirical cumulative distribution functions (CDFs). These correspond to Table 2 and Figures 8 and 9 of the paper.

`results_GDs.py` produces the classification-based comparisons Tables 3 and 4 (SHEPHERD50) and Tables 6 and 7 (Solve-RD and Phenopacket Store).

The comparison with the other tools (Phenomizer, Phrank, RunSolveRD, PubCaseFinder; Table 1 and Figures 5 to 7) is not part of this pipeline: these tools were run separately with their own software, on the same patient files.

[⬆ Back to top](#top)


	
## Computational considerations

Building full disease-disease similarity matrices can be computationally expensive. For development, testing, or exploratory analyses (including reproducibility review), it is strongly recommended to restrict computations using the `mini_rd` parameter.

Example restriction:
```yaml
mini_rd: "ORPHA:100985,ORPHA:100991,ORPHA:1465,ORPHA:329284,ORPHA:34516,ORPHA:412057,ORPHA:663,ORPHA:79445,ORPHA:99949"
```

The Parquet output format is used because it is recommended for large-scale analyses due to improved I/O performance and reduced disk usage.

[⬆ Back to top](#top)


## Troubleshooting

**"HPO files not found"**  
Verify that `input/hpo/` contains the three required HPO files from the official HPO website.

**"No patients found"**  
Ensure Phenopacket files are present in `input/patient/`.

**"YAML parsing error"**  
Check YAML indentation - use spaces, never tabs.

**Memory issues on large datasets**  
Use `mini_rd` and/or `mini_patient` to restrict computation to a smaller subset, and increase available RAM.

[⬆ Back to top](#top)

	
## How to cite

If you use this pipeline in an academic or scientific context, please cite:

> Chahdil, M. et al. *Combining phenotypic similarity and network propagation to improve performance and clinical consistency of rare disease diagnosis.* (Publication details to be added upon journal release.)

For the Dataset :
- Alsentzer E, Li MM, Kobren SN, Noori A, Kohane IS, Zitnik M. Few shot learning for phenotype-driven diagnosis of patients with rare genetic diseases. *npj Digit Med.* 2025;8:380. doi:10.1038/s41746-025-01749-1
- Alsentzer E, Li MM, Kobren SN, Undiagnosed Diseases Network, Kohane IS, Zitnik M. Deep Learning for diagnosing patients with rare genetic diseases. Harvard Dataverse; 2024. doi:10.7910/DVN/TZTPFL
- *Danis D, Bamshad MJ, Bridges Y, et al. A corpus of GA4GH phenopackets: Case-level phenotyping for genomic diagnostics and discovery. *HGG Adv.* 2024;6:100371. doi:10.1016/j.xhgg.2024.100371
[ SOLVE RD]
	

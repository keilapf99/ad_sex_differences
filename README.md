# Sex-specific transcriptional regulatory networks in Alzheimer's disease

Keila A. Pérez-Fuentes, Alejandra P. Pérez González, Joel Ruíz-Hernández, Gabriela, Enrique Lemus, Hugo Tovar

Analysis code for the manuscript. Bulk RNA-seq from two AMP-AD cohorts (ROSMAP: DLPFC, PCC, HCN; Mayo: CRB, TC) is stratified by sex and diagnosis (CERAD 1 = AD, CERAD 4 = control). We compare differential expression, co-expression network topology, gene connectivity and master regulators between females and males.

## Repository layout

Each analysis module matches a section of the Methods. A module holds its code and a README covering who owns it, what it reads, what it writes and how to run it.

| Module | Methods | Owner | Status |
|---|---|---|---|
| [`01_sample_stratification`](analysis/01_sample_stratification) | 3.1–3.2 Datasets and stratification (Table 1) | Paulina (A. P. Pérez González), Keila A. Pérez-Fuentes | results present, code pending |
| [`02_preprocessing_qc`](analysis/02_preprocessing_qc) | 3.3 PCA, NOISeq QC, CPM filter, EDASeq GC correction, TMM, ComBat-seq | Paulina (A. P. Pérez González) | code pending |
| [`03_differential_expression`](analysis/03_differential_expression) | 3.4 DESeq2 (by sex, and sex × diagnosis) | Paulina (A. P. Pérez González) | reference scripts |
| [`04_network_inference`](analysis/04_network_inference) | 3.5 ARACNe-AP co-expression networks | Hugo Tovar | code pending |
| [`05_network_topology`](analysis/05_network_topology) | 3.6 Network topology, AD vs control | Paulina (A. P. Pérez González) | code pending |
| [`06_connectivity`](analysis/06_connectivity) | 3.7–3.8 Degree correlation and sex-specific connectivity changes | Joel Ruíz-Hernández | code pending |
| [`07_master_regulators`](analysis/07_master_regulators) | 3.9 Master regulator analysis (VIPER/msviper) | Keila A. Pérez-Fuentes | done |

```
.
├── analysis/<module>/        code + README per Methods section
│   └── scripts/              numbered in execution order (00_, 01_, …)
├── config/                   paths.R / paths.sh: the only place paths are defined
├── data/                     NOT in git (controlled-access data); see data/README.md
├── results/<module>/         NOT in git; each module writes only to its own folder
├── environment.yml           conda environment
└── CONTRIBUTING.md           how to add your part
```

## Quick start

```bash
git clone git@github.com:keilapf99/ad_sex_differences.git
cd ad_sex_differences
conda env create -f environment.yml && conda activate ad-sex

# Point at the data. Either copy/symlink it into ./data, or use the shared copy:
export ADSEX_DATA=/path/to/shared/data       # optional
export ADSEX_RESULTS=/path/to/shared/results # optional

# Always run scripts from the repository root
Rscript analysis/07_master_regulators/scripts/03_MRA_by_sex.R
```

Scripts never hardcode paths. They `source("config/paths.R")` and build every path from `DATA_DIR`, `RESULTS_DIR` or `module_results("<module>")`. These default to `./data` and `./results`, and the `ADSEX_DATA`, `ADSEX_RESULTS` and `ADSEX_ROOT` environment variables override them.

## Data access

The raw data come from the [AMP-AD Knowledge Portal](https://www.synapse.org/): ROSMAP (syn8456629) and Mayo RNA-seq (syn8466812). They are covered by data use agreements and **must never be committed or pushed**. `data/` and `results/` are git-ignored for that reason. See [`data/README.md`](data/README.md) for the expected contents.

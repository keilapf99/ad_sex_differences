# 07_master_regulators

**Methods section:** 3.9 Master regulator analysis
**Owner:** Keila A. Pérez-Fuentes

Master regulator analysis (VIPER/msviper) of AD vs control, run separately in females and males for each of the 5 regions. The analysis then compares the master regulators between sexes.

## Inputs
- `data/coexpr_networks/` (from `04_network_inference`) and `data/regulators.txt`
- `data/counts_tmm/expression_alzheimer/` (from `02_preprocessing_qc`)
- `data/full_diff_expr/` (from `03_differential_expression`), using the Wald `stat` column as the signature
- `data/raw_counts/` (to build the permutation null model)
- `data/gene_id_to_symbol.tsv`

## Pipeline (run from the repository root, in order)

| Script | What it does | Output (`results/07_master_regulators/`) |
|---|---|---|
| `00_trn_filter.sh` | Keep only TF → target edges from the ARACNe networks | `data/trn_networks/` |
| `01_preprocess_{mayo,rosmap}.R` | Format counts + metadata the same way as the DESeq2 base scripts | `0_preprocess/`, `0_preprocess_rosmap/` |
| `02_null_model_dea.R` | 1000 label permutations × DESeq2, keeping the Wald stat (seed 10) | `null_model_dea/` (~9 GB) |
| `03_MRA_by_sex.R` | msviper per cohort × region × sex | `MRA_by_sex/*_MRA_results.txt` |
| `04_nes_quartile_filter.R` | Keep the top quartile of \|NES\| per combination | `MRA_by_sex/top_quartile_NES/` |
| `05_sex_specific_and_concordance.R` | Sex-specific MRs; for shared MRs, top quartile + concordant/discordant | `MRA_by_sex/sex_specific_MRs/`, `MRA_by_sex/top_quartile_NES_shared/` |
| `06_*.R` | Figures (volcano, discordant heatmap, overview, three-panel summary) | `MRA_by_sex/*.pdf/png` |

`02` and `03` read `NULL_MODEL_NPERM`, `NULL_MODEL_CORES` and `MRA_CORES` from the environment.

`legacy/mra_signature.R` is the original single-combination prototype. It is kept only for reference and still contains absolute paths.

## Summary of current results
- 873 sex-specific master regulators (all regions)
- 1002 shared top-quartile master regulators: 940 concordant, 62 discordant

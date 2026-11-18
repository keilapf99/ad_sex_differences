# data/

Git-ignored. Only this README is tracked. AMP-AD data are controlled-access: request them through Synapse and never commit them.

| Path | Contents | Produced by |
|---|---|---|
| `raw_counts/ROSMAP_counts/` | ROSMAP raw counts + metadata (DLPFC, PCC, HCN) | Synapse syn8456629 |
| `raw_counts/MAYO_rawcounts/` | Mayo raw counts + metadata (CRB, TC) | Synapse syn8466812 |
| `counts_tmm/expression_alzheimer/` | QC'd, GC-corrected, TMM-normalized, batch-corrected matrices per sex/condition/region | `02_preprocessing_qc` |
| `diff_expression/` | DEGs AD vs control per sex/region (filtered) | `03_differential_expression` |
| `full_diff_expr/` | Full DESeq2 tables AD vs control per sex/region (signature for MRA) | `03_differential_expression` |
| `coexpr_networks/` | ARACNe-AP consensus networks per sex/condition/region | `04_network_inference` |
| `trn_networks/` | Co-expression networks filtered to TF → target | `07_master_regulators/00_trn_filter.sh` |
| `regulators.txt` | List of transcription factors | external reference |
| `gene_id_to_symbol.tsv` | Ensembl ID → gene symbol | external reference |

As each module's code is added to the repo, its outputs should move from `data/` to `results/<module>/`. Only raw inputs and external references should stay in `data/`.

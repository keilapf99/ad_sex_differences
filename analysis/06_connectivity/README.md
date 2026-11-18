# 06_connectivity

**Methods section:** 3.7–3.8 Gene connectivity and sex-specific connectivity changes
**Owner:** Joel Ruíz-Hernández

## Inputs
Networks from `04_network_inference`

## Outputs
Degree tables, Spearman correlations AD vs control, female-/male-specific/shared connectivity-change genes

## Method / parameters
Degree per gene. Spearman AD vs control within region × sex. Δdegree standardized as a robust z-score; |z| ≥ 1.96 = outlier. Genes that are outliers in one sex only are sex-specific; outliers in both sexes with opposite direction are "shared"; same direction are excluded.

## How to run
_TBD: add commands here, run from the repository root._

# 03_differential_expression

**Methods section:** 3.4 Differential expression analysis
**Owner:** Paulina (Alejandra P. Pérez González)

## Inputs
Raw counts + metadata

## Outputs
DESeq2 tables (currently in `data/diff_expression/` and `data/full_diff_expr/`)

## Method / parameters
(1) AD vs control within each sex, `~ sequencingBatch + pmi + RIN + condition`. (2) All samples with sex × diagnosis interaction. DEG: padj < 0.05, |log2FC| > 0.5.

`scripts/expr_diff_base_{mayo,rosmap}.R` are the reference DESeq2 scripts. Paulina should confirm they are the final version and update their paths to use `config/paths.R`. They currently take paths as command-line arguments.

## How to run
_TBD: add commands here, run from the repository root._

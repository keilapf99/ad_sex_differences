# 02_preprocessing_qc

**Methods section:** 3.3 RNA-seq preprocessing and quality control
**Owner:** Paulina (Alejandra P. Pérez González)

## Inputs
Raw counts per sex/condition/region

## Outputs
Normalized, batch-corrected matrices (currently in `data/counts_tmm/`), PCA and NOISeq QC plots

## Method / parameters
PCA, NOISeq bias diagnostics (composition, GC, length), CPM > 1 condition-aware filter, EDASeq GC correction, TMM (NOISeq), ComBat-seq (sva) for sequencing batch, PCA before/after.

## How to run
_TBD: add commands here, run from the repository root._

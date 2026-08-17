#!/bin/bash

####################################################################
# Coexpression networks filter to transcriptional regulatory network
# By: Keila A. Perez-Fuentes
# Date: 14-08-2026
# Objective: The objective of this script is to filter the complete
# coexpression network inferenced by aracne into a trn with just
# tfs as regulator for downstream analysis

# Input:
# - Coexpression networks
# - List of regulators (transcription factors)

# Output:
# - Transcriptional regulatory networks

####################################################################
REGULATORS="$HOME/ad_sex_differences/data/regulators.txt"
INPUT_DIR="$HOME/ad_sex_differences/data/coexpr_networks"
OUTPUT_DIR="$HOME/ad_sex_differences/data/trn_networks"

# Regulator column
REG_COL=1

####################################################################

mkdir -p "$OUTPUT_DIR"

for f in "$INPUT_DIR"/*; do

    base=$(basename "$f")

    # Extract phenotype and brain region
    # ROSMAP_RNAseq_female_AD_counts_CRB.network.txt
    #              ↓
    # female_AD_CRB_trn.txt

    name=$(echo "$base" | sed -E 's/^ROSMAP_RNAseq_//; s/_counts//; s/\.network.*$/_trn.txt/')

    echo "Procesando $base -> $name ..."

    (
        head -n 1 "$f"
        awk -F'\t' -v col="$REG_COL" '
            NR==FNR {a[$1]; next}
            ($col in a)
        ' "$REGULATORS" "$f"
    ) > "$OUTPUT_DIR/$name"

done

echo "Filtrado terminado"

######################## END #######################################

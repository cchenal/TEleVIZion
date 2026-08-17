#!/bin/bash
# ------------------------------------------------------------------- [bsub] ---
#BSUB -o analyses/farm_logs/televizion-%J.o
#BSUB -e analyses/farm_logs/televizion-%J.e
#BSUB -J televizion
#BSUB -q long
#BSUB -n 1
#BSUB -M 256000
#BSUB -R "select[mem>256000] rusage[mem=256000] span[hosts=1]"

# ----------------------------------------------------------------- [module] ---

set -euo pipefail

# ------------------------------------------------------------------ [conda] ---

conda activate televizion

# ------------------------------------------------------------------- [code] ---

project="/lustre/scratch126/tol/teams/lawniczak/users/cc54/projects/TEleVIZion/"

# name="qqMarElec17"
# name="qqMarMung10"
name="qqMarSpei1"

genome=${project}"data/peacock_spiders/"${name}"/"${name}".hap1.chroms.tsv"
repeatmasker=${project}"data/peacock_spiders/"${name}"/"${name}".hap1.RM.out"

windowsize=500000
layout="vertical"
figsize="14,10"


if [[ ! -r "$genome" ]]; then
    echo "Error: genome file not found: $genome" >&2
    exit 1
fi

if [[ ! -r "$repeatmasker" ]]; then
    echo "Error: RepeatMasker file not found: $repeatmasker" >&2
    exit 1
fi

echo "name         = $name"
echo "genome       = $genome"
echo "repeatmasker = $repeatmasker"
echo "windowsize   = $windowsize"
echo "layout       = $layout"
echo "figsize      = $figsize"



python3 scripts/televizion_cli.py \
    --name "$name" \
    --genome "$genome" \
    --repeatmasker "$repeatmasker" \
    --windowsize "$windowsize" \
    --layout "$layout" \
    --figsize "$figsize"

conda deactivate
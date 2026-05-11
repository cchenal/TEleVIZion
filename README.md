# TEleVIZion

TEleVIZion is a toolkit for generating insightful visualisations to explore and understand the distribution and dynamics of transposable elements (TEs) and repetitive sequences across genomic data.

# Table of contents

# Introduction 

## Installation

TEleVIZion uses both Python packages and R/Bioconductor packages. The recommended setup is the Conda environment provided in `environment.yml`.

```bash
conda env create -f environment.yml
conda activate televizion
```

If you use Mamba or Micromamba, the same file works:

```bash
mamba env create -f environment.yml
mamba activate televizion
```

Check the install with:

```bash
python3 scripts/televizion_cli.py --help
Rscript scripts/televizion/plot_chromosomes_kimura_optimized.R --help
```

RepeatMasker and EDTA are upstream annotation tools. They are not required inside this environment unless you want to generate those annotation files yourself; TEleVIZion only needs their output files as inputs.

# Pipeline overview

# Input files 

## Genome 

## GC content

## Accessibility 

## Repeats annotation

### RepeatMasker

#### Kimura distance 

```
<Command line to generate kimura file>
```

### EDTA

# Output files 

## Whole genome statistics

## Chromosome specific statistics

## Whole genome repeat landascape 

## Order specific repeat landscape

# Parameters and options

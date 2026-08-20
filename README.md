# TEleVIZion

TEleVIZion is a research workflow for visualising transposable element (TE)
annotations across genomes. It converts RepeatMasker or EDTA output into
windowed genome summaries, whole-genome bar plots, and chromosome-scale
karyotype plots that show where repeat classes, repeat families, divergence
bins, identity bins, GC content, and optional gene-content tracks occur.

The project is designed for exploratory comparative genomics: start from an
annotated genome, decide which chromosomes or scaffolds you want to inspect,
choose a window size, and generate a compact set of figures and tabular
intermediate files under `analyses/<run_name>/`.

## Contents

- [What TEleVIZion does](#what-televizion-does)
- [How the workflow works](#how-the-workflow-works)
- [Installation](#installation)
- [Repository layout](#repository-layout)
- [Input files](#input-files)
- [Quick start](#quick-start)
- [Command-line use cases](#command-line-use-cases)
- [Output files](#output-files)
- [Command-line options](#command-line-options)
- [Troubleshooting](#troubleshooting)
- [License](#license)

## What TEleVIZion Does

TEleVIZion answers questions such as:

- Which repeat classes dominate a genome by insertion count or base-pair span?
- Which repeat families make up each class?
- Are repeats concentrated on particular chromosomes or chromosome arms?
- Are some regions enriched for young or old repeats?
- Do repeat-rich windows coincide with GC-rich/GC-poor regions or an external
  gene-content track?
- How does an EDTA annotation compare with a RepeatMasker annotation when both
  are summarised into the same windowed genome view?

TEleVIZion currently supports two annotation sources:

- `RepeatMasker`: standard `.out` annotation files, optionally paired with a
  RepeatMasker Kimura divergence summary.
- `EDTA`: GFF3 files produced by EDTA, using `classification`, `Name`, and
  `identity` attributes in column 9.

It does not run RepeatMasker or EDTA for you. Those tools are upstream
annotation steps; TEleVIZion consumes their output.

## How The Workflow Works

The main entry point is:

```bash
python3 scripts/televizion_cli.py --help
```

At a high level, the CLI performs these steps:

1. Builds genome windows from a genome metadata table.
   Each chromosome/scaffold is split into fixed-size non-overlapping windows.

2. Parses repeat annotations.
   RepeatMasker records are read from `.out` files. EDTA records are read from
   GFF3 files. For each annotation, TEleVIZion extracts chromosome, start, end,
   repeat class, repeat family, repeat element name, and divergence or identity
   information when available.

3. Handles overlapping annotations inside each window.
   The parser splits overlapping intervals with an interval tree so that
   base-pair coverage is not naively double-counted. If overlapping segments
   have the same class and family, they are counted together. If they represent
   conflicting repeat types, the segment is assigned to `Ambiguous`.

4. Aggregates repeat signal.
   For every window and repeat class/family, TEleVIZion records:

   - insertion count
   - base-pair span
   - fraction of the window covered by that class

5. Adds divergence or identity summaries.
   RepeatMasker Kimura values are grouped into bins:
   `0-10`, `10-20`, `20-30`, `30-40`, and `40-70`.

   EDTA or RepeatMasker identity values are grouped into bins:
   `1-0.9`, `0.9-0.8`, `0.8-0.7`, `0.7-0.6`, and `0.6-0`.

6. Generates plots.
   Python and Matplotlib produce whole-genome and optional per-chromosome bar
   plots. R and `karyoploteR` produce chromosome-scale karyotype plots.

7. Writes reusable intermediate tables.
   The karyotype input tables are saved in
   `analyses/<run_name>/karyoplot_tables/`.

## Installation

TEleVIZion uses Python packages plus R/Bioconductor packages. The recommended
setup is the Conda environment provided in `environment.yml`.

### Option 1: Recommended — Conda

```bash
conda env create -f environment.yml
conda activate televizion

# Mamba or Micromamba can use the same environment file:
# mamba env create -f environment.yml
# mamba activate televizion
```

### Option 2: Without Conda

You can also install TEleVIZion without Conda, but this requires managing the Python and R dependencies separately.

Create a Python virtual environment and install the Python packages:

```bash
python3 -m venv .venv_televizion
source .venv_televizion/bin/activate
python -m pip install --upgrade pip
python -m pip install numpy pandas matplotlib seaborn intervaltree biopython
```

Then install the R dependencies in your R environment:

* `r-base`
* `optparse`
* `karyoploteR`
* `GenomicRanges`
* `IRanges`

You will also need a working R installation (`r-base`) on your system before running these commands.

```r
install.packages("optparse")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "karyoploteR",
    "GenomicRanges",
    "IRanges"
))
```

### Check installation

Check that both the Python CLI and the R plotting script are available:

```bash
python3 scripts/televizion_cli.py --help
Rscript scripts/televizion/plotting_landscape.R --help
```

Run commands from the repository root. Several scripts currently use paths such
as `scripts/utils/create_gc_content.py` and
`scripts/televizion/plotting_landscape.R` relative to the current
working directory.

## Repository Layout

```text
.
|-- README.md
|-- environment.yml
|-- scripts/
|   |-- televizion_cli.py
|   |-- utils/
|   |   |-- create_chroms.py
|   |   |-- create_gc_content.py
|   |   `-- create_gene_content.py
|   `-- televizion/
|       |-- aggregation.py
|       |-- io.py
|       |-- plotting_stats.py
|       `-- plotting_landscape.R
|-- data/
|   `-- example input files, if present locally
`-- analyses/
    `-- generated output directories
```

The main code paths are:

- `scripts/televizion_cli.py`: command-line orchestration.
- `scripts/televizion/io.py`: RepeatMasker, Kimura, and EDTA input parsing.
- `scripts/televizion/aggregation.py`: genome windows, overlap splitting, and
  table export.
- `scripts/televizion/plotting_stats.py`: Matplotlib bar plots and calls into R.
- `scripts/televizion/plotting_landscape.R`: karyotype plots.
- `scripts/utils/create_chroms.py`: helper to make `chroms.tsv` from FASTA.
- `scripts/utils/create_gc_content.py`: helper to make GC-content windows from FASTA.
- `scripts/utils/create_gene_content.py`: helper to make gene-content windows from GFF.

## Input Files

### Genome Metadata

The genome metadata file defines the chromosomes or scaffolds to plot. It is a
tab-delimited file with this shape:

```text
chr	start	end	name	gieStain
OX030893.1	1	28269272	X	chalk
OX030891.1	1	126027569	2RL	chalk
OX030892.1	1	95248607	3RL	chalk
```

Columns:

- `chr`: chromosome or scaffold ID used in the annotation file.
- `start`: usually `1`.
- `end`: sequence length in bp.
- `name`: display label used on the karyotype plot.
- `gieStain`: cytoband-like colour label consumed by `karyoploteR`; `chalk` or
  `gneg` are common simple values.

The `chr` values must match the sequence IDs in your RepeatMasker `.out` or EDTA
GFF3 file. The `name` values can be prettier labels, for example using
`OX030891.1` as the true sequence ID while displaying `2RL` on the plot.

You can generate a starter metadata file from a FASTA file:

```bash
python3 scripts/utils/create_chroms.py \
    -o data/AcolN3/AcolN3.chroms.tsv \
    -f data/AcolN3/AcolN3.fasta \
    --min-length 10000000
# Wrote: data/AcolN3/AcolN3.chroms.tsv
# Sequences kept: 3, skipped (< 10000000 bp): 126
```

Only use the FASTA-derived file directly if the FASTA IDs match the annotation
IDs. If your annotation uses renamed chromosomes, edit the `chr` column or build
the metadata from an assembly report.

### RepeatMasker Annotation

Use a standard RepeatMasker `.out` file:

```text
SW   perc perc perc  query       position in query ... repeat class/family ...
18   30.6 1.4  1.4   OX030925.1  7305 7378 ... (ATTTT)n Simple_repeat ...
```

TEleVIZion reads:

- query sequence as chromosome/scaffold
- query begin/end as coordinates
- repeat name as the repeat element
- class/family as repeat class and repeat family
- percent divergence as a fallback identity-like value in the parser

For current RepeatMasker runs, pass a Kimura summary with `--kimura` so the
karyotype can include the K2P/Kimura panel.

### RepeatMasker Kimura Summary

The Kimura file is the RepeatMasker divergence summary with columns like:

```text
Class	Repeat	absLen	wellCharLen	Kimura%
DNA/	Acol_otherMITEs_Eles1	5705	5623	11.70
DNA/CACTA	Acol_CACTA_Ele1	58882	57414	39.95
```

One common way to generate it from a RepeatMasker alignment file is:

```bash
calcDivergenceFromAlign.pl -s data/my_genome/repeats.divsum \
  data/my_genome/repeats.align
```

Depending on your RepeatMasker installation, the utility may live under the
RepeatMasker `util/` directory and may need to be called with `perl`.

### EDTA Annotation

Use an EDTA GFF3 annotation file:

```text
2RL	EDTA	LTR_retrotransposon	6037	6194	698	-	.	ID=TE_homo_4;Name=TE_00001204_LTR;classification=LTR/unknown;identity=0.784;method=homology
```

TEleVIZion reads these column 9 attributes:

- `classification`: split into class/family, such as `LTR/unknown`.
- `Name`: repeat element name.
- `identity`: match identity, used for identity-bin karyotype panels.

Records with `identity=NA` are still counted for class/family abundance, but do
not contribute to identity-bin summaries.

### GC Content Track

Create the GC-content track before running TEleVIZion, then pass its path with `--gc`:

```bash
python3 scripts/utils/create_gc_content.py data/my_genome.fasta \
  --window 10000 \
  --out data/my_genome/my_genome_gc_windows_10000.tsv
```

`--gc` and `--gene-content` are mutually exclusive; choose one contextual track per run.
The GC table is tab-delimited:

```text
chr	start	end	name	itemRgb
chr1	0	10000	gc1-0.42	#66c2a5
```

### Gene Content Track

Create the gene-content track with `scripts/utils/create_gene_content.py`, then pass its path with `--gene-content`. The file is a five-column, tab-delimited coloured track:

```text
chr	start	end	name	itemRgb
2RL	0	100000	genes1-12	#6aaed6
2RL	100000	200000	genes2-7	#bad6eb
```

The `chr` values must match the genome metadata. Gene content is displayed with a Blues colour bar from 0 to the global maximum genes per window.

### Palette File

By default, TEleVIZion assigns repeat-class colours from Matplotlib's `turbo`
colour map. You can override colours with a TSV palette:

```text
blue	65	133	190	#4185be	LTR
navy	30	58	95	#1e3a5f	DNA
light_blue	60	200	203	#3cc8cb	LINE
```

The bundled `data/palettes/personalised_palette_example.tsv` follows this
six-column format without a header. The last column can contain comma-separated
class aliases, such as `unknown,Undetermined,Unknown`. If you use `--palette`,
every class in `--classesorder` must be present in the palette's category column.

## Quick Start

The following command uses the bundled `Acol_lib_GCA_943734845.1` RepeatMasker
example paths, if they are present in your local checkout:

```bash
python3 scripts/televizion_cli.py \
  --name MyGenome \
  --genome data/Acol_lib_GCA_943734845.1/chroms.tsv \
  --gc data/Acol_lib_GCA_943734845.1/GCA_943734845.1_gc_windows_10000.tsv \
  --repeatmasker data/Acol_lib_GCA_943734845.1/GCA_943734845.1.out \
  --kimura data/Acol_lib_GCA_943734845.1/GCA_943734845.1.kimura \
  --windowsize 5000000 \
  --chromtoplot OX030923.1,OX030924.1,OX030925.1
```

```
bsub -M50000 -R"select[mem>50000] rusage[mem=50000] span[hosts=1]" -n 1 -q small -Is bash -lc "conda activate televizion && python scripts/televizion_cli.py --name AcolN3 --genome data/AcolN3/AcolN3.chroms_names.tsv --gc data/AcolN3/AcolN3_gc_windows_10000.tsv --repeatmasker data/AcolN3/AcolN3.RM.out --windowsize 500000 --layout horizontal --chromtoplot OX030891.1,OX030892.1,OX030893.1 --figsize 11,8"
```

With the default `--output-formats pdf`, this creates:

```text
analyses/MyGenome/
|-- MyGenome_whole_genome_stacked_counts_by_class.pdf
|-- MyGenome_whole_genome_stacked_lengths_by_class.pdf
|-- MyGenome_whole_genome_contiguous_counts_by_class.pdf
|-- MyGenome_whole_genome_contiguous_lengths_by_class.pdf
|-- MyGenome_5000000_karyoplot_stacked_counts_by_class.pdf
|-- MyGenome_5000000_karyoplot_stacked_percentage_by_class.pdf
`-- karyoplot_tables/
    |-- MyGenome_5000000_repeat_classes.bed
    `-- MyGenome_5000000_kimura.bed
```

## Command-Line Use Cases

```
bsub -M150000 -R"select[mem>150000] rusage[mem=150000] span[hosts=1]" -n 1 -q small -Is conda activate televizion; pwd; python3 scripts/televizion_cli.py --name AcolN3 --genome data/AcolN3/AcolN3.chroms_names.tsv --gc data/AcolN3/AcolN3_gc_windows_10000.tsv --repeatmasker data/AcolN3/AcolN3.RM.out --kimura data/AcolN3/AcolN3.RM.kimura --windowsize 5000000 --layout vertical; conda deactivate
```

### 1. RepeatMasker With Kimura Divergence

Use this when you have RepeatMasker `.out` annotations and a Kimura divergence
summary:

```bash
python3 scripts/televizion_cli.py \
  --name An_funestus_RM \
  --genome data/Acol_lib_GCA_943734845.1/chroms.tsv \
  --repeatmasker data/Acol_lib_GCA_943734845.1/GCA_943734845.1.out \
  --kimura data/Acol_lib_GCA_943734845.1/GCA_943734845.1.kimura \
  --windowsize 500000 \
  --chromtoplot OX030923.1,OX030924.1,OX030925.1
```

This produces whole-genome bar plots, a class/family karyotype, and a Kimura
K2P panel in the karyotype figures.

### 2. EDTA Annotation With Identity Tracks

Use this when you have an EDTA GFF3 file and matching genome metadata:

```bash
python3 scripts/televizion_cli.py \
  --name MyGenome_EDTA \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --windowsize 500000 \
  --chromtoplot 2RL,3RL,X
```

EDTA runs create an identity-bin table:

```text
analyses/MyGenome_EDTA/karyoplot_tables/MyGenome_EDTA_500000_identity.bed
```

The karyotype plots then include an identity panel instead of a Kimura panel.

### 3. Plot All Chromosomes

`--chromtoplot all` is the default. It plots every sequence listed in the genome
metadata file:

```bash
python3 scripts/televizion_cli.py \
  --name WholeAssembly \
  --genome data/my_genome/chroms.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --kimura data/my_genome/repeats.divsum \
  --windowsize 1000000 \
  --chromtoplot all
```

For assemblies with many scaffolds, consider filtering the genome metadata file
or using `scripts/utils/create_chroms.py --min-length` before plotting.

### 4. Generate Per-Chromosome Summary Bar Plots

Add `--perchromosome` to produce the four summary bar plots for each chromosome
or scaffold:

```bash
python3 scripts/televizion_cli.py \
  --name PerChromosome \
  --genome data/my_genome/chroms.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --kimura data/my_genome/repeats.divsum \
  --windowsize 500000 \
  --chromtoplot chr2,chr3,chrX \
  --perchromosome
```

Per-chromosome plots are named with the `name` column from the genome metadata,
for example with the default PDF format:

```text
analyses/PerChromosome/per_chromosome/PerChromosome_chr2_stacked_counts_by_class.pdf
```

### 5. Generate Per-Class Karyotype Figures

Add `--perclass` to create one karyotype figure per repeat class in addition to the
stacked all-class plots:

```bash
python3 scripts/televizion_cli.py \
  --name PerClass \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --windowsize 250000 \
  --chromtoplot 2RL,3RL,X \
  --perclass
```

Example per-class outputs with the default PDF format:

```text
analyses/PerClass/per_class/PerClass_250000_karyoplot_percentage_LTR.pdf
analyses/PerClass/per_class/PerClass_250000_karyoplot_counts_LTR.pdf
```

### 6. Include GC Content

Pass a precomputed GC-content track with `--gc`:

```bash
python3 scripts/televizion_cli.py \
  --name WithGC \
  --genome data/my_genome/chroms.tsv \
  --gc data/my_genome/my_genome_gc_windows_10000.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --kimura data/my_genome/repeats.divsum \
  --windowsize 500000 \
  --chromtoplot 2RL,3RL,X
```

The GC track is drawn on the ideogram and a GC colour bar is added to the
karyotype legend.

### 7. Include Gene Content

Pass a precomputed gene-content track with `--gene-content`:

```bash
python3 scripts/televizion_cli.py \
  --name WithGeneContent \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --gene-content data/my_genome/gene_content_windows.tsv \
  --windowsize 500000 \
  --chromtoplot 2RL,3RL,X
```

Generate this track with `scripts/utils/create_gene_content.py`; its chromosome
IDs must match the genome metadata.

### 8. Change Class Order

By default, repeat classes are ordered by total base-pair span. Override the
order with `--classesorder`:

```bash
python3 scripts/televizion_cli.py \
  --name OrderedClasses \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --windowsize 500000 \
  --chromtoplot 2RL,3RL,X \
  --classesorder unknown,DNA,LTR,LINE,SINE,MITE
```

Only classes in `--classesorder` are carried through to the ordered plots and
karyotype tables, so include every class you want to display.

### 9. Use A Custom Palette

Combine `--classesorder` with `--palette`:

```bash
python3 scripts/televizion_cli.py \
  --name CustomPalette \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --windowsize 500000 \
  --chromtoplot 2RL,3RL,X \
  --classesorder unknown,DNA,LTR,LINE,SINE,MITE \
  --palette scripts/palette.tsv
```

Make sure the palette includes each listed class. If you add `Ambiguous` to the
class order, add an `Ambiguous` category to the palette as well.

### 10. Change Figure Size For General Statistics Plots

The `--figsize` option controls only the Python/Matplotlib general statistics
bar plots:

- whole-genome stacked counts and lengths
- whole-genome contiguous counts and lengths
- per-chromosome bar plots generated with `--perchromosome`

It does not affect any R/karyoploteR karyotype figures, including the
`*_karyoplot_*.<format>` outputs.

```bash
python3 scripts/televizion_cli.py \
  --name WideBars \
  --genome data/my_genome/chroms.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --kimura data/my_genome/repeats.divsum \
  --windowsize 500000 \
  --figsize 14,8
```

Use `W,H` in inches.

### 11. Choose Output Formats And DPI

By default, TEleVIZion writes figure outputs as PDF. Use `--output-formats` to
request one or more formats among `pdf`, `png`, and `jpg`:

```bash
python3 scripts/televizion_cli.py \
  --name MultiFormat \
  --genome data/my_genome/chroms.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --kimura data/my_genome/repeats.divsum \
  --windowsize 500000 \
  --output-formats pdf,png,jpg \
  --dpi 300
```

The `--dpi` option controls PNG/JPG raster resolution and defaults to `300`.
Values below `300` are rejected. PDF output remains vector-based where possible.

### 12. Vertical Karyotype Layout

Use `--layout vertical` to stack chromosomes vertically:

```bash
python3 scripts/televizion_cli.py \
  --name VerticalView \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --windowsize 500000 \
  --chromtoplot 2RL,3RL,X \
  --layout vertical
```

The default is `horizontal`.

### 13. Zoom Into A Region

Use `--zoom` with `chrom:start-end`:

```bash
python3 scripts/televizion_cli.py \
  --name ZoomedRegion \
  --genome data/my_genome/chroms.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --kimura data/my_genome/repeats.divsum \
  --windowsize 50000 \
  --chromtoplot 2RL,3RL,X \
  --zoom 2RL:10000000-15000000
```

Zoom coordinates must be whole-number coordinates inside the chromosome bounds
defined in the genome metadata. When zooming, choose a window size small enough
that the region contains plotted windows.

## Output Files

All outputs are written below:

```text
analyses/<name>/
```

### Output Formats

Figure outputs are PDF by default. If `--output-formats` contains multiple
formats, each figure is written once per requested extension. For example,
`--output-formats pdf,png,jpg` creates matching `.pdf`, `.png`, and `.jpg`
files. PNG and JPG are written at `--dpi` resolution, with `300` dpi as the
default minimum.

### Whole-Genome Summary Plots

These figures are always produced:

```text
<name>_whole_genome_stacked_counts_by_class.<format>
<name>_whole_genome_stacked_lengths_by_class.<format>
<name>_whole_genome_contiguous_counts_by_class.<format>
<name>_whole_genome_contiguous_lengths_by_class.<format>
```

Meaning:

- `stacked_counts`: insertion counts stacked by family within each class.
- `stacked_lengths`: total base-pair span stacked by family within each class.
- `contiguous_counts`: horizontal family-level count ranking.
- `contiguous_lengths`: horizontal family-level length ranking.

### Per-Chromosome Summary Plots

Produced only with `--perchromosome` and written under `per_chromosome/`:

```text
per_chromosome/<name>_<chromosome_display_name>_stacked_counts_by_class.<format>
per_chromosome/<name>_<chromosome_display_name>_stacked_lengths_by_class.<format>
per_chromosome/<name>_<chromosome_display_name>_contiguous_counts_by_class.<format>
per_chromosome/<name>_<chromosome_display_name>_contiguous_lengths_by_class.<format>
```

### Karyotype Plots

These figures are always produced:

```text
<name>_<windowsize>_karyoplot_stacked_percentage_by_class.<format>
<name>_<windowsize>_karyoplot_stacked_counts_by_class.<format>
```

With `--perclass`, additional figures are written under `per_class/`:

```text
per_class/<name>_<windowsize>_karyoplot_percentage_<class>.<format>
per_class/<name>_<windowsize>_karyoplot_counts_<class>.<format>
```

With `--zoom`, the output filename includes a zoom suffix:

```text
<name>_<windowsize>_zoom_<chrom>_<start>_<end>_karyoplot_stacked_percentage_by_class.<format>
```

### Karyoplot Tables

Intermediate tables are written to:

```text
analyses/<name>/karyoplot_tables/
```

Class abundance table:

```text
<name>_<windowsize>_repeat_classes.bed
```

This table contains one row per window with:

- `chrom`, `start`, `end`, `barycenter`
- `<class>_count`
- `<class>_count_stacked`
- `<class>_length`
- `<class>_pct`
- `<class>_pct_stacked`

RepeatMasker Kimura table, when `--kimura` is supplied:

```text
<name>_<windowsize>_kimura.bed
```

EDTA identity table:

```text
<name>_<windowsize>_identity.bed
```

The tables are BED-like TSV files intended for plotting. The exported
window `start` is 0-based, while `end` is the window end coordinate.

## Command-Line Options

```text
--name NAME
    Output prefix. Outputs go to analyses/NAME/.

--genome GENOME
    Required genome metadata TSV.

--gc GC
    Optional precomputed GC-content track from scripts/utils/create_gc_content.py.
    Mutually exclusive with --gene-content.

--repeatmasker REPEATMASKER
    RepeatMasker .out input. Mutually exclusive with --edta. One of
    --repeatmasker or --edta is required.

--edta EDTA
    EDTA GFF3 input. Mutually exclusive with --repeatmasker. One of
    --repeatmasker or --edta is required.

--kimura KIMURA
    RepeatMasker Kimura divergence summary. Use with --repeatmasker.

--windowsize WINDOWSIZE
    Window size in bp. Default: 10000.

--chromtoplot CHROMTOPLOT
    Comma-separated chromosome/scaffold IDs, or all. Default: all.

--perchromosome
    Generate per-chromosome bar plots.

--classesorder CLASSESORDER
    Comma-separated repeat class order override.

--perclass
    Generate per-class karyotype plots.

--gene-content GENE_CONTENT
    Optional precomputed gene-content track from scripts/utils/create_gene_content.py.
    Mutually exclusive with --gc. The legend spans 0 to the global maximum.

--figsize W,H
    Python/Matplotlib general statistics bar plot size in inches. This does
    not affect R/karyoploteR karyotype figures.

--output-formats FORMAT[,FORMAT...]
    Figure output formats. Choose one or more of pdf, png, jpg. Default: pdf.

--dpi DPI
    Raster output resolution for png/jpg. Must be at least 300. Default: 300.

--palette PALETTE
    TSV palette file for repeat class colours.

--layout horizontal|vertical
    Karyotype layout. Default: horizontal.

--zoom CHR:START-END
    Plot one region only.
```

## Choosing A Window Size

The best `--windowsize` depends on genome size and the question:

- `50000` to `100000`: fine-scale regional inspection.
- `250000` to `500000`: chromosome-arm level summaries.
- `1000000` to `5000000`: fast whole-genome overview.

Smaller windows create larger tables and slower karyotype plots. If you zoom
into a region, reduce the window size so the interval contains enough windows to
be informative.

## Troubleshooting

### The Karyotype Is Empty Or Missing Chromosomes

Check that the chromosome IDs in the annotation file match the `chr` column of
the genome metadata. `--chromtoplot` also expects `chr` IDs, not the display
names unless both columns use the same values.

### A Palette Run Fails With A Missing Class

If `--palette` is used, every class in `--classesorder` must exist in the
palette category column. Add missing classes such as `Ambiguous`, or remove
`--palette`.

### A Zoom Run Fails

Zoom coordinates must be inside the chromosome bounds in the genome metadata.
The zoomed interval also needs to overlap at least one generated window. Use a
smaller `--windowsize` or a wider zoom interval.

### RepeatMasker Without Kimura

RepeatMasker runs work with or without `--kimura`. If `--kimura` is omitted,
TEleVIZion builds an identity-like panel from the per-insertion percent
divergence in the `.out` file. If you want the K2P/Kimura panel instead,
generate the RepeatMasker divergence summary first and pass it with `--kimura`.

### Matplotlib Cache Warnings

If Matplotlib reports that its cache directory is not writable, set
`MPLCONFIGDIR` to a writable directory before running the CLI:

```bash
mkdir -p .cache/matplotlib
MPLCONFIGDIR=.cache/matplotlib python3 scripts/televizion_cli.py --help
```

### R Package Errors

The karyotype step requires `karyoploteR`, `GenomicRanges`, `IRanges`, and
`optparse`. Recreate or update the Conda environment if R reports missing
packages:

```bash
conda env update -f environment.yml --prune
```

## Development Notes

Useful smoke checks:

```bash
python3 scripts/televizion_cli.py --help
Rscript scripts/televizion/plotting_landscape.R --help
```

The generated files in `analyses/` are analysis outputs. Use a unique `--name`
for each run if you want to keep multiple parameter sets.

## License

This project is distributed under the GNU General Public License v3. See
`LICENSE.md` for the full license text.

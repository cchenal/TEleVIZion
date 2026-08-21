# TEleVIZion

TEleVIZion is a research workflow for visualising transposable element (TE) and repeat annotations across genomes. This version  converts **RepeatMasker**, **EDTA** or **TRASH** output into windowed genome summaries, whole-genome bar plots, and chromosome-scale karyotype plots.

TEleVIZion is designed for exploratory and comparative genomics: choose the chromosomes or scaffolds you want to inspect, select an appropriate window size, and generate a consistent set of figures and reusable intermediate tables. TEleVIZion can help answer questions such as:

- Which repeat types dominate a genome by insertion count or base-pair span?
- Are repeats concentrated on particular chromosomes, chromosome arms, or genomic regions?
- Are some regions enriched for relatively young or old repeats?

> **In short:** TEleVIZion helps you move from TE and repeat annotations to interpretable genome-wide and chromosome-scale visualisations.

The workflow supports:

- **RepeatMasker** `.out` annotations, optionally with a `.divsum` Kimura divergence summary.
- **EDTA** `.gff3` annotations containing `classification`, `Name`, and `identity` attributes.
- **TRASH** (beta) summary annotations containing `name`, `ave.score` and `most.freq.value.N` attributes.
- Optional **GC-content** or **gene-content** tracks.
- Whole-genome and per-chromosome summary plots.
- Chromosome-scale karyotype plots.
- Per-class karyotype plots.
- Regional zooming.
- PDF, PNG, and JPG output.

TEleVIZion does **not** run annotation tools; TEleVIZion consumes their output.

---

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [How the workflow works](#how-the-workflow-works)
- [Inputs](#inputs)
- [Command-line options](#command-line-options)
- [Common use cases](#common-use-cases)
- [Output files](#output-files)
- [Repository layout](#repository-layout)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Tutorial](TUTORIAL.md)

---

## Installation

TEleVIZion uses Python packages plus R/Bioconductor packages. The recommended setup is the Conda environment provided in `environment.yml`.

### Recommended: Conda

```bash
conda env create -f environment.yml
conda activate televizion

# Mamba or Micromamba can use the same environment file:
# mamba env create -f environment.yml
# mamba activate televizion
```

### Without Conda

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

### Check the installation

Check both the Python CLI and the R plotting script:

```bash
python3 scripts/televizion_cli.py --help
Rscript scripts/televizion/plotting_landscape.R --help
```

**NB**: Run TEleVIZion from the **repository root**. Some scripts use repository-relative paths.

---

## Quick start

First create and activate the recommended Conda environment:

```bash
conda env create -f environment.yml
conda activate televizion
```

A minimal RepeatMasker run looks like:

```bash
python3 scripts/televizion_cli.py \
  --name MyGenome \
  --genome data/my_genome/chroms.tsv \
  --repeatmasker data/my_genome/repeats.out \
  --windowsize 500000
```

A minimal EDTA run looks like:

```bash
python3 scripts/televizion_cli.py \
  --name MyGenome_EDTA \
  --genome data/my_genome/chroms.tsv \
  --edta data/my_genome/my_genome.fa.mod.EDTA.TEanno.gff3 \
  --windowsize 500000
```

Outputs are written under:

```text
analyses/<name>/
```

**Important**: chromosome IDs must agree

The chromosome or scaffold IDs in the `chr` column of the genome metadata must match the sequence IDs used in the annotation file. `--chromtoplot` also expects these true sequence IDs, rather than the prettier display labels in the metadata `name` column. This is one of the most important checks to make before running TEleVIZion.

For a guided example that starts from genome files and builds the required inputs, see the [tutorial](analyses/README.md).

---

## How the workflow works

At a high level, TEleVIZion performs the following steps.

### 1. Build genome windows

The chromosomes or scaffolds listed in the genome metadata are split into fixed-size, non-overlapping windows.

### 2. Parse repeat annotations

TEleVIZion reads either:

- RepeatMasker `.out` records,
- EDTA GFF3 records, or
- TRASH summary file.

For each annotation it extracts the genomic interval and classification information, including class, family, element name, and divergence or identity information when available.

### 3. Resolve overlapping annotations

Overlapping intervals are split so that base-pair coverage is not naively double-counted.

When overlapping segments have the same repeat class and family, they are counted together. Conflicting repeat types are assigned to `Ambiguous`.

### 4. Aggregate signal

For each genomic window and repeat class/family, TEleVIZion summarises:

- insertion count;
- base-pair span;
- fraction of the window covered by that repeat class.

### 5. Summarise divergence or identity

RepeatMasker Kimura values are grouped into:

```text
0-10, 10-20, 20-30, 30-40, 40-70
```

EDTA identity values, or RepeatMasker-derived identity-like values when no Kimura file is supplied, are grouped into:

```text
1-0.9, 0.9-0.8, 0.8-0.7, 0.7-0.6, 0.6-0
```

TRASH identity values are the average score within a repeated region for a given consensus sequence.

### 6. Generate plots

Python/Matplotlib generates whole-genome and optional per-chromosome summary plots.

R and [karyoploteR](https://academic.oup.com/bioinformatics/article/33/19/3088/3857734?login=true) generate chromosome-scale karyotype plots.

### 7. Write reusable tables

Window-level tables used for the karyotype plots are retained under:

```text
analyses/<name>/karyoplot_tables/
```

These can also be useful for downstream analysis.

---

## Inputs 

### Genome metadata

The [genome metadata file](data/examples/chroms_metadata_GCA_943734735.2.tsv) defines the chromosomes or scaffolds to plot. It is a tab-delimited file.

Columns:

- `chr`: chromosome or scaffold ID used in the annotation file.
- `start`: usually `1`.
- `end`: sequence length in bp.
- `name`: display label used on the karyotype plot.
- `gieStain`: cytoband-like colour label consumed by `karyoploteR`; `chalk` or `gneg` are common simple values.

The `chr` values must match the sequence IDs in your annotation file. The `name` values can be prettier labels, for example using
`OX030907.1` as the true sequence ID while displaying `2RL` on the plot. 

You can [generate a starter metadata file](scripts/utils/create_chroms.py) from the `.fasta` file used to create the annotation file (see the [tutorial](analyses/README.md)). If an NCBI sequence report is available, it can also be supplied to the helper where appropriate.


### Repeat annotation 

#### RepeatMasker 

Use a standard RepeatMasker `.out` file, more complete that the potentially generated `.gff`:

```text
SW   perc perc perc  query       position in query ... repeat class/family ...
18   30.6 1.4  1.4   OX030925.1  7305 7378 ... (ATTTT)n Simple_repeat ...
```

TEleVIZion reads the query sequence, coordinates, repeat name, repeat class/family and percent divergence information.

#### Optional: RepeatMasker Kimura summary

A RepeatMasker Kimura divergence summary can be supplied with `--kimura`.

Example structure:

```text
Class      Repeat                 absLen  wellCharLen  Kimura%
DNA/       Acol_otherMITEs_Eles1  5705    5623         11.70
DNA/CACTA  Acol_CACTA_Ele1        58882   57414        39.95
```

One common way to generate the summary from a RepeatMasker alignment file (initial run of RepeatMasker done with `-a` option) is:

```bash
calcDivergenceFromAlign.pl \
  -s data/my_genome/repeats.divsum \
  data/my_genome/repeats.align
```

If `--kimura` is omitted, RepeatMasker runs still work: TEleVIZion uses per-insertion percent divergence to construct an identity-like summary.

#### EDTA 

Use an EDTA GFF3 annotation. For example:

```text
2RL EDTA LTR_retrotransposon 6037 6194 698 - . ID=TE_homo_4;Name=TE_00001204_LTR;classification=LTR/unknown;identity=0.784;method=homology
```

TEleVIZion uses the following attributes from column 9:

- `classification`: split into class/family;
- `Name`: repeat element name;
- `identity`: match identity used for identity-bin tracks.

Records with `identity=NA` still contribute to class/family abundance but do not contribute to identity-bin summaries.

#### TRASH

TO DO: complete 

### Optional: Contextual tracks 

You may additionally provide **one** contextual track:

- GC content, or
- Gene content.

`--gc` and `--gene-content` are mutually exclusive.

#### GC-content track

You can generate the [GC-content file](scripts/utils/create_gc_content.py) from the `.fasta` file used to create the annotation file (see the [tutorial](analyses/README.md)). 

The table is a five-column tab-delimited coloured track, for example:

```text
chr   start  end    name      itemRgb
chr1  0      10000  gc1-0.42  #66c2a5
```

The chromosome IDs must match the genome metadata. `name` values should be unique; per default: `gc<window>-<gc_ratio>`.

#### Gene-content track

The [gene content file](scripts/utils/create_gene_content.py) from a `.gff` (see the [tutorial](analyses/README.md)).

Example output:

```text
chr  start   end     name       itemRgb
2RL  0       100000  genes1-12  #6aaed6
2RL  100000  200000  genes2-7   #bad6eb
```

The chromosome IDs must match the genome metadata. `name` values should be unique; per default: `genes<window>-<gene_count>`.

### Custom palette

By default, repeat-class colours are assigned from Matplotlib's `turbo` colour map. A [custom TSV palette](data/palettes/personalised_palette_example.tsv) can be supplied with `--palette`.

Columns: 

- Colour name,
- Red (R) value,
- Green (G) value,
- Blue (B) value,
- HEX colour code,
- Repeat type. It can contain comma-separated class aliases (e.g. Unknown,Unclassified,UNKN).

When `--palette` is used together with `--classesorder`, ensure that every class you want to display is represented in the palette. Include `Ambiguous` if it is part of your selected class set.

### Choosing a window size

The best `--windowsize` depends on genome size and the biological question. Smaller windows provide finer spatial resolution, but generate larger intermediate tables and can make karyotype plotting slower.

When using `--zoom`, choose a window size small enough that the selected region contains multiple windows.

---

## Command-line options

The CLI options fall naturally into four groups:

| Category | Options |
|---|---|
| **Core input** | `--name`, `--genome`, `--repeatmasker`, `--edta`, `--kimura` |
| **Genome view** | `--windowsize`, `--chromtoplot`, `--zoom` |
| **Additional tracks** | `--gc`, `--gene-content` |
| **Plotting and output** | `--perchromosome`, `--perclass`, `--classesorder`, `--palette`, `--figsize`, `--layout`, `--output-formats`, `--dpi` |

### Full option reference

| Option | Required? | Value / format | Default | Purpose | Notes / constraints |
|---|---|---|---|---|---|
| `--name` | No | string | `output` | Sets the run name and output prefix. | Outputs are grouped under the corresponding analysis directory. |
| `--genome` | **Yes** | path to TSV | — | Genome metadata defining chromosome/scaffold coordinates and display names. | Chromosome IDs must correspond to those in the TE annotation. |
| `--repeatmasker` | **One of RM / EDTA** | path to `.out` | — | Uses RepeatMasker annotation as TE input. | Mutually exclusive with `--edta`; required when `--kimura` is used. |
| `--edta` | **One of RM / EDTA** | path to GFF3 | — | Uses EDTA annotation as TE input. | Mutually exclusive with `--repeatmasker`; cannot be combined with `--kimura`. |
| `--kimura` | No | path to divergence summary | — | Adds RepeatMasker Kimura divergence information. | Only valid with `--repeatmasker`. |
| `--windowsize` | No | integer, bp | `10000` | Sets the genomic aggregation window size. | Smaller windows increase spatial resolution and table size. |
| `--chromtoplot` | No | comma-separated IDs or `all` | `all` | Selects chromosomes/scaffolds for karyotype plotting. | Uses `chr` IDs, not display names unless they are identical. |
| `--gc` | No | path to TSV | — | Adds a precomputed GC-content track. | Mutually exclusive with `--gene-content`. |
| `--gene-content` | No | path to TSV | — | Adds a precomputed gene-content track. | Mutually exclusive with `--gc`. |
| `--perchromosome` | No | flag | off | Generates general-statistics plots for each chromosome/scaffold. | Add the flag without a value to enable it. |
| `--perclass` | No | flag | off | Generates additional karyotype figures for individual repeat classes. | Add the flag without a value to enable it. |
| `--classesorder` | No | comma-separated classes | automatic | Overrides repeat-class ordering. | Include every class you want retained in ordered plots/tables. |
| `--palette` | No | path to TSV | built-in colours | Overrides the default repeat-class colours. | Useful for consistent colours across analyses. |
| `--figsize` | No | `W,H` | `10,8` | Sets Python general-statistics plot size in inches. | Values are parsed as integers. Does not affect karyoploteR figures. |
| `--output-formats` | No | comma-separated `pdf`, `png`, `jpg` | `pdf` | Selects one or more figure formats. | Values are case-insensitive and duplicates are removed. |
| `--dpi` | No | integer >= 300 | `300` | Sets raster output resolution. | Applies to PNG/JPG; values below 300 are rejected. |
| `--layout` | No | `horizontal` or `vertical` | `horizontal` | Sets karyotype chromosome arrangement. | Use `vertical` to stack chromosomes vertically. |
| `--zoom` | No | `chrom:start-end` | — | Restricts karyotype plotting to a genomic interval. | Example: `Chr1:1000000-2000000`. |

The authoritative CLI help can always be inspected with:

```bash
python3 scripts/televizion_cli.py --help
```

---

## Output files

All outputs are written below:

```text
analyses/<name>/
```

### Whole-genome summary plots

A repeat family may have many short insertions but relatively little genomic coverage, while another may have fewer, much longer insertions that occupy a large fraction of the genome. Reporting both provides a more complete picture of repeat composition and expansion. Both metrics are important because they capture different aspects of **repeat abundance**.

> **In short**: Insertion counts show how frequently a repeat class or family occurs, regardless of the size of individual insertions. Base-pair (bp) span shows how much of the genome is occupied by those repeats.

These are generated for each run:

```text
<name>_whole_genome_stacked_counts_by_class.<format>
<name>_whole_genome_stacked_lengths_by_class.<format>
<name>_whole_genome_contiguous_counts_by_class.<format>
<name>_whole_genome_contiguous_lengths_by_class.<format>
```

- `stacked_counts`: total insertion counts for each repeat class, with bars subdivided (stacked) by repeat family;
- `stacked_lengths`: total base-pair coverage for each repeat class, with bars subdivided (stacked) by repeat family;
- `contiguous_counts`: repeat families ranked horizontally according to their total insertion counts across the whole genome;
- `contiguous_lengths`: repeat families ranked horizontally according to their total base-pair coverage across the whole genome.

### Per-chromosome summary plots

With `--perchromosome`, equivalent plots are generated under:

```text
analyses/<name>/per_chromosome/
```

### Karyotype plots

Whole-genome summaries can hide where repeats are concentrated or depleted. Window-based plots reveal the spatial organisation of repeat content. They allow to:

- Identify repeat-rich and repeat-poor regions, rather than only measuring genome-wide abundance.
- Detect local expansions or clusters of particular repeat classes or families.
- Compare genomic regions, chromosomes, scaffolds, or haplotypes for differences in repeat composition.
- Relate repeats to genomic features, such as genes, centromeres, telomeres, or regions of structural variation.
- Distinguish different patterns of accumulation: high bp coverage might result from a few long insertions or many short ones, which becomes particularly informative when counts and bp span are visualised together along the genome.

> **In short**: Whole-genome summaries tell you how much repeat sequence is present, while genome-wide profiles show where it is and how it is organised.

The standard karyotype outputs include:

```text
<name>_<windowsize>_karyoplot_stacked_percentage_by_class.<format>
<name>_<windowsize>_karyoplot_stacked_counts_by_class.<format>
```

With `--perclass`, additional plots are written under:

```text
analyses/<name>/per_class/
```

With `--zoom`, the genomic interval is incorporated into the output filename.

### Intermediate karyotype tables

Tables are written to:

```text
analyses/<name>/karyoplot_tables/
```

The repeat-class abundance table is:

```text
<name>_<windowsize>_repeat_classes.bed
```

It contains window coordinates and class-level quantities such as:

```text
<class>_count
<class>_count_stacked
<class>_length
<class>_pct
<class>_pct_stacked
```

When Kimura information is supplied:

```text
<name>_<windowsize>_kimura.bed
```

For EDTA, or RepeatMasker runs using the identity-like summary:

```text
<name>_<windowsize>_identity.bed
```

The exported tables are BED-like TSV files intended primarily for plotting. Exported window starts are 0-based.

---

## Common use cases & tutorial

See the [tutorial](analyses/README.md).

---

## Repository layout

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
|   |-- palettes/
|   |   `-- personalised_palette_example.tsv
|   `-- examples/
`-- analyses/
    `-- README.md
```

Important code paths:

- `scripts/televizion_cli.py`: command-line orchestration.
- `scripts/televizion/io.py`: RepeatMasker, Kimura and EDTA parsing.
- `scripts/televizion/aggregation.py`: genome windows, overlap handling and table export.
- `scripts/televizion/plotting_stats.py`: Matplotlib plots and calls into R.
- `scripts/televizion/plotting_landscape.R`: karyotype plotting.
- `scripts/utils/create_chroms.py`: genome metadata helper.
- `scripts/utils/create_gc_content.py`: GC-content helper.
- `scripts/utils/create_gene_content.py`: gene-content helper.

---

## Troubleshooting

### The karyotype is empty or chromosomes are missing

Check that chromosome IDs in the annotation match the `chr` column of the genome metadata.

Also remember that `--chromtoplot` expects these IDs, not the display labels in `name`, unless the two are identical.

### A palette run reports a missing class

When a custom palette is used, ensure that every selected class is represented by the palette.

If `Ambiguous` is included in your class order, make sure it has a palette entry too.

### A zoom run fails or appears empty

Check that:

1. the chromosome ID exists in the genome metadata;
2. the coordinates lie inside the chromosome bounds;
3. the interval overlaps generated windows;
4. the selected `--windowsize` is appropriate for the region.

Try a smaller window size or a wider zoom interval.

### RepeatMasker without Kimura

A Kimura file is optional.

Without `--kimura`, TEleVIZion builds an identity-like summary from the per-insertion percent divergence in the RepeatMasker `.out` file.

Supply a RepeatMasker divergence summary with `--kimura` when you specifically want the Kimura/K2P panel.

### Matplotlib cache warnings

If the Matplotlib cache directory is not writable:

```bash
mkdir -p .cache/matplotlib
MPLCONFIGDIR=.cache/matplotlib python3 scripts/televizion_cli.py --help
```

### Missing R packages

The karyotype stage requires `karyoploteR`, `GenomicRanges`, `IRanges`, and `optparse`.

If you are using the Conda environment, try:

```bash
conda env update -f environment.yml --prune
```
---

## License

TEleVIZion is distributed under the GNU General Public License v3. See [LICENSE.md](LICENSE.md) for the full licence text.

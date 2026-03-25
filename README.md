# ColocBoost Pipeline

This repository contains WDL workflows and R utility scripts for performing colocalization analysis across multiple molecular phenotype modalities (transcriptomics, proteomics, and splicing) using the [colocboost](https://github.com/StatFunGen/colocboost) R package. The pipeline is designed for use within the [All of Us Research Program](https://allofus.nih.gov/) and supports both individual-level genotype data and GWAS summary statistics.

## Table of Contents

- [Workflows](#workflows)
  - [colocboost\_cross\_ome](#colocboost_cross_ome-workflowscolocboostwdl)
  - [colocboost\_individual\_level](#colocboost_individual_level-workflowscolocboost_individual_levelwdl)
  - [colocboost\_summarystats](#colocboost_summarystats-workflowscolocboost_summarystatswdl)
- [R Utility Scripts](#r-utility-scripts)
  - [run\_colocboost.R](#run_colocboostr)
  - [colocboost\_v2.R](#colocboost_v2r)
  - [colocboost\_summarystats.R](#colocboost_summarystatsr)
  - [colocboost\_utils.R](#colocboost_utilsr)
  - [shard\_VCF.R](#shard_vcfr)
- [Data Preparation](#data-preparation)
  - [BED Files (Molecular Phenotype Matrices)](#bed-files-molecular-phenotype-matrices)
  - [VCF Files (Individual-Level Genotypes)](#vcf-files-individual-level-genotypes)
  - [PLINK Binary Genotype Files](#plink-binary-genotype-files)
  - [Covariate Files](#covariate-files)
  - [GWAS Summary Statistics](#gwas-summary-statistics)
  - [Phenotype Table](#phenotype-table)
  - [Sample List](#sample-list)
- [Docker Container](#docker-container)
- [Running the Workflows](#running-the-workflows)

---

## Workflows

The workflows are written in [WDL (Workflow Definition Language)](https://openwdl.org/) and are registered on [Dockstore](https://dockstore.org/). They can be executed on any WDL-compatible platform such as [Terra](https://terra.bio/) or a local [Cromwell](https://cromwell.readthedocs.io/) server.

### `colocboost_cross_ome` (`workflows/colocboost.wdl`)

Runs transcriptome–proteome colocalization analysis using **individual-level genotype data** (VCF format). This is the main cross-omics pipeline and is designed to scale across an entire genome by scattering over chromosomes and individual gene regions.

**Workflow steps:**

1. **`split_vcf_by_chromosome`** – Splits the input genome-wide VCF into per-chromosome VCFs using `bcftools`.
2. **`split_vcf`** – For each chromosome, extracts per-gene VCF regions using coordinates from the trimmed proteome BED file with ±1 Mb padding around each gene.
3. **`colocboost`** – Runs `run_colocboost.R` for each gene region to perform transcriptome–proteome colocalization.

**Inputs:**

| Parameter | Type | Description |
|---|---|---|
| `VCF_workflow` | File | Genome-wide VCF file (bgzipped) |
| `VCF_workflow_index` | File | Tabix index for the VCF (`.tbi`) |
| `transcriptome_bed` | File | BED file containing gene expression values |
| `trimmed_proteome_bed` | File | BED file with gene coordinates used to define regions (columns: chr, start, end, gene\_id) |
| `proteome_bed` | File | BED file containing protein expression values |
| `transcriptome_covars` | File | Covariate file for the transcriptomics assay |
| `proteome_covars` | File | Covariate file for the proteomics assay |
| `docker_image` | String | Docker image to use (e.g., `ghcr.io/aou-multiomics-analysis/colocboost:main`) |
| `memory` | Int | Memory in GB for the colocboost task |
| `disk_space` | Int | Disk space in GB |
| `num_threads` | Int | Number of CPU threads |

**Outputs:**

- One `_colocboost_res.RDS` file per gene containing colocalization results.

---

### `colocboost_individual_level` (`workflows/colocboost_individual_level.wdl`)

Runs multi-task colocalization across **three molecular phenotype modalities** (expression, splicing, and protein) using **individual-level PLINK binary genotype data**. This workflow processes all genes in a single call using `colocboost_v2.R`.

**Workflow steps:**

1. Parses the phenotype table to map genes to their expression, splicing, and protein phenotype IDs across BED files.
2. Loads multi-task regional genotype and phenotype data using `pecotmr::load_multitask_regional_data`.
3. Runs the full `colocboost_analysis_pipeline` including QC filters (MAF, MAC, missingness), xQTL colocalization, and joint GWAS testing.

**Inputs:**

| Parameter | Type | Description |
|---|---|---|
| `PhenotypeTable` | File | Table mapping genes to phenotype IDs (see [Phenotype Table](#phenotype-table)) |
| `ExpressionMatrix` | File | BED file with gene expression values |
| `SplicingMatrix` | File | BED file with splicing quantification values |
| `ProteinMatrix` | File | BED file with protein abundance values |
| `ExpressionCovars` | File | Covariate file for the expression assay |
| `SplicingCovars` | File | Covariate file for the splicing assay |
| `ProteinCovars` | File | Covariate file for the protein assay |
| `PlinkBedGenotypes` | File | PLINK binary genotype file (`.bed`; `.bim` and `.fam` must be co-located) |
| `SampleList` | File | Single-column file listing sample IDs to include |
| `OutputPrefix` | String | Prefix for output file names |

**Outputs:**

- `{OutputPrefix}_colocboost_res.rds` – RDS file containing per-gene colocalization results.

---

### `colocboost_summarystats` (`workflows/colocboost_summarystats.wdl`)

Runs colocalization analysis using **GWAS summary statistics** alongside individual-level molecular phenotype data and dosage genotypes. Supports optional partialization to remove covariate effects from genotypes.

**Workflow steps:**

1. Loads the phenotype BED file and preprocesses covariate data.
2. Extracts genotype dosages and computes the LD (linkage disequilibrium) matrix for the cis-window.
3. Loads and harmonizes GWAS summary statistics for the gene region.
4. Runs `colocboost` with the molecular phenotype residuals and GWAS summary statistics.

**Inputs:**

| Parameter | Type | Description |
|---|---|---|
| `GenotypeDosage` | File | Dosage-encoded VCF file (bgzipped) |
| `GenotypeDosageIndex` | File | Tabix index for the dosage VCF (`.tbi`) |
| `BedFile` | File | BED file containing molecular phenotype values |
| `Covars` | File | Covariate file for the molecular phenotype assay |
| `SumstatsGWAS` | Array[File] | One or more GWAS summary statistics files (tabix-indexed) |
| `SumstatsGWASIndex` | Array[File] | Tabix indices for each GWAS summary statistics file (`.tbi`) |
| `PhenotypeID` | String | Phenotype identifier (gene or protein ID) used for output naming |
| `Partialize` | Boolean | Whether to residualize genotypes against covariates (`true`) or only residualize phenotypes (`false`) |
| `NumPrempt` | Int | Number of preemptible retries |
| `memory` | Int | Memory in GB |
| `disk_space` | Int | Disk space in GB |
| `num_threads` | Int | Number of CPU threads |

**Outputs:**

- `{PhenotypeID}_colocboost_res.RDS` – Full colocalization results object.
- `{PhenotypeID}_colocboost_summary.RDS` – Summary of colocalized signals.

---

## R Utility Scripts

All scripts are located in the `utils/` directory and are copied into the Docker container at build time.

### `run_colocboost.R`

Entry-point script for the `colocboost_cross_ome` workflow. Performs two-phenotype (transcriptome and proteome) colocalization for a single gene region.

**What it does:**
1. Parses command-line arguments.
2. Loads transcriptome and proteome BED files and covariate files.
3. Subsets samples to those present in both assays.
4. Calls `proteome_transcriptome_coloc()` from `colocboost_utils.R`.
5. Saves the result as an RDS file.

**Arguments:**

| Argument | Description |
|---|---|
| `--vcf` | Path to the (region-level) VCF file |
| `--transcriptome_bed` | Path to the transcriptome BED file |
| `--proteome_bed` | Path to the proteome BED file |
| `--transcriptome_covars` | Path to transcriptomics covariate file |
| `--proteome_covars` | Path to proteomics covariate file |
| `--phenotype_id` | Phenotype identifier string (used for output naming) |

---

### `colocboost_v2.R`

Entry-point script for the `colocboost_individual_level` workflow. Performs multi-task colocalization across expression, splicing, and protein modalities for all genes defined in a phenotype table.

**What it does:**
1. Parses command-line arguments and constructs data-frame lookups for BED files and covariate files by modality.
2. Loads the phenotype table and iterates over each row (gene/region).
3. For each gene, calls `wrap_colocboost()`, which:
   - Maps phenotype IDs to BED files via `extract_parameters()`.
   - Loads multi-task regional genotype and phenotype data via `extract_regional_data()` (which calls `pecotmr::load_multitask_regional_data`).
   - Runs `colocboost_analysis_pipeline()` with QC filters, xQTL colocalization, and cross-trait analysis.
4. Saves all results as a single RDS file.

**Arguments:**

| Argument | Description |
|---|---|
| `--PhenotypeTable` | Path to the phenotype table (see [Phenotype Table](#phenotype-table)) |
| `--ExpressionMatrix` | Path to the expression BED file |
| `--SplicingMatrix` | Path to the splicing BED file |
| `--ProteinMatrix` | Path to the protein BED file |
| `--ExpressionCovars` | Path to the expression covariate file |
| `--SplicingCovars` | Path to the splicing covariate file |
| `--ProteinCovars` | Path to the protein covariate file |
| `--PlinkBedGenotypes` | Path to the PLINK `.bed` file (prefix without extension) |
| `--SampleList` | Path to the sample list file |
| `--OutputPrefix` | Output file name prefix |

---

### `colocboost_summarystats.R`

Entry-point script for the `colocboost_summarystats` workflow. Performs colocalization between a molecular phenotype and one or more GWAS summary statistics datasets.

**What it does:**
1. Parses command-line arguments.
2. Loads the phenotype BED file and covariate file.
3. Calls `preprocess_gene_coloc_boost()` to extract genotypes, residualize the phenotype, and compute the LD matrix.
4. Loads and harmonizes GWAS summary statistics using `clean_GWAS_data()`.
5. Runs `colocboost()` using either partialized genotypes (covariates projected out) or raw genotypes depending on the `--Partialize` flag.
6. Saves both the full result and a summary object as RDS files.

**Arguments:**

| Argument | Description |
|---|---|
| `--GenotypeDosage` | Path to the dosage VCF file |
| `--BedFile` | Path to the molecular phenotype BED file |
| `--Covars` | Path to the covariate file |
| `--SumstatsGWAS` | Space-separated paths to one or more GWAS summary statistics files |
| `--PhenotypeID` | Phenotype identifier string (used for output naming) |
| `--Partialize` | `TRUE` to residualize genotypes on covariates; `FALSE` to residualize phenotype only |

---

### `colocboost_utils.R`

Shared library of helper functions sourced by `run_colocboost.R` and `colocboost_summarystats.R`.

| Function | Description |
|---|---|
| `load_gwas_data()` | Loads a tabix-indexed GWAS file and filters to a genomic region. Handles missing `SE`, `FRQ`, or `BETA`/`OR` columns. |
| `prep_covar_data()` | Reads a covariate file, removes batch covariates, transposes to samples-as-rows format. |
| `partialize_data()` | Projects covariates out of both genotype and phenotype matrices using the hat matrix. |
| `residualize_molecular_phenotype_data()` | Residualizes a phenotype vector against all covariates using linear regression. |
| `residualize_genotype()` | Residualizes a genotype vector against covariates. |
| `extract_phenotype_vector()` | Extracts the expression/protein values for a specific gene from a BED data frame. |
| `extract_genotype_vector()` | Queries a tabix-indexed VCF for all variants in the ±1 Mb cis-window of a gene. |
| `clean_genotype_data_dosage()` | Processes dosage-format genotype output: creates variant IDs, mean-imputes missing values, mean-centers. |
| `clean_genotype_data()` | Processes diploid genotype calls (0/0, 0/1, 1/1): converts to 0/1/2 dosage, scales. |
| `extract_variant_metadata()` | Parses variant IDs from the genotype matrix columns into chromosome, position, ref, and alt fields. |
| `extract_GWAS_data()` | Joins variant metadata with GWAS summary statistics by chromosome and position. |
| `preprocess_gene_coloc_boost()` | Master preprocessing function: extracts genotypes and phenotype, residualizes, computes LD matrix, returns a structured list. |
| `extract_GWAS_data_list()` | Applies `extract_GWAS_data()` across a list of GWAS datasets. |
| `wrap_colocboost_list()` | Runs colocboost for a preprocessed object plus a list of summary statistics. |
| `proteome_transcriptome_coloc()` | Two-phenotype colocalization wrapper for a transcriptome–proteome pair at a single locus. |
| `get_cormat()` | Computes the LD correlation matrix from a genotype matrix. |

---

### `shard_VCF.R`

Utility script to generate tabix query ranges from a phenotype BED file with configurable padding. Intended for pre-processing VCF sharding steps outside of WDL.

**Arguments:**

| Argument | Description |
|---|---|
| `--vcf` | Path to the VCF file |
| `--phenotype_bed` | Path to the phenotype BED file |
| `--padding` | Number of base pairs to pad around each gene region |

---

## Data Preparation

### BED Files (Molecular Phenotype Matrices)

BED files store quantified molecular phenotype values (expression, splicing, or protein abundances) for all samples. They must be tab-separated with the following column structure:

| Column | Name | Description |
|---|---|---|
| 1 | `#chr` | Chromosome (e.g., `chr1`) |
| 2 | `start` | Gene/feature start coordinate (0-based) |
| 3 | `end` | Gene/feature end coordinate |
| 4 | `gene_id` | Unique phenotype identifier (gene ID or protein ID) |
| 5+ | Sample IDs | One column per sample containing the quantified phenotype value |

**Example:**
```
#chr	start	end	gene_id	SAMPLE_001	SAMPLE_002	...
chr1	65418	71585	ENSG00000186092.4	1.23	-0.45	...
```

> **Note:** For the `colocboost_cross_ome` workflow, the `proteome_bed` and `transcriptome_bed` files must share overlapping sample ID column names, as the workflow subsets to the intersection of samples present in both.

> **Note:** The `trimmed_proteome_bed` used in `colocboost_cross_ome` is a minimal 4-column BED file (chr, start, end, gene\_id) used only to define genomic regions for VCF sharding; it does not need to contain sample columns.

---

### VCF Files (Individual-Level Genotypes)

VCF files used by the `colocboost_cross_ome` workflow must be:
- Compressed with `bgzip` (`.vcf.bgz` or `.vcf.gz`)
- Indexed with `tabix` (`.tbi`)

For the `colocboost_summarystats` workflow, the genotype file must be in **dosage format**: a VCF-like file where genotype fields contain numeric dosage values (0–2) rather than allele calls.

Variant IDs are constructed internally as `{CHROM}-{POS}-{REF}-{ALT}`.

---

### PLINK Binary Genotype Files

The `colocboost_individual_level` workflow requires genotypes in [PLINK binary format](https://www.cog-genomics.org/plink/1.9/formats#bed) (`.bed`, `.bim`, `.fam`). All three files must be co-located and share the same base filename prefix. Pass the path to the `.bed` file; the script will strip the extension to derive the prefix.

---

### Covariate Files

Covariate files are used to regress technical and biological sources of variation out of the phenotype data (and optionally genotypes). They must be tab-separated with the following structure:

- **Row 1 (header):** `ID` followed by sample IDs across columns.
- **Column 1:** Covariate names (e.g., `PC1`, `PC2`, `sex`, `age`).
- Remaining cells contain covariate values.

**Example:**
```
ID	SAMPLE_001	SAMPLE_002	...
PC1	0.012	-0.034	...
PC2	-0.005	0.021	...
age	45	52	...
```

> **Note:** Any covariate whose name contains the string `batch` will be automatically excluded by `prep_covar_data()`.

---

### GWAS Summary Statistics

GWAS summary statistics files are used by the `colocboost_summarystats` workflow. Each file must be:
- Tab-separated
- **Tabix-indexed** for efficient region-based querying (bgzip + tabix with `-s`/`-b`/`-e` set to the chromosome and position columns)

Required columns (column names are case-sensitive):

| Column | Description |
|---|---|
| `CHR` | Chromosome (numeric, e.g., `1`) |
| `BP` | Base pair position |
| `BETA` or `OR` | Effect size (log odds ratio or beta coefficient) |
| `SE` | Standard error (if absent, computed from `P`) |
| `P` | P-value |
| `N` | Sample size |
| `FRQ` | Effect allele frequency (optional; set to `NA` if missing) |

> **Note:** The file name is used to derive the GWAS trait name: everything after the first `_munged_` suffix is stripped from the basename (e.g., `T2D_munged_20240101.tsv` → `T2D`).

---

### Phenotype Table

Used by the `colocboost_individual_level` workflow (`colocboost_v2.R`). A tab-separated file where each row represents one gene/region to be analyzed.

Required columns:

| Column | Description |
|---|---|
| `gene_id` | Gene identifier |
| `region` | Genomic region string (e.g., `chr1:100000-200000`) |
| `Expression` | Phenotype ID in the expression BED file (or `NA` if not available) |
| `Splicing` | Phenotype ID in the splicing BED file (or `NA` if not available) |
| `Protein` | Phenotype ID in the protein BED file (or `NA` if not available) |

**Example:**
```
gene_id	region	Expression	Splicing	Protein
ENSG00000186092	chr1:64418-72585	ENSG00000186092.4	clu_12345_chr1	PROT_ENSG00000186092
```

> **Note:** Rows with `NA` values in a modality column are automatically excluded for that modality; a gene only needs to be present in at least one modality.

---

### Sample List

Used by the `colocboost_individual_level` workflow. A single-column file (no header, or any header will be treated as the first sample ID) listing the sample IDs to include in the analysis.

**Example:**
```
SAMPLE_001
SAMPLE_002
SAMPLE_003
```

---

## Docker Container

The pipeline uses a Docker container built with [Pixi](https://pixi.sh/), available at:

```
ghcr.io/aou-multiomics-analysis/colocboost:main
```

The container includes:
- R 4.4 with the following packages: `colocboost`, `pecotmr`, `susier`, `tidyverse`, `data.table`, `bedr`, `argparse`, `optparse`, `janitor`, `gprofiler2`, `enrichr`
- Bioinformatics tools: `bcftools`, `bedtools`, `htslib`

The image is automatically built and pushed to the GitHub Container Registry via GitHub Actions on every push to the `main` branch.

To pull the image locally:
```bash
docker pull ghcr.io/aou-multiomics-analysis/colocboost:main
```

---

## Running the Workflows

The WDL workflows can be run on any WDL-compatible platform. Below are examples using [Cromwell](https://cromwell.readthedocs.io/).

### Example: `colocboost_cross_ome`

Create an inputs JSON file (e.g., `inputs_cross_ome.json`):

```json
{
  "colocboost_wdl.VCF_workflow": "gs://bucket/genotypes.vcf.bgz",
  "colocboost_wdl.VCF_workflow_index": "gs://bucket/genotypes.vcf.bgz.tbi",
  "colocboost_wdl.transcriptome_bed": "gs://bucket/expression.bed",
  "colocboost_wdl.trimmed_proteome_bed": "gs://bucket/proteome_regions.bed",
  "colocboost_wdl.proteome_bed": "gs://bucket/protein.bed",
  "colocboost_wdl.transcriptome_covars": "gs://bucket/expression_covars.tsv",
  "colocboost_wdl.proteome_covars": "gs://bucket/protein_covars.tsv",
  "colocboost_wdl.docker_image": "ghcr.io/aou-multiomics-analysis/colocboost:main",
  "colocboost_wdl.memory": 32,
  "colocboost_wdl.disk_space": 100,
  "colocboost_wdl.num_threads": 4
}
```

```bash
cromwell run workflows/colocboost.wdl -i inputs_cross_ome.json
```

### Example: `colocboost_individual_level`

```json
{
  "ColocboostWorkflow.PhenotypeTable": "gs://bucket/phenotype_table.tsv",
  "ColocboostWorkflow.ExpressionMatrix": "gs://bucket/expression.bed",
  "ColocboostWorkflow.SplicingMatrix": "gs://bucket/splicing.bed",
  "ColocboostWorkflow.ProteinMatrix": "gs://bucket/protein.bed",
  "ColocboostWorkflow.ExpressionCovars": "gs://bucket/expression_covars.tsv",
  "ColocboostWorkflow.SplicingCovars": "gs://bucket/splicing_covars.tsv",
  "ColocboostWorkflow.ProteinCovars": "gs://bucket/protein_covars.tsv",
  "ColocboostWorkflow.PlinkBedGenotypes": "gs://bucket/genotypes.bed",
  "ColocboostWorkflow.SampleList": "gs://bucket/samples.txt",
  "ColocboostWorkflow.OutputPrefix": "my_analysis"
}
```

```bash
cromwell run workflows/colocboost_individual_level.wdl -i inputs_individual_level.json
```

### Example: `colocboost_summarystats`

```json
{
  "colocboost_wdl.GenotypeDosage": "gs://bucket/genotypes_dosage.vcf.bgz",
  "colocboost_wdl.GenotypeDosageIndex": "gs://bucket/genotypes_dosage.vcf.bgz.tbi",
  "colocboost_wdl.BedFile": "gs://bucket/expression.bed",
  "colocboost_wdl.Covars": "gs://bucket/covars.tsv",
  "colocboost_wdl.SumstatsGWAS": ["gs://bucket/T2D_munged.tsv.bgz", "gs://bucket/BMI_munged.tsv.bgz"],
  "colocboost_wdl.SumstatsGWASIndex": ["gs://bucket/T2D_munged.tsv.bgz.tbi", "gs://bucket/BMI_munged.tsv.bgz.tbi"],
  "colocboost_wdl.PhenotypeID": "ENSG00000186092",
  "colocboost_wdl.Partialize": false,
  "colocboost_wdl.NumPrempt": 3,
  "colocboost_wdl.memory": 16,
  "colocboost_wdl.disk_space": 50,
  "colocboost_wdl.num_threads": 2
}
```

```bash
cromwell run workflows/colocboost_summarystats.wdl -i inputs_summarystats.json
```

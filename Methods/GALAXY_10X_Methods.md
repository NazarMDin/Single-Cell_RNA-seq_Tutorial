# 🌌 Pre-processing of 10X Single-Cell RNA Datasets (Galaxy)

**Source:** [Galaxy Training Network — Pre-processing of 10X scRNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)  
**Platform:** [Galaxy Project](https://usegalaxy.org) (web-based, no coding required)  
**Time Estimate:** ~1 hour  
**Dataset:** Available on [Zenodo (record 3457880)](https://zenodo.org/record/3457880)

---

## Overview

This tutorial covers the **upstream / pre-processing** stage of single-cell RNA-seq analysis — taking raw sequencing output (BCL/FASTQ files) and producing the count matrix that tools like Scanpy then consume. It uses the **Galaxy platform**, a free, browser-based bioinformatics environment, so no local installation is required.

### Key Questions Addressed

- What is 10X Genomics and how does it work?
- What are BCL and MTX files?
- What is STARsolo and how does it compare to Cell Ranger?
- What is an HDF5 file and why does it matter for scRNA-seq?
- How do we distinguish high-quality cells from low-quality ones?

---

## Learning Objectives

1. Demultiplex single-cell FASTQ data from 10X Genomics
2. Understand transparent matrix file formats (MTX, HDF5)
3. Appreciate the importance of cell quality filtering before downstream analysis

---

## Background: The 10X Genomics Workflow

```
Tissue/Blood Sample
      │
      ▼
[10X Chromium Controller]
Cells are encapsulated in droplets, each with a unique
cell barcode (CB) + unique molecular identifier (UMI)
      │
      ▼
[Library Preparation & Illumina Sequencing]
Produces paired-end FASTQ files:
  Read 1 → Cell barcode + UMI (28 bp)
  Read 2 → cDNA sequence (biological read)
      │
      ▼
[Alignment & Demultiplexing]
STARsolo or Cell Ranger maps reads to genome,
assigns each read to a cell and gene
      │
      ▼
[Count Matrix]
Rows = cells (barcodes), Columns = genes
Values = UMI counts (integers)
      │
      ▼
[Quality Filtering]
Remove empty droplets and dying/damaged cells
      │
      ▼
[Output: MTX or HDF5 format]
Ready for Scanpy / Seurat downstream analysis
```

---

## Key File Formats

| Format | Extension | Description |
|--------|-----------|-------------|
| BCL | `.bcl` | Raw base-call files from Illumina sequencer |
| FASTQ | `.fastq.gz` | Sequencing reads with quality scores |
| MTX | `matrix.mtx` + `barcodes.tsv` + `genes.tsv` | Sparse count matrix (Market Exchange Format) |
| HDF5 | `.h5` / `.h5ad` | Hierarchical Data Format — stores all data in one file |

### Why HDF5?

- **Self-describing**: stores metadata alongside data
- **Portable**: single file contains everything
- **Efficient**: supports random access and compression
- Used by Cell Ranger (`.h5`), Scanpy (`.h5ad`), and scvi-tools

---

## STARsolo vs Cell Ranger

| Feature | STARsolo | Cell Ranger |
|---------|----------|-------------|
| Source | Open-source (STAR aligner) | Proprietary (10X Genomics) |
| Speed | Fast | Slower |
| Cost | Free | Free but closed-source |
| Platform | Galaxy / command line | Standalone |
| Output | MTX, HDF5 | MTX, HDF5 |
| Chemistry support | Broad | Optimised for 10X |

---

## Galaxy Workflow Steps

### Step 1: Upload Data

Upload the FASTQ files (Read 1: barcodes/UMI; Read 2: cDNA) to Galaxy via the **Upload Data** button or from Zenodo.

> Zenodo dataset: https://zenodo.org/record/3457880

### Step 2: Run STARsolo

**Tool:** `RNA STARsolo`

Key parameters to set:
- **Single-cell library type**: 10X Chromium v2/v3
- **Cell barcode whitelist**: provide the 10X barcode whitelist
- **Read 1**: barcode+UMI FASTQ
- **Read 2**: cDNA FASTQ
- **Reference genome**: hg38 (human) or mm10 (mouse)

STARsolo simultaneously aligns reads to the genome and demultiplexes by cell barcode.

### Step 3: Inspect the Output Matrix

STARsolo produces three files in MTX format:
```
matrix.mtx     ← sparse count matrix (genes × cells)
barcodes.tsv   ← one cell barcode per line
genes.tsv      ← gene ID and symbol per line
```

### Step 4: Convert to HDF5 (optional)

**Tool:** `Column Join` or dedicated conversion tools in Galaxy

This bundles the three MTX files into a single `.h5` file, which can be directly read by Scanpy:

```python
import scanpy as sc
adata = sc.read_10x_h5("filtered_feature_bc_matrix.h5")
```

### Step 5: Quality Assessment

Inspect summary statistics to distinguish:

| Cell type | Characteristics |
|-----------|----------------|
| ✅ Valid cell | High UMI count, moderate gene count |
| ❌ Empty droplet | Very low UMI count |
| ❌ Doublet | Abnormally high UMI + gene count |
| ❌ Dying cell | High mitochondrial gene % |

Use the knee plot / barcode rank plot to set filtering thresholds.

---

## Connecting to the Scanpy Pipeline (This Repo)

The output of this Galaxy pre-processing workflow feeds directly into the Scanpy analysis documented in `Notebook/pbmc3k_analysis.ipynb`:

```
Galaxy STARsolo Output (MTX / HDF5)
         │
         ▼
sc.read_10x_mtx() or sc.read_10x_h5()
         │
         ▼
AnnData object (adata)
         │
         ▼
QC → Normalisation → HVG → PCA → UMAP → Clustering → Annotation
(see Notebook/pbmc3k_analysis.ipynb)
```

---

## Prerequisites

- Basic familiarity with Galaxy (see [Galaxy Introduction Tutorial](https://training.galaxyproject.org/training-material/topics/introduction))
- Recommended pre-reading: [Pre-processing of Single-Cell RNA Data (slides)](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing/slides.html)

---

## Answer Histories (Completed Examples)

Pre-run Galaxy histories are available for verification:

| Server | Date |
|--------|------|
| [UseGalaxy.org](https://usegalaxy.org/u/videmp/h/gtn-preprocessing-of-10x-scrna-seq-data-april-2025) | April 2025 |
| [UseGalaxy.eu](https://usegalaxy.eu/u/videmp/h/gtn-preprocessing-of-10x-scrna-seq-data-april-2025) | April 2025 |
| [UseGalaxy.org.au](https://usegalaxy.org.au/u/videmp/h/gtn-preprocessing-of-10x-scrna-seq-data-april-2025) | April 2025 |

---

## References

- Galaxy Training Network. Pre-processing of 10X Single-Cell RNA Datasets. https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html
- Dobin, A. et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21.
- 10X Genomics Chromium Platform: https://www.10xgenomics.com/instruments/chromium-x-series
- Zenodo dataset: https://zenodo.org/record/3457880

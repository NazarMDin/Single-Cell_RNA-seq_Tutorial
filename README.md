# 🔬 Single-Cell RNA-seq Analysis with Scanpy

**Dataset:** 3k PBMCs from a Healthy Donor (10x Genomics)  
**Tool:** Scanpy (Python-based single-cell analysis)  
**Tutorial source:** [scverse/scanpy-tutorials — basic-scrna-tutorial.ipynb](https://github.com/scverse/scanpy-tutorials/blob/main/basic-scrna-tutorial.ipynb)  
**Reference:** Wolf et al., *Genome Biology*, 2018

---

## 📋 Table of Contents

- [Overview](#overview)
- [Pipeline Summary](#pipeline-summary)
- [Dataset](#dataset)
- [Repository Structure](#repository-structure)
- [Step-by-Step Workflow](#step-by-step-workflow)
  - [1. Data Loading](#1-data-loading)
  - [2. Quality Control](#2-quality-control)
  - [3. Normalisation](#3-normalisation)
  - [4. Feature Selection](#4-feature-selection)
  - [5. Dimensionality Reduction](#5-dimensionality-reduction)
  - [6. Clustering](#6-clustering)
  - [7. Marker Gene Identification](#7-marker-gene-identification)
  - [8. Cell Type Annotation](#8-cell-type-annotation)
- [Key Results](#key-results)
- [Installation](#installation)
- [Usage](#usage)
- [Discussion](#discussion)
- [References](#references)

---

## Overview

This repository documents the completion of the Scanpy basic scRNA-seq tutorial using the canonical **PBMC3k dataset** — 2,700 peripheral blood mononuclear cells sequenced on the 10x Genomics Chromium platform. The analysis covers the full single-cell RNA-seq downstream analysis workflow: from raw count matrix to annotated cell type clusters.

Single-cell RNA sequencing (scRNA-seq) allows the measurement of gene expression in individual cells rather than bulk tissue, revealing cellular heterogeneity that would otherwise be masked. This tutorial demonstrates how to process, cluster, and annotate single-cell data using the Scanpy ecosystem.

---

## Pipeline Summary

```
Raw count matrix (10x MTX format)
         │
         ▼
  [AnnData object]
  Load into .h5ad format
         │
         ▼
  [Quality Control]
  Filter low-quality cells & genes
  MT gene % · n_genes · n_counts
         │
         ▼
  [Normalisation + Log transform]
  Total-count normalise → log1p
  Freeze raw counts
         │
         ▼
  [Highly Variable Genes]
  Select top HVGs
         │
         ▼
  [Regress out confounders]
  n_counts · percent_mito
         │
         ▼
  [PCA]
  50 principal components
         │
         ▼
  [Neighbourhood graph]
  k-NN graph (k=10, n_pcs=40)
         │
         ▼
  [UMAP + t-SNE]
  2D embeddings for visualisation
         │
         ▼
  [Louvain clustering]
  Community detection on k-NN graph
         │
         ▼
  [Marker genes]
  rank_genes_groups (Wilcoxon)
         │
         ▼
  [Cell type annotation]
  Cluster labels → known cell types
```

---

## Dataset

**PBMC3k** — 2,700 Peripheral Blood Mononuclear Cells from a Healthy Donor  
Freely available from 10x Genomics.

| Property | Value |
|---|---|
| Cells | 2,700 (after QC: ~2,638) |
| Genes | 32,738 (after filtering: HVGs retained) |
| Sequencer | Illumina NextSeq 500 |
| Reads per cell | ~69,000 |
| Chemistry | 10x Genomics v1 |
| Reference genome | hg19 (GRCh37) |

### Download the raw data

```bash
mkdir -p data
cd data
curl -O https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
```

This creates `data/filtered_gene_bc_matrices/hg19/` containing:
- `matrix.mtx` — sparse count matrix
- `genes.tsv` — gene names
- `barcodes.tsv` — cell barcodes

---

## Repository Structure

```
scanpy-scrna-pbmc3k/
├── README.md                              ← This file
├── environment.yml                        ← Conda environment (recommended)
├── requirements.txt                       ← pip install alternative
├── LICENSE                                ← MIT licence
├── CITATION.cff                           ← Citation metadata
├── notebook/
│   └── pbmc3k_analysis.ipynb             ← Complete analysis notebook
├── methods/
│   └── METHODS.md                        ← Detailed methods for each step
├── discussion/
│   └── DISCUSSION.md                     ← Results interpretation
├── results/
│   ├── pbmc3k.h5ad                       ← Final AnnData object (regenerate by running notebook)
│   ├── FIGURES.md                        ← Figure descriptions and how to regenerate
│   └── figures/
│       ├── qc_violin.png                 ← QC metrics violin plot
│       ├── qc_scatter.png                ← Total counts vs genes scatter
│       ├── hvg_plot.png                  ← Highly variable genes plot
│       ├── pca_variance.png              ← PCA variance ratio (elbow plot)
│       ├── umap_louvain.png              ← UMAP coloured by cluster
│       ├── umap_celltype.png             ← UMAP coloured by annotated cell type
│       ├── umap_marker_expression.png    ← UMAP coloured by known marker genes
│       ├── tsne_louvain.png              ← t-SNE coloured by cluster
│       ├── marker_genes_dotplot.png      ← Dot plot of canonical marker genes
│       ├── marker_genes_heatmap.png      ← Heatmap of top marker genes per cluster
│       ├── stacked_violin_markers.png    ← Stacked violin of marker genes
│       └── celltype_proportions.png      ← Bar chart of cell type proportions
└── data/
    └── .gitkeep                          ← Folder tracked; raw data downloaded separately
```

> **Note:** Raw data files and `pbmc3k.h5ad` are not tracked in git (see `.gitignore`). Run the notebook end-to-end to regenerate all outputs.

---

## Step-by-Step Workflow

### 1. Data Loading

The raw count matrix was read from the 10x MTX format into an **AnnData** object using Scanpy.

```python
import scanpy as sc

adata = sc.read_10x_mtx(
    'data/filtered_gene_bc_matrices/hg19/',
    var_names='gene_symbols',
    cache=True,
)
adata.var_names_make_unique()
```

**Initial object:** 2,700 cells × 32,738 genes

| AnnData slot | Contents |
|---|---|
| `adata.X` | Count matrix (cells × genes) |
| `adata.obs` | Cell metadata (QC metrics, cluster labels) |
| `adata.var` | Gene metadata (HVG flags, means, dispersions) |
| `adata.uns` | Unstructured data (PCA, neighbours, colours) |
| `adata.obsm` | Embeddings (PCA, UMAP, t-SNE coordinates) |

---

### 2. Quality Control

Low-quality cells were filtered based on three metrics:

| Metric | Threshold | Rationale |
|---|---|---|
| `n_genes_by_counts` | < 200 | Likely empty droplet |
| `n_genes_by_counts` | > 2,500 | Likely doublet (two cells captured together) |
| `pct_counts_mt` | > 5% | Dying or damaged cell |

```python
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 2500, :]
adata = adata[adata.obs.pct_counts_mt < 5, :]
```

**After QC:** 2,638 cells retained

| Plot | Description |
|---|---|
| ![QC violin](results/figures/qc_violin.png) | Distribution of QC metrics per cell |
| ![QC scatter](results/figures/qc_scatter.png) | Total counts vs number of genes |

---

### 3. Normalisation

```python
sc.pp.normalize_total(adata, target_sum=1e4)   # scale to 10,000 counts per cell
sc.pp.log1p(adata)                              # log(x + 1) transform
adata.raw = adata                               # freeze for later DE testing
```

---

### 4. Feature Selection

```python
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
```

![HVG plot](results/figures/hvg_plot.png)

---

### 5. Dimensionality Reduction

```python
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.tsne(adata)
```

![PCA variance](results/figures/pca_variance.png)

---

### 6. Clustering

```python
sc.tl.louvain(adata)
```

![UMAP Louvain](results/figures/umap_louvain.png)
![t-SNE Louvain](results/figures/tsne_louvain.png)

**Result:** 8 distinct Louvain clusters identified

---

### 7. Marker Gene Identification

```python
sc.tl.rank_genes_groups(adata, 'louvain', method='wilcoxon')
```

![Marker heatmap](results/figures/marker_genes_heatmap.png)
![Marker dotplot](results/figures/marker_genes_dotplot.png)

---

### 8. Cell Type Annotation

Clusters were annotated using known PBMC marker genes:

| Cluster | Cell type | Key markers |
|---|---|---|
| 0 | CD4+ T cells | IL7R, CCR7 |
| 1 | CD14+ Monocytes | CD14, LYZ |
| 2 | CD4+ T cells (memory) | IL7R, S100A4 |
| 3 | B cells | MS4A1, CD79A |
| 4 | CD8+ T cells | CD8A, CD8B |
| 5 | NK cells | GNLY, NKG7 |
| 6 | FCGR3A+ Monocytes | FCGR3A, MS4A7 |
| 7 | DC / Platelets | FCER1A, PPBP |

![UMAP cell types](results/figures/umap_celltype.png)
![Stacked violin](results/figures/stacked_violin_markers.png)
![Cell type proportions](results/figures/celltype_proportions.png)

---

## Key Results

| Step | Before | After |
|---|---|---|
| Cells | 2,700 | 2,638 |
| Genes | 32,738 | HVGs selected |
| PCs used | — | 40 |
| Clusters | — | 8 Louvain clusters |
| Cell types | — | 8 PBMC subtypes |

**Cell type proportions (approximate):**

| Cell type | ~% of cells |
|---|---|
| CD4+ T cells (all) | ~50% |
| CD14+ Monocytes | ~20% |
| B cells | ~10% |
| CD8+ T cells | ~8% |
| NK cells | ~6% |
| FCGR3A+ Monocytes | ~4% |
| DC / Platelets | ~2% |

---

## Installation

### Option A — Conda (recommended)

```bash
conda env create -f environment.yml
conda activate scanpy-pbmc
```

### Option B — pip

```bash
pip install -r requirements.txt
```

### Verify

```python
import scanpy as sc
sc.logging.print_header()
```

---

## Usage

```bash
# 1. Clone the repo
git clone https://github.com/yourusername/scanpy-scrna-pbmc3k.git
cd scanpy-scrna-pbmc3k

# 2. Set up environment
conda env create -f environment.yml
conda activate scanpy-pbmc

# 3. Download data
mkdir -p data && cd data
curl -O https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
cd ..

# 4. Launch notebook
jupyter notebook notebook/pbmc3k_analysis.ipynb
```

---

## Discussion

The PBMC3k dataset is a well-characterised benchmark, and the results closely match the expected biology of peripheral blood. Eight distinct cell populations were recovered corresponding to the major PBMC subtypes. The clean separation of two monocyte populations (classical CD14+ and non-classical FCGR3A+) and three T cell subsets (CD4+ naive, CD4+ memory, CD8+) demonstrates the resolution advantage of scRNA-seq over bulk sequencing approaches.

See [`discussion/DISCUSSION.md`](discussion/DISCUSSION.md) for a full step-by-step interpretation of every result, comparison to expected cell type proportions, and discussion of limitations.

---

## References

- Wolf, F.A. et al. (2018). Scanpy: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- Satija, R. et al. (2015). Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology*, 33, 495–502.
- Blondel, V.D. et al. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics*, P10008.
- McInnes, L. et al. (2018). UMAP: Uniform Manifold Approximation and Projection. *arXiv*, 1802.03426.
- van der Maaten, L. & Hinton, G. (2008). Visualizing data using t-SNE. *JMLR*, 9, 2579–2605.
- 10x Genomics PBMC3k dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
- Scanpy documentation: https://scanpy.readthedocs.io

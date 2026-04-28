# Single-Cell RNA-seq Tutorial Collection

This repository documents a progressive, hands-on journey through single-cell RNA-seq analysis — from raw sequencing data to annotated cell type clusters. It covers three complementary tutorials spanning upstream pre-processing, core data structures, and full downstream analysis.

---

## Tutorials in This Repository

| # | Tutorial | Platform | Focus |
|---|----------|----------|-------|
| 1 | [10X Pre-processing (Galaxy)](#tutorial-1-pre-processing-of-10x-scrna-seq-data) | Galaxy (GUI) | BCL → count matrix |
| 2 | [AnnData Getting Started](#tutorial-2-getting-started-with-anndata) | Python / Jupyter | Core data structure |
| 3 | [PBMC3k Scanpy Analysis](#tutorial-3-single-cell-rna-seq-analysis-with-scanpy) | Python / Jupyter | Full downstream pipeline |

---

## Tutorial 1: Pre-processing of 10X scRNA-seq Data

**Source:** [Galaxy Training Network](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)  
**Platform:** [Galaxy Project](https://usegalaxy.org) — no coding required  
**Notes:** [`Methods/GALAXY_10X_Methods.md`](Methods/GALAXY_10X_Methods.md)

This tutorial covers the **upstream** stage: taking raw 10X Genomics FASTQ files and producing a count matrix using **STARsolo** on the Galaxy platform.

### What you'll learn
- How the 10X Chromium droplet-based sequencing works
- File formats: BCL, FASTQ, MTX, HDF5
- Using STARsolo as a free, open-source alternative to Cell Ranger
- Distinguishing valid cells from empty droplets and doublets

### Workflow Summary
```
FASTQ files (Read1: barcode+UMI, Read2: cDNA)
       │
       ▼
  [STARsolo on Galaxy]
  Align + demultiplex
       │
       ▼
  Count Matrix (MTX or HDF5)
  Rows = cells, Columns = genes
       │
       ▼
  Quality assessment → ready for Scanpy
```

---

## Tutorial 2: Getting Started with AnnData

**Source:** [scverse-tutorials (readthedocs)](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)  
**Platform:** Python / Jupyter Notebook  
**Notes:** [`Notebook/Anndata.ipynb`](Notebook/Anndata.ipynb)

This tutorial introduces **AnnData** — the fundamental data structure underlying Scanpy and the entire scverse ecosystem. Understanding AnnData is essential for working with any single-cell Python tool.

### What you'll learn
- The anatomy of an `AnnData` object and all its slots
- How to load, inspect, and subset `.h5ad` files
- Working with sparse matrices for memory efficiency
- Difference between views and copies

### AnnData Slot Quick Reference

| Slot | Contents |
|------|----------|
| `adata.X` | Active data matrix (cells × genes) |
| `adata.obs` | Cell-level metadata (QC, cluster labels) |
| `adata.var` | Gene-level metadata (HVG flags, stats) |
| `adata.obsm` | Multi-dim embeddings (PCA, UMAP, t-SNE) |
| `adata.obsp` | Cell-cell pairwise matrices (distances, adjacency) |
| `adata.layers` | Alternative matrices (raw counts, CPM, etc.) |
| `adata.uns` | Unstructured metadata (colours, parameters) |

```python
import anndata
adata = anndata.read_h5ad('path/to/file.h5ad')
print(adata)   # inspect the object
```

---

## Tutorial 3: Single-Cell RNA-seq Analysis with Scanpy

**Source:** [scverse/scanpy-tutorials — basic-scrna-tutorial.ipynb](https://github.com/scverse/scanpy-tutorials/blob/main/basic-scrna-tutorial.ipynb)  
**Dataset:** 3k PBMCs from a Healthy Donor (10x Genomics)  
**Notebook:** [`Notebook/pbmc3k_analysis.ipynb`](Notebook/pbmc3k_analysis.ipynb)  
**Reference:** Wolf et al., *Genome Biology*, 2018

This is the full **downstream analysis** pipeline — starting from a count matrix and ending with annotated cell type clusters.

### Pipeline Summary

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

### Dataset

**PBMC3k** — 2,700 Peripheral Blood Mononuclear Cells from a Healthy Donor

| Property | Value |
|----------|-------|
| Cells | 2,700 (after QC: ~2,638) |
| Genes | 32,738 (after filtering: HVGs retained) |
| Sequencer | Illumina NextSeq 500 |
| Chemistry | 10x Genomics v1 |
| Reference genome | hg19 (GRCh37) |

### Key Results

| Step | Before | After |
|------|--------|-------|
| Cells | 2,700 | 2,638 |
| Genes | 32,738 | HVGs selected |
| PCs used | — | 40 |
| Clusters | — | 8 Louvain clusters |
| Cell types | — | 8 PBMC subtypes |

### Cell Type Annotations

| Cluster | Cell Type | Key Markers |
|---------|-----------|-------------|
| 0 | CD4+ T cells | IL7R, CCR7 |
| 1 | CD14+ Monocytes | CD14, LYZ |
| 2 | CD4+ T cells (memory) | IL7R, S100A4 |
| 3 | B cells | MS4A1, CD79A |
| 4 | CD8+ T cells | CD8A, CD8B |
| 5 | NK cells | GNLY, NKG7 |
| 6 | FCGR3A+ Monocytes | FCGR3A, MS4A7 |
| 7 | DC / Platelets | FCER1A, PPBP |

---

## Repository Structure

```
Single-Cell_RNA-seq_Tutorial/
├── README.md                                  ← This file
├── environment.yml                            ← Conda environment
├── requirements.txt                           ← pip alternative
├── CITATION.cff                               ← Citation metadata
│
├── Notebook/
│   └── pbmc3k_analysis.ipynb                 ← Tutorial 3: Full Scanpy pipeline
│   └── Anndata.ipynb                         ← Tutorial 2: Anndata Workflow
│
├── Methods/
│   └── sc-RNA_Methods.md                     ← Detailed methods (Tutorial 3)
│   └── ANNDATA_Methods.md                    ← Detailed methods (Tutorial 2)
│   └── GALAXY_10X_Methods.md                 ← Detailed methods (Tutorial 1)
│
├── Discussion/
│   └── DISCUSSION.md                         ← Results interpretation (Tutorial 3)
│
├── Results/
│   └── Figures_sc                            ← Generated plots from Tutorial 3
│   └── Figures_anndata                       ← Generated plots from Tutorial 2
│   └── Figures_preprocessing                 ← Generated plots from Tutorial 1
│
└── Data/
    └── .gitkeep                              ← Raw data downloaded separately
```

> **Note:** Raw data files and `.h5ad` outputs are not tracked in git. Run the notebook to regenerate them.

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
git clone https://github.com/NazarMDin/Single-Cell_RNA-seq_Tutorial.git
cd Single-Cell_RNA-seq_Tutorial

# 2. Set up environment
conda env create -f environment.yml
conda activate scanpy-pbmc

# 3. Download PBMC3k data (for Tutorial 3)
mkdir -p Data && cd Data
curl -O https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
tar -xzf pbmc3k_filtered_gene_bc_matrices.tar.gz
cd ..

# 4. Launch notebook (Tutorial 3)
jupyter notebook Notebook/pbmc3k_analysis.ipynb
```

For Tutorial 1 (Galaxy), visit [usegalaxy.org](https://usegalaxy.org) — no local installation needed.  
For Tutorial 2 (AnnData), read [`Notebook/Anndata.ipynb`](Notebook/Anndata.ipynb).

---

## Discussion

See [`Discussion/DISCUSSION.md`](Discussion/DISCUSSION.md) for a detailed interpretation of the PBMC3k clustering results, comparison to expected cell type proportions, and discussion of limitations.

---

## References

- Wolf, F.A. et al. (2018). Scanpy: large-scale single-cell gene expression data analysis. *Genome Biology*, 19, 15.
- Satija, R. et al. (2015). Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology*, 33, 495–502.
- Dobin, A. et al. (2013). STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*, 29(1), 15–21.
- Blondel, V.D. et al. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics*, P10008.
- McInnes, L. et al. (2018). UMAP: Uniform Manifold Approximation and Projection. *arXiv*, 1802.03426.
- van der Maaten, L. & Hinton, G. (2008). Visualizing data using t-SNE. *JMLR*, 9, 2579–2605.
- 10x Genomics PBMC3k dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.1.0/pbmc3k
- Scanpy documentation: https://scanpy.readthedocs.io
- AnnData documentation: https://anndata.readthedocs.io
- Galaxy Training Network: https://training.galaxyproject.org

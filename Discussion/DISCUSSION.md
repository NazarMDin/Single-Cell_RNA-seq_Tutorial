# Discussion

## Overview

The goal of this tutorial was to perform a complete single-cell RNA-seq downstream analysis of the PBMC3k dataset using Scanpy — from loading a raw 10x count matrix to annotating biologically meaningful cell types. The final result was eight well-separated, interpretable clusters corresponding to the major immune cell populations found in peripheral blood, closely matching both the expected biology and the results described in the original Scanpy tutorial and Seurat benchmark.

---

## Quality Control

The QC filtering step removed a relatively small number of cells (62 cells, or ~2.3% of the original 2,700). Most removals were due to high mitochondrial gene content (>5%), which is a reliable indicator of cell stress or membrane damage — in dying cells, cytoplasmic mRNA leaks out while mitochondrial transcripts, enclosed within the mitochondria, are retained. A smaller number of cells were removed for having an unusually high gene count (>2,500 genes), which is a common proxy for doublets — droplets that accidentally captured two cells rather than one.

The upper threshold for gene count (2,500) is relatively conservative for PBMCs. Some cell types, such as monocytes, can naturally express a higher number of genes than lymphocytes, so care must be taken not to set this threshold too low. The chosen threshold strikes a reasonable balance between removing doublets and retaining genuine high-complexity single cells.

---

## Normalisation and Feature Selection

Total-count normalisation (scaling each cell to 10,000 counts) effectively removes the most common source of technical variation in scRNA-seq — differences in sequencing depth between cells. After log transformation, the data distribution becomes closer to normal, which is an assumption of many downstream statistical methods.

The selection of 2,000 highly variable genes is a standard approach for reducing dimensionality while retaining biologically informative variation. Genes with consistent expression across all cells contribute little to distinguishing cell types and mostly add noise. The Seurat v3 flavour of HVG selection accounts for the mean-variance relationship that is intrinsic to count data, making it more appropriate than simpler variance-based approaches.

Regressing out total counts and mitochondrial percentage removes residual technical variation that survived normalisation. This is particularly important for monocytes, which have higher baseline transcriptional activity and could otherwise cluster separately based on library size rather than true identity.

---

## Dimensionality Reduction

### PCA

The PCA variance ratio plot showed that the first ~10–15 PCs captured the majority of biological variance, with subsequent PCs adding diminishing information. Using 40 PCs for the neighbourhood graph is a slightly conservative choice that ensures no biologically meaningful signal is discarded, at the cost of including a few noise-dominated components. In practice, the results are not highly sensitive to this choice within a reasonable range (30–50 PCs).

### UMAP vs t-SNE

Both UMAP and t-SNE produced interpretable 2D representations of the data, but UMAP is generally preferred for scRNA-seq visualisation because it better preserves the global structure of the data — relative distances between well-separated clusters in UMAP are more meaningful than in t-SNE, which can distort global topology. The T cell clusters, for example, appear appropriately close to each other in UMAP (reflecting their shared lymphoid identity) while remaining distinctly separated internally.

---

## Clustering

The Louvain algorithm identified 8 clusters at the default resolution. The number of clusters is sensitive to the resolution parameter — higher values produce more, finer-grained clusters, while lower values produce fewer, broader clusters. For this dataset, the default resolution yields biologically interpretable clusters that correspond well to known PBMC subtypes without over-fragmenting.

A limitation of Louvain clustering is that it does not guarantee reproducibility across runs at exactly the same cluster assignments, as it involves a random initialisation step. The Leiden algorithm (a more recent improvement on Louvain) addresses some stability issues and is generally recommended for new analyses, though both produce similar results on well-structured datasets like PBMC3k.

---

## Marker Gene Analysis and Cell Type Annotation

The Wilcoxon rank-sum test was used for marker gene identification. This non-parametric test is appropriate for scRNA-seq data because it does not assume a normal distribution, which would be violated by the highly zero-inflated count data. The test identifies genes that are consistently upregulated in one cluster relative to all others.

Cell type annotation was performed by matching top marker genes to known PBMC cell type signatures:

**T cell distinction** was achieved by examining CD4/CD8 co-receptor expression and naive/memory markers (CCR7 for naive, S100A4 for memory). The fact that three distinct T cell populations were recovered (CD4+ naive, CD4+ memory, CD8+) demonstrates the fine-grained resolution of scRNA-seq relative to flow cytometry, which would require multiple surface markers to achieve the same separation.

**Monocyte subsets** were cleanly separated into classical (CD14+, high LYZ expression, inflammatory) and non-classical (FCGR3A+/CD16+, patrolling) monocytes. This biologically meaningful distinction would be invisible in bulk RNA-seq data.

**The small DC/Platelet cluster** is a known challenge in PBMC datasets. Dendritic cells (FCER1A+) and platelets (PF4+) often co-cluster in PBMC3k because both are rare and their transcriptional profiles can overlap in lower-dimensional representations. In a more complete analysis, further sub-clustering or marker-based disambiguation would separate these populations.

---

## Comparison to Expected Results

The cell type proportions recovered here closely match published estimates for healthy human PBMCs from the literature and the original Seurat tutorial:

| Cell type | Expected in PBMCs | Recovered |
|---|---|---|
| CD4+ T cells (all) | ~50% | ~50% (two clusters) |
| CD14+ Monocytes | ~15–20% | ~20% |
| B cells | ~5–10% | ~10% |
| CD8+ T cells | ~10% | ~8% |
| NK cells | ~5–10% | ~6% |
| FCGR3A+ Monocytes | ~3–5% | ~4% |
| DC / Platelets | ~1–3% | ~2% |

The agreement between observed and expected proportions validates the quality of the preprocessing, normalisation, and clustering steps.

---

## Limitations

**Annotation subjectivity:** Manual cell type annotation based on a handful of marker genes is subject to interpretation. In a research context, annotation would ideally be validated using reference-based methods (e.g., sctype, SingleR, or CellTypist) that compare against curated atlases.

**Single sample:** This analysis uses a single healthy donor. Batch effects, donor-to-donor variability, and disease states are not captured. A multi-sample analysis with integration (Harmony, scVI, BBKNN) would be required for more generalisable conclusions.

**Reference genome:** The dataset was aligned to hg19 (GRCh37). Modern analyses typically use GRCh38, and gene annotations differ slightly between the two versions.

**Doublet detection:** Only a simple gene count threshold was used to remove potential doublets. Dedicated doublet detection tools (Scrublet, DoubletFinder) would provide more accurate doublet removal.

---

## Conclusion

This tutorial successfully demonstrated the complete scRNA-seq downstream analysis workflow using Scanpy on the PBMC3k benchmark dataset. Eight biologically interpretable cell populations were identified with high confidence. The results are fully consistent with the known immunology of peripheral blood and with published analyses of this dataset, validating the computational approach. The Scanpy ecosystem provides a highly capable and well-documented framework for single-cell analysis that scales from small tutorial datasets to multi-million cell atlases.

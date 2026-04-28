# Methods

Detailed description of every step in the scRNA-seq analysis pipeline.

---

## 1. Data Loading

Raw count data was downloaded from the 10x Genomics website as a compressed MTX archive and loaded into an **AnnData** object using `sc.read_10x_mtx()`. Gene names were used as variable identifiers, and duplicate gene names were made unique with `adata.var_names_make_unique()`.

**Initial object:** 2,700 cells × 32,738 genes

---

## 2. Quality Control

Three QC metrics were computed per cell using `sc.pp.calculate_qc_metrics()`:

| Metric | Description |
|---|---|
| `n_genes_by_counts` | Number of genes with ≥1 count |
| `total_counts` | Total UMI count per cell |
| `pct_counts_mt` | Percentage of counts from mitochondrial genes |

Mitochondrial genes were identified by the `MT-` prefix in gene names.

**Filtering thresholds applied:**

| Filter | Threshold | Cells removed |
|---|---|---|
| Minimum genes | < 200 | ~0 (pre-filtered by 10x) |
| Maximum genes (doublet proxy) | > 2,500 | ~18 |
| Mitochondrial content | > 5% | ~57 |
| Minimum cells per gene | < 3 | ~18,000 genes |

**After QC:** 2,638 cells × 13,714 genes

---

## 3. Normalisation

**Total-count normalisation** (`sc.pp.normalize_total`, target = 10,000 counts per cell) was applied to correct for differences in sequencing depth. Each cell's count vector was scaled so that the total counts equal 10,000.

**Log transformation** (`sc.pp.log1p`) was then applied: expression values were transformed as log(x + 1). This stabilises variance and makes the distribution more suitable for PCA and linear models.

The normalised, log-transformed values were frozen in `adata.raw` before further processing to preserve them for differential expression testing.

---

## 4. Highly Variable Gene Selection

Highly variable genes (HVGs) were identified using `sc.pp.highly_variable_genes()` with the following parameters:

| Parameter | Value |
|---|---|
| `min_mean` | 0.0125 |
| `max_mean` | 3 |
| `min_disp` | 0.5 |

These parameters select genes with intermediate-to-high mean expression and high normalised dispersion — genes that vary meaningfully across cells rather than being uniformly expressed or expressing only noise.

**Result:** ~1,838–2,000 HVGs selected

**Confounder regression** (`sc.pp.regress_out`) was applied for `total_counts` and `pct_counts_mt` to remove residual technical variation. Expression values were then **scaled** to zero mean and unit variance, with values clipped at 10 standard deviations (`sc.pp.scale(max_value=10)`).

---

## 5. Principal Component Analysis

PCA was computed on the scaled HVG matrix using `sc.tl.pca(svd_solver='arpack', n_comps=50)`. The ARPACK solver was chosen for computational efficiency on sparse matrices.

The **variance ratio plot** (elbow plot) was inspected to determine the appropriate number of PCs. Based on the inflection point in the cumulative variance curve, **40 PCs** were selected for downstream neighbourhood graph construction.

---

## 6. Neighbourhood Graph Construction

A **k-nearest neighbour (k-NN) graph** was constructed in PCA space using `sc.pp.neighbors(n_neighbors=10, n_pcs=40)`. Each cell is connected to its 10 nearest neighbours in 40-dimensional PCA space. This graph encodes transcriptional similarity between cells and forms the basis of both clustering and UMAP embedding.

---

## 7. Dimensionality Reduction for Visualisation

**UMAP** (`sc.tl.umap`) was computed from the k-NN graph using the UMAP algorithm (McInnes et al. 2018). UMAP preserves both local and global structure and is the primary visualisation used throughout.

**t-SNE** (`sc.tl.tsne`) was computed for comparison. t-SNE preserves local structure but may distort global relationships between clusters.

---

## 8. Clustering

**Louvain community detection** (`sc.tl.louvain`) was applied to the k-NN graph at the default resolution (1.0). The algorithm iteratively merges nodes to maximise graph modularity, identifying communities of transcriptionally similar cells without requiring the number of clusters to be specified in advance.

**Result:** 8 Louvain clusters

---

## 9. Marker Gene Identification

Marker genes for each cluster were identified using `sc.tl.rank_genes_groups()` with the **Wilcoxon rank-sum test** (`method='wilcoxon'`). Each cluster was compared against all remaining cells. The Wilcoxon test was chosen because it is non-parametric, does not assume a normal distribution, and is robust to the zero-inflation characteristic of scRNA-seq count data.

Results were visualised using rank gene plots, dot plots, stacked violin plots, and heatmaps.

---

## 10. Cell Type Annotation

Cluster identities were assigned by comparing the top marker genes for each cluster against canonical PBMC cell type markers from the literature (Satija et al. 2015, 10x Genomics Human PBMC reference). Annotation was performed manually.

| Cluster | Cell type | Key markers used |
|---|---|---|
| 0 | CD4+ T cells | IL7R, CCR7 |
| 1 | CD14+ Monocytes | CD14, LYZ |
| 2 | CD4+ T cells (memory) | IL7R, S100A4 |
| 3 | B cells | MS4A1, CD79A |
| 4 | CD8+ T cells | CD8A, CD8B |
| 5 | NK cells | GNLY, NKG7 |
| 6 | FCGR3A+ Monocytes | FCGR3A, MS4A7 |
| 7 | DC / Platelets | FCER1A, PPBP |

---

## Software Versions

| Package | Version |
|---|---|
| scanpy | ≥ 1.9.0 |
| anndata | ≥ 0.9.0 |
| numpy | ≥ 1.23.0 |
| pandas | ≥ 1.5.0 |
| matplotlib | ≥ 3.6.0 |
| umap-learn | ≥ 0.5.3 |
| leidenalg | ≥ 0.9.0 |
| louvain | ≥ 0.8.0 |

Exact versions can be reproduced using the provided `environment.yml` or `requirements.txt`.

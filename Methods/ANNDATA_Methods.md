# 📦 Getting Started with AnnData

**Source:** [scverse-tutorials — anndata_getting_started](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)  
**Package:** `anndata` (core data structure of the scverse ecosystem)  
**Dataset:** Preprocessed PBMC3k (2,638 cells × 11,505 genes, `.h5ad` format)

---

## What is AnnData?

`AnnData` (Annotated Data) is the central data structure used by Scanpy and the broader scverse ecosystem. It stores a data matrix alongside rich metadata for both cells (observations) and genes (variables) in a single, self-contained object.

```
AnnData object with n_obs × n_vars = 2638 × 11505
    obs:  'n_genes', 'percent_mito', 'n_counts', 'louvain_cell_types'
    var:  'gene_names', 'n_cells', 'gene_ids'
    uns:  'louvain', 'louvain_colors', 'pca'
    obsm: 'X_pca', 'X_tsne', 'X_umap'
    layers: 'raw'
    obsp: 'distances_all'
```

---

## AnnData Slot Reference

| Slot | Shape / Type | Description |
|------|-------------|-------------|
| `adata.X` | `(n_obs, n_vars)` sparse | **Active data matrix** — normalized/log-transformed counts used for analysis |
| `adata.obs` | `(n_obs,)` DataFrame | **Cell-level metadata** — QC metrics, cluster labels, etc. |
| `adata.var` | `(n_vars,)` DataFrame | **Gene-level metadata** — HVG flags, mean expression, dispersion |
| `adata.uns` | dict | **Unstructured metadata** — color palettes, PCA results, neighbour graphs |
| `adata.obsm` | dict of arrays | **Multi-dim cell annotations** — PCA, UMAP, t-SNE coordinates |
| `adata.varm` | dict of arrays | **Multi-dim gene annotations** — loadings, etc. |
| `adata.obsp` | dict of sparse matrices | **Cell-cell pair annotations** — distance/adjacency matrices |
| `adata.varp` | dict of sparse matrices | **Gene-gene pair annotations** |
| `adata.layers` | dict of matrices | **Alternative data matrices** (same shape as X) — e.g. raw counts |
| `adata.obs_names` | Index | Cell barcodes |
| `adata.var_names` | Index | Gene symbols |

---

## Loading a `.h5ad` File

```python
import anndata

# Load from disk
adata = anndata.read_h5ad('path/to/file.h5ad')

# Or download using pooch
import pooch
datapath = pooch.retrieve(
    path=pooch.os_cache("scverse_tutorials"),
    url="https://exampledata.scverse.org/tutorials/scverse-getting-started-anndata-pbmc3k_processed.h5ad",
    known_hash="md5:b80deb0997f96b45d06f19c694e46243",
)
adata = anndata.read_h5ad(datapath)
```

---

## Key Concepts Covered

### 1. The Active Data Matrix (`adata.X`)

- Stored as a **scipy sparse matrix** (CSR format) for memory efficiency
- Rows = cells, Columns = genes
- Contains normalized + log1p-transformed counts by default

```python
# Fraction of non-zero entries (sparsity check)
print(adata.X.nnz / (adata.X.shape[0] * adata.X.shape[1]))
# ≈ 0.068 → only ~7% of entries are non-zero
```

### 2. Layers — Alternative Data Versions

Layers hold additional matrices with the **same shape as `adata.X`**:

```python
# Access raw integer counts
adata.layers["raw"]

# Add a new layer (e.g. CPM-normalised)
import numpy as np
adata.layers["cpm"] = adata.X / adata.X.sum(axis=1) * 1e6
```

### 3. Cell and Gene Annotations (`obs` / `var`)

```python
# Cell metadata
adata.obs.head()          # DataFrame with rows = cells
adata.obs["louvain_cell_types"].value_counts()

# Gene metadata
adata.var.head()          # DataFrame with rows = genes
```

### 4. Indexing with `obs_names` and `var_names`

```python
# Boolean subsetting
adata[adata.obs["louvain_cell_types"] == "CD4 T cells", :]

# Named subsetting
adata[["AAACATACAACCAC-1", "AAACATTGAGCTAC-1"], ["CD3D", "MS4A1"]]
```

### 5. Multi-dimensional Annotations (`obsm` / `varm`)

```python
# UMAP coordinates (shape: n_cells × 2)
adata.obsm["X_umap"]

# PCA coordinates (shape: n_cells × n_PCs)
adata.obsm["X_pca"]
```

### 6. Pairwise Annotations (`obsp` / `varp`)

```python
# Cell-cell distance matrix (sparse)
adata.obsp["distances_all"]
```

### 7. Unstructured Metadata (`uns`)

```python
# Cluster color palette
adata.uns["louvain_colors"]

# PCA variance explained
adata.uns["pca"]["variance_ratio"]
```

### 8. Views vs. Copies

```python
# Slicing creates a VIEW (memory-efficient, but read-only)
view = adata[:100, :]
print(view.is_view)   # True

# Convert to an independent copy
copy = view.copy()
print(copy.is_view)   # False
```

---

## Saving AnnData

```python
adata.write_h5ad("output/my_results.h5ad")
```

---

## Key Takeaways

1. **AnnData is the universal container** for single-cell data in the scverse ecosystem — used by Scanpy, scvi-tools, squidpy, and others.
2. **`adata.X`** holds the primary matrix; **`layers`** preserves alternative representations (raw counts, CPM, etc.).
3. **Sparse matrices** are used throughout to handle the high dimensionality and sparsity of scRNA-seq data efficiently.
4. **Slicing returns views** — always `.copy()` if you intend to modify a subset.
5. The `.h5ad` file format is an HDF5-based format that stores all slots in a single portable file.

---

## References

- AnnData documentation: https://anndata.readthedocs.io
- scverse tutorials: https://scverse-tutorials.readthedocs.io
- Tutorial notebook (download): https://scverse-tutorials.readthedocs.io/en/latest/_sources/notebooks/anndata_getting_started.ipynb

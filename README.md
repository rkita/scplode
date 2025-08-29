# scplode

**"Scplode" your disk space once and get faster, memory-efficient access to your single cell data.**

## Why scplode?

- **Disk space is cheap. Time and memory are not.** 
- A one-time memory-map creation pays dividends in faster data access (~0-3x compared to anndata backed).  
- Fast performance for ML data loading and manual data exploration.
- Familiar `read_h5ad()` and indexing commands for seamless data exploration.
- Simple, lightweight, and familiar API enables easy transition from your current data access.

## How it compares

| Tool | Pros | Cons |
|------|------|------|
| **Backed AnnData** | Classic API | Slower access times |
| **Bionemo-scdl** | Similar memory mapping | Loss of obs/var , ML-focused |
| **ScDataset** | No format conversion necessary | Loss of obs/var, ML-focused |
| **scplode** | Simple, fast, memory-efficient, familiar API | Requires disk space for mmap. |

## Installation

```bash
pip install scplode
```

## Quick Start

```python
import scplode as sp

#On the first time, memory maps are created (be patient!)
adata = sp.read_h5ad('your_data.h5ad')

#Access your data like usual, for example:
adata[0:10]
#or adata[cell_barcodes] 

#Subsequent calls, memory maps are identified and quickly accessed. 
adata = sp.read_h5ad('your_data.h5ad')

```
```
#OUTPUT
[INFO] Creating index
[INFO] Creating index: reading adata file
[INFO] Creating index: writing mmap dat file
100%
 1/1 [00:00<00:00, 84.30it/s]
[INFO] Creating index: packing obs
[INFO] Creating index: packing var
[INFO] Loading index: obs
[INFO] Loading index: var
[INFO] Loading index: dat (implicitly)
```

```python
#Accessed data return AnnData object, as usual:
adata[0:10]
```

```
#OUTPUT
View of AnnData object with n_obs × n_vars = 10 × 50
    obs: 'cell_type'
    var: 'gene_name'
```
```python
#ML data loaders can use .get to skip AnnData object creation
adata.get([indices])

X = adata.get([barcodes])
type(X)
```
```
#OUTPUT
numpy.ndarray
```

## Examples

Located in examples directory:
- `00_example`: Demonstrates use of scplode, and tests for equivalent results
- `01_benchmark`: Compares scplode and anndata for random and contiguous indexing
- `02_state_benchmark`: Compares scplode and anndata when using scplode with the Arc Institute State Model's data loader. This demonstrates the performance improvement in ML applications. See https://github.com/rkita/cell-load/tree/scplode-integration for an example of transitioning to scplode-based access.

## Requirements

- Python 3.8+
- AnnData
- Pandas
- Numpy

## License

MIT
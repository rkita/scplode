# scplode

**"Scplode" your disk space once and get faster access to your single cell data with less memory.**

## Why scplode?

- **Disk space is cheap. Time and memory are not.** 
- A one-time index pays dividends in faster data access
- Familiar `read_h5ad()` and indexing commands for seamless data exploration

## How it compares

| Tool | Pros | Cons |
|------|------|------|
| **Backed AnnData** | Full-featured | Slower access times |
| **BioNeMo** | Similar memory mapping | Heavier weight, unfamiliar API |
| **scplode** | Fast, memory-efficient, familiar API | Requires disk space for index |

## Installation

```bash
pip install scplode
```

## Quick Start

```python
import scplode as sp

# Your familiar workflow, now faster
adata = sp.read_h5ad('your_data.h5ad')
```

## Examples

Located in examples directory:
- `00_example`: Demonstrates use of scplode, and tests for equivalent results
- `01_benchmark`: Compares scplode and anndata for random and contiguous indexing
- `02_state_benchmark`: Compares scplode and anndata when using scplode with the Arc Institute State Model

## Requirements

- Python 3.8+
- AnnData
- Pandas
- Numpy

## License

MIT
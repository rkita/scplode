#!/usr/bin/env python3

import logging
from pathlib import Path
from typing import Union, List

import anndata as ad
import numpy as np
import pandas as pd
from tqdm.auto import tqdm

logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("[%(levelname)s] %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.INFO)


def read_h5ad(adata_path: Union[str, Path]) -> "Scplode":
    """
    Reads a `.h5ad` file and returns a `Scplode` object. If the index files
    (.obs, .var, .dat) do not exist, they are created.

    Args:
        adata_path (str | Path): Path to `.h5ad` file.

    Returns:
        Scplode: Wrapper object providing indexed access to AnnData.
    """
    adata_path = Path(adata_path)
    assert adata_path.exists(), f"File does not exist: {adata_path}"
    scdata = Scplode(adata_path)
    if not scdata.has_index():
        scdata.create_index()
    scdata.load_index()
    return scdata



class Scplode:
    def __init__(self, adata_path: Union[str, Path]) -> None:
        """
        Wrapper around AnnData for fast, chunked access using memory-mapped data.

        Args:
            adata_path (str | Path): Path to `.h5ad` file.
        """
        self.adata_path = Path(adata_path)

        # Other paths
        self.obs_path = self.adata_path.with_suffix(".obs")
        self.dat_path = self.adata_path.with_suffix(".dat")
        self.var_path = self.adata_path.with_suffix(".var")
        self.mtime_path = self.adata_path.with_suffix(".mtime")  # store file mtime


        # Variables updated later
        self.obs: Union[pd.DataFrame, None] = None
        self.var: Union[pd.DataFrame, None] = None
        self.n_obs: int = 0
        self.n_vars: int = 0
        self.adata_mtime: float = 0.0  # stored mtime

    def has_index(self) -> bool:
        """
        Check whether index files (.obs, .var, .dat) exist.

        Returns:
            bool: True if index files exist, False otherwise.
        """
        return self.obs_path.exists() and self.var_path.exists() and self.dat_path.exists()

    def load_index(self) -> bool:
        """
        Load precomputed index files into memory.

        Returns:
            bool: True if successful.
        """
        assert self.has_index(), "Index files are missing. Run create_index() first."
        self._check_mtime()  # verify file consistency
        logger.info("Loading index: obs")
        self.obs = pd.read_pickle(self.obs_path)

        logger.info("Loading index: var")
        self.var = pd.read_pickle(self.var_path)

        logger.info("Loading index: dat (implicitly)")
        self.n_vars = self.var.shape[0]
        self.n_obs = self.obs.shape[0]

        return True

    def _save_mtime(self) -> None:
        """Save current .h5ad modification time to a file."""
        self.adata_mtime = self.adata_path.stat().st_mtime
        self.mtime_path.write_text(str(self.adata_mtime))

    def _load_mtime(self) -> float:
        """Load stored modification time."""
        if self.mtime_path.exists():
            return float(self.mtime_path.read_text())
        return 0.0

    def _check_mtime(self) -> None:
        """Check that the current file mtime matches the stored mtime."""
        stored = self._load_mtime()
        current = self.adata_path.stat().st_mtime
        if stored != current:
            raise RuntimeError(
                f"The .h5ad file has changed since the index was created.\n"
                f"Stored mtime: {stored}, current mtime: {current}"
            )

    def create_index(self, chunk_size: int = 1000) -> None:
        """
        Create index files by processing the `.h5ad` file in chunks.

        Args:
            chunk_size (int, optional): Number of rows per chunk when writing to memory-mapped file. Default = 1000.
        """
        logger.info("Creating index")
        logger.info("Creating index: reading adata file")
        adata = ad.read_h5ad(self.adata_path, backed="r")
        n_obs, n_vars = adata.shape

        logger.info("Creating index: writing mmap dat file")
        with open(self.dat_path, "wb") as f:
            for start_idx in tqdm(range(0, n_obs, chunk_size)):
                end_idx = min(start_idx + chunk_size, n_obs)
                if type(adata.X).__name__ == "Dataset":
                    chunk = np.array(adata.X[start_idx:end_idx, :]).astype("float32")
                elif type(adata.X).__name__ == "_CSRDataset":
                    chunk = adata.X[start_idx:end_idx, :].to_memory().toarray().astype("float32")
                else:
                    chunk = adata.X[start_idx:end_idx, :].toarray().astype("float32")
                chunk.tofile(f)
                del chunk

        logger.info("Creating index: packing obs")
        adata.obs.to_pickle(self.obs_path)

        logger.info("Creating index: packing var")
        adata.var.to_pickle(self.var_path)

        self._save_mtime()

    def delete_index(self) -> None:
        """
        Delete the `.obs`, `.var`, and `.dat` files created by Scplode for a given `.h5ad`.
    
        Args:
            adata_path (str | Path): Path to the original `.h5ad` file.
    
        Notes:
            - Will not delete the original `.h5ad`.
            - Logs which files are deleted or missing.
        """
        for path in [self.obs_path, self.var_path, self.dat_path, self.mtime_path]:
            if path.exists():
                path.unlink()
                logger.info(f"Deleted: {path}")
            else:
                logger.warning(f"File not found, skipping: {path}")

    def get(self, indices: Union[List[str], pd.Index, np.ndarray]) -> np.ndarray:
        """
        Retrieve data for given indices as a NumPy array.

        Args:
            indices (list[str] | pd.Index | np.ndarray): List or array of observation IDs.

        Returns:
            np.ndarray: Array of shape (len(indices), n_vars).
        """
        row_positions = self.obs.index.get_indexer(indices)
        data = np.memmap(self.dat_path, dtype="float32", mode="r").reshape(self.n_obs, -1)
        return data[row_positions, :]

    def get_adata(self, indices: Union[List[str], pd.Index, np.ndarray]) -> ad.AnnData:
        """
        Retrieve a subset of the data as an AnnData object.

        Args:
            indices (list[str] | pd.Index | np.ndarray): Observation IDs.

        Returns:
            AnnData: Subset of the AnnData object.
        """
        data = self.get(indices)
        obs = self.obs.loc[indices].copy()
        obs.index = obs.index.astype(str)

        # Normalize categorical columns
        for col in obs.select_dtypes(["category"]).columns:
            obs[col] = obs[col].astype(str).astype("category")

        return ad.AnnData(X=data, obs=obs, var=self.var)

    def __getitem__(self, key: Union[int, slice, list, np.ndarray, pd.Index, tuple]) -> ad.AnnData:
        """
        Retrieve subset of the data using indexing.
    
        Supports:
            - scdata[row]
            - scdata[row_slice]
            - scdata[row_list_or_array]
            - scdata[row_labels]
            - scdata[row_slice, col_slice_or_names]
    
        Args:
            key: int, slice, list[int], np.ndarray[int], pd.Index, or tuple
                If tuple, interpreted as (rows, columns), columns can be:
                    - int / slice / list of ints
                    - list of column names (from self.var.index)
    
        Returns:
            AnnData: Subset of the dataset.
        """
        # Row-only indexing
        if not isinstance(key, tuple):
            row_key = key
            col_key = slice(None)  # all columns
        else:
            if len(key) != 2:
                raise IndexError("Only 2D indexing supported: (rows, columns)")
            row_key, col_key = key
    
        # Convert rows to observation labels
        if isinstance(row_key, int):
            row_indices = self.obs.index[row_key:row_key + 1]
        elif isinstance(row_key, slice):
            row_indices = self.obs.index[row_key]
        elif isinstance(row_key, (list, np.ndarray)) and np.issubdtype(np.array(row_key).dtype, np.integer):
            row_indices = self.obs.index[row_key]
        elif isinstance(row_key, (pd.Index, list, np.ndarray)):
            row_indices = row_key
        else:
            raise TypeError(f"Unsupported row index type: {type(row_key)}")
    
        # Retrieve AnnData for selected rows
        adata_subset = self.get_adata(row_indices)
    
        # Handle columns
        if isinstance(col_key, int):
            col_indices = [col_key]
        elif isinstance(col_key, slice):
            col_indices = range(self.var.shape[0])[col_key]
        elif isinstance(col_key, list) and all(isinstance(c, int) for c in col_key):
            col_indices = col_key
        elif isinstance(col_key, list) and all(isinstance(c, str) for c in col_key):
            # Convert column names to positions
            col_indices = [self.var.index.get_loc(c) for c in col_key]
        elif isinstance(col_key, np.ndarray):
            if np.issubdtype(col_key.dtype, np.integer):
                col_indices = col_key.tolist()
            else:
                raise TypeError("Unsupported numpy array type for columns")
        else:
            raise TypeError(f"Unsupported column index type: {type(col_key)}")
    
        # Slice columns of AnnData
        return adata_subset[:, col_indices]

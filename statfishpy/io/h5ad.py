#statfish/io/h5ad.py

""" NOT COMPLETED/VERIFIED """

import anndata as ad
import numpy as np


def read_h5ad(
    path,
    counts_layer=None,
    spatial_key="spatial",
    x_key="x",
    y_key="y",
    volume_key="cell_volume",
    batch_key=None,
    make_raw=True,
):
    """
    Read an h5ad file and standardize it for MERFISH analysis.

    Parameters
    ----------
    path : str
        Path to .h5ad file.
    counts_layer : str or None
        Layer containing raw counts. If None, adata.X is assumed to be raw counts.
    spatial_key : str
        Key in adata.obsm for spatial coordinates.
    x_key, y_key : str
        Column names in adata.obs if spatial coords are not in obsm.
    volume_key : str
        Column in adata.obs containing cell volume.
    batch_key : str or None
        Optional batch column in adata.obs.
    make_raw : bool
        Whether to store raw counts in adata.raw.

    Returns
    -------
    adata : AnnData
        Standardized AnnData object.
    """
    adata = ad.read_h5ad(path)

    _setup_counts(adata, counts_layer, make_raw)
    _setup_spatial(adata, spatial_key, x_key, y_key)
    _setup_volume(adata, volume_key)
    _validate_batch(adata, batch_key)

    _record_provenance(
        adata,
        source="h5ad",
        path=path,
        counts_layer=counts_layer,
        spatial_key=spatial_key,
        volume_key=volume_key,
        batch_key=batch_key,
    )

    return adata


# ------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------

def _setup_counts(adata, counts_layer, make_raw):
    """
    Ensure raw MERFISH counts live in adata.layers['counts'].
    """
    if counts_layer is None:
        # Assume adata.X contains raw counts
        if "counts" not in adata.layers:
            adata.layers["counts"] = adata.X.copy()
    else:
        if counts_layer not in adata.layers:
            raise KeyError(
                f"Layer '{counts_layer}' not found in adata.layers. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        adata.layers["counts"] = adata.layers[counts_layer].copy()

    if make_raw and adata.raw is None:
        adata.raw = adata.copy()


def _setup_spatial(adata, spatial_key, x_key, y_key):
    """
    Ensure spatial coordinates are stored in adata.obsm['spatial'].
    """
    if spatial_key in adata.obsm:
        coords = adata.obsm[spatial_key]
    else:
        if x_key not in adata.obs or y_key not in adata.obs:
            raise KeyError(
                f"Spatial coordinates not found. Expected either "
                f"adata.obsm['{spatial_key}'] or obs columns "
                f"'{x_key}' and '{y_key}'."
            )
        coords = adata.obs[[x_key, y_key]].to_numpy()

    coords = np.asarray(coords)

    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(
            "Spatial coordinates must be a 2D array with at least two columns."
        )

    adata.obsm["spatial"] = coords


def _setup_volume(adata, volume_key):
    """
    Ensure cell volume exists and is numeric.
    """
    if volume_key not in adata.obs:
        raise KeyError(
            f"Cell volume '{volume_key}' not found in adata.obs. "
            "Volume normalization requires this field."
        )

    if not np.issubdtype(adata.obs[volume_key].dtype, np.number):
        raise TypeError(
            f"Cell volume '{volume_key}' must be numeric, "
            f"found dtype {adata.obs[volume_key].dtype}."
        )


def _validate_batch(adata, batch_key):
    """
    Validate batch column if provided.
    """
    if batch_key is None:
        return

    if batch_key not in adata.obs:
        raise KeyError(
            f"Batch key '{batch_key}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )


def _record_provenance(adata, **kwargs):
    """
    Store I/O provenance in adata.uns.
    """
    adata.uns.setdefault("merfishpy_io", {})
    adata.uns["merfishpy_io"].update(kwargs)

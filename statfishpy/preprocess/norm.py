import numpy as np
import scipy.sparse as sp

def n_by_volume(
    #normalize counts per cell volume
    adata,
    counts_layer="counts",
    volume_key="volume",
    output_layer="volume_norm",
    scale_factor=1.0,
        ):
    """
    Normalize MERFISH counts by cell volume.

    normalized_ij = counts_ij / volume_i * scale_factor

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    counts_layer : str
        Layer containing raw counts.
    volume_key : str
        Column in adata.obs containing cell volumes.
    output_layer : str
        Name of layer to store normalized counts.
    scale_factor : float
        Optional scaling factor (e.g. median volume or 1e3).
    """

    if counts_layer not in adata.layers:
        raise KeyError(f"Layer '{counts_layer}' not found in adata.layers.")

    if volume_key not in adata.obs:
        raise KeyError(f"'{volume_key}' not found in adata.obs.")

    X = adata.layers[counts_layer]
    volumes = adata.obs[volume_key].to_numpy()

    if np.any(volumes <= 0):
        raise ValueError("Cell volumes must be strictly positive.")

    inv_vol = scale_factor / volumes

    if sp.issparse(X):
        X_norm = X.multiply(inv_vol[:, None])
    else:
        X_norm = X * inv_vol[:, None]

    adata.layers[output_layer] = X_norm

    return adata
    

def n_by_area(
    adata,
    counts_layer="counts",
    area_key="area",
    output_layer="area_norm",
    scale_factor=1.0,
):
    """
    Normalize MERFISH counts by cell area.

    normalized_ij = counts_ij / area_i * scale_factor
    """

    if counts_layer not in adata.layers:
        raise KeyError(f"Layer '{counts_layer}' not found in adata.layers.")

    if area_key not in adata.obs:
        raise KeyError(f"'{area_key}' not found in adata.obs.")

    X = adata.layers[counts_layer]
    areas = adata.obs[area_key].to_numpy()

    if np.any(areas <= 0):
        raise ValueError("Cell areas must be strictly positive.")

    inv_area = scale_factor / areas

    if sp.issparse(X):
        X_norm = X.multiply(inv_area[:, None])
    else:
        X_norm = X * inv_area[:, None]

    adata.layers[output_layer] = X_norm
    return adata


def log1p(
    adata,
    counts_layer="counts",
    output_layer="log1p",
    target_sum=1e4,
):
    """
    Standard log normalization:
        1) library size normalize
        2) log1p transform
    """

    if counts_layer not in adata.layers:
        raise KeyError(f"Layer '{counts_layer}' not found in adata.layers.")

    X = adata.layers[counts_layer]

    if sp.issparse(X):
        n_counts = np.asarray(X.sum(axis=1)).ravel()
        scale = target_sum / n_counts
        X_norm = X.multiply(scale[:, None])
        X_log = X_norm.log1p()
    else:
        n_counts = X.sum(axis=1)
        scale = target_sum / n_counts
        X_norm = X * scale[:, None]
        X_log = np.log1p(X_norm)

    adata.layers[output_layer] = X_log
    return adata




    
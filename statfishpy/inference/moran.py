import numpy as np
import scipy.sparse as sp

def morans_i(
    adata,
    input_layer="gp_smooth",
    graph_key="spatial_connectivities",
    output_key="morans_i",
):
    """
    Compute Moran's I for all genes.

    Stores results in:
        adata.var[output_key]
    """

    if input_layer not in adata.layers:
        raise KeyError(f"{input_layer} not found in adata.layers")

    if graph_key not in adata.obsp:
        raise KeyError(f"{graph_key} not found in adata.obsp")

    X = adata.layers[input_layer]  # (cells x genes)
    W = adata.obsp[graph_key]      # (cells x cells)

    if sp.issparse(X):
        X = X.toarray()

    # Center gene expression
    mean = X.mean(axis=0)
    X_centered = X - mean

    # denominator: variance per gene
    denom = (X_centered ** 2).sum(axis=0)

    # numerator: spatial covariance
    WX = W @ X_centered
    num = (X_centered * WX).sum(axis=0)

    N = X.shape[0]
    W_sum = W.sum()

    I = (N / W_sum) * (num / denom)

    adata.var[output_key] = I

    return adata

def get_spatially_correlated_genes(
    adata,
    moran_key="morans_i",
    threshold=0.0,
    sort=True,
):
    """
    Return genes with positive spatial autocorrelation.
    """

    if moran_key not in adata.var:
        raise KeyError(f"{moran_key} not found in adata.var")

    vals = adata.var[moran_key].values

    mask = vals > threshold
    genes = adata.var_names[mask]

    if sort:
        order = np.argsort(vals[mask])[::-1]
        genes = genes[order]

    return list(genes)
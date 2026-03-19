import numpy as np
from sklearn.decomposition import PCA
import scipy.sparse as sp

def pca_embed(
    adata,
    input_layer="gp_smooth",
    n_components=50,
    output_key="X_pca",
):

    if input_layer not in adata.layers:
        raise KeyError(f"{input_layer} not found in adata.layers")

    X = adata.layers[input_layer]

    if sp.issparse(X):
        X = X.toarray()

    model = PCA(n_components=n_components)

    embedding = model.fit_transform(X)

    adata.obsm[output_key] = embedding
    adata.uns["pca_variance_ratio"] = model.explained_variance_ratio_

    return adata

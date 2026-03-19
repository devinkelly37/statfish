import numpy as np
import scipy.sparse as sp
from sklearn.neighbors import NearestNeighbors

def build_spatial_graph(
    adata,
    x_key="center_x",
    y_key="center_y",
    k=30
):

    coords = np.column_stack([
        adata.obs[x_key].values,
        adata.obs[y_key].values
    ])

    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(coords)

    dist, idx = nn.kneighbors(coords)

    n = coords.shape[0]

    rows = np.repeat(np.arange(n), k)
    cols = idx.flatten()
    vals = dist.flatten()

    D = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))

    adata.obsp["spatial_distances"] = D

    return adata

def build_spatial_kernel(
    adata,
    bandwidth=50
):

    D = adata.obsp["spatial_distances"].copy()

    # Gaussian kernel
    D.data = np.exp(-(D.data**2) / (2 * bandwidth**2))

    # normalize rows
    row_sum = np.array(D.sum(axis=1)).flatten()
    row_sum[row_sum == 0] = 1

    W = sp.diags(1 / row_sum) @ D

    adata.obsp["spatial_kernel"] = W

    return adata
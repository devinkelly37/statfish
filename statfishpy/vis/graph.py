import matplotlib.pyplot as plt
import numpy as np

def plot_embedding(
    adata,
    embedding="X_umap",
    color=None,
    ax=None,
    s=8,
    alpha=0.8,
    cmap="tab20",
    title=None,
):

    if embedding not in adata.obsm:
        raise KeyError(f"{embedding} not found in adata.obsm")

    coords = adata.obsm[embedding]

    if ax is None:
        fig, ax = plt.subplots(figsize=(5,5))
    else:
        fig = ax.figure

    x = coords[:,0]
    y = coords[:,1]

    if color is None:

        ax.scatter(x, y, s=s, alpha=alpha)

    else:

        if color in adata.obs:
            c = adata.obs[color].values
        else:
            raise KeyError(f"{color} not found in adata.obs")

        sc = ax.scatter(
            x,
            y,
            c=c,
            s=s,
            alpha=alpha,
            cmap=cmap
        )

        fig.colorbar(sc, ax=ax)

    ax.set_xlabel(f"{embedding} 1")
    ax.set_ylabel(f"{embedding} 2")

    if title is None:
        title = embedding

    ax.set_title(title)

    return fig, ax

def plot_spatial_clusters(
    adata,
    cluster_key="cluster",
    x_key="center_x",
    y_key="center_y",
    ax=None,
    s=10,
    alpha=0.9,
    cmap="tab20",
    title="Spatial clusters",
):

    if cluster_key not in adata.obs:
        raise KeyError(f"{cluster_key} not found in adata.obs")

    x = adata.obs[x_key].values
    y = adata.obs[y_key].values
    clusters = adata.obs[cluster_key].values

    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))
    else:
        fig = ax.figure

    sc = ax.scatter(
        x,
        y,
        c=clusters,
        s=s,
        cmap=cmap,
        alpha=alpha
    )

    ax.set_xlabel(x_key)
    ax.set_ylabel(y_key)

    ax.set_title(title)

    fig.colorbar(sc, ax=ax)

    ax.set_aspect("equal")

    return fig, ax
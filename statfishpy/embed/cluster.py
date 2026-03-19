import umap
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering

def umap_embed(
    adata,
    input_key="X_pca",
    n_neighbors=15,
    min_dist=0.3,
    output_key="X_umap",
):

    if input_key not in adata.obsm:
        raise KeyError(f"{input_key} not found in adata.obsm")

    X = adata.obsm[input_key]

    model = umap.UMAP(
        n_neighbors=n_neighbors,
        min_dist=min_dist,
    )

    embedding = model.fit_transform(X)

    adata.obsm[output_key] = embedding

    return adata


def tsne_embed(
    adata,
    input_key="X_pca",
    perplexity=30,
    output_key="X_tsne",
):

    if input_key not in adata.obsm:
        raise KeyError(f"{input_key} not found in adata.obsm")

    X = adata.obsm[input_key]

    model = TSNE(
        perplexity=perplexity,
        n_components=2,
    )

    embedding = model.fit_transform(X)

    adata.obsm[output_key] = embedding

    return adata


def spatial_cluster(
    adata,
    input_key="X_pca",
    n_clusters=10,
    output_key="cluster",
):

    if input_key not in adata.obsm:
        raise KeyError(f"{input_key} not found")

    X = adata.obsm[input_key]

    model = SpectralClustering(
        n_clusters=n_clusters,
        affinity="nearest_neighbors",
    )

    labels = model.fit_predict(X)

    adata.obs[output_key] = labels

    return adata
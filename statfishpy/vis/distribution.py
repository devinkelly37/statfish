import numpy as np
import scipy.sparse as sp
import matplotlib.pyplot as plt


def plot_mean_variance(
    X,
    log=False,
    ax=None,
    s=6,
    alpha=0.5,
    title="Mean vs Variance (per gene)",
):
    """
    Plot gene-wise mean vs variance.

    Parameters
    ----------
    X : ndarray or sparse matrix
        Cell x gene matrix.
    log : bool
        Whether to log-scale both axes.
    ax : matplotlib Axes or None
        Axis to plot on.
    s : int
        Marker size.
    alpha : float
        Marker transparency.
    title : str
        Plot title.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """

    if sp.issparse(X):
        mean = np.asarray(X.mean(axis=0)).ravel()
        var = np.asarray(X.power(2).mean(axis=0)).ravel() - mean**2
    else:
        mean = X.mean(axis=0)
        var = X.var(axis=0)

    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 5))
    else:
        fig = ax.figure

    ax.scatter(mean, var, s=s, alpha=alpha)

    ax.set_xlabel("Mean expression")
    ax.set_ylabel("Variance")

    if log:
        ax.set_xscale("log")
        ax.set_yscale("log")

    ax.set_title(title)

    return fig, ax


def plot_gene_histogram(
    X,
    gene_index,
    bins=50,
    log=False,
    ax=None,
    alpha=0.7,
    title=None,
):
    """
    Plot histogram of expression values across cells for a single gene.

    Parameters
    ----------
    X : ndarray or sparse matrix
        Cell x gene matrix.
    gene_index : int
        Index of the gene to plot.
    bins : int
        Number of histogram bins.
    log : bool
        Log scale on x-axis.
    ax : matplotlib Axes or None
        Axis to plot on.
    alpha : float
        Histogram transparency.
    title : str or None
        Plot title.

    Returns
    -------
    fig, ax : matplotlib Figure and Axes
    """

    # extract gene vector
    if sp.issparse(X):
        gene_vals = X[:, gene_index].toarray().ravel()
    else:
        gene_vals = X[:, gene_index]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6,5))
    else:
        fig = ax.figure

    ax.hist(gene_vals, bins=bins, alpha=alpha)

    ax.set_xlabel("Expression value across cells")
    ax.set_ylabel("Number of cells")

    if log:
        ax.set_xscale("log")

    if title is None:
        title = f"Gene {gene_index} expression histogram"

    ax.set_title(title)

    return fig, ax
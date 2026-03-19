import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.optim as optim

def pearson_resid(
    adata,
    counts_layer="counts",
    output_layer="pearson_residuals",
    theta=100.0,
    clip=10.0,
):

    if counts_layer not in adata.layers:
        raise KeyError(f"Layer '{counts_layer}' not found.")

    X = adata.layers[counts_layer]

    # convert to CSR for efficient row operations
    if sp.issparse(X):
        X = X.tocsr()

    # library size per cell
    n_counts = np.asarray(X.sum(axis=1)).ravel()

    if np.any(n_counts == 0):
        raise ValueError("Cells with zero counts detected.")

    # gene frequencies
    gene_totals = np.asarray(X.sum(axis=0)).ravel()
    p = gene_totals / gene_totals.sum()

    # expected counts μ = n_c * p_g
    mu = n_counts[:, None] * p[None, :]

    # variance under NB model
    var = mu + (mu ** 2) / theta

    # compute residuals
    if sp.issparse(X):
        X = X.toarray()

    resid = (X - mu) / np.sqrt(var + 1e-8)

    # clip extreme residuals
    resid = np.clip(resid, -clip, clip)

    adata.layers[output_layer] = resid.astype(np.float32)

    return adata

def neg_bin_glm(
    adata,
    counts_layer="counts",
    batch_key=None,
    covariates=None,
    offset_key="volume",
    output_layer="nb_residuals",
    lr=1e-2,
    max_epochs=100,
    clip=10.0,
    verbose=True,
):
    """
    This performs:
        - Negative binomial regression
        - Library size normalization
        - Optional batch correction
        - Optional covariate regression
        - Variance stabilization

    Parameters
    ----------
    adata : AnnData
    counts_layer : str
        Raw counts layer.
    batch_key : str or None
        Column in adata.obs defining batches.
    covariates : list[str] or None
        Continuous covariates (volume, area, etc).
    offset_key : str or None
    output_layer : str
        Layer to store residuals.
    n_latent : int
        Latent dimension (0 = pure GLM).
    max_epochs : int
        Training epochs.
    """

    X = adata.layers[counts_layer]
    if sp.issparse(X):
        X = X.toarray()
    X = torch.tensor(X, dtype=torch.float32)  # might still be huge

    n_cells, n_genes = X.shape

    # Build design matrix (drop first batch level)
    design = [torch.ones(n_cells, 1)]
    if covariates:
        for cov in covariates:
            if cov == offset_key:
                continue
            col = torch.tensor(adata.obs[cov].values, dtype=torch.float32).view(-1, 1)
            design.append(col)
    if batch_key:
        # Use pandas get_dummies with drop_first=True
        import pandas as pd
        dummies = pd.get_dummies(adata.obs[batch_key], drop_first=True)
        design.append(torch.tensor(dummies.values, dtype=torch.float32))
    D = torch.cat(design, dim=1)  # (cells, p)
    p = D.shape[1]

    # Offset
    if offset_key:
        offset = torch.log(torch.tensor(adata.obs[offset_key].values, dtype=torch.float32))
        offset = offset.view(-1, 1)
    else:
        offset = torch.zeros(n_cells, 1)

    # Gene-wise fitting (parallel loop recommended)
    beta_all = torch.zeros(n_genes, p)
    log_theta_all = torch.zeros(n_genes)

    for g in range(n_genes):
        y = X[:, g]
        beta = torch.zeros(p, requires_grad=True)
        log_theta = torch.tensor(0.0, requires_grad=True)
        optimizer = optim.Adam([beta, log_theta], lr=lr)

        for epoch in range(max_epochs):
            optimizer.zero_grad()
            eta = D @ beta + offset.squeeze()
            mu = torch.exp(eta.clamp(max=15))  # prevent overflow
            theta = torch.exp(log_theta) + 1e-4
            nll = -(torch.lgamma(y + theta) - torch.lgamma(theta) -
                    torch.lgamma(y + 1) + theta * torch.log(theta / (theta + mu)) +
                    y * torch.log(mu / (theta + mu)))
            loss = nll.mean()
            loss.backward()
            optimizer.step()

        beta_all[g] = beta.detach()
        log_theta_all[g] = log_theta.detach()

    # Compute residuals (using fitted parameters)
    with torch.no_grad():
        eta = D @ beta_all.T + offset
        mu = torch.exp(eta.clamp(max=15))
        theta = torch.exp(log_theta_all) + 1e-4
        var = mu + mu**2 / theta
        resid = (X - mu) / torch.sqrt(var + 1e-8)

    resid = resid.clamp(-clip, clip).numpy()
    adata.layers[output_layer] = resid
    return adata


def spatial_smooth(
    adata,
    input_layer="pearson_residuals",
    output_layer="gp_smooth"
):

    if "spatial_kernel" not in adata.obsp:
        raise KeyError("Run build_spatial_kernel first")

    X = adata.layers[input_layer]
    W = adata.obsp["spatial_kernel"]

    smoothed = W @ X

    adata.layers[output_layer] = smoothed

    return adata
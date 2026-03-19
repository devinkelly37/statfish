import statfishpy as sf
import matplotlib.pyplot as plt

# --------------------------------------------------
# Load data
# --------------------------------------------------

adata = sf.io.read_h5ad(
    "raw_counts_mouse4_sagittal.h5ad",
    x_key="center_x",
    y_key="center_y",
    volume_key="volume",
)

print("Dataset shape:", adata.shape)

# --------------------------------------------------
# Preprocessing
# --------------------------------------------------

sf.preprocess.pearson_resid(adata)

# --------------------------------------------------
# Build spatial graph + kernel
# --------------------------------------------------

sf.spatial.build_spatial_graph(
    adata,
    x_key="center_x",
    y_key="center_y",
    k=30
)

sf.spatial.build_spatial_kernel(
    adata,
    bandwidth=50
)

# --------------------------------------------------
# Spatial smoothing
# --------------------------------------------------

sf.preprocess.spatial_smooth(
    adata,
    input_layer="pearson_residuals",
    output_layer="gp_smooth"
)

# --------------------------------------------------
# Compute Moran's I
# --------------------------------------------------

sf.inference.morans_i(
    adata,
    input_layer="gp_smooth",
    graph_key="spatial_kernel",
    output_key="morans_i"
)

# --------------------------------------------------
# Print correlated genes
# --------------------------------------------------

# Get all genes with positive spatial autocorrelation
genes = sf.inference.get_spatially_correlated_genes(
    adata,
    threshold=0.9
)

print(f"\nNumber of spatially correlated genes: {len(genes)}")

# Print top 20 genes
print("\nTop 20 spatially correlated genes:")
for g in genes[:20]:
    print(g)


from pathlib import Path

import scanpy as sc


INPUT_FILE = Path("output/processed_dataset.h5ad")
OUTPUT_FILE = Path("output/preliminary_marker_genes.csv")

print("Loading processed AnnData file...")
adata = sc.read_h5ad(INPUT_FILE)

# Use existing biological cell-type labels when available.
if "cell_type" in adata.obs.columns:
    GROUP_COLUMN = "cell_type"
elif "Annotation" in adata.obs.columns:
    GROUP_COLUMN = "Annotation"
else:
    raise ValueError(
        "No cell-type or annotation column was found."
    )

print(f"Using labels from: {GROUP_COLUMN}")

# Remove cells without a label.
adata = adata[
    adata.obs[GROUP_COLUMN].notna()
].copy()

# Use highly variable genes when available to reduce memory and runtime.
if "vst.variable" in adata.var.columns:
    variable_genes = adata.var["vst.variable"].fillna(False).astype(bool)

    if variable_genes.sum() >= 50:
        adata = adata[:, variable_genes].copy()
        print(f"Using {adata.n_vars} variable genes.")

# Make sure labels are categorical.
adata.obs[GROUP_COLUMN] = adata.obs[GROUP_COLUMN].astype("category")

print(
    f"Running marker analysis for "
    f"{adata.obs[GROUP_COLUMN].nunique()} groups..."
)

# Use the unscaled expression matrix.
sc.tl.rank_genes_groups(
    adata,
    groupby=GROUP_COLUMN,
    method="wilcoxon",
    use_raw=False,
    n_genes=20,
)

markers = sc.get.rank_genes_groups_df(
    adata,
    group=None,
)

markers.to_csv(OUTPUT_FILE, index=False)

print(f"Saved preliminary markers to: {OUTPUT_FILE}")
print(f"Rows written: {len(markers):,}")
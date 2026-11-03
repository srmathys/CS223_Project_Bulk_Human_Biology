# The linkage file that runs as part of the main
# It will take a saved object in output
import anndata as ad
import numpy as np
import pandas as pd

from rpy2.robjects import numpy2ri
FILE_PATH = './data/human_bone_marrow_hematopoieticd.h5ad'

def scBubbleTreeToNetworkX(file, pca_val):
# From the processed output file that is a .h5ad pull the requisite values
# Feed these values into scBubbleTree

    adata = ad.read_h5ad(file)
# Should I do backend loading if it's rather large?
    ad_pca_matrix = adata.obsm['X_pca'][:,pca_val]
    ad_metadata = adata.obs['cell_type']




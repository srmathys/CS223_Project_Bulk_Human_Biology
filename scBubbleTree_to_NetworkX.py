# The linkage file that runs as part of the main
# It will take a saved object in output
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from rpy2.robjects as ro


from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
FILE_PATH = './data/human_bone_marrow_hematopoieticd.h5ad'



def scBubbleTreeToNetworkX(file, pca_val):
# From the processed output file that is a .h5ad pull the requisite values
# Feed these values into scBubbleTree

    adata = ad.read_h5ad(file)
# Should I do backend loading if it's rather large? probably
    pca_matrix = adata.obsm['X_pca'][:,pca_val]
    metadata = adata.obs['cell_type']
    marker_genes = sc.get.rank_genes_groups_df( adata, group=None)








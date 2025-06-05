"""Class for defining methods to package pipeline outputs into AnnData objects
"""

import os
import pandas as pd
import anndata
import scanpy as sc
import seaborn as sns
import numpy as np

## Attempt to import celltypist, and prompt installation if not found
import celltypist
from celltypist import models

from .pipeline import ENACT


class CellTypistPipeline(ENACT):
    """Class for running CellAssign algorithm"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run_cell_typist(self):
        """Runs CellTypist"""
        bin_assign_results = self.merge_files_sparse(self.bin_assign_dir)
        cell_lookup_df = self.merge_files(self.cell_ix_lookup_dir, save=False)

        spatial_cols = ["cell_x", "cell_y"]
        stat_columns = ["num_shared_bins", "num_unique_bins", "num_transcripts"]
        cell_lookup_df.loc[:, "id"] = cell_lookup_df["id"].astype(str)
        cell_lookup_df = cell_lookup_df.set_index("id")
        cell_lookup_df["num_transcripts"] = cell_lookup_df["num_transcripts"].fillna(0)

        bin_assign_result_sparse, gene_columns = bin_assign_results
        adata = anndata.AnnData(X=bin_assign_result_sparse, obs=cell_lookup_df.copy())
        adata.var_names = gene_columns

        adata.obsm["spatial"] = cell_lookup_df[spatial_cols].astype(int)
        adata.obsm["stats"] = cell_lookup_df[stat_columns].astype(int)

        lib_size = adata.X.sum(1)
        adata.obs["size_factor"] = lib_size / np.mean(lib_size)
        adata.obs["lib_size"] = lib_size

        #  normalize adata to the log1p normalised format (to 10,000 counts per cell)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Check if the model file exists on disk
        if os.path.isfile(self.cell_typist_model):
            pass
            # The model file already exists on disk; no further action is required.
        # Append '.pkl' if not already present and recheck
        elif not self.cell_typist_model.endswith(".pkl"):
            self.cell_typist_model = self.cell_typist_model + ".pkl"
            if not os.path.isfile(self.cell_typist_model):
                models.download_models(model=self.cell_typist_model)
        predictions = celltypist.annotate(adata, model=self.cell_typist_model)
        adata = predictions.to_adata(
            insert_labels=True, insert_conf=True, insert_prob=True
        )

        adata.obs.rename(columns={"predicted_labels": "cell_type"}, inplace=True)
        adata.obs[adata.obsm["spatial"].columns] = adata.obsm["spatial"]
        adata.obs[adata.obsm["stats"].columns] = adata.obsm["stats"]
        adata.obs["chunk_name"] = cell_lookup_df["chunk_name"]
        results_df = adata.obs.drop(columns=adata.obs["cell_type"].unique().tolist())
        results_df.to_csv(
            os.path.join(self.cellannotation_results_dir, "merged_results.csv")
        )


if __name__ == "__main__":
    # Creating CellAssignPipeline object
    cell_typist = CellTypistPipeline(configs_path="config/configs.yaml")
    cell_typist.run_cell_typist()

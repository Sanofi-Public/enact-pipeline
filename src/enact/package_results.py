"""Class for defining methods to package pipeline outputs into AnnData objects
"""

import os
import yaml
import json
import shutil
import anndata
import pandas as pd
from PIL import Image
import numpy as np
from scipy.sparse import csr_matrix

# import squidpy as sq

from .pipeline import ENACT


class PackageResults(ENACT):
    """Class for packaging ENACT pipeline outputs"""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.files_to_ignore = [
            "merged_results.csv",
            "merged_results_old.csv",
            "cells_adata.h5",
            ".ipynb_checkpoints",
        ]

    def merge_cellassign_output_files(self):
        """Merges the CellAssign results with gene counts

        Returns:
            _type_: _description_
        """
        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(self.bin_assign_dir)
        cell_by_gene_list = []
        for chunk_name in chunk_list:
            if chunk_name in self.files_to_ignore:
                continue
            index_lookup = pd.read_csv(
                os.path.join(self.cell_ix_lookup_dir, chunk_name)
            )
            trancript_counts = pd.read_csv(
                os.path.join(self.bin_assign_dir, chunk_name)
            ).drop(columns=["Unnamed: 0"])
            cell_by_gene_chunk = pd.concat(
                [index_lookup["id"], trancript_counts], axis=1
            )
            cell_by_gene_list.append(cell_by_gene_chunk)
        cell_by_gene_df = pd.concat(cell_by_gene_list, axis=0)
        return cell_by_gene_df

    def merge_sargent_output_files(self):
        """Merges the Sargent chunk results into a single results file

        Returns:
            _type_: _description_
        """
        os.makedirs(self.cellannotation_results_dir, exist_ok=True)
        # Merge the sargent_results_chunks data and gene_to_cell_assignment_chunks_ix_lookup
        chunks = os.listdir(self.cellannotation_results_dir)
        sargent_results_list = []
        cell_by_gene_list = []
        for chunk_name in chunks:
            if chunk_name in self.files_to_ignore:
                continue
            cell_labels = pd.read_csv(
                os.path.join(self.cellannotation_results_dir, chunk_name)
            )
            index_lookup = pd.read_csv(
                os.path.join(self.cell_ix_lookup_dir, chunk_name)
            )
            trancript_counts = pd.read_csv(
                os.path.join(self.bin_assign_dir, chunk_name)
            ).drop(columns=["Unnamed: 0"])

            sargent_result_chunk = pd.concat([index_lookup, cell_labels["x"]], axis=1)
            cell_by_gene_chunk = pd.concat(
                [index_lookup["id"], trancript_counts], axis=1
            )
            sargent_result_chunk.drop("Unnamed: 0", axis=1, inplace=True)
            sargent_results_list.append(sargent_result_chunk)
            cell_by_gene_list.append(cell_by_gene_chunk)
        sargent_results_df = pd.concat(sargent_results_list, axis=0)
        sargent_results_df = sargent_results_df.rename(columns={"x": "cell_type"})
        cell_by_gene_df = pd.concat(cell_by_gene_list, axis=0)
        sargent_results_df.to_csv(
            os.path.join(self.cellannotation_results_dir, "merged_results.csv"), index=False
        )
        return sargent_results_df, cell_by_gene_df

    def df_to_adata(self, results_df, cell_by_gene_df):
        """Converts pd.DataFrame object with pipeline results to AnnData

        Args:
            results_df (_type_): _description_

        Returns:
            anndata.AnnData: Anndata with pipeline outputs
        """
        file_columns = results_df.columns
        spatial_cols = ["cell_x", "cell_y"]
        stat_columns = ["num_shared_bins", "num_unique_bins", "num_transcripts"]
        results_df.loc[:, "id"] = results_df["id"].astype(str)
        results_df = results_df.set_index("id")
        results_df["num_transcripts"] = results_df["num_transcripts"].fillna(0)
        results_df["cell_type"] = results_df["cell_type"].str.lower()
        adata = anndata.AnnData(cell_by_gene_df.set_index("id"))
        adata.obs = adata.obs.merge(results_df, on="id").drop_duplicates(keep='first')

        adata.obsm["spatial"] = adata.obs[spatial_cols].astype(int)
        adata.obsm["stats"] = adata.obs[stat_columns].astype(int)
        
        # This column is the output of cell type inference pipeline
        adata.obs["cell_type"] = adata.obs[["cell_type"]].astype("category")
        adata.obs["patch_id"] = adata.obs[["chunk_name"]]
        adata.obs = adata.obs[["cell_type", "patch_id"]]

        # Converting the Anndata cell transcript counts to sparse format for more efficient storage
        adata.X = csr_matrix(adata.X).astype(np.float32)
        return adata

    def create_tmap_file(self):
        """Creates a tmap file for the sample being run on ENACT
        """
        # The following three files need to be in the same directory:
        # cells_adata.h5, wsi file, experiment_tmap.tmap
        tmap_template_path = "./templates/tmap_template.tmap"
        with open(tmap_template_path, "r") as stream:
            tmap_template = yaml.safe_load(stream)
            tmap_template["filename"] = self.configs["analysis_name"]
            bin_to_cell_method = self.configs["params"]["bin_to_cell_method"]
            cell_annotation_method = self.configs["params"]["cell_annotation_method"]
            wsi_src_path = self.configs["paths"]["wsi_path"]
            wsi_fname = "wsi.tif"
            run_name = f"{bin_to_cell_method}|{cell_annotation_method}"
            tmap_template["markerFiles"][0]["title"] = f"ENACT Results: {run_name.replace('|', ' | ')}"
            tmap_template["markerFiles"][0]["expectedHeader"].update(
                {
                    "X": "/obsm/spatial/cell_x",
                    "Y": "/obsm/spatial/cell_y",
                    "gb_col": "/obs/cell_type/",
                }
            )
            tmap_template["layers"][0].update(
                    {"name": wsi_fname, "tileSource": f"{wsi_fname}.dzi"}
                )
            tmap_template["markerFiles"][0]["path"] = f"{run_name}_cells_adata.h5"

            # save tmap file at a separate directory "tmap"
            tmap_output_dir = os.path.join(self.cache_dir, "tmap")
            os.makedirs(tmap_output_dir, exist_ok=True)
            tmap_file_path = os.path.join(tmap_output_dir, f"{run_name}_tmap.tmap")
            with open(tmap_file_path, "w") as outfile:
                outfile.write(json.dumps(tmap_template, indent=4))

            # Copy the anndata file to the "tmap" directory
            adata_src_path = os.path.join(
                self.cellannotation_results_dir, "cells_adata.h5"
            )
            adata_dst_path = os.path.join(tmap_output_dir, f"{run_name}_cells_adata.h5")
            shutil.copy(adata_src_path, adata_dst_path)

            # Copy the cells_layer.png file to the "tmap" directory
            layer_src_path = os.path.join(
                self.cache_dir, "cells_layer.png"
            )
            layer_dst_path = os.path.join(tmap_output_dir, "cells_layer.png")
            if os.path.exists(layer_src_path):
                shutil.copy(layer_src_path, layer_dst_path)

            # Saving a cropped version (lite version) of the image file to the "tmap" directory
            wsi_dst_path = os.path.join(tmap_output_dir, wsi_fname)
            cropped_image, _ = self.load_image()
            cropped_image = Image.fromarray(cropped_image)
            cropped_image.save(wsi_dst_path)

            message = f"""
            Sample ready to visualize on TissUUmaps. To install TissUUmaps, follow the instructions at:\n
            https://tissuumaps.github.io/TissUUmaps-docs/docs/intro/installation.html#. 
            
            To view the the sample, follow the instructions at:\n
            https://tissuumaps.github.io/TissUUmaps-docs/docs/starting/projects.html#loading-projects
            
            TissUUmaps project file is located here:\n
            {tmap_file_path}
            """
            print (message)

    # def run_neighborhood_enrichment(self, adata):
    #     """Sample function to run Squidpy operations on AnnData object

    #     Args:
    #         adata (_type_): _description_

    #     Returns:
    #         _type_: _description_
    #     """
    #     sq.gr.spatial_neighbors(adata)
    #     sq.gr.nhood_enrichment(adata, cluster_key="cell_type")
    #     return adata

    def save_adata(self, adata):
        """Save the anndata object to disk

        Args:
            adata (_type_): _description_
        """
        adata.write(
            os.path.join(self.cellannotation_results_dir, "cells_adata.h5"),
            compression="gzip",
        )


if __name__ == "__main__":
    # Creating ENACT object
    so_hd = PackageResults(configs_path="config/configs.yaml")
    results_df, cell_by_gene_df = so_hd.merge_sargent_output_files()
    adata = so_hd.df_to_adata(results_df, cell_by_gene_df)
    # adata = so_hd.run_neighborhood_enrichment(adata) # Example integration with SquiPy
    so_hd.save_adata(adata)
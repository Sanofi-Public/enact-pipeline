"""Class for defining methods for VisiumHD pipeline
"""

import os
from csbdeep.utils import normalize
import geopandas as gpd
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
import scanpy as sc
from scipy import sparse
from scipy.spatial import Voronoi
import shapely
from shapely.geometry import Polygon, Point
from shapely import wkt
from stardist.models import StarDist2D
import tifffile as tifi
from tqdm import tqdm
import yaml
import anndata
import scvi
import seaborn as sns
from scvi.external import CellAssign
import logging
import ssl
import argparse
import ast

Image.MAX_IMAGE_PIXELS = None
from .utils.logging import get_logger
from .assignment_methods.naive import naive_assignment
from .assignment_methods.weight_by_area import weight_by_area_assignment
from .assignment_methods.weight_by_gene import (
    weight_by_gene_assignment,
    weight_by_cluster_assignment,
)


class ENACT:
    """Class for methods for the ENACT pipeline"""

    def __init__(
        self,
        cache_dir="",
        wsi_path="",
        visiumhd_h5_path="",
        tissue_positions_path="",
        analysis_name="enact_demo",
        seg_method="stardist",
        image_type="he",
        nucleus_expansion=True,
        expand_by_nbins=2,
        patch_size=4000,
        use_hvg=True,
        n_hvg=1000,
        destripe_norm=False,
        n_clusters=4,
        n_pcs=250,
        bin_representation="polygon",
        bin_to_cell_method="weighted_by_cluster",
        cell_annotation_method="celltypist",
        cell_typist_model="",
        run_synthetic=False,
        segmentation=True,
        bin_to_geodataframes=True,
        bin_to_cell_assignment=True,
        cell_type_annotation=True,
        block_size=4096,
        prob_thresh=.005,
        overlap_thresh=.001,
        min_overlap=28,
        context=128,
        n_tiles=(4,4,1),
        stardist_modelname="2D_versatile_he",
        channel_to_segment=2,
        cell_markers={},
        chunks_to_run=[],
        configs_dict={},
    ):
        """
        Initialize the class with the following parameters:

        Args:
            cache_dir (str): Directory to cache ENACT results. This must be populated 
                by the user.
            wsi_path (str): Path to the Whole Slide Image (WSI) file. This must be 
                populated by the user.
            visiumhd_h5_path (str): Path to the Visium HD h5 file containing spatial 
                transcriptomics data. This must be populated by the user.
            tissue_positions_path (str): Path to the tissue positions file that 
                contains spatial locations of barcodes. This must be populated by the 
                user.
            analysis_name (str): Name of the analysis, used for output directories and 
                results. Default is "enact_demo".
            seg_method (str): Cell segmentation method. Default is "stardist". 
                Options: ["stardist"].
            image_type (str): Specify if image is H&E (he) or IF(if) image. Default is "he". 
                Options: ["he", "if"].
            patch_size (int): Size of patches (in pixels) to process the image. Use a 
                smaller patch size to reduce memory requirements. Default is 4000.
            use_hvg (bool): Whether to use highly variable genes (HVG) during the 
                analysis. Default is True. Options: [True].
            n_hvg (int): Number of highly variable genes to use if `use_hvg` is True. 
                Default is 1000.
            n_clusters (int): Number of clusters. Only used if `bin_to_cell_method` is 
                "weighted_by_cluster". Default is 4.
            bin_representation (str): Representation type for VisiumHD bins. Default is 
                "polygon". Options: ["polygon"].
            bin_to_cell_method (str): Method to assign bins to cells. Default is 
                "weighted_by_cluster". Options: ["naive", "weighted_by_area", 
                "weighted_by_gene", "weighted_by_cluster"].
            cell_annotation_method (str): Method for annotating cell types. Default is 
                "celltypist". Options: ["celltypist", "sargent" (if installed), 
                "cellassign"].
            cell_typist_model (str): Path to the pre-trained CellTypist model for cell 
                type annotation. Only used if `cell_annotation_method` is 
                "celltypist". Refer to https://www.celltypist.org/models for a full 
                list of models. Default is an empty string.
            run_synthetic (bool): Whether to run synthetic data generation for testing 
                purposes. Default is False.
            segmentation (bool): Flag to run the image segmentation step. Default is 
                True.
            bin_to_geodataframes (bool): Flag to convert the bins to GeoDataFrames. 
                Default is True.
            bin_to_cell_assignment (bool): Flag to run bin-to-cell assignment. Default 
                is True.
            cell_type_annotation (bool): Flag to run cell type annotation. Default is 
                True.
            block_size (int): stardist parameter, the size of image blocks the model 
                processes at a time.
            prob_thresh (float): stardist parameter, value between 0 and 1, higher values 
                lead to fewer segmented objects, but will likely avoid false positives.
            overlap_thresh (float): stardist parameter,value between 0 and 1, higher 
                values allow segmented objects to overlap substantially.
            min_overlap (int): stardist parameter, overlap between blocks, should it 
                be larger than the size of a cell
            context (int): stardist parameter, context pixels around the blocks to be 
                included during prediction
            n_tiles (iterable): stardist parameter, This parameter denotes a tuple of 
                the number of tiles for every image axis
            stardist_modelname(str): Name or Path to the pre-trained Stardist model for image segmentation. 
                Refer to https://github.com/stardist/stardist?tab=readme-ov-file#pretrained-models-for-2d for a full 
                list of models. Default is "2D_versatile_he".
            channel_to_segment(int): Only applicable for IF images. This is the image channel to segment 
                (usually the DAPI channel). Default is 2.
            cell_markers (dict): A dictionary of cell markers used for annotation. Only 
                used if `cell_annotation_method` is one of ["sargent", "cellassign"].
            chunks_to_run (list): Specific chunks of data to run the analysis on for 
                debugging purposes. Default is an empty list (runs all chunks).
            configs_dict (dict): Dictionary containing ENACT configuration parameters. 
                If provided, the values in `configs_dict` will override any 
                corresponding parameters passed directly to the class constructor. This 
                is useful for running ENACT with a predefined configuration for 
                convenience and consistency. Default is an empty dictionary (i.e., 
                using the parameters defined in the class constructor).
        """

        # Todo: add class documentation
        user_configs = {
            "analysis_name": analysis_name,
            "run_synthetic": run_synthetic,
            "cache_dir": cache_dir,
            "paths": {
                "wsi_path": wsi_path,
                "visiumhd_h5_path": visiumhd_h5_path,
                "tissue_positions_path": tissue_positions_path,
            },
            "steps": {
                "segmentation": segmentation,
                "bin_to_geodataframes": bin_to_geodataframes,
                "bin_to_cell_assignment": bin_to_cell_assignment,
                "cell_type_annotation": cell_type_annotation,
            },
            "params": {
                "seg_method": seg_method,
                "image_type": image_type,
                "nucleus_expansion": nucleus_expansion,
                "expand_by_nbins": expand_by_nbins,
                "patch_size": patch_size,
                "bin_representation": bin_representation,
                "bin_to_cell_method": bin_to_cell_method,
                "cell_annotation_method": cell_annotation_method,
                "cell_typist_model": cell_typist_model,
                "use_hvg": use_hvg,
                "n_hvg": n_hvg,
                "destripe_norm": destripe_norm,
                "n_clusters": n_clusters,
                "n_pcs": n_pcs,
                "chunks_to_run": chunks_to_run,
            },
            "stardist": {
                "block_size": block_size,
                "prob_thresh": prob_thresh,
                "overlap_thresh": overlap_thresh,
                "min_overlap": min_overlap,
                "context": context,
                "n_tiles": n_tiles,
                "stardist_modelname": stardist_modelname,
                "channel_to_segment": channel_to_segment
            },
            
            "cell_markers": cell_markers,
        }
        self.configs = user_configs
        if configs_dict != {}:
            # If user specifies a configs_dict -> it overwrites all other user-specified parameters
            self.overwrite_configs(configs_dict)

        if self.configs["cache_dir"] == "":
            raise ValueError(f"Error: Please provide a value for 'cache_dir'.")

        if self.configs["steps"]["cell_type_annotation"]:
            if self.configs["params"]["cell_annotation_method"] == "celltypist":
                if self.configs["params"]["cell_typist_model"] == "":
                    raise ValueError(
                        f"Error: Please provide a value for 'cell_typist_model'. "
                        "Refer to https://www.celltypist.org/models for a full list of models."
                    )

            if self.configs["params"]["cell_annotation_method"] in [
                "sargent",
                "cellassign",
            ]:
                if self.configs["cell_markers"] == {}:
                    raise ValueError(f"Error: Please provide a value for 'cell_markers'.")

        # Load input files
        core_paths = ["wsi_path", "visiumhd_h5_path", "tissue_positions_path"]
        for core_path in core_paths:
            if self.configs["paths"][core_path] == "":
                raise ValueError(f"Error: Please provide a value for '{core_path}'. ")
        self.initiate_instance_variables()
        self.load_configs()

    def overwrite_configs(self, configs_dict):
        """
        Function overwrites the configurations with the content in the configs_dict
        """
        for key, value in configs_dict.items():
            if key in ["cell_markers"] or not isinstance(value, dict):
                self.configs[key] = value
            else:
                for sub_key, sub_value in value.items():
                    self.configs[key][sub_key] = sub_value

    def initiate_instance_variables(self):
        """
        Creates instance variables for all the keys defined in the config file
        """
        run_details = "ENACT running with the following configurations: \n"
        kwargs = {}
        for key, value in self.configs.items():
            if key in ["cell_markers"] or not isinstance(value, dict):
                setattr(self, key, value)
                run_details += f" {key}: {value}\n"
                kwargs[key] = value
            else:
                for sub_key, sub_value in value.items():
                    setattr(self, sub_key, sub_value)
                    run_details += f" {sub_key}: {sub_value}\n"
                    kwargs[sub_key] = sub_value
        self.run_details = run_details
        self.kwargs = kwargs

    def load_configs(self):
        """Loading the configuations and parameters"""
        # Generating paths
        self.cache_dir = os.path.join(self.cache_dir, self.analysis_name)
        self.nuclei_df_path = os.path.join(self.cache_dir, "nuclei_df.csv")
        self.cells_df_path = os.path.join(self.cache_dir, "cells_df.csv")
        self.cells_layer_path = os.path.join(self.cache_dir, "cells_layer.png")
        self.cell_chunks_dir = os.path.join(self.cache_dir, "chunks", "cells_gdf")
        self.bin_chunks_dir = os.path.join(self.cache_dir, "chunks", "bins_gdf")
        self.bin_assign_dir = os.path.join(
            self.cache_dir, "chunks", self.bin_to_cell_method, "bin_to_cell_assign"
        )
        self.cell_ix_lookup_dir = os.path.join(
            self.cache_dir, "chunks", self.bin_to_cell_method, "cell_ix_lookup"
        )
        os.makedirs(self.cache_dir, exist_ok=True)
        self.cellannotation_results_dir = os.path.join(
            self.cache_dir,
            "chunks",
            self.bin_to_cell_method,
            f"{self.cell_annotation_method}_results",
        )
        os.makedirs(self.cellannotation_results_dir, exist_ok=True)
        os.makedirs(self.cell_chunks_dir, exist_ok=True)
        self.logger = get_logger("ENACT", self.cache_dir)
        self.logger.info(f"<initiate_instance_variables> {self.run_details}")

    def load_image(self, file_path=None):
        """Load image from given file path
        Arguments:
            file_path {string} : path to the file that we are trying to load
        Returns:
            np.array -- loaded image as numpy array
        """
        if file_path == None:
            file_path = self.wsi_path
        img_arr = tifi.imread(file_path)
        crop_bounds = self.get_image_crop_bounds()
        x_min, y_min, x_max, y_max = crop_bounds
        if self.image_type == "if":
            # IF images have a different shape: (# channels, height, width)
            img_arr = img_arr[self.channel_to_segment, y_min:y_max, x_min:x_max]
        else:
            # H&E images have a shape: (height, width, # channels)
            img_arr = img_arr[y_min:y_max, x_min:x_max, :]
        self.logger.info("<load_image> Successfully loaded image!")
        return img_arr, crop_bounds

    def get_image_crop_bounds(self):
        """Get the crop location of the image to adjust the coordinates accordingly

        Returns:
            _type_: _description_
        """
        tissue_pos_list = pd.read_parquet(self.tissue_positions_path)

        # Cleaning up, removing negative coords,removing out of tissue bins
        tissue_pos_list_filt = tissue_pos_list[tissue_pos_list.in_tissue == 1]
        tissue_pos_list_filt = tissue_pos_list_filt.copy()
        tissue_pos_list_filt["pxl_row_in_fullres"] = tissue_pos_list_filt[
            "pxl_row_in_fullres"
        ].astype(int)
        tissue_pos_list_filt["pxl_col_in_fullres"] = tissue_pos_list_filt[
            "pxl_col_in_fullres"
        ].astype(int)
        tissue_pos_list_filt = tissue_pos_list_filt.loc[
            (tissue_pos_list_filt.pxl_row_in_fullres >= 0)
            & (tissue_pos_list_filt.pxl_col_in_fullres >= 0)
        ]
        x_min = tissue_pos_list_filt["pxl_col_in_fullres"].min()
        y_min = tissue_pos_list_filt["pxl_row_in_fullres"].min()
        x_max = tissue_pos_list_filt["pxl_col_in_fullres"].max()
        y_max = tissue_pos_list_filt["pxl_row_in_fullres"].max()
        return (x_min, y_min, x_max, y_max)

    def normalize_image(self, image, min_percentile=5, max_percentile=95):
        """_summary_

        Args:
            image (_type_): _description_
            min_percentile (int, optional): _description_. Defaults to 5.
            max_percentile (int, optional): _description_. Defaults to 95.

        Returns:
            _type_: _description_
        """
        # Adjust min_percentile and max_percentile as needed
        image_norm = normalize(image, min_percentile, max_percentile)
        self.logger.info("<normalize_image> Successfully normalized image!")
        return image_norm

    def segment_cells(self, image,):
        """_summary_

        Args:
            image (_type_): _description_
            prob_thresh (float, optional): _description_. Defaults to 0.005.

        Returns:
            _type_: _description_
        """
        if self.seg_method == "stardist":
            # Adjust nms_thresh and prob_thresh as needed
            # ssl._create_default_https_context = ssl._create_unverified_context
            self.stardist_model = StarDist2D.from_pretrained(self.stardist_modelname)
            if self.stardist_modelname=="2D_versatile_he" and self.image_type=="if":
                self.logger.warning(
                    f"<segment_cells> User requested using {self.stardist_modelname} on an IF image."
                    " This may lead to poor segmentation performance. Please consider setting "
                    "'stardist_modelname' to '2D_versatile_fluo' for IF images."
                )
            if isinstance(self.n_tiles, str):
                n_tiles = ast.literal_eval(self.n_tiles)  # Evaluate if it's a string
            else:
                n_tiles = self.n_tiles  # Use it as-is if it's not a string
            if self.image_type == "if":
                # IF using an IF image, pick a channel for segmentation. Usually the DAPI channel.
                n_tiles = n_tiles[:2]
                axes = "YX"
            else:
                axes = "YXC"
            labels, polys = self.stardist_model.predict_instances_big(
                image,
                axes=axes,
                block_size=self.block_size,
                prob_thresh=self.prob_thresh,
                nms_thresh=self.overlap_thresh,
                min_overlap=self.min_overlap,
                context=self.context,
                n_tiles=n_tiles,
                normalizer=None
            )
            #store resulting labels as sparse matrix NPZ - super efficient space wise
            labels_sparse = sparse.csr_matrix(labels)
            labels_npz_path = self.cache_dir + "stardist_labels.npz"
            sparse.save_npz(labels_npz_path, labels_sparse)
            nr_of_labels = len(np.unique(labels_sparse.data))
            self.logger.info(f"<run_segmentation> Found {nr_of_labels} objects")
            self.logger.info("<run_segmentation> Successfully segmented cells!")
            return labels, polys
        else:
            self.logger.warning("<run_segmentation> Invalid cell segmentation model!")
            return None, None

    def expand_nuclei_with_voronoi(self, gdf, expansion_size):
        """
        Expands the nuclei polygons within their corresponding Voronoi cells.

        Parameters:
        - gdf: GeoDataFrame with a geometry column containing nuclei polygons.
        - expansion_size: Distance by which to expand the nuclei (in same units as the GeoDataFrame's CRS).

        Returns:
        - GeoDataFrame with a new column 'expanded_geometry' containing the expanded polygons.
        """
        # Step 1: Compute centroids of the nuclei
        nuclei_polygons = gdf.geometry
        centroids = [poly.centroid for poly in nuclei_polygons]
        points = np.array([(pt.x, pt.y) for pt in centroids])

        # Step 2: Generate Voronoi tessellation
        vor = Voronoi(points)
        voronoi_cells = []
        voronoi_cell_mapping = dict(
            zip(range(len(nuclei_polygons)), [(None, None)] * len(nuclei_polygons))
        )

        num_infinite_regions = 0
        for i, region_index in enumerate(vor.point_region):
            vertices = vor.regions[region_index]
            if -1 not in vertices and len(vertices) > 0:  # Exclude infinite regions
                cell_polygon = Polygon([vor.vertices[v] for v in vertices])
                voronoi_cells.append(cell_polygon)
                voronoi_cell_mapping[i] = (
                    nuclei_polygons[i],
                    cell_polygon,
                )  # Map region index to original nucleus
            else:
                num_infinite_regions += 1

        nuclei_voronoi_pair_list = list(voronoi_cell_mapping.values())

        # Step 3: Expand nuclei within their Voronoi cells
        expanded_polygons = []
        num_unexpanded_cells = 0
        for nuclei_voronoi_pair in nuclei_voronoi_pair_list:
            nucleus, voronoi_cell = nuclei_voronoi_pair
            if voronoi_cell is None:
                expanded_polygons.append(nucleus)  # Append None if no Voronoi cell
                num_unexpanded_cells += 1
                continue
            expanded_nucleus = nucleus.buffer(expansion_size)  # Expand the nucleus
            clipped_expansion = expanded_nucleus.intersection(
                voronoi_cell
            )  # Clip to Voronoi cell
            expanded_polygons.append(clipped_expansion)

        # Step 4: Create new GeoDataFrame with expanded polygons
        gdf["geometry"] = expanded_polygons
        self.logger.info(
            f"<expand_nuclei_with_voronoi> Number of unexpanded cells: {num_unexpanded_cells}"
        )
        return gdf

    def convert_stardist_output_to_gdf(self, cell_polys, save_path=None):
        """Convert stardist output to geopandas dataframe

        Args:
            cell_polys (_type_): _description_
            save_path (_type_, optional): _description_. Defaults to None.
        """
        if save_path == None:
            save_path = self.nuclei_df_path
        # Creating a list to store Polygon geometries
        geometries = []
        centroids = []
        cell_x, cell_y = [], []

        # Iterating through each nuclei in the 'polys' DataFrame
        for nuclei in range(len(cell_polys["coord"])):
            # Extracting coordinates for the current nuclei and converting them to (y, x) format
            coords = [
                (y, x)
                for x, y in zip(
                    cell_polys["coord"][nuclei][0], cell_polys["coord"][nuclei][1]
                )
            ]
            # Creating a Polygon geometry from the coordinates
            poly = Polygon(coords)
            centroid = list(poly.centroid.coords)[0]
            centroids.append(centroid)
            geometries.append(poly)
            cell_x.append(centroid[0])
            cell_y.append(centroid[1])

        # Creating a GeoDataFrame using the Polygon geometries
        gdf = gpd.GeoDataFrame(geometry=geometries)
        gdf["id"] = [f"ID_{i+1}" for i, _ in enumerate(gdf.index)]
        gdf["cell_x"] = cell_x
        gdf["cell_y"] = cell_y
        gdf["centroid"] = centroids
        gdf.to_csv(save_path)
        self.logger.info(
            f"<convert_stardist_output_to_gdf> Mean nuclei area: {gdf.geometry.area.mean()}"
        )
        return gdf

    def split_df_to_chunks(self, df, x_col, y_col, output_dir):
        """
        Break the cells df into files, of size patch_size x patch_size

        Args:
            df (_type_): _description_
        """
        os.makedirs(output_dir, exist_ok=True)
        i = 0
        # Need to divide into chunks of patch_size pixels by patch_size pixels
        df[["patch_x", "patch_y"]] = (df[[x_col, y_col]] / self.patch_size).astype(int)
        df["patch_id"] = df["patch_x"].astype(str) + "_" + df["patch_y"].astype(str)
        self.logger.info(
            f"<split_df_to_chunks> Splitting into chunks. output_dir: {output_dir}"
        )
        unique_patches = df.patch_id.unique()
        for patch_id in tqdm(unique_patches, total=len(unique_patches)):
            patch_cells = df[df.patch_id == patch_id]
            if len(patch_cells) == 0:
                continue
            patch_cells.to_csv(os.path.join(output_dir, f"patch_{patch_id}.csv"))

    def destripe(self, adata, quantile=0.99):
        """Adaptation of the destripe method from Bin2cell:
        https://github.com/Teichlab/bin2cell/blob/main/bin2cell/bin2cell.py
        All credit goes to authors of Bin2cell (https://academic.oup.com/bioinformatics/article/40/9/btae546/7754061)
        Args:
            adata (_type_): _description_
            quantile (float, optional): _description_. Defaults to 0.99.

        Returns:
            _type_: _description_
        """
        self.logger.info(f"<destripe> Running destripe normalization")

        counts_key = "n_counts"
        factor_key = "destripe_factor"
        adjusted_counts_key = "n_counts_adjusted"
        adata.obs[counts_key] = adata.X.sum(axis=1)
        quant = adata.obs.groupby("array_row")[counts_key].quantile(quantile)
        # divide each row by its quantile (order of obs[counts_key] and obs[array_row] match)
        adata.obs[factor_key] = adata.obs[counts_key] / adata.obs["array_row"].map(
            quant
        )

        # repeat on columns
        quant = adata.obs.groupby("array_col")[factor_key].quantile(quantile)
        adata.obs[factor_key] /= adata.obs["array_col"].map(quant)

        # propose adjusted counts as the global quantile multipled by the destripe factor
        adata.obs[adjusted_counts_key] = adata.obs[factor_key] * np.quantile(
            adata.obs[counts_key], quantile
        )

        sc._utils.view_to_actual(adata)
        # adjust the count matrix to have n_counts_adjusted sum per bin (row)
        # premultiplying by a diagonal matrix multiplies each row by a value: https://solitaryroad.com/c108.html
        bin_scaling = sparse.diags(
            adata.obs[adjusted_counts_key] / adata.obs[counts_key]
        )
        adata_scaled = adata.copy()
        adata_scaled.X = bin_scaling.dot(adata_scaled.X)

        # adata_scaled.write_h5ad(
        #     os.path.join(self.cache_dir, "destriped_adata.h5"),
        #     compression="gzip",
        # )
        self.logger.info(f"<destripe> Successfully ran destripe normalization")
        return adata_scaled

    def get_bin_size(self):
        """Gets the bin size
        Returns:
            int: bin size in pixels
        """
        # Load the Spatial Coordinates
        df_tissue_positions = pd.read_parquet(self.tissue_positions_path)

        first_row = df_tissue_positions[
            (df_tissue_positions["array_row"] == 0)
            & (df_tissue_positions["array_col"] == 0)
        ]["pxl_col_in_fullres"]
        second_row = df_tissue_positions[
            (df_tissue_positions["array_row"] == 0)
            & (df_tissue_positions["array_col"] == 1)
        ]["pxl_col_in_fullres"]
        bin_size = np.abs(second_row.iloc[0] - first_row.iloc[0])
        self.logger.info(
            f"<get_bin_size> Bin size computed: {bin_size} pixels"
        )
        return bin_size

    def load_visiumhd_dataset(self, crop_bounds, destripe=False):
        """Loads the VisiumHD dataset and adjusts the
        coordinates to the cropped image

        Args:
            crop_bounds (tuple): crop bounds

        Returns:
            AnnData: AnnData object with the VisiumHD data
            int: bin size in pixels
        """

        # Accounting for crop bounds
        if crop_bounds is not None:
            x1, y1, _, _ = crop_bounds
        else:
            x1, y1 = (0, 0)
        # Load Visium HD data
        adata = sc.read_10x_h5(self.visiumhd_h5_path)

        # Load the Spatial Coordinates
        df_tissue_positions = pd.read_parquet(self.tissue_positions_path)

        # Set the index of the dataframe to the barcodes
        df_tissue_positions = df_tissue_positions.set_index("barcode")

        # Create an index in the dataframe to check joins
        df_tissue_positions["index"] = df_tissue_positions.index

        # *Important step*: Representing coords in the cropped WSI frame
        df_tissue_positions["pxl_row_in_fullres"] = (
            df_tissue_positions["pxl_row_in_fullres"] - y1
        )
        df_tissue_positions["pxl_col_in_fullres"] = (
            df_tissue_positions["pxl_col_in_fullres"] - x1
        )
        # Adding the tissue positions to the meta data
        adata.obs = pd.merge(
            adata.obs, df_tissue_positions, left_index=True, right_index=True
        )

        if destripe:
            adata = self.destripe(adata)

        first_row = df_tissue_positions[
            (df_tissue_positions["array_row"] == 0)
            & (df_tissue_positions["array_col"] == 0)
        ]["pxl_col_in_fullres"]
        second_row = df_tissue_positions[
            (df_tissue_positions["array_row"] == 0)
            & (df_tissue_positions["array_col"] == 1)
        ]["pxl_col_in_fullres"]
        bin_size = np.abs(second_row.iloc[0] - first_row.iloc[0])
        if self.configs["params"]["use_hvg"]:
            # Keeping the top n highly variable genes + the user requested cell markers
            n_genes = self.configs["params"]["n_hvg"]
            # Normalizing to median total counts
            adata_norm = adata.copy()
            sc.pp.normalize_total(adata_norm)
            # Logarithmize the data
            sc.pp.log1p(adata_norm)
            sc.pp.highly_variable_genes(adata_norm, n_top_genes=n_genes)

            hvg_mask = adata_norm.var["highly_variable"]
            cell_markers = [
                item
                for sublist in self.configs["cell_markers"].values()
                for item in sublist
            ]
            hv_genes = set(adata_norm[:,hvg_mask].var.index)
            all_genes = set(adata_norm.var.index)
            missing_markers = set(cell_markers) - all_genes # Genes missing from .h5 file
            self.logger.info(
                f"<load_visiumhd_dataset> Missing the following markers: {missing_markers}"
            )
            # Removing missing cell markers from consideration
            cell_markers = set(cell_markers) - missing_markers
            gene_list = list(hv_genes | cell_markers)
            hvg_mask = hvg_mask.copy()
            hvg_mask.loc[gene_list] = True
            adata = adata[:, hvg_mask]
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        return adata, bin_size

    def generate_bin_polys(self, bins_df, x_col, y_col, bin_size):
        """Represents the bins as Shapely polygons

        Args:
            bins_df (pd.DataFrame): bins dataframe
            x_col (str): column with the bin centre x-coordinate
            y_col (str): column with the bin centre y-coordinate
            bin_size (int): bin size in pixels

        Returns:
            list: list of Shapely polygons
        """
        geometry = []
        # Generates Shapely polygons to represent each bin
        if self.bin_representation == "point":
            # Geometry column is just the centre (x, y) for a VisiumHD bin
            geometry = [Point(xy) for xy in zip(bins_df[x_col], bins_df[y_col])]
        elif self.bin_representation == "polygon":
            self.logger.info(
                f"<generate_bin_polys> Generating bin polygons. num_bins: {len(bins_df)}"
            )
            half_bin_size = bin_size / 2
            bbox_coords = pd.DataFrame(
                {
                    "min_x": bins_df[x_col] - half_bin_size,
                    "min_y": bins_df[y_col] - half_bin_size,
                    "max_x": bins_df[x_col] + half_bin_size,
                    "max_y": bins_df[y_col] + half_bin_size,
                }
            )
            geometry = [
                shapely.geometry.box(min_x, min_y, max_x, max_y)
                for min_x, min_y, max_x, max_y in tqdm(
                    zip(
                        bbox_coords["min_x"],
                        bbox_coords["min_y"],
                        bbox_coords["max_x"],
                        bbox_coords["max_y"],
                    ),
                    total=len(bins_df),
                )
            ]
        else:
            self.logger.warning("<generate_bin_polys> Select a valid mode!")
        return geometry

    def convert_adata_to_cell_by_gene(self, adata):
        """Converts the AnnData object from bin-by-gene to
        cell-by-gene AnnData object.

        Args:
            adata (AnnData): bin-by-gene AnnData

        Returns:
            AnnData: cell-by-gene AnnData
        """
        # Group the data by unique cell IDs
        groupby_object = adata.obs.groupby(["id"], observed=True)

        # Extract the gene expression counts from the AnnData object
        counts = adata.X.tocsr()

        # Obtain the number of unique nuclei and the number of genes in the expression data
        N_groups = groupby_object.ngroups
        N_genes = counts.shape[1]

        # Initialize a sparse matrix to store the summed gene counts for each nucleus
        summed_counts = sparse.lil_matrix((N_groups, N_genes))

        # Lists to store the IDs of polygons and the current row index
        polygon_id = []
        row = 0
        # Iterate over each unique polygon to calculate the sum of gene counts.
        for polygons, idx_ in groupby_object.indices.items():
            summed_counts[row] = counts[idx_].sum(0)
            row += 1
            polygon_id.append(polygons)
        # Create an AnnData object from the summed count matrix
        summed_counts = summed_counts.tocsr()
        cell_by_gene_adata = anndata.AnnData(
            X=summed_counts,
            obs=pd.DataFrame(polygon_id, columns=["id"], index=polygon_id),
            var=adata.var,
        )
        return cell_by_gene_adata

    def generate_bins_gdf(self, adata, bin_size):
        """Convert the bins Anndata object to a geodataframe

        Args:
            adata (_type_): _description_

        Returns:
            _type_: _description_
        """
        bin_coords_df = adata.obs.copy()
        geometry = self.generate_bin_polys(
            bins_df=bin_coords_df,
            x_col="pxl_col_in_fullres",
            y_col="pxl_row_in_fullres",
            bin_size=bin_size,
        )
        bins_gdf = gpd.GeoDataFrame(bin_coords_df, geometry=geometry)
        return bins_gdf

    def assign_bins_to_cells(self, adata, crop_bounds):
        """Assigns bins to cells based on method requested by the user

        Args:
            adata (_type_): _description_
        """
        os.makedirs(self.bin_assign_dir, exist_ok=True)
        os.makedirs(self.cell_ix_lookup_dir, exist_ok=True)
        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(self.cell_chunks_dir)
        self.logger.info(
            f"<assign_bins_to_cells> Assigning bins to cells using {self.bin_to_cell_method} method"
        )
        for chunk in tqdm(chunk_list, total=len(chunk_list)):
            bin_to_cell_chunk_path = os.path.join(self.bin_assign_dir, chunk)
            if os.path.exists(bin_to_cell_chunk_path) or chunk in {".ipynb_checkpoints"}:
                continue

            bin_gdf_chunk_path = os.path.join(self.bin_chunks_dir, chunk)
            cell_gdf_chunk_path = os.path.join(self.cell_chunks_dir, chunk)

            if not (os.path.exists(bin_gdf_chunk_path) and os.path.exists(cell_gdf_chunk_path)):
                continue

            # Loading the cells geodataframe
            cell_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(cell_gdf_chunk_path))
            cell_gdf_chunk = cell_gdf_chunk[~cell_gdf_chunk["geometry"].isna()]
            cell_gdf_chunk["geometry"] = cell_gdf_chunk["geometry"].apply(wkt.loads)
            cell_gdf_chunk = gpd.GeoDataFrame(cell_gdf_chunk, geometry="geometry")

            # Loading the bins geodataframe
            bin_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(bin_gdf_chunk_path))
            bin_gdf_chunk["geometry"] = bin_gdf_chunk["geometry"].apply(wkt.loads)
            bin_gdf_chunk.set_geometry("geometry", inplace=True)
            # Perform a spatial join to check which coordinates are in a cell nucleus
            result_spatial_join = gpd.sjoin(
                bin_gdf_chunk,
                cell_gdf_chunk[["geometry", "id", "cell_x", "cell_y"]],
                how="left",
                predicate="intersects",
            )

            # Only keeping the bins that overlap with a cell. index_right = index of the cell
            result_spatial_join = result_spatial_join[
                ~result_spatial_join["index_right"].isna()
            ]

            # Getting unique bins and overlapping bins. bins are repeated for the cells that share them
            barcodes_in_overlaping_polygons = pd.unique(
                result_spatial_join[result_spatial_join.duplicated(subset=["index"])][
                    "index"
                ]
            )
            result_spatial_join["unique_bin"] = ~result_spatial_join["index"].isin(
                barcodes_in_overlaping_polygons
            )
            # Filter the adata object to contain only the barcodes in result_spatial_join
            # shape: (#bins_overlap x #genes)
            expanded_adata = adata[result_spatial_join["index"]]
            # Adding the cell ids to the anndata object (the cell that the bin is assigned to)
            # Can have duplicate bins (i.e. "expanded") if a bin is assigned to more than one cell
            expanded_adata = expanded_adata.copy() 
            expanded_adata.obs["id"] = result_spatial_join["id"].tolist()

            # Reshape the anndata object to (#cells x #genes)
            filtered_result_spatial_join = result_spatial_join[
                result_spatial_join["unique_bin"]
            ]
            filtered_adata = adata[filtered_result_spatial_join["index"]]
            filtered_adata = filtered_adata.copy()
            filtered_adata.obs["id"] = filtered_result_spatial_join["id"].tolist()

            unfilt_result_spatial_join = result_spatial_join.copy()
            self.logger.info("<assign_bins_to_cells> done spatial join")

            if result_spatial_join.empty:
                self.logger.info(
                    "result_spatial_join is empty, skipping bin-to-cell assignment."
                )
                continue

            elif self.bin_to_cell_method == "naive":
                result_spatial_join = naive_assignment(result_spatial_join)
                expanded_adata = filtered_adata.copy()

            elif self.bin_to_cell_method == "weighted_by_area":
                result_spatial_join, expanded_adata = weight_by_area_assignment(
                    result_spatial_join, expanded_adata, cell_gdf_chunk
                )

            elif self.bin_to_cell_method == "weighted_by_gene":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_gene_assignment(
                    result_spatial_join, expanded_adata, unique_cell_by_gene_adata
                )

            elif self.bin_to_cell_method == "weighted_by_cluster":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_cluster_assignment(
                    result_spatial_join,
                    expanded_adata,
                    unique_cell_by_gene_adata,
                    n_clusters=self.configs["params"]["n_clusters"],
                    n_pcs=self.configs["params"]["n_pcs"],
                )
            else:
                self.logger.info("ERROR", self.bin_to_cell_method)
            self.logger.info("<assign_bins_to_cells> convert_adata_to_cell_by_gene")
            cell_by_gene_adata = self.convert_adata_to_cell_by_gene(expanded_adata)
            # TODO: CK save AnnData object
            # TODO: This save has to be moved to a place outside this for loop
            # self.logger.info("<assign_bins_to_cells> save expanded_adata to '.h5ad' file")
            expanded_adata.write_h5ad(
                os.path.join(self.bin_assign_dir, chunk + ".h5ad")
            )
            del expanded_adata

            # Save the gene to cell assignment results to a .csv file
            chunk_gene_to_cell_assign_df = pd.DataFrame(
                cell_by_gene_adata.X.toarray(),
                columns=cell_by_gene_adata.var_names,
            )

            chunk_gene_to_cell_assign_df = chunk_gene_to_cell_assign_df.loc[
                :, ~chunk_gene_to_cell_assign_df.columns.duplicated()
            ].copy()
            # Saving counts to cache
            chunk_gene_to_cell_assign_df.to_csv(
                os.path.join(self.bin_assign_dir, chunk)
            )

            # Getting number of bins shared between cells
            overlaps_df = (
                unfilt_result_spatial_join.groupby(["id", "unique_bin"])
                .count()["in_tissue"]
                .reset_index()
            )
            overlaps_df = overlaps_df.pivot(
                index="id", columns="unique_bin", values="in_tissue"
            ).fillna(0)
            try:
                overlaps_df.columns = ["num_shared_bins", "num_unique_bins"]
            except:
                overlaps_df.columns = ["num_unique_bins"]
                overlaps_df["num_shared_bins"] = 0
            cell_gdf_chunk = cell_gdf_chunk.merge(
                overlaps_df, how="left", left_on="id", right_index=True
            )
            cell_gdf_chunk[["num_shared_bins", "num_unique_bins"]] = cell_gdf_chunk[
                ["num_shared_bins", "num_unique_bins"]
            ].fillna(0)
            # Adjusting cell location based on crop boundaries
            if crop_bounds is not None:
                x1, y1, _, _ = crop_bounds
            else:
                x1, y1 = (0, 0)
            # Save index lookup to store x and y values and cell index
            index_lookup_df = cell_by_gene_adata.obs.merge(
                cell_gdf_chunk, how="left", left_index=True, right_on="id"
            )[
                ["cell_x", "cell_y", "num_shared_bins", "num_unique_bins", "id"]
            ].reset_index(
                drop=True
            )
            index_lookup_df["num_transcripts"] = chunk_gene_to_cell_assign_df.sum(
                axis=1
            )
            index_lookup_df["chunk_name"] = chunk
            index_lookup_df.to_csv(os.path.join(self.cell_ix_lookup_dir, chunk))
            self.logger.info(
                f"<assign_bins_to_cells> Processed {chunk} using {self.bin_to_cell_method}. Mean count per cell: {chunk_gene_to_cell_assign_df.sum(axis=1).mean()}"
            )
        
        self.logger.info(
            f"<assign_bins_to_cells> Successfully assigned bins to cells!"
        )

    def assign_bins_to_cells_synthetic(self):
        """Assigns bins to cells based on method requested by the user

        Args:
            adata (_type_): _description_
        """
        os.makedirs(self.bin_assign_dir, exist_ok=True)
        os.makedirs(self.cell_ix_lookup_dir, exist_ok=True)
        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(self.cell_chunks_dir)

        self.logger.info(
            f"<assign_bins_to_cells_synthetic> Assigning bins to cells using {self.bin_to_cell_method} method"
        )
        for chunk in tqdm(chunk_list, total=len(chunk_list)):
            if os.path.exists(os.path.join(self.cell_ix_lookup_dir, chunk)):
                continue
            if chunk in [".ipynb_checkpoints"]:
                continue
            # Loading the cells geodataframe
            cell_gdf_chunk_path = os.path.join(self.cell_chunks_dir, chunk)
            cell_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(cell_gdf_chunk_path))
            cell_gdf_chunk["geometry"] = cell_gdf_chunk["geometry"].apply(wkt.loads)
            cell_gdf_chunk.set_geometry("geometry", inplace=True)
            cell_gdf_chunk["geometry"] = cell_gdf_chunk["geometry"].buffer(0)
            # Loading the bins geodataframe
            bin_gdf_chunk_path = os.path.join(self.bin_chunks_dir, chunk)
            bin_gdf_chunk = gpd.GeoDataFrame(pd.read_csv(bin_gdf_chunk_path))
            bin_gdf_chunk["geometry"] = bin_gdf_chunk["geometry"].apply(wkt.loads)
            bin_gdf_chunk.set_geometry("geometry", inplace=True)

            # Perform a spatial join to check which coordinates are in a cell nucleus
            result_spatial_join = gpd.sjoin(
                bin_gdf_chunk[["geometry", "assigned_bin_id", "row", "column"]],
                cell_gdf_chunk[["geometry", "cell_id"]],
                how="left",
                predicate="intersects",
            )

            # Only keeping the bins that overlap with a cell
            result_spatial_join = result_spatial_join[
                ~result_spatial_join["index_right"].isna()
            ]

            # Getting unique bins and overlapping bins
            barcodes_in_overlaping_polygons = pd.unique(
                result_spatial_join[
                    result_spatial_join.duplicated(subset=["assigned_bin_id"])
                ]["assigned_bin_id"]
            )
            result_spatial_join["unique_bin"] = ~result_spatial_join[
                "assigned_bin_id"
            ].isin(barcodes_in_overlaping_polygons)
            bin_gdf_chunk = bin_gdf_chunk.set_index("assigned_bin_id")
            adata = anndata.AnnData(
                bin_gdf_chunk.drop(columns=["geometry", "row", "column"])
            )

            # Filter the adata object to contain only the barcodes in result_spatial_join
            # shape: (#bins_overlap x #genes)
            adata.obs_names_make_unique()
            expanded_adata = adata[result_spatial_join["assigned_bin_id"]]
            # Adding the cell ids to the anndata object (the cell that the bin is assigned to)
            # Can have duplicate bins (i.e. "expanded") if a bin is assigned to more than one cell
            expanded_adata.obs["id"] = result_spatial_join["cell_id"].tolist()
            expanded_adata.obs["index"] = result_spatial_join[
                "assigned_bin_id"
            ].tolist()

            # Reshape the anndata object to (#cells x #genes)
            filtered_result_spatial_join = result_spatial_join[
                result_spatial_join["unique_bin"]
            ]
            filtered_adata = adata[filtered_result_spatial_join["assigned_bin_id"]]
            filtered_adata.obs["id"] = filtered_result_spatial_join["cell_id"].tolist()
            if not sparse.issparse(filtered_adata.X):
                filtered_adata.X = sparse.csr_matrix(filtered_adata.X)
            if not sparse.issparse(expanded_adata.X):
                expanded_adata.X = sparse.csr_matrix(expanded_adata.X)
            self.logger.info("<assign_bins_to_cells> done spatial join")

            cell_gdf_chunk.rename(columns={"cell_id": "id"}, inplace=True)
            result_spatial_join.rename(
                columns={"assigned_bin_id": "index"}, inplace=True
            )
            result_spatial_join.rename(columns={"cell_id": "id"}, inplace=True)

            unfilt_result_spatial_join = result_spatial_join.copy()
            self.logger.info("<assign_bins_to_cells> done spatial join")
            if self.bin_to_cell_method == "naive":
                result_spatial_join = naive_assignment(result_spatial_join)
                expanded_adata = filtered_adata.copy()

            elif self.bin_to_cell_method == "weighted_by_area":
                # cell_gdf_chunk = cell_gdf_chunk.set_index('index_right')
                result_spatial_join, expanded_adata = weight_by_area_assignment(
                    result_spatial_join, expanded_adata, cell_gdf_chunk
                )

            elif self.bin_to_cell_method == "weighted_by_gene":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_gene_assignment(
                    result_spatial_join, expanded_adata, unique_cell_by_gene_adata
                )

            elif self.bin_to_cell_method == "weighted_by_cluster":
                unique_cell_by_gene_adata = self.convert_adata_to_cell_by_gene(
                    filtered_adata
                )
                result_spatial_join, expanded_adata = weight_by_cluster_assignment(
                    result_spatial_join,
                    expanded_adata,
                    unique_cell_by_gene_adata,
                    n_clusters=self.configs["params"]["n_clusters"],
                    n_pcs=self.configs["params"]["n_pcs"],
                )

            self.logger.info("<assign_bins_to_cells> convert_adata_to_cell_by_gene")

            if not sparse.issparse(expanded_adata.X):
                expanded_adata.X = sparse.csr_matrix(expanded_adata.X)
            cell_by_gene_adata = self.convert_adata_to_cell_by_gene(expanded_adata)
            del expanded_adata

            # Save the gene to cell assignment results to a .csv file
            chunk_gene_to_cell_assign_df = pd.DataFrame(
                cell_by_gene_adata.X.toarray(),
                columns=cell_by_gene_adata.var_names,
            )
            # Saving counts to cach
            chunk_gene_to_cell_assign_df.insert(
                0, "id", cell_by_gene_adata.obs["id"].values
            )

            chunk_gene_to_cell_assign_df = chunk_gene_to_cell_assign_df.loc[
                :, ~chunk_gene_to_cell_assign_df.columns.duplicated()
            ].copy()

            # Saving counts to cache
            chunk_gene_to_cell_assign_df.to_csv(
                os.path.join(self.bin_assign_dir, chunk)
            )

            # Getting number of bins shared between cells
            overlaps_df = (
                unfilt_result_spatial_join.groupby(["id", "unique_bin"])
                .count()["index"]
                .reset_index()
            )
            overlaps_df = overlaps_df.pivot(
                index="id", columns="unique_bin", values="index"
            ).fillna(0)

            try:
                overlaps_df.columns = ["num_shared_bins", "num_unique_bins"]
            except:
                overlaps_df.columns = ["num_unique_bins"]
                overlaps_df["num_shared_bins"] = 0
            cell_gdf_chunk = cell_gdf_chunk.merge(
                overlaps_df, how="left", left_on="id", right_index=True
            )
            cell_gdf_chunk[["num_shared_bins", "num_unique_bins"]] = cell_gdf_chunk[
                ["num_shared_bins", "num_unique_bins"]
            ].fillna(0)
            if "cell_x" not in cell_gdf_chunk.columns:
                cell_gdf_chunk["cell_x"] = cell_gdf_chunk["geometry"].centroid.x
            if "cell_y" not in cell_gdf_chunk.columns:
                cell_gdf_chunk["cell_y"] = cell_gdf_chunk["geometry"].centroid.y
            index_lookup_df = cell_by_gene_adata.obs.merge(
                cell_gdf_chunk, how="left", left_index=True, right_on="id"
            )[
                ["cell_x", "cell_y", "num_shared_bins", "num_unique_bins", "id"]
            ].reset_index(
                drop=True
            )

            index_lookup_df["num_transcripts"] = chunk_gene_to_cell_assign_df.drop(
                columns=["id"], axis=1
            ).sum(axis=1)
            index_lookup_df.to_csv(os.path.join(self.cell_ix_lookup_dir, chunk))
            self.logger.info(
                f"Number of shared bins {overlaps_df['num_shared_bins'].sum()}"
            )
            self.logger.info(f"{chunk} finished")
            self.logger.info(
                f"{self.bin_to_cell_method} mean count per cell: {index_lookup_df['num_transcripts'].mean()}"
            )

    def merge_files_sparse(self, input_folder):
        """
        
        For the large bin assignment directory (self.bin_assign_dir),
        this method loads each CSV as a sparse matrix to avoid creating a giant dense DataFrame.
        
        Args:
            input_folder (str): Directory path containing the CSV files.
        
        Returns:
            If directory equals self.bin_assign_dir:
                tuple: (X_sparse, gene_columns)
                    - X_sparse: a scipy.sparse.csr_matrix with stacked count data.
                    - gene_columns: list of gene names corresponding to the matrix columns.
            Otherwise:
                pd.DataFrame: the merged DataFrame.
        """
        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(input_folder)
        # For the large count files, load each as sparse.
        sparse_list = []
        gene_columns = []
        for filename in chunk_list:
            if filename in ["annotated.csv", ".ipynb_checkpoints"]:
                continue
            # Read each .csv file and append it to the list
            file_path = os.path.join(input_folder, filename)
            df = pd.read_csv(file_path)
            if "Unnamed: 0" in df.columns:
                df = df.drop(columns=["Unnamed: 0"])
            if gene_columns == []:
                gene_columns = df.columns.tolist()
            df = df.fillna(0).astype(int)
            # Convert to a sparse matrix and append.
            sparse_list.append(sparse.csr_matrix(df.values).astype(np.float32))
        # Vertically stack the sparse matrices.
        X_sparse = sparse.vstack(sparse_list)
        return (X_sparse, gene_columns)

    def merge_files(
        self, input_folder, output_file_name="merged_results.csv", save=True
    ):
        """Merges all files in a specified input folder into a single output file.

        Args:
            input_folder (str): The path to the folder containing the input files to be merged.
            output_file_name (str): The name of the output file.
        """
        # List to store the DataFrames
        csv_list = []
        output_file = os.path.join(input_folder, output_file_name)

        if self.configs["params"]["chunks_to_run"]:
            chunk_list = self.configs["params"]["chunks_to_run"]
        else:
            chunk_list = os.listdir(input_folder)
        # Loop through all files in the directory
        for filename in chunk_list:
            if filename in ["annotated.csv", ".ipynb_checkpoints"]:
                continue
            if "merged" in filename:
                continue

            # Read each .csv file and append it to the list
            file_path = os.path.join(input_folder, filename)
            df = pd.read_csv(file_path)
            csv_list.append(df)

        # Concatenate all DataFrames in the list
        concatenated_df = pd.concat(csv_list, ignore_index=True)

        if save:
            # Save the concatenated DataFrame to the output file
            concatenated_df.to_csv(output_file, index=False)
            self.logger.info(
                f"<merge_files> files have been merged and saved to {output_file}"
            )
        return concatenated_df

    def run_cell_type_annotation(self):
        """Runs cell type annotation"""
        ann_method = self.configs["params"]["cell_annotation_method"]
        if ann_method == "sargent":
            self.logger.info(
                f"<run_cell_type_annotation> Will launch Sargent separately. "
                "Please ensure Sargent is installed."
            )
        elif ann_method == "cellassign":
            from .cellassign import CellAssignPipeline

            cellassign_obj = CellAssignPipeline(**self.kwargs)
            cellassign_obj.format_markers_to_df()
            cellassign_obj.run_cell_assign()
            self.logger.info(
                f"<run_cell_type_annotation> Successfully ran CellAssign on Data."
            )

        elif ann_method == "celltypist":
            from .celltypist import CellTypistPipeline

            celltypist_obj = CellTypistPipeline(**self.kwargs)

            celltypist_obj.run_cell_typist()
            self.logger.info(
                f"<run_cell_type_annotation> Successfully ran CellTypist on Data."
            )
        else:
            self.logger.info(
                "<run_cell_type_annotation> Please select a valid cell annotation "
                "method. options=['cellassign', 'sargent']"
            )

    def package_results(self):
        """Packages the results of the pipeline"""
        from .package_results import PackageResults

        pack_obj = PackageResults(**self.kwargs)
        ann_method = self.configs["params"]["cell_annotation_method"]
        if ann_method == "sargent":
            results_df, cell_by_gene_df = pack_obj.merge_sargent_output_files()
            adata = pack_obj.df_to_adata(results_df, cell_by_gene_df)
            pack_obj.save_adata(adata)
            pack_obj.create_tmap_file()
            self.logger.info("<package_results> Packaged Sargent results")
        elif ann_method == "cellassign":
            cell_by_gene_df = pack_obj.merge_cellassign_output_files()
            results_df = pd.read_csv(
                os.path.join(self.cellannotation_results_dir, "merged_results.csv")
            )
            adata = pack_obj.df_to_adata(results_df, cell_by_gene_df)
            pack_obj.save_adata(adata)
            pack_obj.create_tmap_file()
            self.logger.info("<package_results> Packaged CellAssign results")
        elif ann_method == "celltypist":
            cell_by_gene_df = pack_obj.merge_cellassign_output_files()
            results_df = pd.read_csv(
                os.path.join(self.cellannotation_results_dir, "merged_results.csv")
            )
            adata = pack_obj.df_to_adata(results_df, cell_by_gene_df)
            pack_obj.save_adata(adata)
            pack_obj.create_tmap_file()
            self.logger.info("<package_results> Packaged CellTypist results")
        else:
            self.logger.info(
                f"<package_results> Please select a valid cell annotation method"
            )

    def convert_stardist_output_to_image(self, wsi_shape, cells_gdf):
        """Converts the stardist segmentation outputs to a .png with
        the cell outlines highlighted.

        Args:
            wsi_shape (tuple): Shape of the output image (height, width).
            cells_gdf (GeoDataFrame): GeoDataFrame containing cell geometries.
        """
        # Initialize an empty mask
        cells_mask = np.zeros(wsi_shape, dtype=np.uint8)

        # Create a PIL image and draw object
        pil_image = Image.fromarray(cells_mask)
        draw = ImageDraw.Draw(pil_image)

        # Iterate over polygons and draw them
        for poly in cells_gdf["geometry"].tolist():
            if poly is None:
                continue
            # Convert polygon coordinates to integer tuples
            poly_coords = [tuple(map(int, coord)) for coord in poly.exterior.coords]
            draw.line(poly_coords, fill=255, width=1)  # Draw polygon outline

        # Save the resulting mask as a PNG
        pil_image.save(self.cells_layer_path, format="PNG")
        return pil_image

    def run_enact(self):
        """Runs ENACT given the user-specified configs"""
        if not self.run_synthetic:
            # Loading image and getting shape and cropping boundaries (if applicable)
            wsi, crop_bounds = self.load_image()

            # Run cell segmentation
            if self.segmentation:
                if not os.path.exists(self.nuclei_df_path):
                    # If nuclei_df file does not exist, run Stardist to get cell nuclei
                    wsi_norm = self.normalize_image(
                        image=wsi, min_percentile=5, max_percentile=95
                    )
                    cell_labels, cell_polys = self.segment_cells(image=wsi_norm)
                    # cells_gdf has the nuclei boundaries
                    nuclei_gdf = self.convert_stardist_output_to_gdf(cell_polys=cell_polys)
                else:
                    # If nuclei_df file is present, load it to save time
                    nuclei_gdf = pd.read_csv(self.nuclei_df_path)
                    nuclei_gdf = nuclei_gdf[~nuclei_gdf['geometry'].isna()]
                    nuclei_gdf["geometry"] = nuclei_gdf["geometry"].apply(
                        wkt.loads
                    )
                    nuclei_gdf = gpd.GeoDataFrame(nuclei_gdf, geometry="geometry")
                    print(f"Mean cell area: {nuclei_gdf.geometry.area.mean()}")

                # If no expansion is requested, use the nuclei as cell boundaries
                cells_gdf = nuclei_gdf.copy()

                if self.nucleus_expansion:
                    # If user requested nuclei expansion, run expand nuclei
                    self.logger.info(f"<run_enact> Expanding nuclei to get cell boundaries")
                    bin_size = self.get_bin_size()
                    # Expanding the nuclei by the size of 2 bins all around
                    cells_gdf = self.expand_nuclei_with_voronoi(
                        nuclei_gdf, expansion_size=self.expand_by_nbins * bin_size
                    )
                    self.logger.info(
                        f"<convert_stardist_output_to_gdf> Mean nuclei area after expansion: {cells_gdf.geometry.area.mean()}"
                    )
                    # Save results to disk
                    cells_gdf.to_csv(self.cells_df_path)

                # Saving segmentation as a .png
                self.convert_stardist_output_to_image(
                    wsi_shape=wsi.shape, cells_gdf=cells_gdf
                )
                self.split_df_to_chunks(
                    df=cells_gdf,
                    x_col="cell_x",
                    y_col="cell_y",
                    output_dir=self.cell_chunks_dir,
                )
                del cells_gdf, nuclei_gdf # Clearing memory
            del wsi  # Clearing memory

            # Loading the VisiumHD reads
            if self.bin_to_geodataframes:
                bins_adata, bin_size = self.load_visiumhd_dataset(
                    crop_bounds, destripe=self.destripe_norm
                )
                # Convert VisiumHD reads to geodataframe objects
                bins_gdf = self.generate_bins_gdf(bins_adata, bin_size)
                # Splitting the bins geodataframe object
                self.split_df_to_chunks(
                    df=bins_gdf,
                    x_col="pxl_col_in_fullres",
                    y_col="pxl_row_in_fullres",
                    output_dir=self.bin_chunks_dir,
                )
                del bins_gdf

            # Run bin-to-cell assignment
            if self.bin_to_cell_assignment:
                bins_adata, bin_size = self.load_visiumhd_dataset(
                    crop_bounds, destripe=self.destripe_norm
                )
                self.assign_bins_to_cells(bins_adata, crop_bounds)

            # Run cell type annotation
            if self.cell_type_annotation:
                self.run_cell_type_annotation()
                self.package_results()

        else:
            # Generating synthetic data
            if self.analysis_name in ["xenium", "xenium_nuclei", "xenium_nuclei_debug"]:
                cells_gdf = pd.read_csv(self.cells_df_path)
                self.split_df_to_chunks(
                    df=cells_gdf,
                    x_col="cell_x",
                    y_col="cell_y",
                    output_dir=self.cell_chunks_dir,
                )
            self.assign_bins_to_cells_synthetic()


if __name__ == "__main__":
    # Creating ENACT object
    parser = argparse.ArgumentParser(description="Specify ENACT config file location.")
    parser.add_argument(
        "--configs_path", type=str, required=False, help="Config file location"
    )
    args = parser.parse_args()
    if not args.configs_path:
        configs_path = "config/configs.yaml"
    else:
        configs_path = args.configs_path
    print(f"<ENACT> Loading configurations from {configs_path}")
    with open(configs_path, "r") as stream:
        configs = yaml.safe_load(stream)
    so_hd = ENACT(configs_dict=configs)
    so_hd.run_enact()
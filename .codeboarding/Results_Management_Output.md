```mermaid

graph LR

    ENACT_Pipeline_Orchestrator["ENACT Pipeline Orchestrator"]

    Results_Formatter_Saver["Results Formatter & Saver"]

    ENACT_Pipeline_Orchestrator -- "orchestrates" --> Results_Formatter_Saver

    ENACT_Pipeline_Orchestrator -- "provides data to" --> Results_Formatter_Saver

    click ENACT_Pipeline_Orchestrator href "https://github.com/Sanofi-Public/enact-pipeline/blob/main/.codeboarding//ENACT_Pipeline_Orchestrator.md" "Details"

```

[![CodeBoarding](https://img.shields.io/badge/Generated%20by-CodeBoarding-9cf?style=flat-square)](https://github.com/CodeBoarding/GeneratedOnBoardings)[![Demo](https://img.shields.io/badge/Try%20our-Demo-blue?style=flat-square)](https://www.codeboarding.org/demo)[![Contact](https://img.shields.io/badge/Contact%20us%20-%20contact@codeboarding.org-lightgrey?style=flat-square)](mailto:contact@codeboarding.org)



## Component Details



This subsystem is critical for transforming raw and intermediate pipeline outputs into standardized, usable, and shareable formats for downstream analysis and visualization. It comprises two fundamental components: the `ENACT Pipeline Orchestrator` and the `Results Formatter & Saver`. These two components are fundamental to the "Results Management & Output" subsystem because they collectively ensure that the ENACT pipeline's efforts culminate in actionable and accessible data. The `ENACT Pipeline Orchestrator` is the driver that initiates and manages the entire result packaging workflow. It acts as the central decision-maker, ensuring that the correct steps are taken based on the pipeline's configuration (e.g., the chosen cell annotation method). Without this component, the complex process of result packaging would not be initiated or coordinated, leading to fragmented and unusable outputs. The `Results Formatter & Saver` is the worker that performs the actual, intricate data transformation, merging, and saving into standardized formats (AnnData) and visualization-ready files (TMAP). This is crucial for the usability and interpretability of the ENACT pipeline's output. Without this component, the raw results would remain in intermediate, often disparate, formats, making downstream analysis, integration with other tools, or visualization extremely difficult or impossible. Together, these components form a cohesive and essential subsystem that bridges the gap between raw computational results and their practical application in biological research.



### ENACT Pipeline Orchestrator

This component, primarily represented by the `ENACT` class, serves as the central coordinator for the entire ENACT pipeline. It manages the overall workflow, including the critical step of packaging results. It determines which cell annotation method was used and then delegates the specific result merging, formatting, and saving tasks to the `Results Formatter & Saver` component. It also directly handles the merging of chunked intermediate files.





**Related Classes/Methods**:



- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/pipeline.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.pipeline.ENACT.package_results` (0:0)</a>

- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/pipeline.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.pipeline.ENACT.merge_files` (0:0)</a>





### Results Formatter & Saver

This component, embodied by the `PackageResults` class, is responsible for the detailed handling of ENACT pipeline outputs. Its core functions include merging specific output files based on the cell annotation method (e.g., Sargent, CellAssign), converting processed DataFrames into AnnData objects, saving these AnnData objects, and generating TMAP files for visualization in TissUUmaps. It ensures that the final results are in standard and usable formats.





**Related Classes/Methods**:



- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/package_results.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.package_results.PackageResults.__init__` (0:0)</a>

- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/package_results.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.package_results.PackageResults.merge_sargent_output_files` (0:0)</a>

- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/package_results.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.package_results.PackageResults.merge_cellassign_output_files` (0:0)</a>

- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/package_results.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.package_results.PackageResults.df_to_adata` (0:0)</a>

- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/package_results.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.package_results.PackageResults.save_adata` (0:0)</a>

- <a href="https://github.com/Sanofi-Public/enact-pipeline/blob/master/src/enact/package_results.py#L0-L0" target="_blank" rel="noopener noreferrer">`enact.package_results.PackageResults.create_tmap_file` (0:0)</a>









### [FAQ](https://github.com/CodeBoarding/GeneratedOnBoardings/tree/main?tab=readme-ov-file#faq)
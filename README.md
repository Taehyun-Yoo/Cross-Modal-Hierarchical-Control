# Cross-Modal-Hierarchical-Control
The study investigated a graded rostro-caudal pattern of activation, representation, and functional connectivity in the prefrontal cortex (PFC).
Here, we provide the code and behavioral data for this study.

## Description

### `code/`

This folder contains scripts related to behavioral/fMRI analysis.

- `behavior/Behavioral_analysis.Rmd`: script for behavioral analysis

- `fMRI/EPI/`: scripts for fMRI analysis
  
  - `univairate_analysis/`: scripts for univariate analysis (each hierarchy)
    - `Preprocessing/`: script for preprocessing
    - `hierarchy_parametric/`: scripts for parametric modulation analysis
    - `ROI/`: scripts for ROI analysis
    - `SVD/`: scripts for SVD analysis
      
  - `multivairate_analysis/`: scripts for MVPA (each hierarchy)
    - `Preprocessing/`: scripts for preprocessing using the Least-Squares Single (LSS) approach
    - `searchlight/`: scripts for searchlight-based MVPA
    - `ROI/`: scripts for ROI-based MVPA
      
  - `FC_analysis/`: script for functional connectivity analysis (each hierarchy)

### `data/`

This folder contains both raw and modified behavioral data.
Raw fMRI data are also available via OpenNeuro (https://openneuro.org/datasets/ds006628)

- `behavior_raw/`: raw behavior data for each task (response: two runs / feature: four runs / dimension: two runs)
- `01~53/`: modified individual behaviora data for analysis

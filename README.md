# Differential Structural Profiles of the Hippocampus and Thalamus in Unilateral MTLE

This repository contains the R scripts and preprocessed datasets used for the analyses reported in our study on hippocampal and thalamic structural abnormalities and their cognitive relevance in unilateral mesial temporal lobe epilepsy (MTLE).

## Overview

The repository is organized to support the main analytical components of the study, including:

1. High-resolution structural abnormality mapping
2. Structure–cognition analyses aligned to epileptogenic laterality
3. Additive models of hippocampal and thalamic volume measures
4. Clinical and demographic data organization

The goal of this repository is to improve transparency and reproducibility by providing the processed data and R scripts required to reproduce the main results, figures, and tables.

## Repository contents

### R scripts
- **High-resolution structural abnormality mapping**  
  R code for case-versus-control structural comparisons across cortical, subcortical, hippocampal subfield, and thalamic nucleus measures.

- **Structure_Cognition_EpileptogenicLaterality_integrated**  
  R code for structure–cognition analyses aligned to epileptogenic laterality, including ipsilateral/contralateral recoding and partial correlation analyses.

- **Additive models of hippocampal and thalamic volume measures**  
  R code for hierarchical/additive modeling to compare the explanatory value of hippocampal and thalamic volumetric measures for cognitive outcomes.

- **clinical information**  
  R code for organizing and summarizing demographic, clinical, and neuropsychological information.

### Preprocessed datasets
- **df_raw_clean**  
  Preprocessed patient dataset used for the main analyses.

- **df_hc_clean**  
  Preprocessed healthy control dataset used as the reference group.

## Data description

The shared datasets are preprocessed versions of the study data and were prepared for statistical analysis in R.  
They include the variables necessary to reproduce the analyses presented in the manuscript, such as:

- demographic information
- clinical variables
- neuropsychological measures
- hippocampal subfield volumes
- thalamic nucleus volumes
- other derived structural metrics used in the study

## Reproducibility

The scripts are organized according to the main analysis modules.  
Users may run the scripts independently or follow the workflow below:

1. Load the preprocessed datasets (`df_raw_clean` and `df_hc_clean`)
2. Run the structural abnormality mapping analyses
3. Run the structure–cognition analyses aligned to epileptogenic laterality
4. Run the additive models
5. Export tables and figures as needed

Because local file paths may differ across systems, users may need to modify input and output directories before running the scripts.

## Software requirements

The analyses were conducted in **R**.  
Required packages are specified within the scripts and may include commonly used packages such as:

- tidyverse
- dplyr
- tidyr
- ggplot2
- readr
- stringr
- purrr
- flextable
- officer

Please install any missing packages before running the scripts.

## Notes

- These files are provided to support transparency and reproducibility of the reported analyses.
- The repository includes **preprocessed data only**.
- Users should ensure that any shared data comply with relevant ethical and privacy requirements.
- If directory structures differ from those used in the original analysis environment, file paths should be updated accordingly.

## Citation

If you use this repository, please cite the associated manuscript:

Differential Structural Profiles of the Hippocampus and Thalamus and Their Cognitive Relevance in Unilateral Mesial Temporal Lobe Epilepsy


**[Your name]**  
**[Your institution]**  
**[Your email]**

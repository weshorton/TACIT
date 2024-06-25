# Spatial Deconvolution of Cell Types and Cell States at Scale Utilizing TACIT

![Overview](https://github.com/huynhkl953/TACITomics/blob/main/image/overview.png)

## Table of Contents
- [Introduction](#introduction)
- [Highlights](#highlights)
- [Installation](#installation)
- [Dependency](#dependency)
- [Input](#input)
- [Output](#output)
- [Usage](#usage)
- [License](#license)
- [Help](#Help)
- [Citation](#citation)
## Introduction
Identifying cell types and states remains a time-consuming and error-prone challenge for spatial biology. While deep learning is increasingly used, it is difficult to generalize due to variability at the level of cells, neighborhoods, and niches in health and disease. To address this, we developed TACIT, an unsupervised algorithm for cell annotation using predefined signatures that operates without training data, using unbiased thresholding to distinguish positive cells from background, focusing on relevant markers to identify ambiguous cells in multiomics assays.

Using five datasets (5,000,000-cells; 50-cell types) from three niches (colon, intestine, gland), TACIT outperformed existing unsupervised methods in accuracy and scalability. Integration of TACIT-identified cells with a novel Shiny app revealed new phenotypes in two inflammatory gland diseases. Finally, using combined spatial transcriptomics and proteomics, we discover under- and overrepresented immune cell types and states in regions of interest, suggesting multimodality is essential for translating spatial biology to clinical applications.

Paper: 

## Highlights
1. Disease-agnostic, tissue-agnostic, assay-agnostic, modality-agnostic, and species-agnostic tool for scalable single-cell spatial biological analyses.
2. Outperformance to best-in-class toolkits for cell type identification and cell state analyses, even when cell segmentation errors are especially challenging in highly inflamed regions of interest.
3. CELLxFEATURE matrix construction of both mRNA and protein of unlinked assays (Akoya’s Phenocycler-Fusion 2.0 and 10x Xenium) using mask transfer and multimodal cell segmentation.
4. Detailed analyses to understand the impact of mRNA only, protein only, and mRNA and protein (mixed) CELLxFEATURE matrices on cell identification and cell state analyses.

## Installation

### Hardware requirement

CPU: i5

Memory: 16G or more

Storage: 10GB or more

### Software requirement

R version: >= 4.0 suggested

### Install

You can install the development version of TACIT:

```R
# install.packages("devtools")
devtools::install_github("huynhkl953/TACIT")
```

## Dependency
TACIT requires dependency on the following R packages:

- Seurat (>= version 5.0.3): for performing microclustering
- class (>= version 7.3-22): for k-Nearest Neighbour Classification
- segmented (>= version 2.0-3): for identifying the breakpoints

## Input
TACIT requires two inputs:
1. CELLxFEATURE matrix: 
Features like probe intensity (protein antibodies) and count values (mRNA probes) are quantified, normalized, and stored in a CELLxFEATURE matrix.![expression](https://github.com/huynhkl953/TACITomics/blob/main/image/expression.png)
2. TYPExMARKER matrix: 
The TYPExMARKER matrix is derived from expert knowledge, with values between 0 and 1, indicating the relevance of markers for defining cell types. ![signature](https://github.com/huynhkl953/TACITomics/blob/main/image/signature.png)

Parameters:
1. **r (resolution)**: Depends on the data but aims to create microcluster cell communities with sizes averaging between 0.1–0.5% of cells per microcluster. The higher the resolution, the greater the number of microclusters.
2. **p (dimension)**: Number of dimensions used for microclusters.

## Output

TACIT Outputs:
- **final_threshold**: Dataframe containing the threshold for each cell type relevant score.
- **Group_threshold_data**: Binary matrix using the threshold values.
- **data_clean_final**: Dataframe of final annotations.
- **clusters**: Index of microclusters for each cell.


## Usage
```R
library(Seurat)
library(class)
library(segmented)
library(readr)

CELLxFEATURE=read_csv("CELLxFEATURE.csv")
TYPExMARKER=read_csv("TYPExMARKER.csv")

TACIT=TACIT(CELLxFEATURE,TYPExMARKER,r=10,p=10)
print(TACIT)
```


## Help
If you encounter a clear bug, please file an issue with a minimal reproducible example on GitHub. For questions and other discussion, please email: huynhk4@vcu.edu.

## Citation
Spatial Deconvolution of Cell Types and Cell States at Scale Utilizing TACIT https://www.biorxiv.org/content/10.1101/2024.05.31.596861v1.abstract.



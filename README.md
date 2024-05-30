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
- [Contributing](#contributing)
- [License](#license)

## Introduction
Identifying cell types and states remains a time-consuming and error-prone challenge for spatial biology. While deep learning is increasingly used, it is difficult to generalize due to variability at the level of cells, neighborhoods, and niches in health and disease. To address this, we developed TACIT, an unsupervised algorithm for cell annotation using predefined signatures that operates without training data, using unbiased thresholding to distinguish positive cells from background, focusing on relevant markers to identify ambiguous cells in multiomics assays.

Using five datasets (4,600,000-cells; 51-cell types) from three niches (colon, intestine, gland), TACIT outperformed existing unsupervised methods in accuracy and scalability. Integration of TACIT-identified cells with a novel Shiny app revealed new phenotypes in two inflammatory gland diseases. Finally, using combined spatial transcriptomics and proteomics, we discover under- and overrepresented immune cell types and states in regions of interest, suggesting multimodality is essential for translating spatial biology to clinical applications.

Paper: 

## Highlights
1. Disease-agnostic, tissue-agnostic, assay-agnostic, modality-agnostic, and species-agnostic tool for scalable single-cell spatial biological analyses.
2. Outperformance to best-in-class toolkits for cell type identification and cell state analyses, even when cell segmentation errors are especially challenging in highly inflamed regions of interest.
3. CELLxFEATURE matrix construction of both mRNA and protein of unlinked assays (Akoya’s Phenocycler-Fusion 2.0 and 10x Xenium) using mask transfer and multimodal cell segmentation.
4. Detailed analyses to understand the impact of mRNA only, protein only, and mRNA and protein (mixed) CELLxFEATURE matrices on cell identification and cell state analyses.

## Installation
You can install the development version of CELESTA:

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
Features like probe intensity (protein antibodies) and count values (mRNA probes) are quantified, normalized, and stored in a CELLxFEATURE matrix.(https://github.com/huynhkl953/TACITomics/blob/main/image/expression.png)
2. TYPExMARKER matrix: 
The TYPExMARKER matrix is derived from expert knowledge, with values between 0 and 1, indicating the relevance of markers for defining cell types. (https://github.com/huynhkl953/TACITomics/blob/main/image/signature.png)


## Usage
```R
library(Seurat)
library(class)
library(segmented)
library(readr)

data_expression=read_csv("expression.csv")
Signature_matrix=read_csv("expression.csv")

TACIT=TACIT(data_expression,Signature_matrix,r=10,p=10)
print(TACIT)
```







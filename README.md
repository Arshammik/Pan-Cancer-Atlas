# Project Overview 
### This is a meta analysis on pan cancer focuses on TCGA and GEO data
this readme file will be expanding...
## Study Design:
#### [Flowchart](https://github.com/Arshammik/Lung-Cancer-TCGA-/blob/main/Flowchart.pdf)
<br/>
In this project we will utlize from the GEO data in order to evaluate the candidates gene networks and subsequently train a machine learning model. The model will be tested with TCGA data in order to construct a pan-cancer tool.

## Prerequisites:

* Python (3.7)
* R (4.2.1)
* ggplot2 (3.3.6)
* limma (3.52.2)
* DESeq2 (1.36.0)
* PyTorch (0.4.1)
* NumPy (1.14.3)
* SciPy (1.0.0)
* pandas (0.22.0)

# Data
## Lung Cancer:
### Lung Cancer from GEO 
| GEO accession | Experiment type | Sample Size  | Age Declaration | Gender Declaration  |
| -----         | -----------     | ------------ | --------------  | ---------------     |
| GSE175616     | Microarray      | 123          | Yes             | Yes                 |
| GSE160769     | Microarray      | 30           | No              | Yes                 |
| GSE29016      | Microarray      | 109          | Yes             | Yes                 |
| GSE19027      | Microarray      | 60           | Yes             | Yes                 |
|     ...       |
### Lung Cancer from TCGA — LUAD Project
DEGs evaluated using *DESeq2 package in R* (find the table in `Results/dif.csv`)


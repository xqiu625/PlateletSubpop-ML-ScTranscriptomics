# PlateletSubpop-ML-ScTranscriptomics

[![DOI](https://zenodo.org/badge/DOI/10.3390/ijms25115941.svg)](https://doi.org/10.3390/ijms25115941)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-%3E%3D4.2.0-blue)](https://www.r-project.org/)

This repository contains the analysis code for the paper **"Deciphering Abnormal Platelet Subpopulations in COVID-19, Sepsis, and Systemic Lupus Erythematosus through Machine Learning and Single-Cell Transcriptomics"** (Qiu et al., 2024, *International Journal of Molecular Sciences*).

The computational core of the study: integrate **12 public scRNA-seq cohorts** into a batch-corrected atlas, model platelet state heterogeneity probabilistically, and cast clinical outcome prediction as a **supervised learning problem** solved with gradient boosting and deep neural networks.

## Overview

![Platelet Subpopulations](PlateletSubpop-ML-ScTranscriptomics-overview.jpeg)

### Key Findings

- **Platelet to T Cell Ratio as a Prognostic Biomarker**: The proportion of platelets to T cells in peripheral blood mononuclear cells (PBMC) was identified as the most potent predictor for distinguishing survivors from fatal patients.

- **Discovery of Distinct Platelet Subpopulations**: Identification of different platelet subgroups, including active coagulation, hypoxic, and quiescent clusters, in fatal COVID-19 patients.

- **Key Observations in Severe and Fatal Conditions**:
  - Increased platelet aggregation with monocytes
  - Platelets amplify endothelial dysfunction
  - Reduction in lymphocyte activation and differentiation

---

## Mathematical Methods

### 1. Multi-Cohort Integration via Harmony

Twelve datasets introduce severe batch effects. Let $\mathbf{z}_i \in \mathbb{R}^{d}$ be the PCA embedding of cell $i$ and $\phi_i$ its one-hot batch indicator. Harmony alternates between **soft $k$-means clustering** and **mixture-of-experts batch correction**. Cells are assigned to centroids $\{Y_k\}$ with responsibility

$$R_{ki} = \frac{\exp\left(-d(\mathbf{z}_i, Y_k)/\sigma\right)}{\sum_{k'}\exp\left(-d(\mathbf{z}_i, Y_{k'})/\sigma\right)},$$

and each centroid learns a batch-specific correction $W_k \in \mathbb{R}^{B \times d}$ by ridge regression of $\mathbf{z}_i - Y_k$ onto $\phi_i$ with penalty $\lambda \lVert W_k\rVert^2$. The diversity objective penalizes batch-poor clusters,

$$\max_{\theta}\ \sum_{i,k} R_{ki}\  \phi_i^{\top} \theta_k \quad \text{s.t. cluster entropy regularization},$$

yielding corrected embeddings $\hat{\mathbf{z}}_i = \mathbf{z}_i - \sum_k R_{ki}\  \phi_i^{\top} W_k$ in which biological, not technical, variation dominates.

### 2. Outcome Prediction with XGBoost

Prediction of patient outcome (survivor vs. fatal) is posed as additive function estimation $\hat{y}_i = \sum_{t=1}^{T} f_t(\mathbf{x}_i)$, $f_t \in \mathcal{F}$ (CART space), minimizing the regularized objective

$$\mathcal{L} = \sum_{i=1}^{n} \ell\left(y_i, \hat{y}_i\right) + \sum_{t=1}^{T} \Omega(f_t), \qquad \Omega(f) = \gamma T + \frac{1}{2}\lambda \lVert \mathbf{w} \rVert^2.$$

Each boosting round optimizes a second-order Taylor approximation with gradient/Hessian statistics $g_i = \partial_{\hat{y}}\ell$, $h_i = \partial^2_{\hat{y}}\ell$:

$$\tilde{\mathcal{L}}^{(t)} \simeq \sum_{i=1}^{n}\left[ g_i f_t(\mathbf{x}_i) + \frac{1}{2} h_i f_t^2(\mathbf{x}_i) \right] + \Omega(f_t).$$

For a leaf $j$ with index set $I_j$, writing $G_j = \sum_{i \in I_j} g_i$ and $H_j = \sum_{i \in I_j} h_i$, the optimal weight and the split gain are closed-form:

$$w_j^{*} = -\frac{G_j}{H_j + \lambda}, \qquad \mathrm{Gain} = \frac{1}{2}\left[ \frac{G_L^2}{H_L + \lambda} + \frac{G_R^2}{H_R + \lambda} - \frac{(G_L + G_R)^2}{H_L + H_R + \lambda} \right] - \gamma.$$

### 3. Deep Neural Network Classifier

For multi-class cell-state classification over $C$ states, an $L$-layer network $\mathbf{h}^{(l)} = \mathrm{ReLU}(\mathbf{W}^{(l)}\mathbf{h}^{(l-1)} + \mathbf{b}^{(l)})$ with softmax output $\hat{p}_{ic} = e^{o_{ic}} / \sum_{c'} e^{o_{ic'}}$ is trained by minimizing the categorical cross-entropy

$$\mathcal{L}_{\mathrm{CE}} = -\frac{1}{n}\sum_{i=1}^{n}\sum_{c=1}^{C} y_{ic} \log \hat{p}_{ic},$$

with gradients $\partial \mathcal{L} / \partial \mathbf{o}_i = \hat{\mathbf{p}}_i - \mathbf{y}_i$ propagated by the chain rule.

### 4. Differential Expression: MAST Hurdle Model

Gene-wise differential expression uses a two-part generalized linear **hurdle model**. For gene $g$ in cell $i$ with (log-transformed) expression $y_{ig}$, let $z_{ig} = \mathbb{1}[y_{ig} > 0]$ indicate detection:

$$\mathrm{logit} P(z_{ig} = 1) = \mathbf{x}_i^{\top}\boldsymbol{\beta}_g^{c}, \qquad y_{ig} \mid (z_{ig} = 1) \sim \mathcal{N}\left(\mathbf{x}_i^{\top}\boldsymbol{\beta}_g^{g},\ \sigma_g^2\right).$$

The discrete part models dropout frequency; the continuous part models positive expression level. Wald tests on $\boldsymbol{\beta}^g$ and $\boldsymbol{\beta}^c$ are combined and FDR-controlled.

### 5. Pathway Module Scoring

For a KEGG gene set $\mathcal{M}$, per-cell activity is scored as the mean expression of the module minus a matched background — control genes $\mathcal{C}_b$ drawn from $B$ expression bins:

$$s_i^{(\mathcal{M})} = \frac{1}{|\mathcal{M}|}\sum_{g \in \mathcal{M}} \tilde{x}_{ig} \ -\  \frac{1}{B}\sum_{b=1}^{B} \frac{1}{|\mathcal{C}_b|}\sum_{g \in \mathcal{C}_b} \tilde{x}_{ig}.$$

### 6. Trajectory Inference

Pseudotime is assigned by embedding cells into a principal graph on the UMAP manifold (reversed graph embedding) and projecting each cell to its nearest graph point; pseudotime is the geodesic distance to the chosen root:

$$\tau(i) = d_{\mathcal{G}}\left(\mathrm{proj}_{\mathcal{G}}(\mathbf{u}_i),\ \text{root}\right).$$

### 7. The Prognostic Ratio

The headline biomarker is the scalar feature $r = n_{\text{platelet}} / n_{\text{T}}$ computed per patient from PBMC composition — a one-dimensional predictor whose ROC analysis ($\mathrm{AUC} = P(r_{\text{fatal}} > r_{\text{survivor}})$) outperformed high-dimensional transcriptomic signatures.

---

## Results Summary

### Platelet Cluster Distribution by Disease Severity

| Cluster | Annotation | Key Characteristics |
|---------|------------|---------------------|
| C0-C2 | Active Coagulation | Elevated in fatal COVID-19 (FT) |
| C3 | Quiescent | Predominant in healthy controls (HC) |
| C4 | Hypoxic/Stress | Highly enriched in fatal cases (38.5%) |
| C5-C6 | Intermediate | Variable across conditions |
| C11 | Severe/Fatal-specific | 78.4% frequency in fatal COVID-19 |

### Machine Learning Model Performance

| Model | Accuracy | AUC | Application |
|-------|----------|-----|-------------|
| XGBoost | ~85% | 0.89 | Outcome prediction (Survivor vs Non-survivor) |
| DNN | ~83% | 0.87 | Multi-class cell state classification |

### Key Visualizations

| Figure | Description | Script |
|--------|-------------|--------|
| Fig 1 | Cell type composition and platelet ratios | `Figure1B-G.R`, `Figure1H.R` |
| Fig 2 | Differential expression and pathway analysis | `Figure2A.R` - `Figure2D,E.R` |
| Fig 3 | UMAP clustering visualization | `Figure3A,B.R` |
| Fig 4 | Platelet subpopulation markers | `Figure4A,B.R` - `Figure4E.R` |
| Fig 5 | GSVA pathway enrichment | `Figure5A-F.R`, `Figure5G,H.R` |
| Fig 6 | Trajectory analysis (Monocle3) | `Figure6A,B.R` - `Figure6I.R` |
| Fig 7 | Cross-disease comparison | `Figure7A,B.R`, `Figure7C,D.R` |
| Fig 8 | Statistical comparisons | `Figure8.R` |
| Fig 9 | Module score analysis | `Figure9A.R` - `Figure9F,G.R` |

## Project Structure

```
PlateletSubpop-ML-ScTranscriptomics/
├── scripts/
│   ├── data/                    # Data processing pipelines
│   │   ├── download_data.sh     # Download raw datasets
│   │   ├── process_PMID*.R      # Individual dataset processing (12 scripts)
│   │   ├── integrate_platelet_data_harmony.R  # Harmony batch correction
│   │   ├── seurat_v4_cell_type.R              # Cell type annotation
│   │   └── corrected_gene_names.R             # Gene name standardization
│   │
│   ├── models/                  # Machine learning models
│   │   ├── DNN_model.R          # Deep Neural Network classifier
│   │   └── XGBoost_model.R      # XGBoost gradient boosting classifier
│   │
│   ├── analysis/                # Statistical analysis
│   │   └── Module_score_analysis.R  # KEGG pathway module scoring
│   │
│   └── visualization/           # Figure generation scripts
│       └── Figure*.R            # Scripts for Figures 1-9
│
├── results/
│   └── table/
│       ├── Supp_Table2.csv      # Platelet cluster frequencies
│       └── Supp_Table3.csv      # Additional statistics
│
├── notebooks/                   # Jupyter notebooks (demo)
│   └── ML_Classification_Demo.ipynb
│
├── requirements.txt             # R package dependencies
├── LICENSE
└── README.md
```

## Installation

### Prerequisites

- R >= 4.2.0
- Python >= 3.8 (for TensorFlow/Keras)

### R Package Installation

```r
# Install CRAN packages
install.packages(c("tidyverse", "Matrix", "xgboost", "caret", "keras",
                   "reticulate", "cowplot", "ggpubr", "pheatmap",
                   "RColorBrewer", "VennDiagram", "viridis", "Cairo",
                   "reshape2", "rstatix", "ROCR", "pROC", "HGNChelper"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Seurat", "SeuratDisk", "harmony", "SingleR",
                       "celldex", "scCATCH", "SingleCellExperiment",
                       "scater", "monocle3", "MAST", "sctransform",
                       "org.Hs.eg.db", "GO.db", "KEGGREST",
                       "clusterProfiler", "enrichplot", "DOSE",
                       "AnnotationHub", "GSVA", "GSEABase", "limma",
                       "ComplexHeatmap", "circlize", "EnhancedVolcano",
                       "ggnewscale", "ggupset", "MeSHDbi", "meshes"))

# TensorFlow/Keras setup
library(reticulate)
install_tensorflow()
library(keras)
install_keras()
```

## Usage

### 1. Data Processing

```r
# Download and process individual datasets
source("scripts/data/download_data.sh")
source("scripts/data/process_PMID32514174_scRNAseq.R")

# Integrate datasets using Harmony
source("scripts/data/integrate_platelet_data_harmony.R")

# Annotate cell types
source("scripts/data/seurat_v4_cell_type.R")
```

### 2. Machine Learning Classification

```r
# Train XGBoost classifier
source("scripts/models/XGBoost_model.R")

# Train Deep Neural Network
source("scripts/models/DNN_model.R")
```

### 3. Generate Figures

```r
# Run visualization scripts
source("scripts/visualization/Figure1B-G.R")
# ... repeat for other figures
```

## Dataset Links and Corresponding Papers

| Dataset Link | Paper Title |
|--------------|-------------|
| [GSE150728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728) | [A single-cell atlas of the peripheral immune response in patients with severe COVID-19](https://www.nature.com/articles/s41591-020-0944-y) |
| [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689) | [Immunophenotyping of COVID-19 and influenza highlights role of type I interferons](https://www.science.org/doi/10.1126/sciimmunol.abd1554) |
| [EGAS00001004571](https://ega-archive.org/studies/EGAS00001004571) | [Severe COVID-19 is marked by a dysregulated myeloid cell compartment](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7405822/) |
| [GSE155673](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155673) | [Systems biological assessment of immunity to mild versus severe COVID-19](https://www.science.org/doi/full/10.1126/science.abc6261) |
| [GSE158055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158055) | [COVID-19 immune features revealed by large-scale single-cell transcriptome atlas](https://www.cell.com/cell/fulltext/S0092-8674(21)00148-3) |
| [COVID-19 Cell Atlas](https://www.covid19cellatlas.org/index.patient.html) | [Time-resolved systems immunology reveals a late juncture linked to fatal COVID-19](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7874909/) |
| [E-MTAB-10026](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10026/) | [Single-cell multi-omics analysis of the immune response in COVID-19](https://www.nature.com/articles/s41591-021-01329-2) |
| [GSE151263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151263) | [Single cell RNA sequencing identifies an early monocyte gene signature in ARDS](https://insight.jci.org/articles/view/135678) |
| [Bernardes_2020_COVID19](https://www.fastgenomics.org/platform/) | [Longitudinal Multi-omics Analyses Identify Responses of Megakaryocytes](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7689306/) |
| [GSE163668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163668) | [Global Absence and Targeting of Protective Immune States in Severe COVID-19](https://www.nature.com/articles/s41586-021-03234-7) |
| [GSE167363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363) | [Dynamic changes in human single-cell transcriptional signatures during fatal sepsis](https://jlb.onlinelibrary.wiley.com/doi/10.1002/JLB.5MA0721-825R) |
| [GSE142016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142016) | [Transcriptomic, epigenetic, and functional analyses implicate neutrophil diversity in SLE](https://www.pnas.org/content/116/50/25222) |

## Citation

If you use this code, please cite:

```bibtex
@article{qiu2024platelet,
  title={Deciphering Abnormal Platelet Subpopulations in COVID-19, Sepsis, and
         Systemic Lupus Erythematosus through Machine Learning and
         Single-Cell Transcriptomics},
  author={Qiu, Xinru and others},
  journal={International Journal of Molecular Sciences},
  year={2024},
  volume={25},
  number={11},
  pages={5941},
  doi={10.3390/ijms25115941}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please open a GitHub issue or contact the authors.

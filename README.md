# PlateletSubpop-ML-ScTranscriptomics
[![DOI](https://zenodo.org/badge/DOI/10.3390/ijms25115941.svg)](https://doi.org/10.3390/ijms25115941)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the analysis code for the paper "Deciphering Abnormal Platelet Subpopulations in COVID-19, Sepsis, and Systemic Lupus Erythematosus through Machine Learning and Single-Cell Transcriptomics" (Qiu et al., 2021, Journal of Leukocyte Biology).


## Overview

![Platelet Subpopulations](https://github.com/xqiu625/PlateletSubpop-ML-ScTranscriptomics/assets/26670165/14ac3d98-7811-4b37-a2a7-f0b1037697de)

### Key Findings:
- **Platelet to T Cell Ratio as a Prognostic Biomarker**:
  - The proportion of platelets to T cells in peripheral blood mononuclear cells (PBMC) was identified as the most potent predictor for distinguishing survivors from fatal patients, highlighting its potential as a prognostic biomarker.

- **Discovery of Distinct Platelet Subpopulations**:
  - The identification of different platelet subgroups, including active coagulation, hypoxic, and quiescent clusters, in fatal COVID-19 patients suggests potential targeted treatment strategies.

- **Key Observations in Severe and Fatal Conditions**:
  - **Platelet Aggregation with Monocytes**: Increased platelet aggregation with monocytes was observed in severe and fatal cases.
  - **Amplification of Endothelial Dysfunction**: Platelets were found to amplify endothelial dysfunction, contributing to the severity of the condition.
  - **Reduction in Lymphocyte Activation**: A decrease in lymphocyte activation and differentiation due to platelet activity was noted, indicating a broader role of platelets in inflammatory and immune responses.

# Dataset Links and Corresponding Papers

| Dataset Link                                                                                      | Paper Title and Link                                                                                                                                          |
|---------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------|
| [GSE150728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728)                         | [A single-cell atlas of the peripheral immune response in patients with severe COVID-19](https://www.nature.com/articles/s41591-020-0944-y)                   |
| [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)                         | [Immunophenotyping of COVID-19 and influenza highlights role of type I interferons in development of severe COVID-19](https://www.science.org/doi/10.1126/sciimmunol.abd1554) |
| [EGAS00001004571](https://ega-archive.org/studies/EGAS00001004571) | [Severe COVID-19 is marked by a dysregulated myeloid cell compartment](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7405822/)                                    |
| [GSE155673](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155673)                         | [Systems biological assessment of immunity to mild versus severe COVID-19 infection in humans](https://www.science.org/doi/full/10.1126/science.abc6261)       |
| [GSE158055](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158055)                         | [COVID-19 immune features revealed by a large-scale single-cell transcriptome atlas](https://www.cell.com/cell/fulltext/S0092-8674(21)00148-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867421001483%3Fshowall%3Dtrue) |
| [COVID-19 Cell Atlas](https://www.covid19cellatlas.org/index.patient.html)                        | [Time-resolved systems immunology reveals a late juncture linked to fatal COVID-19](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7874909/)                    |
| [E-MTAB-10026](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-10026/)                      | [Single-cell multi-omics analysis of the immune response in COVID-19](https://www.nature.com/articles/s41591-021-01329-2)                                      |
| [GSE151263](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151263)                         | [Single cell RNA sequencing identifies an early monocyte gene signature in acute respiratory distress syndrome](https://insight.jci.org/articles/view/135678)  |
| [Bernardes_2020_COVID19](https://www.fastgenomics.org/platform/) | [Longitudinal Multi-omics Analyses Identify Responses of Megakaryocytes, Erythroid Cells, and Plasmablasts as Hallmarks of Severe COVID-19](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7689306/) |
| [GSE163668](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163668)                         | [Global Absence and Targeting of Protective Immune States in Severe COVID-19](https://www.nature.com/articles/s41586-021-03234-7)                              |
| [GSE167363](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167363)                         | [Dynamic changes in human single-cell transcriptional signatures during fatal sepsis](https://jlb.onlinelibrary.wiley.com/doi/10.1002/JLB.5MA0721-825R)       |
| [GSE142016](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE142016)                         | [Transcriptomic, epigenetic, and functional analyses implicate neutrophil diversity in the pathogenesis of systemic lupus erythematosus](https://www.pnas.org/content/116/50/25222) |

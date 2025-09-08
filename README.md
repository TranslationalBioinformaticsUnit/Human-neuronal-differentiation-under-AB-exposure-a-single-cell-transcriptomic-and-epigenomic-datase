# 🧠 Adult Hippocampal Neurogenesis Multi-Omics Dataset

This repository provides a time-resolved single-cell multi-omics dataset profiling adult hippocampal neurogenesis (AHN) in vitro, with and without Alzheimer’s disease (AD)-related pathology.

Human iPSC-derived neural progenitor cells (NPCs) expressing XCL1-DCXpGFP were differentiated over four key time points (Day 0, 7, 13, 20). In parallel, a subset of cultures was exposed weekly to physiologically relevant amyloid-β (Aβ) peptide 1-42 to model AD conditions. Single-cell RNA sequencing (scRNA-seq) and ATAC sequencing (scATAC-seq) were performed across all time points and conditions.

Key findings include:

Cell fate trajectories: Clear neuronal and astrocytic differentiation with stage-specific gene signatures.

Aβ-induced effects: Altered gene expression during NPC-to-neuroblast transition and synaptogenesis, without changes in lineage proportions.

Regulatory disruption: scATAC-seq revealed perturbed networks, notably involving transcription factors POU2F2 and ONECUT2, linked to impaired neuronal maturation.

Clinical relevance: Aβ-affected genes overlapped with dysregulation observed in human hippocampal AD samples.

📊 This dataset serves as a reference resource for studying human neuronal differentiation and AD-associated dysregulation, supporting both mechanistic and translational research.


# 📂 Data Availability

All raw and processed data files are publicly available in the Gene Expression Omnibus (GEO):

Accession: [GSE307094](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE307094)

Contents: Raw and processed scRNA-seq and scATAC-seq data


# 🔎 Interactive Exploration

To support rapid and user-friendly data exploration, we provide a Shiny web application for interactive visualization and inspection of the scRNA-seq dataset:

👉 [Launch Shiny App](https://translationalbio.shinyapps.io/neurogenesis/)


# 💻 Code & Reproducibility

All analysis code used to generate the results and figures from the study is available in this repository. This includes:

Preprocessing of scRNA-seq and scATAC-seq data

Integration and trajectory analysis

Differential expression and enhancer gene regulatory network analysis

Visualization scripts for reproducing key figures

🔧 The code is organized into Python and R scripts, enabling full reproducibility of the study’s findings. Users can adapt these workflows for their own single-cell datasets.

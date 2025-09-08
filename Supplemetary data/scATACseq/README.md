# 🧬 scATAC-seq Fragments & Peaks

This repository includes scATAC-seq fragment and peak files capturing chromatin accessibility dynamics during human adult hippocampal neurogenesis (AHN) across developmental stages, for both control and amyloid-β (Aβ)-treated conditions. These files provide the foundation for analyzing regulatory landscapes and integrating with scRNA-seq data.


# 📌 Experimental Conditions

Control: Untreated neural progenitor cells (NPCs)

Aβ-treated: NPCs exposed weekly to physiologically relevant amyloid-β 1-42


# 🗓 Developmental Time Points

Day 0: iPSC/NPC baseline

Day 7: Early differentiation

Day 13: NPC-to-neuroblast transition

Day 20: Neuronal maturation


# 📂 File Contents

[Fragment files](https://figshare.com/s/576414a1251ad58742de) (.tsv.gz):

Chromosome, start and end positions of fragments

Cell barcodes linking fragments to individual cells

Fragment counts per cell

[Peak files](https://figshare.com/s/8e8fd4fc4e11341a3518) (.bed or .tsv):

Called accessible regions across cells

Can be used for downstream analyses like gene regulatory network inference, trajectory mapping, or integration with scRNA-seq


# ⚡ Usage

These files are compatible with popular scATAC-seq analysis tools, including ArchR, Signac, SnapATAC, and Cicero, enabling reproducible exploration of chromatin dynamics during neuronal differentiation and AD-associated perturbations.

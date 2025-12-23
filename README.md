# pcSecYeastSpecies: Cross-species proteome-constrained modeling of yeast protein secretion

This repository contains the code and models used in the study:

**Cross-species proteome-constrained modeling reveals trade-offs in yeast protein secretion under temperature and glycosylation stress**

We built on the pcSecYeast model for *Saccharomyces cerevisiae* and developed proteome-constrained protein secretion models for *Komagataella phaffii* and *Kluyveromyces marxianus*. By integrating genome-scale metabolism, detailed secretory pathway representations, proteome allocation constraints, temperature-dependent enzyme kinetics, and humanized glycosylation modules, this framework enables quantitative, cross-species analysis of secretion capacity, metabolic trade-offs, and stress responses under industrially relevant conditions.

---

## Model Construction

Proteome-constrained secretion models were constructed from a common template model using species-specific build scripts.  
Models for *S. cerevisiae*, *K. phaffii*, and *K. marxianus* are generated using `buildModel_pcSecYeast.m`, `buildModel_pcSecPichia.m`, and `buildModel_pcSecKmarx.m`, respectively.

---

## Simulation and Analysis

All simulation scripts required to reproduce the analyses in the manuscript are included in the corresponding species folders under `Code/`.  
Scripts to reproduce all manuscript figures are located in `Code/Figures/`. Each figure script loads the processed results and generates the corresponding plots.  
Simulation outputs are saved in the `Results/` directory.

---

## Installation

### Required software

- MATLAB (R2020b or later)
- COBRA Toolbox for MATLAB
- RAVEN Toolbox
- SoPlex solver

---

## Contact

**Lizheng Liu** (GitHub: @Zephyr-112), Institute of Biopharmaceutical and Health Engineering, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China  

**Feiran Li** (GitHub: @feiranl), Institute of Biopharmaceutical and Health Engineering, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China


# pcSecYeastSpecies: Cross-species proteome-constrained modeling of yeast protein secretion

This repository contains the code and models used in the study:

**Cross-species proteome-constrained modeling reveals trade-offs in yeast protein secretion under temperature and glycosylation stress**

We built on the pcSecYeast model for *Saccharomyces cerevisiae* and developed proteome-constrained protein secretion models for *Komagataella phaffii* and *Kluyveromyces marxianus*. By integrating genome-scale metabolism, detailed secretory pathway representations, proteome allocation constraints, temperature-dependent enzyme kinetics, and humanized glycosylation modules, this framework enables quantitative, cross-species analysis of secretion capacity, metabolic trade-offs, and stress responses under industrially relevant conditions.

---

## Model Construction

Proteome-constrained secretion models were constructed from a common template model using species-specific build scripts.  
Models for *S. cerevisiae*, *K. phaffii*, and *K. marxianus* are generated using  [`buildModel_pcSecYeast.m`](Code/pcSecYeast/buildModel_pcSecYeast.m),  [`buildModel_pcSecPichia.m`](Code/pcSecPichia/buildModel_pcSecPichia.m), and  [`buildModel_pcSecKmarx.m`](Code/pcSecKmarx/buildModel_pcSecKmarx.m), respectively.

---

## Simulation and Analysis

All simulation scripts required to reproduce the analyses in the manuscript are provided in the corresponding species-specific folders under [`Code/`](Code/).  
Scripts to reproduce all manuscript figures are located in [`Code/Figures/`](Code/Figures/); each figure script loads the processed results and generates the corresponding plots.  
Simulation outputs are saved in the [`Results/`](Results/) directory.

---

## Installation

### Required software

- MATLAB (R2020b or later)
- [COBRA Toolbox for MATLAB](https://github.com/opencobra/cobratoolbox)
- [RAVEN Toolbox](https://github.com/SysBioChalmers/RAVEN)
- solver [SoPlex](https://soplex.zib.de/)

Please ensure that all required toolboxes and solvers are properly installed and added to the MATLAB path before running the scripts.
---

## Contact

**Lizheng Liu** ([GitHub: @Zephyr-112](https://github.com/Zephyr-112)), Institute of Biopharmaceutical and Health Engineering, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China  

**Feiran Li** ([GitHub: @feiranl](https://github.com/feiranl)), Institute of Biopharmaceutical and Health Engineering, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China


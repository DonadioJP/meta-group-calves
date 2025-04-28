# Meta-Analysis: Effects of Early Group Housing on Dairy Calves

<div align="center">
  <img src="figures/groupcalves.jpg" alt="Calves" width="400"/>
</div>


## Introduction

This repository contains data, scripts, and supplementary materials for the meta-analysis:  
**"A meta-analysis approach to evaluate the effects of early group housing on calf performance, health, and behavior during the preweaning period"**  
Published in *Journal of Dairy Science* (2025). [DOI: 10.3168/jds.2024-25159](https://doi.org/10.3168/jds.2024-25159)

---

## 📋 Study Overview

### Background
Dairy calves are traditionally separated from dams shortly after birth and raised individually during the preweaning period. While this practice is common, it contrasts with natural social dynamics and raises welfare concerns. Social isolation can lead to stress, fearfulness, and reduced cognitive development. Group housing (pairs or small groups) has emerged as a potential alternative, but its effects on performance, health, and behavior remain inconsistently reported across studies. This meta-analysis synthesizes evidence from 51 articles (85 studies) to evaluate the impact of early group housing on dairy calves.

### Key Findings
- **Performance**:  
  🟢 Group-housed calves had **higher average daily gain (ADG)** (+0.06 kg/d, *P* = 0.001) and **weaning weight** (+1.44 kg, *P* = 0.037).  
  🟢 Increased **concentrate intake** in group housing (*P* = 0.021).  
- **Behavior**:  
  🟢 More **active behaviors** (feeding, playing) and fewer **stress-related behaviors** (self-grooming, pen interactions) in group-housed calves.  
  🟢 Group-housed calves vocalized more in novel environments but interacted less with humans.  
- **Health**:  
  ➖ No significant differences in blood parameters (e.g., glucose, TNF-α).  
  ⚠️ Limited data on robust health outcomes (e.g., disease incidence).  

---

## 📂 Repository Structure
├── data/

├── scripts/

├── figures/

├── manuscript/

├── README.md

└── LICENSE

## 🛠️ Tools used
1. **Data**: Raw data from 51 included studies is in `data/meta_analysis_data.xlsx`.
2. **Scripts**:
   - Run `meta_analysis.R` to reproduce Tables 2-4 and Figures 2-4.
   - Use `forest_plots.R` for visualization (STATA scripts provided for compatibility).
3. **Dependencies**: ![R](https://img.shields.io/badge/R-276DC3?logo=r&logoColor=white) ![metafor](https://img.shields.io/badge/metafor-4.2--0-8B9DC3) ![dplyr](https://img.shields.io/badge/dplyr-1.1.0-1D6F42) ![ggplot2](https://img.shields.io/badge/ggplot2-3.4.0-3A7CB8) ![Excel](https://img.shields.io/badge/Excel-217346?logo=microsoftexcel&logoColor=white)

## ✉️ Contact
**Corresponding Author**:  
Prof. Matheus Deniz
Email: [m.deniz@unesp.br](mailto:m.deniz@unesp.br)  
Affiliation: Universidade Estadual Paulista, Brazil

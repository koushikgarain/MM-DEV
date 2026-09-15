# Data for MM-DEV

This directory contains the simulated, empirical, and experimental datasets used in the analyses presented in the manuscript:

**“Early warning signals to anticipate critical transitions in high-dimensional systems”**

The datasets are grouped into four categories:

1. Ricker model simulations
2. Lotka–Volterra model simulations
3. Lake Zurich empirical ecological data
4. Experimental methane-producing bioreactor data

The corresponding R scripts for analyzing these datasets and reproducing the figures are provided in the main directory of this repository.

---

## 1. Ricker Model Data

### Files

- `500RemoveSeq1_RK21.csv`
- `Ricker21_1000.csv`

These files contain simulated multivariate time-series data generated from a high-dimensional Ricker model. The simulations used in the present study were generated in **R**.

The model formulation and simulation framework follow the approach described in:

> Grziwotz, F., Chang, C.-W., Dakos, V., van Nes, E. H., Schwarzländer, M., Kamps, O., Heßler, M., Tokuda, I. T., Telschow, A., & Hsieh, C.-h. (2023). **Anticipating the occurrence and type of critical transitions.** *Science Advances*, 9(1), eabq4558.  
> DOI: 10.1126/sciadv.abq4558

The data and R analytical procedures associated with the original study are publicly archived on Zenodo:

**Zenodo DOI:** 10.5281/zenodo.7379358

The Ricker datasets included in the present repository were generated for the MM-DEV analysis and are used to evaluate the ability of MM-DEV to anticipate critical transitions in high-dimensional model-generated systems.

---

## 2. Lotka–Volterra Model Data

### Files

- `lot_vol_6_1_data.mat`
- `lot_vol_6_2_data.mat`
- `lot_vol_6_1_Isomap.csv`
- `lot_vol_6_2_Isomap.csv`

These files contain simulated data generated from Lotka–Volterra-type ecological systems using **MATLAB**.

The model construction and simulation approach used in the present study follow the framework described in:

> Lever, J. J., Van Nes, E. H., Scheffer, M., & Bascompte, J. (2023). **Five fundamental ways in which complex food webs may spiral out of control.** *Ecology Letters*, 26(10), 1765–1779.  
> DOI: 10.1111/ele.14293

The code associated with the original study is publicly archived on Zenodo:

**Zenodo DOI:** 10.5281/zenodo.8132025

The `.mat` files contain model-generated data used in the present study. The corresponding `*_Isomap.csv` files contain the datasets used for the ISOMAP/manifold-learning component of the MM-DEV analysis.

These simulations were generated for the present analysis following the modeling framework of the above study rather than being experimental observations.

---

## 3. Lake Zurich Data

### File

- `lake_zurich.csv`

This file contains long-term monthly observations of the plankton community and environmental conditions in **Lake Zurich, Switzerland**.

The dataset contains time series of **13 plankton functional groups**, together with **phosphate concentration (PO4)** and **water temperature**, covering the period from **January 1978 to December 2019**.

The Lake Zurich dataset used in the present study is associated with the data and analyses reported in the following studies.

### Original ecological dataset

> Merz, E., Saberski, E., Gilarranz, L. J., Isles, P. D. F., Sugihara, G., Berger, C., & Pomati, F. (2023). **Disruption of ecological networks in lakes by climate change and nutrient fluctuations.** *Nature Climate Change*, 13, 389–396.  
> DOI: 10.1038/s41558-023-01615-6

The original study investigated long-term plankton community dynamics and environmental changes across Swiss lakes, including Lake Zurich.

The associated ecological data are publicly available through an open-access data repository:

**Data repository DOI:** 10.25678/0007VX

The Lake Zurich observations were collected as part of long-term ecological monitoring, and the datasets analyzed in the above study were assembled by the authors.

### Previous dynamical analysis of Lake Zurich

The Lake Zurich time series was also analyzed in:

> Medeiros, L. P., Sorenson, D. K., Johnson, B. J., Palkovacs, E. P., & Munch, S. B. (2025). **Revealing unseen dynamical regimes of ecosystems from population time-series data.** *Proceedings of the National Academy of Sciences*, 122(24), e2416637122.  
> DOI: 10.1073/pnas.2416637122

This study used the long-term Lake Zurich time series to investigate dynamical regimes and ecological transitions from population time-series data.

### Use in the present study

In the present MM-DEV study, `lake_zurich.csv` is used as a high-dimensional empirical ecological case study for evaluating changes in dynamical stability and the ability of MM-DEV to provide early-warning information associated with ecological regime shifts.

The use and interpretation of the Lake Zurich dataset in the present study are informed by **both Merz et al. (2023) and Medeiros et al. (2025)**.

Users of the Lake Zurich data should cite the original data source and the relevant previous studies, as appropriate.

---

## 4. Experimental Bioreactor Data

### Files

- `KTU_biorec_quan_D4.csv`
- `UASB_daughter4_reactors.xlsx`
- `UASB_daughter_reactors_environmental_parameter_table_first110days_20241230update_Expertise.xlsx`
- `UASB_daughter_reactors_environmental_parameter_table_first110days_20241230update_Navie.xlsx`

These files contain experimental data obtained from methane-producing anaerobic bioreactor systems.

The experimental data were generated by members of the research team and collaborators. The researchers responsible for the experimental work are included as co-authors of the associated MM-DEV manuscript.

The related experimental system, microbial community measurements, reactor-performance measurements, and methodology are described in:

> Chang, C.-J. et al. (2025). **Bioenergetically constrained dynamical microbial interactions govern the performance and stability of methane-producing bioreactors.** *npj Biofilms and Microbiomes*, 11, Article 31.  
> DOI: 10.1038/s41522-025-00668-z

The 16S rRNA gene amplicon sequencing data associated with the experimental study are publicly available through the **NCBI Sequence Read Archive (SRA)** under:

**BioProject:** PRJNA1047765

The reactor-performance metadata and analytical resources associated with the original study are also publicly available as described in the corresponding publication.

The bioreactor datasets included in this repository are used to evaluate MM-DEV on experimental high-dimensional microbial systems.

---

## Data Provenance Summary

| Dataset | Type | Source/Generation | Primary References |
|---|---|---|---|
| Ricker model | Simulated | Generated in R for the present analysis using a previously published Ricker-model framework | Grziwotz et al. (2023) |
| Lotka–Volterra model | Simulated | Generated in MATLAB for the present analysis using a previously published food-web modeling framework | Lever et al. (2023) |
| Lake Zurich | Empirical | Long-term publicly available ecological observations | Merz et al. (2023); Medeiros et al. (2025) |
| Bioreactor | Experimental | Generated by members of the research team and collaborators | Chang et al. (2025) |

---

## Data Reuse and Citation

The simulated datasets in this repository are provided to facilitate reproduction of the analyses and figures presented in the associated MM-DEV manuscript.

For empirical datasets originating from or associated with previously published studies, users should cite the corresponding original publications and data sources in addition to citing the MM-DEV repository.

For the Lake Zurich dataset, users should refer to **Merz et al. (2023)** and **Medeiros et al. (2025)** and cite the original data repository where appropriate.

For the bioreactor datasets, users should cite the corresponding experimental study by **Chang et al. (2025)** and the MM-DEV repository when using the processed datasets provided here.

---

## Reproducibility

The analysis scripts corresponding to these datasets are provided in the main directory of the MM-DEV GitHub repository:

`https://github.com/koushikgarain/MM-DEV`

The scripts reproduce the model analyses and figures presented in the associated manuscript.

For additional information about running the analyses, software dependencies, and citation of the MM-DEV method, please see the main repository `README.md` and `CITATION.cff`.
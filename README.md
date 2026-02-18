# Blood Metabolic Ages Are Associated with the Progression of Central Biomarkers of Alzheimer's Disease

Tianchuan Gao<sup>1†</sup>, Pradeep Varathan Pugalenthi<sup>1†</sup>, Yen-Ning Huang<sup>2†</sup>, Matthias Arnold<sup>3,4</sup>, Rima Kaddurah-Daouk<sup>3,5,6</sup>, Andrew J Saykin<sup>2</sup>, Jingwen Yan<sup>1,2\*</sup>, Kwangsik Nho<sup>2\*</sup>

<sup>1</sup> Department of Biomedical Engineering and Informatics, Indiana University Indianapolis, Indianapolis, IN, USA
<sup>2</sup> Department of Radiology and Imaging Sciences, Indiana University School of Medicine, Indianapolis, IN, USA
<sup>3</sup> Department of Psychiatry and Behavioral Sciences, Duke University, Durham, NC, USA
<sup>4</sup> Institute of Computational Biology, Helmholtz Zentrum München, German Research Center for Environmental Health, Neuherberg, Germany
<sup>5</sup> Duke Institute of Brain Sciences, Duke University, Durham, NC, USA
<sup>6</sup> Department of Medicine, Duke University, Durham, NC, USA

<sup>†</sup> Equal contribution
<sup>\*</sup> Corresponding authors: jingyan@iu.edu; knho@iu.edu

---

## Overview

This repository contains the analysis code for MetaboAge, a blood-based metabolic age prediction model trained on NMR Nightingale metabolomics data from 67,964 healthy UK Biobank participants. The trained model is applied to 1,263 ADNI participants (CN, MCI, AD) to compute MetaboAge delta (predicted metabolic age − chronological age). Three downstream analyses examine the association of MetaboAge delta with A/T/N biomarkers (CSF, FDG-PET, Amyloid PET, MRI) and cognitive composite scores, with sex-stratified analyses throughout.

---

## Code

**`age_prediction.ipynb`**
Trains sex-stratified MLP models (whole, male, female) on UK Biobank data with polynomial bias correction (degree = 3), applies them to ADNI, and reports prediction performance. Outputs predicted age CSV files used by all downstream analyses.

**`cross_sectional_association.R`**
Tests cross-sectional associations between MetaboAge delta and A/T/N biomarkers and cognitive scores at baseline using linear regression (AssociatoR package), stratified by sex. All p-values are FDR-corrected. Outputs a summary results table.

**`cross_sectional_analysis.R`**
Generates scatter plots of MetaboAge delta vs. each significant biomarker with regression overlay, using mixed-effects model population-level predictions stratified by MetaboAge delta sign (Figure 3).

**`longitudinal_analysis.R`**
Fits linear mixed-effects models to test whether baseline MetaboAge delta moderates the longitudinal rate of change in biomarkers and cognitive scores over up to 8 years (delta_age × year interaction). Generates trajectory plots (Figure 4).

**`survival_analysis.R`**
Kaplan-Meier survival analysis testing whether MetaboAge delta (above vs. at-or-below median) predicts time to cognitive progression under three scenarios (CN/EMCI→AD, CN→MCI/AD, MCI→AD), stratified by APOE ε4 carrier status (Figure 6).

---

## Dependencies

**Python** (`age_prediction.ipynb`): `numpy`, `pandas`, `matplotlib`, `seaborn`, `scipy`, `scikit-learn`

**R** (all `.R` scripts):

```r
install.packages(c("survival", "survminer", "ggsurvfit", "lme4",
                   "ggpubr", "ggtext", "hrbrthemes", "dplyr", "rempsyc", "readxl"))
devtools::install_github("PradoVarathan/AssociatoR")
```

---

## Acknowledgements

Data used in preparation of this article were obtained from the Alzheimer's Disease Metabolomics Consortium (ADMC; https://sites.duke.edu/adnimetab/) and the Alzheimer's Disease Neuroimaging Initiative (ADNI; adni.loni.usc.edu). This research was supported by NIH grants R01 AG081951, U01AG068057, U19AG074879, R01 LM012535, and NSF CAREER award 1942394, among others.

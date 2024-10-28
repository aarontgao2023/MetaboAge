# Blood metabolic ages are associated with the progression of central biomarkers of Alzheimer’s disease
Pradeep Varathan Pugalenthi<sup>1</sup>, Tianchuan Gao<sup>1</sup>, Yen-ning Huang<sup>2</sup>, Matthias Arnold3<sup>3,4</sup>, Rima Kaddurah-Daouk<sup>3,5,6</sup>, Andrew J Saykin<sup>2</sup>, Jingwen Yan<sup>1,2* </sup>, Kwangsik Nho<sup>2* </sup>
<br>
<sup>1 </sup> Department of Biomedical Engineering and Informatics, Indiana University Indianapolis, Indianapolis, IN, USA<br>
<sup>2 </sup> Department of Radiology and Imaging Sciences, Indiana University School of Medicine, Indianapolis, IN, USA <br>
<sup>3 </sup> Department of Psychiatry and Behavioral Sciences, Duke University, Durham, NC, USA<br>
<sup>4 </sup> Institute of Computational Biology, Helmholtz Zentrum München, German Research Center for Environmental Health, Neuherberg, Germany<br>
<sup>5 </sup> Duke Institute of Brain Sciences, Duke University, Durham, NC, USA<br>
<sup>6 </sup> Department of Medicine, Duke University, Durham, NC, USA<br>

## Abstract
![](Model_Desc1.png?raw=true)
Alzheimer's disease (AD) is characterized by the accumulation of amyloid plaques and tau tangles, resulting in gradual memory loss, with aging being a key contributing factor. Biological age, which reflects an individual's physiological function level, can be estimated using -omics markers. It has been suggested that AD patients might display signs of accelerated aging, leading to a notable gap between their chronological and biological ages. In this study, we aim to develop a blood-based metabolic age prediction model (MetaboAge) to assess metabolic changes associated with normal aging and to investigate the potential of MetaboAge delta (the difference between biological and chronological ages) as an indicator of accelerated metabolic aging related to AD pathology and neurodegeneration. We identified a cluster of metabolic parameters linked to aging, with variation observed between males and females. Notably, metabolite levels in females showed stronger associations with aging. Furthermore, we found a significant association between MetaboAge delta and longitudinal changes in memory and brain structure, highlighting its potential for monitoring disease progression or treatment response in AD.

## Code Inventory
1. [Age_Prediction.ipynb](Age_Prediction.ipynb) --> File descirption of the original code for deployment of MetaboAge and AgeBias training on UKBB <br>

2. [CrossSectional_Analysis.R](CrossSectional_Analysis.R) --> Association analysis using the original package AssociatoR (https://github.com/PradoVarathan/AssociatoR)

3. [Longitudinal_analysis.R](Longitudinal_analysis.R) --> Association with longitudinal changes in AD biomarkers and cognitive performances.

4. [SurvivalAnalysis.R](SurvivalAnalysis.R) --> Survival analysis between multiple groups to test metaboAge,  genetic effects in AD and the progression of AD.


## Contact
Tianchuan Gao(tg11@iu.edu)<br>
Jingwen Yan(jingyan@iu.edu)<br>
Kwangsik Nho(knho@iu.edu)

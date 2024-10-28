# This file is to run cross-sectional analysis

library(AssociatoR)
library(rempsyc)
library(dplyr)
#Cross Sectional
prep_factors = function(data, add_on){
  if(add_on == ""){
    data$Sex = as.factor(data$Sex)
  }
  data$APOE4 = as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == 'VISCODE2')] = 'timepoint'
  data$delta_age = data$Corrected_Predicted_Age - data$Actual_Age
  return(data)

}


add_on = ""
all_p_vals = data.frame('Phenotype' = character(),'Population' = character(),"Delta_Age_pvalue" = numeric(),"Beta" = numeric(),"Std_Error" = numeric())
for(add_on in c("","__FEMALE","__MALE")){
  mri_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/MRI_ALL_ADNI',add_on,'.csv'))
  csf_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/CSF_ALL_ADNI',add_on,'.csv'))
  plasma_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/Plasma_ALL_ADNI',add_on,'.csv'))
  fdg_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/FDG_ALL_ADNI',add_on,'.csv'))
  cogn_cross = readRDS(paste0('Data/CrossSectionalDataPrepared/Cognitivie_Scores',add_on,'_ALL.RDS'))



  mri_cross = prep_factors(mri_cross, add_on)
  if(add_on == ""){
    formula = 'delta_age + Sex + APOE4 + Education + scale( ICV ) + Actual_Age + BMI'
  }else{
    formula = 'delta_age + APOE4 + Education + scale( ICV ) + Actual_Age + BMI'
  }
  mri_p_vals = run_cross_sectional_timepoints(data = mri_cross,
                                              targets = c('HippVol','entorhinal_thickness'),
                                              timepoints = c('bl','m12','m24'),
                                              formula = formula)
  mri_p_vals = mri_p_vals[which(mri_p_vals$timepoint == 'bl'),]
  all_p_vals = rbind(all_p_vals,
                     data.frame('Phenotype' = rownames(mri_p_vals),
                                'Population' = rep(add_on,nrow(mri_p_vals)),
                                "Delta_Age_pvalue" = mri_p_vals$delta_age_p_val,"Beta" = mri_p_vals$delta_age_Coef,"Std_Error" =  mri_p_vals$delta_age_Var_Exp))
                     
  #write.csv(mri_p_vals,paste0('CrossSectionalResults/MRI_ROI',add_on,'.csv'), quote = F)

  csf_cross = prep_factors(csf_cross, add_on)
  csf_cross$ABETA = log(csf_cross$ABETA)
  csf_cross$TAU = log(csf_cross$TAU)
  csf_cross$PTAU = log(csf_cross$PTAU)
  if(add_on == ""){
    formula = 'delta_age + Sex + APOE4 + Actual_Age + BMI'
  }else{
    formula = 'delta_age + APOE4 + Actual_Age + BMI'
  }  
  csf_p_vals = run_cross_sectional_timepoints(data = csf_cross,
                                              targets = c('ABETA','PTAU'),
                                              timepoints = c('bl','m12','m24'),
                                              formula = formula)
  
  csf_p_vals = csf_p_vals[which(csf_p_vals$timepoint == 'bl'),]
  all_p_vals = rbind(all_p_vals,
                     data.frame('Phenotype' = rownames(csf_p_vals),
                                'Population' = rep(add_on,nrow(csf_p_vals)),
                                "Delta_Age_pvalue" = csf_p_vals$delta_age_p_val,"Beta" = csf_p_vals$delta_age_Coef,"Std_Error" =  csf_p_vals$delta_age_Var_Exp))
  #write.csv(csf_p_vals,paste0('CrossSectionalResults/CSF',add_on,'.csv'), quote = F)



  # formula = 'delta_age + APOE4 + Actual_Age'
  # plasma_cross = prep_factors(plasma_cross, add_on)
  # plasma_p_vals = run_cross_sectional_timepoints(data = plasma_cross,
  #                                             targets = c('AB40','AB42'),
  #                                             timepoints = c('bl','m12','m24'),
  #                                             formula = formula)
  # write.csv(plasma_p_vals,paste0('CrossSectionalResults/Plasma',add_on,'.csv'), quote = F)


  if(add_on == ""){
    formula = 'delta_age + Sex + APOE4 + Actual_Age + BMI'
  }else{
    formula = 'delta_age + APOE4 + Actual_Age + BMI'
  }
  fdg_cross = prep_factors(fdg_cross, add_on)
  fdg_p_vals = run_cross_sectional_timepoints(data = fdg_cross,
                                              targets = c('Ang_Mean'),
                                              timepoints = c('bl','m12','m24'),
                                              formula = formula)
  fdg_p_vals = fdg_p_vals[which(fdg_p_vals$timepoint == 'bl'),]
  all_p_vals = rbind(all_p_vals,
                     data.frame('Phenotype' = rownames(fdg_p_vals),
                                'Population' = rep(add_on,nrow(fdg_p_vals)),
                                "Delta_Age_pvalue" = fdg_p_vals$delta_age_p_val,"Beta" = fdg_p_vals$delta_age_Coef,"Std_Error" =  fdg_p_vals$delta_age_Var_Exp))
  # write.csv(fdg_p_vals,paste0('CrossSectionalResults/FDG',add_on,'.csv'), quote = F)


  cogn_cross_main = merge(prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[,c('RID','timepoint','ADNI_MEM')],
                          prep_factors(cogn_cross$`_ADNI_EF`, add_on)[,c('RID','timepoint','ADNI_EF')],
                          by = c('RID','timepoint'))
  cogn_cross_main = merge(cogn_cross_main, prep_factors(cogn_cross$`_ADNI_LAN`, add_on),by = c('RID','timepoint'))
  if(add_on == ""){
    formula = 'delta_age + Sex + Education + Actual_Age + BMI'
  }else{
    formula = 'delta_age + Education + Actual_Age + BMI'
  }
  cogn_p_vals = run_cross_sectional_timepoints(data = cogn_cross_main,
                                              targets = c('ADNI_MEM','ADNI_EF'),
                                              timepoints = c('bl','m12','m24'),
                                              formula = formula)
  
  cogn_p_vals = cogn_p_vals[which(cogn_p_vals$timepoint == 'bl'),]
  all_p_vals = rbind(all_p_vals,
                     data.frame('Phenotype' = rownames(cogn_p_vals),
                                'Population' = rep(add_on,nrow(cogn_p_vals)),
                                "Delta_Age_pvalue" = cogn_p_vals$delta_age_p_val,"Beta" = cogn_p_vals$delta_age_Coef,"Std_Error" =  cogn_p_vals$delta_age_Var_Exp))
  
  # write.csv(cogn_p_vals,paste0('CrossSectionalResults/Cognitive',add_on,'.csv'), quote = F)
}

library(ggpubr)

annotate_figure(ggarrange(
making_data_counts(mri_cross,target_cols = c('HippVol','entorhinal_thickness'), title = 'MRI',thresh = 8),
making_data_counts(csf_cross,target_cols = c('ABETA','TAU','PTAU'), title = 'CSF',thresh = 8),
making_data_counts(plasma_cross,target_cols = c('AB40','AB42'), title = 'Plasma',thresh = 8),
making_data_counts(fdg_cross,target_cols = c('Cing_Mean','LTemp_Mean','Ang_Mean'), title = 'FDG',thresh = 8),
making_data_counts(cogn_cross_main,target_cols = c('ADNI_EF','ADNI_MEM','ADNI_LAN'), title = 'Cognitive Scores',thresh = 8),nrow = 3,ncol =2, top = text_grob(paste0('Female Subjects'),
                                                                                                                                                          color = "red", face = "bold", size = 14)))



for(add_on in c("","__FEMALE","__MALE")){
  amyloid_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/Amyloid_ALL_ADNI',add_on,'.csv'))
  amyloid_cross = prep_factors(amyloid_cross, add_on)
  if(add_on == ""){
    formula = 'delta_age + Sex + APOE4 + Actual_Age + BMI'
  }else{
    formula = 'delta_age + APOE4 + Actual_Age + BMI'
  }
  amyloid_p_vals = run_cross_sectional_timepoints(data = amyloid_cross,
                                              targets = c('SUMMARYSUVR_WHOLECEREBNORM'),
                                              timepoints = c('bl','m24','m48'),
                                              formula = formula)
  
  amyloid_p_vals = amyloid_p_vals[which(amyloid_p_vals$timepoint == 'bl'),]
  all_p_vals = rbind(all_p_vals,
                     data.frame('Phenotype' = rownames(amyloid_p_vals),
                                'Population' = rep(add_on,nrow(amyloid_p_vals)),
                                "Delta_Age_pvalue" = amyloid_p_vals$delta_age_p_val,"Beta" = amyloid_p_vals$delta_age_Coef,"Std_Error" =  amyloid_p_vals$delta_age_Var_Exp))
  print(amyloid_p_vals$delta_age_p_val)
  # write.csv(amyloid_p_vals,paste0('CrossSectionalResults/Amyloid',add_on,'.csv'), quote = F)
}


write.csv(all_p_vals,paste0('CrossSectionalResults/With_BMI_P_Vals.csv'), quote = F)

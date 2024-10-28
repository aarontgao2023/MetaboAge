# This file is to run longitudinal analysis

library(AssociatoR)
library(rempsyc)
library(dplyr)
library(ggpubr)
prep_factors = function(data, add_on){
  
  if(add_on == ""){
    data$Sex = as.factor(data$Sex)
  }
  data$APOE4 = as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == 'VISCODE2')] = 'timepoint'
  #data$delta_age = data$Corrected_Predicted_Age - data$Actual_Age
  return(data)
  
}

longitudinal_p_values = read.csv('LongitudinalPlots/FDR_Corrected_Longitudinal_pvalues.csv')

make_plots_longitudinal = function(data, formula, target, thresh_year, target_label, include_scatter = F,p_value_given){
  
  
  data = data %>% filter(year <= thresh_year)
  
  columns_needed = unique(strsplit(formula,' ')[[1]][sapply(strsplit(formula,' ')[[1]],function(str){!grepl("[^A-Za-z0-9_ ]", str)})])
  columns_needed =  columns_needed[suppressWarnings(is.na(as.numeric(columns_needed)))]
  data = data[,c(columns_needed,target)]
  target_factor = strsplit(formula,' ')[[1]][c(which(strsplit(formula,' ')[[1]] == '*')-1,which(strsplit(formula,' ')[[1]] == '*')+1)]
  data = na.omit(data)
  ## Fitting line
  start_prep = T
  effect_plots = list()
  residual_plots = list()
  out_data_frame  = list()
  beta_data_frame = list()
  
  target_formula = paste(target,formula, sep = ' ~ ')
  
  og = lmer(as.formula(target_formula),data, REML = FALSE)
  og.n = lmer(as.formula(chartr("*","+",target_formula)),data, REML = FALSE)
  temp_og = summary(og)
  k = anova(og,og.n)
  p_val = k$`Pr(>Chisq)`[2]

  
  Quartiles <- quantile(data[,colnames(data) == "delta_age"], probs = c(0, 0.25, 0.5, 0.75, 1))
  data$Quartiles <- cut(data[,colnames(data) == "delta_age"], breaks = Quartiles, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  
  pred1=predict(og,re.form=NA,se.fit=TRUE,nsim=500)
  # save prediction results
  temp_model_data = data
  temp_model_data$boot_fit = pred1$fit
  
  data$Quartiles[data$Quartiles == "Q2"] = "Q1"
  data$Quartiles[data$Quartiles == "Q3"] = "Q4"
  
  # Only plot Q1 and Q4
  plot_Q1Q4 = temp_model_data %>% filter(Quartiles=="Q1"|Quartiles=="Q4")
  if(add_on == ""){
    title_sub = "Whole Population"
  }else if(add_on == "__MALE"){
    title_sub = "Male Population"
  }else if(add_on == "__FEMALE"){
    title_sub = "Female Population"}
  
  
  p_val = p_value_given
  if(p_val < 0.001){
    p_val_out = "< 0.001"
  }else{
    p_val_out = as.character(round(p_val,3))
  }
  
  if(p_val < 0.05){
    p_val_out = as.formula(paste0("'P-value '~ bold('",p_val_out,"')"))
  }else{
    p_val_out = paste0("P-value ",p_val_out)
  }
  
  if(!include_scatter){
    w = ggplot(plot_Q1Q4,aes(x=year,y=boot_fit,col=Quartiles))+geom_smooth(aes(year,boot_fit),linewidth=1.5,method=lm,se=TRUE)+
      ylab(target_label)+xlab("Year")+
      labs(title = title_sub,
           subtitle = (p_val_out))+
      scale_color_manual(values=c("skyblue","coral1"),name = "",labels=c("sub-median of metaboAge delta", "super-median of metaboAge delta")) 
  }else{
    w = ggplot(plot_Q1Q4,aes(x=year,y=!! sym(target),col=Quartiles)) + geom_point()+geom_smooth(aes(year,boot_fit),linewidth=1.5,method=lm,se=TRUE)+
      ylab(target_label)+xlab("Year")+
      labs(title = title_sub,
           subtitle = (p_val_out))+
      scale_color_manual(values=c("skyblue","coral1"),name = "",labels=c("sub-median of metaboAge delta", "super-median of metaboAge delta")) 
  }
  

  
  return(w)
}


paths = list(
  
  'Data/LongitudianlDataPrepared/Amyloid_Longitudinal_ADNI',
  'Data/LongitudianlDataPrepared/CSF_Longitudinal_ADNI',
  'Data/LongitudianlDataPrepared/MRI_Longitudinal_ADNI',
  'Data/LongitudianlDataPrepared/Cognitivie_Scores',
  'Data/LongitudianlDataPrepared/Cognitivie_Scores'
  
)

preprocess = list(
  
  function(x){
    amyloid_cross = read.csv(paste0(x,add_on,'.csv'))
    amyloid_cross = amyloid_cross[which(amyloid_cross$year%%1 != 0.5),]
    amyloid_cross = amyloid_cross[which(!amyloid_cross$year %in% c(1,3,10)),]
    amyloid_cross = prep_factors(amyloid_cross, add_on)
    amyloid_bl_delta_age_subjects_more_than_one = names(table(amyloid_cross$RID)[table(amyloid_cross$RID) > 2 ])
    amyloid_bl_delta_age_subjects_more_than_one = amyloid_bl_delta_age_subjects_more_than_one[!amyloid_bl_delta_age_subjects_more_than_one %in% c(668,677)]
    amyloid_cross = amyloid_cross %>% filter(RID %in% amyloid_bl_delta_age_subjects_more_than_one)
    return(amyloid_cross)
  },
  function(x){
    csf_cross = read.csv(paste0(x,add_on,'.csv'))
    if(add_on != ""){
      csf_cross = csf_cross %>% filter(BATCH == "MEDIAN")
      
    }
    csf_cross = prep_factors(csf_cross, add_on)
    csf_cross$ABETA = log(csf_cross$ABETA)
    csf_cross$TAU = log(csf_cross$TAU)
    csf_cross$PTAU = log(csf_cross$PTAU)
    #removing subject with only baseline timepoint
    csf_bl_delta_age_subjects_more_than_one = names(table(csf_cross$RID)[table(csf_cross$RID) != 1])
    csf_cross = csf_cross %>% filter(RID %in% csf_bl_delta_age_subjects_more_than_one)
    return(csf_cross)
  },
  function(x){
    mri_cross = read.csv(paste0(x,add_on,'.csv'))
    mri_cross = prep_factors(mri_cross, add_on)
    return(mri_cross)
  },
  function(x){
    cogn_cross = readRDS(paste0(x,add_on,'.RDS'))
    cogn_cross_main = merge(prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[,c('RID','timepoint','ADNI_MEM')],
                            prep_factors(cogn_cross$`_ADNI_EF`, add_on),
                            by = c('RID','timepoint'))
    cogn_cross_main = na.omit(cogn_cross_main)
    return(cogn_cross_main)
  },
  function(x){
    cogn_cross = readRDS(paste0(x,add_on,'.RDS'))
    cogn_cross_main = merge(prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[,c('RID','timepoint','ADNI_MEM')],
                            prep_factors(cogn_cross$`_ADNI_EF`, add_on),
                            by = c('RID','timepoint'))
    cogn_cross_main = na.omit(cogn_cross_main)
    return(cogn_cross_main)
  }
  
)


formula_list_whole = list(
  
  '1 + APOE4 + delta_age * year + Sex + Actual_Age + ( 1 + year | RID )' ,
  '1 + Actual_Age + APOE4 + delta_age * year + Sex + ( 1 + year | RID )',
  '1 + Actual_Age + Sex + APOE4 + delta_age * year + scale( ICV ) + Education + ( 1 + year | RID )',
  '1 + Sex + Actual_Age + Education + delta_age * year + ( 1 + year | RID )',
  '1 + Sex + Actual_Age + Education + delta_age * year + ( 1 + year | RID )'
  
)

formula_list_male_female = list(
  
  '1 + APOE4 + delta_age * year + Actual_Age + ( 1 + year | RID )' ,
  '1 + Actual_Age + APOE4 + delta_age * year + ( 1 + year | RID )',
  '1 + Actual_Age + APOE4 + delta_age * year + scale( ICV ) + Education + ( 1 + year | RID )',
  '1 + Actual_Age + Education + delta_age * year + ( 1 + year | RID )',
  '1 + Actual_Age + Education + delta_age * year + ( 1 + year | RID )'
  
)

years_threshold = c(4,3,8,8,8)

targets = c('SUMMARYSUVR_COMPOSITE_REFNORM','PTAU','entorhinal_thickness','ADNI_EF', 'ADNI_MEM')
target_names = c("Amyloid SUVR","CSF PTAU","Entorhinal Thickness","Executive Function","Memory")
rownames(longitudinal_p_values) = longitudinal_p_values$Phenotype
fdr_corrected_p_value_whole = longitudinal_p_values[c('Amyloid SUVR Composite Score',"CSF PTAU","MRI Entorhinal Thickness","Cognition Composite Executive Function","Cognition Composite Memory"),2]
fdr_corrected_p_value_male = longitudinal_p_values[c('Amyloid SUVR Composite Score',"CSF PTAU","MRI Entorhinal Thickness","Cognition Composite Executive Function","Cognition Composite Memory"),3]
fdr_corrected_p_value_female = longitudinal_p_values[c('Amyloid SUVR Composite Score',"CSF PTAU","MRI Entorhinal Thickness","Cognition Composite Executive Function","Cognition Composite Memory"),4]

long_effect_plots = list()
long_effect_plots_w_scatter = list()



for(i in 1:length(targets)){
  
  
  for(add_on in c("","__MALE","__FEMALE")){
    
    data = lapply(paths[i], preprocess[[i]])[[1]]
    if (add_on == ""){
      formula = formula_list_whole[[i]]
    }else{
      formula = formula_list_male_female[[i]]
    }
    
    if (add_on == ""){
      p_value_corr = fdr_corrected_p_value_whole[i]
    }else if(add_on == "__MALE"){
      p_value_corr = fdr_corrected_p_value_male[i]
    }else{
      p_value_corr = fdr_corrected_p_value_female[i]
    }
    
    
    target = targets[[i]]
    thresh_year = years_threshold[[i]]
    
    long_effect_plots[[paste0(target_names[i],add_on)]] = make_plots_longitudinal(data = data,
                                formula = formula,
                                target = target,
                                thresh_year = thresh_year,
                                target_label = target_names[i],
                                p_value_given = p_value_corr)
    
    long_effect_plots_w_scatter[[paste0(target_names[i],add_on)]] = make_plots_longitudinal(data = data,
                                                                          formula = formula,
                                                                          target = target,
                                                                          thresh_year = thresh_year,
                                                                          target_label = target_names[i],
                                                                          include_scatter = T,
                                                                          p_value_given = p_value_corr)
  }
  
  
  
}


annotate_figure(ggarrange(#long_effect_plots[[1]],long_effect_plots[[2]],long_effect_plots[[3]],
                          #long_effect_plots[[4]],long_effect_plots[[5]],long_effect_plots[[6]],
                          long_effect_plots[[7]],long_effect_plots[[8]],long_effect_plots[[9]],
                          long_effect_plots[[10]],long_effect_plots[[11]],long_effect_plots[[12]],
                          long_effect_plots[[13]],long_effect_plots[[14]],long_effect_plots[[15]],
                          labels = c("A", "","","B", "","","C", "",""), common.legend = TRUE,
                          ncol = 3, nrow = 3))   #
annotate_figure(ggarrange(#long_effect_plots_w_scatter[[1]],long_effect_plots_w_scatter[[2]],long_effect_plots_w_scatter[[3]],
                          long_effect_plots_w_scatter[[4]],long_effect_plots_w_scatter[[5]],long_effect_plots_w_scatter[[6]],
                          long_effect_plots_w_scatter[[7]],long_effect_plots_w_scatter[[8]],long_effect_plots_w_scatter[[9]],
                          long_effect_plots_w_scatter[[10]],long_effect_plots_w_scatter[[11]],long_effect_plots_w_scatter[[12]],
                          long_effect_plots_w_scatter[[13]],long_effect_plots_w_scatter[[14]],long_effect_plots_w_scatter[[15]],
                          labels = c("A", "","","B", "","","C", "","","D", "",""), common.legend = TRUE,
                          ncol = 3, nrow = 4))









# CROSS-SECTIONAL DATA SUPPLEMENTARY --------------------------------------


#install.packages("ggstatsplot")

library(ggstatsplot)
mri_cross = read.csv(paste0('Data/CrossSectionalDataPrepared/MRI_ALL_ADNI',add_on,'.csv'))
mri_cross$delta_age = mri_cross$Corrected_Predicted_Age - mri_cross$Actual_Age
grouped_ggscatterstats(
  data             = dplyr::filter(mri_cross, VISCODE2 %in% c("bl", "m12", "m24")),
  x                = Corrected_Predicted_Age,
  y                = HippVol,
  grouping.var     = VISCODE2,
  xlab             = "Corrected_Predicted_Age",
  ggtheme          = ggplot2::theme_grey(),
  ggplot.component = list(ggplot2::scale_x_continuous(breaks = seq(40, 100, 10), limits = (c(40, 100)))),
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Relationship between Corrected_Predicted_Age and Hippocaml Volume")
)

grouped_ggscatterstats(
  data             = dplyr::filter(mri_cross, VISCODE2 %in% c("bl", "m12", "m24")),
  x                = delta_age,
  y                = HippVol,
  grouping.var     = VISCODE2,
  xlab             = "Delta Age",
  ggtheme          = ggplot2::theme_grey(),
  ggplot.component = list(ggplot2::scale_x_continuous(breaks = seq(-24, 57, 10), limits = (c(-25, 60)))),
  plotgrid.args    = list(nrow = 1),
  annotation.args  = list(title = "Relationship between Delta Age and Hippocaml Volume")
)



### Cross Sectional plots
library(hrbrthemes)

prep_factors = function(data, add_on){
  if(add_on == ""){
    data$Sex = as.factor(data$Sex)
  }
  data$APOE4 = as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == 'VISCODE2')] = 'timepoint'
  data$delta_age = data$Corrected_Predicted_Age - data$Actual_Age
  data = data %>% filter(timepoint %in% c('bl','m12','m24'))
  return(data)
  
}
make_plots_crossSectional = function(data, formula, target, thresh_year, target_label, include_scatter = F,p_value_given){
  
  
  data = data %>% filter(year <= thresh_year)
  
  columns_needed = unique(strsplit(formula,' ')[[1]][sapply(strsplit(formula,' ')[[1]],function(str){!grepl("[^A-Za-z0-9_ ]", str)})])
  columns_needed =  columns_needed[suppressWarnings(is.na(as.numeric(columns_needed)))]
  data = data[,c(columns_needed,target)]
  target_factor = strsplit(formula,' ')[[1]][c(which(strsplit(formula,' ')[[1]] == '*')-1,which(strsplit(formula,' ')[[1]] == '*')+1)]
  data = na.omit(data)
  ## Fitting line
  start_prep = T
  effect_plots = list()
  residual_plots = list()
  out_data_frame  = list()
  beta_data_frame = list()
  
  target_formula = paste(target,formula, sep = ' ~ ')
  
  og = lmer(as.formula(target_formula),data, REML = FALSE)
  og.n = lmer(as.formula(chartr("*","+",target_formula)),data, REML = FALSE)
  temp_og = summary(og)
  k = anova(og,og.n)
  p_val = k$`Pr(>Chisq)`[2]
  
  
  Quartiles <- quantile(data[,colnames(data) == "delta_age"], probs = c(0, 0.25, 0.5, 0.75, 1))
  data$Quartiles <- cut(data[,colnames(data) == "delta_age"], breaks = Quartiles, labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  
  pred1=predict(og,re.form=NA,se.fit=TRUE,nsim=500)
  # save prediction results
  temp_model_data = data
  temp_model_data$boot_fit = pred1$fit
  
  data$Quartiles[data$delta_age > 0] = "Q4"
  data$Quartiles[data$delta_age < 0] = "Q1"
  
  
  # Only plot Q1 and Q4
  plot_Q1Q4 = temp_model_data %>% filter(Quartiles=="Q1"|Quartiles=="Q4")
  if(add_on == ""){
    title_sub = "Whole Population"
  }else if(add_on == "__MALE"){
    title_sub = "Male Population"
  }else if(add_on == "__FEMALE"){
    title_sub = "Female Population"}
  
  
  p_val = p_value_given
  if(p_val < 0.001){
    p_val_out = "< 0.001"
  }else{
    p_val_out = as.character(round(p_val,4))
  }
  
  if(p_val < 0.05){
    p_val_out = as.formula(paste0("'P-value '~ bold('",p_val_out,"')"))
  }else{
    p_val_out = paste0("P-value ",p_val_out)
  }
  
  if(!include_scatter){
    w = ggplot(plot_Q1Q4,aes(x=year,y=boot_fit,col=Quartiles))+geom_smooth(aes(year,boot_fit),linewidth=1.5,method=lm,se=TRUE)+
      ylab(target_label)+xlab("Year")+
      labs(title = title_sub,
           subtitle = (p_val_out))+
      scale_color_manual(values=c("skyblue","coral1"),name = "",labels=c("metaboAge < 0", "metaboAge > 0")) 
  }else{
    w = ggplot(plot_Q1Q4,aes(x=year,y=!! sym(target),col=Quartiles)) + geom_point()+geom_smooth(aes(year,boot_fit),linewidth=1.5,method=lm,se=TRUE)+
      ylab(target_label)+xlab("Year")+
      labs(title = title_sub,
           subtitle = (p_val_out))+
      scale_color_manual(values=c("skyblue","coral1"),name = "",labels=c("metaboAge < 0", "metaboAge > 0")) 
  }
  
  
  
  return(w)
}

paths = list(
  
  'Data/CrossSectionalDataPrepared/CSF_ALL_ADNI',
  'Data/CrossSectionalDataPrepared/FDG_ALL_ADNI',
  'Data/CrossSectionalDataPrepared/MRI_ALL_ADNI',
  'Data/CrossSectionalDataPrepared/MRI_ALL_ADNI',
  'Data/CrossSectionalDataPrepared/Cognitivie_Scores'
  
)
  
preprocess = list(
  
  function(x){
    csf_cross = read.csv(paste0(x,add_on,'.csv'))
    csf_cross = prep_factors(csf_cross, add_on)
    csf_cross$ABETA = log(csf_cross$ABETA)
    csf_cross$TAU = log(csf_cross$TAU)
    csf_cross$PTAU = log(csf_cross$PTAU)
    csf_bl_delta_age_subjects_more_than_one = names(table(csf_cross$RID)[table(csf_cross$RID) != 1])
    csf_cross = csf_cross %>% filter(RID %in% csf_bl_delta_age_subjects_more_than_one)
    return(csf_cross)
  },
  function(x){
    fdg_cross = read.csv(paste0(x,add_on,'.csv'))
    fdg_cross = prep_factors(fdg_cross, add_on)
    return(fdg_cross)
  },
  function(x){
    mri_cross = read.csv(paste0(x,add_on,'.csv'))
    mri_cross = prep_factors(mri_cross, add_on)
    return(mri_cross)
  },
  function(x){
    mri_cross = read.csv(paste0(x,add_on,'.csv'))
    mri_cross = prep_factors(mri_cross, add_on)
    return(mri_cross)
  },
  function(x){
    cogn_cross = readRDS(paste0(x,add_on,'_ALL.RDS'))
    cogn_cross_main = merge(prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[,c('RID','timepoint','ADNI_MEM')],
                            prep_factors(cogn_cross$`_ADNI_EF`, add_on),
                            by = c('RID','timepoint'))
    cogn_cross_main = na.omit(cogn_cross_main)
    return(cogn_cross_main)
  }
  
)

cross_sectional_p_values = read.csv('CrossSectionalResults/FDR_Corrected_All_pvalues.csv')
#timepoint_cs = read.csv('CrossSectionalResults/All_CrossSectional_Pvalues.csv')
colnames(cross_sectional_p_values) = c('Phenotype','p_value_whole','p_value_male','p_value_female')
#cross_sectional_p_values$Timepoint = timepoint_cs$X.1[2:nrow(timepoint_cs)]
#cross_sectional_p_values = cross_sectional_p_values %>% filter(Timepoint == 'bl') %>% 
cross_sectional_p_values = cross_sectional_p_values %>% filter(Phenotype %in% c('CSF PTAU','PET Angular Mean','MRI Entorhinal Thickness','MRI Hippocampal Volume','Cognition Composite Memory'))
cross_sectional_p_values$Timepoint = NULL



targets = c('CSF PTAU','PET Angular Mean','MRI Entorhinal Thickness','MRI Hippocampal Volume','Cognition Composite Memory')
targets_name = c("PTAU","Ang_Mean","entorhinal_thickness","HippVol","ADNI_MEM")
rownames(cross_sectional_p_values) = cross_sectional_p_values$Phenotype
fdr_corrected_p_value_whole = cross_sectional_p_values[targets,2]
fdr_corrected_p_value_male = cross_sectional_p_values[targets,3]
fdr_corrected_p_value_female = cross_sectional_p_values[targets,4]


cross_plots = list()

for(i in 1:length(targets)){
  
  
  for(add_on in c("","__MALE","__FEMALE")){
    
    data = lapply(paths[i], preprocess[[i]])[[1]]
    #boxplot.stats(data$delta_age)$out
    data = data %>% filter(delta_age < 30)
    
    if (add_on == ""){
      p_value_corr = fdr_corrected_p_value_whole[i]
    }else if(add_on == "__MALE"){
      p_value_corr = fdr_corrected_p_value_male[i]
    }else{
      p_value_corr = fdr_corrected_p_value_female[i]
    }
    
    target = targets[[i]]
    temp_y = targets_name[i]
    colnames(data)[which(colnames(data)==temp_y)] = 'Temp_Y'
    if(add_on == ""){
      title_sub = "Whole Population"
    }else if(add_on == "__MALE"){
      title_sub = "Male Population"
    }else if(add_on == "__FEMALE"){
      title_sub = "Female Population"}
    
    p_val = p_value_corr
    if(p_val < 0.001){
      p_val_out = "< 0.001"
    }else{
      p_val_out = as.character(round(p_val,4))
    }
    
    if(p_val < 0.05){
      p_val_out = as.formula(paste0("'P-value '~ bold('",p_val_out,"')"))
    }else{
      p_val_out = paste0("P-value ",p_val_out)
    }
    

    
      cross_plots[[paste0(targets[i],add_on)]] = ggplot(data, aes(x=delta_age, y=Temp_Y)) +
        geom_point() +
        labs(title = title_sub,
             subtitle = (p_val_out)) + 
        geom_smooth(method=lm , color="coral1", fill="skyblue", se=T) + xlab('Delta Age') + ylab(target) +
      theme_ipsum()
    
     }
  
  
  
}


annotate_figure(ggarrange(cross_plots[[1]],cross_plots[[2]],cross_plots[[3]],
                          cross_plots[[4]],cross_plots[[5]],cross_plots[[6]],
                          cross_plots[[7]],cross_plots[[8]],cross_plots[[9]],
                          cross_plots[[10]],cross_plots[[11]],cross_plots[[12]],
                          cross_plots[[13]],cross_plots[[14]],cross_plots[[15]],
  labels = c("A", "","","B", "","","C", "","","D", "","","E","",""), common.legend = TRUE,
  ncol = 3, nrow = 5))   #



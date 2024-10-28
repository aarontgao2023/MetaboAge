# This file is to run survival analysis
#########################################################################
## Data Preparation -
getwd()
library(dplyr)
library("survival")
library("survminer")
library(ggsurvfit)
setwd('C:/Users/ppugale/OneDrive - Indiana University/Documents/AgePredictionProject/')
#install.packages(c("survival", "survminer"))
#install.packages("ggsurvfit")
# Obtaining demographic information ---------------------------------------


demo_adni = read.csv('ADNI_Data/ADNI_Long_Nightingale_Demographics_at_BL.csv')
demo_adni$Sex[which(demo_adni$Sex == 'Male')] = 1
demo_adni$Sex[which(demo_adni$Sex == 'Female')] = 0
demo_adni$Sex = as.factor(demo_adni$Sex)
demo_adni$APOE4 = as.factor(demo_adni$APOE4)

# Obtaining longitudinal diagnosis information
demo_adni_D = readxl::read_xlsx('Data/ADNI1GO23_selectedvars_12242020.xlsx', sheet = 'DX_Demo')
diangosis_long = data.frame('RID' = numeric(), 'VISCODE2' = character(), 'Diagnosis' = numeric())
for(r in demo_adni_D$RID){

  temp = demo_adni_D[which(demo_adni_D$RID == r), grep('_DXGrp',colnames(demo_adni_D))]
  diangosis_long = rbind(diangosis_long,
                         data.frame('RID' = rep(r,17), 'VISCODE2' = c(0,seq(6,24,6),seq(36,168,12)), 'Diagnosis' = unlist(as.vector(temp[1,-1]))))

}
diangosis_long = na.omit(diangosis_long)
#checking if diagnosis changes for any subjects over time
change_in_Diag = c()
for(r in diangosis_long$RID){
  if(length(unique(diangosis_long$Diagnosis[which(diangosis_long$RID == r)])) == 1){
    cat('-')
  }else{
    change_in_Diag = c(change_in_Diag,r)
    cat('~')
  }
}

#########################################################################
# Importing metabolite Age file -------------------------------------------

predicted_adni_brain_age = read.csv('Oct2023_Results/Predicted_Age_Whole_ADNI_Population_PolyBias.csv')
predicted_adni_brain_age_male = read.csv('Oct2023_Results/Predicted_Age_Male_ADNI_Population_PolyBias.csv')
predicted_adni_brain_age_female = read.csv('Oct2023_Results/Predicted_Age_Female_ADNI_Population_PolyBias.csv')

predicted_adni_brain_age$delta_age = 0
predicted_adni_brain_age_male$delta_age = 0
predicted_adni_brain_age_female$delta_age = 0

# logging a delta_age_bl column


update_data_to_survival_data = function(predicted_adni_brain_age, demo_adni, demo_adni_D,type = 3){
  #type 4 --> NO CN.. MCI-->AD
  #type 5 --> Starting as CN --> MCI or AD
  #type 3 --> CN and EMCI as starting point
  #type 2 --> CN only as starting point
  #type 1 --> NO CN, EMCI as starting point
  VISCODE_Converter = function(viscode_list){
    new_e = c()
    for(e in viscode_list){
      if(e == 'bl'){
        new_e = c(new_e, 0)
      }else{
        new_e = c(new_e, strsplit(e,'m')[[1]][2])
      }
    }
    return(new_e)
  }
  predicted_adni_brain_age$delta_age = 0
  for(i in unique(predicted_adni_brain_age$RID)){
    temp_idx = which(predicted_adni_brain_age$RID == i & predicted_adni_brain_age$VISCODE2 == 'bl')
    predicted_adni_brain_age$delta_age[which(predicted_adni_brain_age$RID == i)] = predicted_adni_brain_age$Corrected_Predicted_Age[temp_idx] - predicted_adni_brain_age$Actual_Age[temp_idx]
  }
  complete_data_w_covariate = merge(predicted_adni_brain_age, demo_adni[,c('RID','Sex','APOE4','BMI','Education')], by = 'RID')
  complete_data_w_covariate$VISCODE2 = sapply(complete_data_w_covariate$VISCODE2,VISCODE_Converter )
  complete_data_w_covariate$Diagnosis = NULL
  complete_data_w_covariate = complete_data_w_covariate[,c('RID','delta_age','Sex','APOE4','BMI','Education')]
  complete_data_w_covariate = merge(complete_data_w_covariate, diangosis_long, by = c('RID'))
  complete_data_w_covariate_w_o_bl_AD = unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 4 & complete_data_w_covariate$VISCODE2 == 0),'RID'])
  complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 5),'RID']))
  
  if(type == 1){
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 1 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 3 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    
    complete_data_w_covariate = complete_data_w_covariate[!(complete_data_w_covariate$RID %in% complete_data_w_covariate_w_o_bl_AD),]
    complete_data_w_covariate = na.omit(complete_data_w_covariate)
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 2)] = 1
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 3)] = 2
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 4)] = 3
  }else if(type == 2){
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 2 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 3 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate = complete_data_w_covariate[!(complete_data_w_covariate$RID %in% complete_data_w_covariate_w_o_bl_AD),]
    complete_data_w_covariate = na.omit(complete_data_w_covariate)
  }else if(type == 3){
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 3 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate = complete_data_w_covariate[!(complete_data_w_covariate$RID %in% complete_data_w_covariate_w_o_bl_AD),]
    complete_data_w_covariate = na.omit(complete_data_w_covariate)
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 2)] = 1
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 3)] = 2
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 4)] = 3
  }else if(type == 4){
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 1 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate = complete_data_w_covariate[!(complete_data_w_covariate$RID %in% complete_data_w_covariate_w_o_bl_AD),]
    complete_data_w_covariate = na.omit(complete_data_w_covariate)
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 2)] = 1
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 3)] = 1
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 4)] = 2
  }else if(type == 5){
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 3 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate_w_o_bl_AD = c(complete_data_w_covariate_w_o_bl_AD,unique(complete_data_w_covariate[which(complete_data_w_covariate$Diagnosis == 2 & complete_data_w_covariate$VISCODE2 == 0),'RID']))
    complete_data_w_covariate = complete_data_w_covariate[!(complete_data_w_covariate$RID %in% complete_data_w_covariate_w_o_bl_AD),]
    complete_data_w_covariate = na.omit(complete_data_w_covariate)
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 3)] = 2
    complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$Diagnosis == 4)] = 2
  }
  
  
  

  complete_data_w_covariate$Time = as.numeric(complete_data_w_covariate$VISCODE2)/12
  complete_data_w_covariate$Diagnosis = as.numeric(complete_data_w_covariate$Diagnosis)
  complete_data_w_covariate$Progression_AD = 1 # right censored
  for(i in unique(complete_data_w_covariate$RID)){
    if(max(unique(complete_data_w_covariate$Diagnosis)) %in% complete_data_w_covariate$Diagnosis[which(complete_data_w_covariate$RID == i)]){ 
      complete_data_w_covariate$Progression_AD[which(complete_data_w_covariate$RID == i)] = 2 #event
      t = min(complete_data_w_covariate$Time[which(complete_data_w_covariate$RID == i & complete_data_w_covariate$Diagnosis == max(unique(complete_data_w_covariate$Diagnosis)) )])
      complete_data_w_covariate$Time[which(complete_data_w_covariate$RID == i)] = t
    }
  }
  out_df = data.frame('RID' = numeric(),
                      #'delta_age_group' = numeric(),
                      'delta_age' = numeric(),
                      'Progression_AD' = numeric(),
                      'APOE4' = numeric(),
                      'Time' = numeric())
  for(i in unique(complete_data_w_covariate$RID)){

    out_df = rbind(out_df,data.frame('RID' = i,
                                    # 'delta_age_group' = unique(complete_data_w_covariate$delta_age_group[which(complete_data_w_covariate$RID == i)]),
                                     'delta_age' = unique(complete_data_w_covariate$delta_age[which(complete_data_w_covariate$RID == i)]),
                                     'Progression_AD' = unique(complete_data_w_covariate$Progression_AD[which(complete_data_w_covariate$RID == i)]),
                                     'APOE4' = unique(complete_data_w_covariate$APOE4[which(complete_data_w_covariate$RID == i)]),
                                     'Time' = max(unique(complete_data_w_covariate$Time[which(complete_data_w_covariate$RID == i)]))))

  }
  q2 = quantile(out_df$delta_age)[3]
  q3 = quantile(out_df$delta_age)[3]
  out_df$delta_age_group = 1
  
  #Change this PARAMETER IF CHANGING TO DELTA AGE > 0  AND NOT Q1 ADN Q4
  out_df$delta_age_group[which(out_df$delta_age > q3)] = 1
  out_df$delta_age_group[which(out_df$delta_age <= q2)] = 0
  # 
  # out_df$delta_age_group[which(out_df$delta_age > 0)] = 1
  # out_df$delta_age_group[which(out_df$delta_age <= 0)] = 0
  # 
  
  return(out_df)
}

#complete_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age, demo_adni, demo_adni_D, type = 5)

#male_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age_male, demo_adni, demo_adni_D)

#female_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age_female, demo_adni, demo_adni_D)


# Survival Analysis Fit ---------------------------------------------------



for(t in c(3,4,5)){
  
  complete_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age, demo_adni, demo_adni_D, type = t)
  
  complete_data_w_covariate = complete_data_w_covariate %>% filter(Time <= 10)
  complete_data_w_covariate$APOE4[which(complete_data_w_covariate$APOE4 == 2)] = 1
  complete_data_w_covariate_w_o_APOE = complete_data_w_covariate[which(complete_data_w_covariate$APOE4 == 0),]
  complete_data_w_covariate_w_APOE = complete_data_w_covariate[which(complete_data_w_covariate$APOE4 == 1),]
  
  if(t == 1){
    title = paste0("Delta Age For Survival \n",
                   "NO CN, EMCI as starting point")
  }else if(t == 2){
    title = paste0("Delta Age For Survival \n",
                   "CN as starting point")
  }else if(t == 3){
    title = paste0("Delta Age For Survival \n",
                   "CN, EMCI as starting point")
  }else if(t == 4){
    title = paste0("Delta Age For Survival \n",
                   "MCI as starting point -->AD")
  }else if(t == 5){
    title = paste0("Delta Age For Survival \n",
                   "CN as starting point -->MCI or AD")
  }
  # Running Survival function
  fit1 <- survminer::surv_fit(Surv(Time,Progression_AD) ~ delta_age_group, data = complete_data_w_covariate)
  w = ggsurvplot(fit1, data = complete_data_w_covariate,
                 legend.title = title,
                 legend.labs = c('delta_age Q1','delta_age Q4'),
                 # Add p-value
                 pval = TRUE,
                 conf.int = TRUE,
                 risk.table = TRUE,
                 tables.height = 0.2,
                 tables.theme = theme_cleantable(), palette = c('skyblue','orange'),
                 xlab="Time(years)",
                 ylab= "Survival probability(t10)")
  

  print(paste0("Type ",t))
  w = survminer:::.build_ggsurvplot(x = w,
                                surv.plot.height = NULL,
                                risk.table.height = NULL,
                                ncensor.plot.height = NULL)

  ggsave(paste0("LongitudinalPlots/SurvivalPlots/Whole_Population_Type",t,".png"), w)
  
}




# Survival Analysis with APOE ---------------------------------------------


#type 3 --> CN and EMCI as starting point
#type 2 --> CN only as starting point
#type 1 --> NO CN, EMCI as starting point
#type 4 --> NO CN, EMCI and LMCI as starting points


for(t in c(3,4,5)){

  complete_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age, demo_adni, demo_adni_D, type = t)
  male_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age_male, demo_adni, demo_adni_D, type = t)
  female_data_w_covariate = update_data_to_survival_data(predicted_adni_brain_age_female, demo_adni, demo_adni_D, type = t)
  time_filter = 10
   
  complete_data_w_covariate = complete_data_w_covariate %>% filter(Time <= time_filter)
  male_data_w_covariate = male_data_w_covariate %>% filter(Time <= time_filter)
  female_data_w_covariate = female_data_w_covariate %>% filter(Time <= time_filter)

  compare_to_apoe = function(complete_data_w_covariate, title = "Whole Population - EMCI,CN --> AD"){
  
    complete_data_w_covariate$APOE4[which(complete_data_w_covariate$APOE4 == 2)] = 1
    data_w_o_APOE = complete_data_w_covariate[which(complete_data_w_covariate$APOE4 == 0),]
    data_w_APOE = complete_data_w_covariate[which(complete_data_w_covariate$APOE4 != 0),]
    
    fit0 <- survminer::surv_fit(Surv(Time,Progression_AD) ~ delta_age_group,data = complete_data_w_covariate)
    ggsurv_w_all <- ggsurvplot(fit = fit0, data = complete_data_w_covariate, pval = round(surv_pvalue(fit0)$pval,2),title = 'All Population',
                                tables.height = 0.2,           risk.table = TRUE,legend.labs = c('delta_age Q1','delta_age Q4'),
                                tables.theme = theme_cleantable(),
                                ggtheme = theme_bw(), palette = c('skyblue','orange'))
    fit1 <- survminer::surv_fit(Surv(Time,Progression_AD) ~ delta_age_group,data = data_w_APOE)
    ggsurv_w_APOE <- ggsurvplot(fit = fit1, data = data_w_APOE, pval = round(surv_pvalue(fit1)$pval,2),title = 'APOE e4 carriers',
                                tables.height = 0.2,           risk.table = TRUE,legend.labs = c('delta_age Q1','delta_age Q4'),
                                tables.theme = theme_cleantable(),
                                ggtheme = theme_bw(), palette = c('skyblue','orange'))
    fit2 <- survminer::surv_fit(Surv(Time,Progression_AD) ~ delta_age_group, data = data_w_o_APOE)
    print(summary(fit1))

    ggsurv_w_o_APOE <- ggsurvplot(fit2, data = data_w_o_APOE, pval = round(surv_pvalue(fit2)$pval,2),title = 'APOE e4 non-carriers',
                                tables.height = 0.2,           risk.table = TRUE,legend.labs = c('delta_age Q1','delta_age Q4'),
                                tables.theme = theme_cleantable(),
                                ggtheme = theme_bw(), palette = c('skyblue','orange'))
  
    outplot = arrange_ggsurvplots(list(ggsurv_w_APOE,ggsurv_w_o_APOE, ggsurv_w_all),print = F,title = title,
                        ncol = 3, nrow = 1, risk.table.height = 0.4,)
  
    return(outplot)
  
  }
  
  if(t == 1){
    title = paste0("Delta Age For Survival \n",
                   "NO CN, EMCI as starting point")
  }else if(t == 2){
    title = paste0("Delta Age For Survival \n",
                   "CN as starting point")
  }else if(t == 3){
    title = paste0("Delta Age For Survival \n",
                   "CN, EMCI as starting point")
  }else if(t == 4){
    title = paste0("Delta Age For Survival \n",
                   "MCI as starting point -->AD")
  }else if(t == 5){
    title = paste0("Delta Age For Survival \n",
                   "CN as starting point -->MCI or AD")
  }
  
  w2 = compare_to_apoe(complete_data_w_covariate, title = title)
  ggsave(paste0("LongitudinalPlots/SurvivalPlots/Whole_Population_APOE_Type",t,".png"), w2)
  m2 = compare_to_apoe(male_data_w_covariate, title = title)
  ggsave(paste0("LongitudinalPlots/SurvivalPlots/Male_Population_APOE_Type",t,".png"), m2)
  f2 = compare_to_apoe(female_data_w_covariate, title = title)
  ggsave(paste0("LongitudinalPlots/SurvivalPlots/Female_Population_APOE_Type",t,".png"), f2)
  
}


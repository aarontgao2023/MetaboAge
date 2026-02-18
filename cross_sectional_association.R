# =============================================================================
# Cross-Sectional Association Analysis: MetaboAge Delta vs. Brain/Cognitive Phenotypes
# =============================================================================
#
# Description: Runs linear regression models (via AssociatoR) to test the
#              cross-sectional association between metabolic age delta
#              (metaboAge - actual age) and brain imaging biomarkers (MRI,
#              CSF, FDG-PET, Amyloid PET) and cognitive composite scores in
#              the ADNI cohort. Analyses are stratified by sex (whole, male,
#              female). Baseline (bl) p-values, beta coefficients, and standard
#              errors are compiled into a single results table.
#
# Dependencies:
#   devtools::install_github("PradoVarathan/AssociatoR")
#
# Input files:
#   - CrossSectionalDataPrepared/MRI_ALL_ADNI[_MALE|_FEMALE].csv
#   - CrossSectionalDataPrepared/CSF_ALL_ADNI[_MALE|_FEMALE].csv
#   - CrossSectionalDataPrepared/FDG_ALL_ADNI[_MALE|_FEMALE].csv
#   - CrossSectionalDataPrepared/Cognitivie_Scores[_MALE|_FEMALE]_ALL.RDS
#   - CrossSectionalDataPrepared/Amyloid_ALL_ADNI[_MALE|_FEMALE].csv
#
# Output:
#   - CrossSectionalResults/With_BMI_P_Vals.csv
# =============================================================================

# Set the working directory to the project root before running this script.
# setwd('/path/to/AgePredictionProject/')

library(AssociatoR)
library(rempsyc)
library(dplyr)
library(ggpubr)


# -----------------------------------------------------------------------------
# Helper Function
# -----------------------------------------------------------------------------

# Prepare common factor variables and compute delta_age.
# For whole-population runs (add_on == ""), Sex is coerced to factor.
prep_factors <- function(data, add_on) {
  if (add_on == "") {
    data$Sex <- as.factor(data$Sex)
  }
  data$APOE4 <- as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == "VISCODE2")] <- "timepoint"
  data$delta_age <- data$Corrected_Predicted_Age - data$Actual_Age
  return(data)
}


# -----------------------------------------------------------------------------
# Cross-Sectional Association Analysis
# (MRI, CSF, FDG-PET, and Cognitive Scores)
# -----------------------------------------------------------------------------

all_p_vals <- data.frame(
  "Phenotype"        = character(),
  "Population"       = character(),
  "Delta_Age_pvalue" = numeric(),
  "Beta"             = numeric(),
  "Std_Error"        = numeric()
)

for (add_on in c("", "__FEMALE", "__MALE")) {

  # Load modality-specific data
  mri_cross  <- read.csv(paste0("CrossSectionalDataPrepared/MRI_ALL_ADNI",  add_on, ".csv"))
  csf_cross  <- read.csv(paste0("CrossSectionalDataPrepared/CSF_ALL_ADNI",  add_on, ".csv"))
  fdg_cross  <- read.csv(paste0("CrossSectionalDataPrepared/FDG_ALL_ADNI",  add_on, ".csv"))
  cogn_cross <- readRDS(paste0("CrossSectionalDataPrepared/Cognitivie_Scores", add_on, "_ALL.RDS"))

  # --- MRI: hippocampal volume and entorhinal cortex thickness ---
  mri_cross <- prep_factors(mri_cross, add_on)
  formula_mri <- if (add_on == "") {
    "delta_age + Sex + APOE4 + Education + scale( ICV ) + Actual_Age + BMI"
  } else {
    "delta_age + APOE4 + Education + scale( ICV ) + Actual_Age + BMI"
  }
  mri_p_vals <- run_cross_sectional_timepoints(
    data       = mri_cross,
    targets    = c("HippVol", "entorhinal_thickness"),
    timepoints = c("bl", "m12", "m24"),
    formula    = formula_mri
  )
  # Retain baseline timepoint only
  mri_p_vals <- mri_p_vals[mri_p_vals$timepoint == "bl", ]
  all_p_vals <- rbind(all_p_vals, data.frame(
    "Phenotype"        = rownames(mri_p_vals),
    "Population"       = rep(add_on, nrow(mri_p_vals)),
    "Delta_Age_pvalue" = mri_p_vals$delta_age_p_val,
    "Beta"             = mri_p_vals$delta_age_Coef,
    "Std_Error"        = mri_p_vals$delta_age_Var_Exp
  ))

  # --- CSF: amyloid-beta and phosphorylated tau (log-transformed) ---
  csf_cross <- prep_factors(csf_cross, add_on)
  csf_cross$ABETA <- log(csf_cross$ABETA)
  csf_cross$TAU   <- log(csf_cross$TAU)
  csf_cross$PTAU  <- log(csf_cross$PTAU)
  formula_csf <- if (add_on == "") {
    "delta_age + Sex + APOE4 + Actual_Age + BMI"
  } else {
    "delta_age + APOE4 + Actual_Age + BMI"
  }
  csf_p_vals <- run_cross_sectional_timepoints(
    data       = csf_cross,
    targets    = c("ABETA", "PTAU"),
    timepoints = c("bl", "m12", "m24"),
    formula    = formula_csf
  )
  csf_p_vals <- csf_p_vals[csf_p_vals$timepoint == "bl", ]
  all_p_vals <- rbind(all_p_vals, data.frame(
    "Phenotype"        = rownames(csf_p_vals),
    "Population"       = rep(add_on, nrow(csf_p_vals)),
    "Delta_Age_pvalue" = csf_p_vals$delta_age_p_val,
    "Beta"             = csf_p_vals$delta_age_Coef,
    "Std_Error"        = csf_p_vals$delta_age_Var_Exp
  ))

  # --- FDG-PET: angular gyrus mean uptake ---
  fdg_cross <- prep_factors(fdg_cross, add_on)
  formula_fdg <- if (add_on == "") {
    "delta_age + Sex + APOE4 + Actual_Age + BMI"
  } else {
    "delta_age + APOE4 + Actual_Age + BMI"
  }
  fdg_p_vals <- run_cross_sectional_timepoints(
    data       = fdg_cross,
    targets    = c("Ang_Mean"),
    timepoints = c("bl", "m12", "m24"),
    formula    = formula_fdg
  )
  fdg_p_vals <- fdg_p_vals[fdg_p_vals$timepoint == "bl", ]
  all_p_vals <- rbind(all_p_vals, data.frame(
    "Phenotype"        = rownames(fdg_p_vals),
    "Population"       = rep(add_on, nrow(fdg_p_vals)),
    "Delta_Age_pvalue" = fdg_p_vals$delta_age_p_val,
    "Beta"             = fdg_p_vals$delta_age_Coef,
    "Std_Error"        = fdg_p_vals$delta_age_Var_Exp
  ))

  # --- Cognitive composites: memory and executive function ---
  cogn_cross_main <- merge(
    prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[, c("RID", "timepoint", "ADNI_MEM")],
    prep_factors(cogn_cross$`_ADNI_EF`,  add_on)[, c("RID", "timepoint", "ADNI_EF")],
    by = c("RID", "timepoint")
  )
  cogn_cross_main <- merge(
    cogn_cross_main,
    prep_factors(cogn_cross$`_ADNI_LAN`, add_on),
    by = c("RID", "timepoint")
  )
  formula_cogn <- if (add_on == "") {
    "delta_age + Sex + Education + Actual_Age + BMI"
  } else {
    "delta_age + Education + Actual_Age + BMI"
  }
  cogn_p_vals <- run_cross_sectional_timepoints(
    data       = cogn_cross_main,
    targets    = c("ADNI_MEM", "ADNI_EF"),
    timepoints = c("bl", "m12", "m24"),
    formula    = formula_cogn
  )
  cogn_p_vals <- cogn_p_vals[cogn_p_vals$timepoint == "bl", ]
  all_p_vals <- rbind(all_p_vals, data.frame(
    "Phenotype"        = rownames(cogn_p_vals),
    "Population"       = rep(add_on, nrow(cogn_p_vals)),
    "Delta_Age_pvalue" = cogn_p_vals$delta_age_p_val,
    "Beta"             = cogn_p_vals$delta_age_Coef,
    "Std_Error"        = cogn_p_vals$delta_age_Var_Exp
  ))
}


# -----------------------------------------------------------------------------
# Amyloid PET: Composite SUVR (Whole Cerebellum Normalised)
# -----------------------------------------------------------------------------

for (add_on in c("", "__FEMALE", "__MALE")) {

  amyloid_cross <- read.csv(paste0("CrossSectionalDataPrepared/Amyloid_ALL_ADNI",
                                   add_on, ".csv"))
  amyloid_cross <- prep_factors(amyloid_cross, add_on)

  formula_amyloid <- if (add_on == "") {
    "delta_age + Sex + APOE4 + Actual_Age + BMI"
  } else {
    "delta_age + APOE4 + Actual_Age + BMI"
  }

  amyloid_p_vals <- run_cross_sectional_timepoints(
    data       = amyloid_cross,
    targets    = c("SUMMARYSUVR_WHOLECEREBNORM"),
    timepoints = c("bl", "m24", "m48"),
    formula    = formula_amyloid
  )
  amyloid_p_vals <- amyloid_p_vals[amyloid_p_vals$timepoint == "bl", ]
  all_p_vals <- rbind(all_p_vals, data.frame(
    "Phenotype"        = rownames(amyloid_p_vals),
    "Population"       = rep(add_on, nrow(amyloid_p_vals)),
    "Delta_Age_pvalue" = amyloid_p_vals$delta_age_p_val,
    "Beta"             = amyloid_p_vals$delta_age_Coef,
    "Std_Error"        = amyloid_p_vals$delta_age_Var_Exp
  ))
}


# -----------------------------------------------------------------------------
# Save Results
# -----------------------------------------------------------------------------

write.csv(all_p_vals, "CrossSectionalResults/With_BMI_P_Vals.csv", quote = FALSE)

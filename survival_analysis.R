# =============================================================================
# Survival Analysis: MetaboAge Delta and Risk of Cognitive Progression in ADNI
# =============================================================================
#
# Description: Kaplan-Meier survival analysis examining whether metabolic age
#              delta (metaboAge - actual age) predicts time to cognitive
#              progression (to MCI or AD) in the ADNI cohort. Three progression
#              scenarios are examined:
#                (A) CN/EMCI --> AD  (type 3)
#                (B) CN      --> MCI/AD  (type 5)
#                (C) MCI     --> AD  (type 4)
#              Results are stratified by APOE e4 carrier status.
#
# Input files:
#   - ADNI_Long_Nightingale_Demographics_at_BL.csv
#   - ADNI1GO23_selectedvars_12242020.xlsx  (sheet: DX_Demo)
#   - Oct2023_Results/Predicted_Age_Whole_ADNI_Population_PolyBias.csv
#   - Oct2023_Results/Predicted_Age_Male_ADNI_Population_PolyBias.csv
#   - Oct2023_Results/Predicted_Age_Female_ADNI_Population_PolyBias.csv
#
# Output:
#   - Combined_Survival_Plot_Final.png
# =============================================================================

library(dplyr)
library(survival)
library(survminer)
library(ggsurvfit)
library(ggplot2)
library(ggpubr)


# -----------------------------------------------------------------------------
# Data Loading and Preprocessing
# -----------------------------------------------------------------------------

# Baseline demographic information
demo_adni <- read.csv("ADNI_Long_Nightingale_Demographics_at_BL.csv")
demo_adni$Sex[which(demo_adni$Sex == "Male")]   <- 1
demo_adni$Sex[which(demo_adni$Sex == "Female")] <- 0
demo_adni$Sex   <- as.factor(demo_adni$Sex)
demo_adni$APOE4 <- as.factor(demo_adni$APOE4)

# Longitudinal diagnosis data (one row per subject x timepoint)
demo_adni_D    <- readxl::read_xlsx("ADNI1GO23_selectedvars_12242020.xlsx",
                                    sheet = "DX_Demo")
diangosis_long <- data.frame("RID"      = numeric(),
                             "VISCODE2" = character(),
                             "Diagnosis" = numeric())

for (r in demo_adni_D$RID) {
  temp <- demo_adni_D[which(demo_adni_D$RID == r),
                      grep("_DXGrp", colnames(demo_adni_D))]
  diangosis_long <- rbind(
    diangosis_long,
    data.frame("RID"       = rep(r, 17),
               "VISCODE2"  = c(0, seq(6, 24, 6), seq(36, 168, 12)),
               "Diagnosis" = unlist(as.vector(temp[1, -1])))
  )
}
diangosis_long <- na.omit(diangosis_long)

# Predicted metabolic age results (whole / male / female models)
predicted_adni_brain_age        <- read.csv("Oct2023_Results/Predicted_Age_Whole_ADNI_Population_PolyBias.csv")
predicted_adni_brain_age_male   <- read.csv("Oct2023_Results/Predicted_Age_Male_ADNI_Population_PolyBias.csv")
predicted_adni_brain_age_female <- read.csv("Oct2023_Results/Predicted_Age_Female_ADNI_Population_PolyBias.csv")

predicted_adni_brain_age$delta_age        <- 0
predicted_adni_brain_age_male$delta_age   <- 0
predicted_adni_brain_age_female$delta_age <- 0


# -----------------------------------------------------------------------------
# Data Preparation: Build Survival Dataset
# -----------------------------------------------------------------------------

# Converts the predicted age data into a survival-ready data frame.
#
# Progression scenarios (type argument):
#   type 3 --> Starting as CN or EMCI; event = AD
#   type 4 --> Starting as MCI; event = AD
#   type 5 --> Starting as CN; event = MCI or AD
#
# Returns one row per subject with:
#   - delta_age       : metabolic age delta at baseline (continuous)
#   - delta_age_group : binary grouping (0 = sub-median, 1 = super-median)
#   - Progression_AD  : 1 = censored, 2 = event
#   - Time            : time to event or last follow-up (years)
#   - APOE4           : APOE e4 carrier status

update_data_to_survival_data <- function(predicted_adni_brain_age,
                                         demo_adni,
                                         demo_adni_D,
                                         type = 3) {

  # Convert visit code strings (e.g. "m24") to numeric months
  VISCODE_Converter <- function(viscode_list) {
    new_e <- c()
    for (e in viscode_list) {
      if (e == "bl") {
        new_e <- c(new_e, 0)
      } else {
        new_e <- c(new_e, strsplit(e, "m")[[1]][2])
      }
    }
    return(new_e)
  }

  # Anchor delta_age to each subject's baseline observation
  predicted_adni_brain_age$delta_age <- 0
  for (i in unique(predicted_adni_brain_age$RID)) {
    bl_idx <- which(predicted_adni_brain_age$RID == i &
                      predicted_adni_brain_age$VISCODE2 == "bl")
    predicted_adni_brain_age$delta_age[
      which(predicted_adni_brain_age$RID == i)] <-
      predicted_adni_brain_age$Corrected_Predicted_Age[bl_idx] -
      predicted_adni_brain_age$Actual_Age[bl_idx]
  }

  # Merge with demographics and longitudinal diagnosis
  complete_data <- merge(
    predicted_adni_brain_age,
    demo_adni[, c("RID", "Sex", "APOE4", "BMI", "Education")],
    by = "RID"
  )
  complete_data$VISCODE2  <- sapply(complete_data$VISCODE2, VISCODE_Converter)
  complete_data$Diagnosis <- NULL
  complete_data <- complete_data[, c("RID", "delta_age", "Sex", "APOE4", "BMI", "Education")]
  complete_data <- merge(complete_data, diangosis_long, by = "RID")

  # Identify subjects who are already at maximum disease severity at baseline
  excl_bl_AD <- unique(complete_data[which(complete_data$Diagnosis == 4 &
                                             complete_data$VISCODE2 == 0), "RID"])
  excl_bl_AD <- c(excl_bl_AD, unique(complete_data[which(complete_data$Diagnosis == 5), "RID"]))

  # Apply inclusion/exclusion criteria and recode diagnosis by progression type
  if (type == 3) {
    # CN/EMCI starting point; exclude LMCI and AD at baseline
    excl_bl_AD <- c(excl_bl_AD,
                    unique(complete_data[which(complete_data$Diagnosis == 3 &
                                                 complete_data$VISCODE2 == 0), "RID"]))
    complete_data <- complete_data[!(complete_data$RID %in% excl_bl_AD), ]
    complete_data <- na.omit(complete_data)
    complete_data$Diagnosis[complete_data$Diagnosis == 2] <- 1
    complete_data$Diagnosis[complete_data$Diagnosis == 3] <- 2
    complete_data$Diagnosis[complete_data$Diagnosis == 4] <- 3

  } else if (type == 4) {
    # MCI starting point; exclude CN at baseline
    excl_bl_AD <- c(excl_bl_AD,
                    unique(complete_data[which(complete_data$Diagnosis == 1 &
                                                 complete_data$VISCODE2 == 0), "RID"]))
    complete_data <- complete_data[!(complete_data$RID %in% excl_bl_AD), ]
    complete_data <- na.omit(complete_data)
    complete_data$Diagnosis[complete_data$Diagnosis == 2] <- 1
    complete_data$Diagnosis[complete_data$Diagnosis == 3] <- 1
    complete_data$Diagnosis[complete_data$Diagnosis == 4] <- 2

  } else if (type == 5) {
    # CN starting point; exclude MCI and AD at baseline
    excl_bl_AD <- c(excl_bl_AD,
                    unique(complete_data[which(complete_data$Diagnosis == 3 &
                                                 complete_data$VISCODE2 == 0), "RID"]),
                    unique(complete_data[which(complete_data$Diagnosis == 2 &
                                                 complete_data$VISCODE2 == 0), "RID"]))
    complete_data <- complete_data[!(complete_data$RID %in% excl_bl_AD), ]
    complete_data <- na.omit(complete_data)
    complete_data$Diagnosis[complete_data$Diagnosis == 3] <- 2
    complete_data$Diagnosis[complete_data$Diagnosis == 4] <- 2
  }

  # Convert visit months to years
  complete_data$Time      <- as.numeric(complete_data$VISCODE2) / 12
  complete_data$Diagnosis <- as.numeric(complete_data$Diagnosis)

  # Set event indicator: 1 = right-censored, 2 = event (progression)
  complete_data$Progression_AD <- 1
  max_dx <- max(unique(complete_data$Diagnosis))
  for (i in unique(complete_data$RID)) {
    if (max_dx %in% complete_data$Diagnosis[complete_data$RID == i]) {
      complete_data$Progression_AD[complete_data$RID == i] <- 2
      t_event <- min(complete_data$Time[complete_data$RID == i &
                                          complete_data$Diagnosis == max_dx])
      complete_data$Time[complete_data$RID == i] <- t_event
    }
  }

  # Collapse to one row per subject (baseline delta_age, max observed time)
  out_df <- data.frame("RID"          = numeric(),
                       "delta_age"    = numeric(),
                       "Progression_AD" = numeric(),
                       "APOE4"        = numeric(),
                       "Time"         = numeric())

  for (i in unique(complete_data$RID)) {
    idx <- which(complete_data$RID == i)
    out_df <- rbind(out_df, data.frame(
      "RID"            = i,
      "delta_age"      = unique(complete_data$delta_age[idx]),
      "Progression_AD" = unique(complete_data$Progression_AD[idx]),
      "APOE4"          = unique(complete_data$APOE4[idx]),
      "Time"           = max(unique(complete_data$Time[idx]))
    ))
  }

  # Binary grouping: above vs. at-or-below the median delta_age
  q_median <- quantile(out_df$delta_age)[3]
  out_df$delta_age_group                                    <- 1
  out_df$delta_age_group[out_df$delta_age >  q_median]     <- 1
  out_df$delta_age_group[out_df$delta_age <= q_median]     <- 0

  return(out_df)
}


# -----------------------------------------------------------------------------
# Kaplan-Meier Plotting Function (with APOE Stratification)
# -----------------------------------------------------------------------------

# Generates a three-panel KM plot for a single progression scenario:
#   - Panel 1: APOE e4 carriers
#   - Panel 2: APOE e4 non-carriers
#   - Panel 3: Combined (all subjects)
#
# Each panel shows KM curves for sub-median vs. super-median metaboAge delta,
# with Hazard Ratio and transition counts in the subtitle.
compare_to_apoe <- function(data, title = "") {

  # Collapse APOE4 = 2 (homozygous carrier) into APOE4 = 1 (carrier)
  data$APOE4 <- as.numeric(as.character(data$APOE4))
  data$APOE4[data$APOE4 == 2] <- 1
  data_w_APOE   <- data[data$APOE4 != 0, ]
  data_w_o_APOE <- data[data$APOE4 == 0, ]

  # Compute Hazard Ratio (Cox regression) and transition counts for subtitle
  get_stats_txt <- function(df) {
    if (nrow(df) == 0) return("")
    cox_fit <- coxph(Surv(Time, Progression_AD) ~ delta_age_group, data = df)
    s       <- summary(cox_fit)
    hr      <- round(s$conf.int[1], 2)
    hr_ci   <- paste0(round(s$conf.int[3], 2), "-", round(s$conf.int[4], 2))
    n_sub   <- sum(df$Progression_AD[df$delta_age_group == 0] == 2)
    n_super <- sum(df$Progression_AD[df$delta_age_group == 1] == 2)
    return(paste0("HR: ", hr, " (", hr_ci, ")\n",
                  "Transitions: sub-median=", n_sub, ", super-median=", n_super))
  }

  # Build a single KM subplot for a given data subset
  make_sub_plot <- function(df, sub_title) {
    if (nrow(df) == 0) return(ggplot() + theme_void())

    fit   <- surv_fit(Surv(Time, Progression_AD) ~ delta_age_group, data = df)
    stats <- get_stats_txt(df)

    p <- ggsurvplot(
      fit,
      data       = df,
      pval       = TRUE,
      title      = sub_title,
      subtitle   = stats,
      legend.labs = c("sub-median of metaboAge delta",
                      "super-median of metaboAge delta"),
      palette    = c("skyblue", "orange"),
      ggtheme    = theme_bw()
    )

    p$plot <- p$plot + theme(
      plot.subtitle = element_text(size = 9, face = "plain"),
      plot.title    = element_text(size = 11, face = "bold"),
      legend.text   = element_text(size = 10)
    )

    return(p$plot)
  }

  # Build the three panels
  p1 <- make_sub_plot(data_w_APOE,
                      paste0("APOE e4 carriers (n=", nrow(data_w_APOE), ")"))
  p2 <- make_sub_plot(data_w_o_APOE,
                      paste0("APOE e4 non-carriers (n=", nrow(data_w_o_APOE), ")"))
  p3 <- make_sub_plot(data,
                      paste0("Combined Population (n=", nrow(data), ")"))

  # Arrange panels horizontally and add a row title
  arranged_row  <- ggarrange(p1, p2, p3, ncol = 3, nrow = 1,
                             common.legend = TRUE, legend = "top")
  annotated_row <- annotate_figure(arranged_row,
                                   top = text_grob(title, face = "bold", size = 14))
  return(annotated_row)
}


# -----------------------------------------------------------------------------
# Main Analysis: Generate Combined Survival Plot (3 Progression Scenarios)
# -----------------------------------------------------------------------------

titles   <- c(
  "(A) Progression of CN/EMCI to AD",
  "(B) Progression of CN to MCI/AD",
  "(C) Progression of MCI to AD"
)
type_map <- c(3, 5, 4)

plots_list <- list()

for (i in 1:3) {
  t_val     <- type_map[i]
  temp_data <- update_data_to_survival_data(
    predicted_adni_brain_age, demo_adni, demo_adni_D, type = t_val
  )
  # Restrict follow-up to 10 years; exclude zero-time observations
  temp_data <- temp_data %>% filter(Time <= 10 & Time > 0)

  plots_list[[i]] <- compare_to_apoe(temp_data, title = titles[i])
}


# -----------------------------------------------------------------------------
# Assemble and Save Final Combined Figure
# -----------------------------------------------------------------------------

final_combined_plot <- ggarrange(plotlist = plots_list, ncol = 1, nrow = 3)

ggsave("Combined_Survival_Plot_Final.png", final_combined_plot,
       width = 15, height = 15, dpi = 300)

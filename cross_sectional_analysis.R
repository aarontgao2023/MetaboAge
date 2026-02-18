# =============================================================================
# Cross-Sectional Analysis: MetaboAge Delta vs. Brain and Cognitive Phenotypes
# =============================================================================
#
# Description: Generates cross-sectional scatter plots examining the association
#              between metabolic age delta (metaboAge - actual age) and brain
#              imaging biomarkers (CSF, FDG-PET, MRI) and cognitive scores in
#              the ADNI cohort. Analyses are stratified by sex.
#
# Input files:
#   - CrossSectionalDataPrepared/[modality]_ALL_ADNI[_MALE|_FEMALE].csv
#   - CrossSectionalDataPrepared/Cognitivie_Scores[_MALE|_FEMALE]_ALL.RDS
#   - FDR_Corrected_All_pvalues.csv
#
# Output:
#   - CrossSectionalPlots_Final.png
# =============================================================================

library(ggpubr)
library(hrbrthemes)
library(lme4)
library(dplyr)


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Prepare common factor variables and compute delta_age; filter to
# timepoints bl, m12, m24.
prep_factors <- function(data, add_on) {
  if (add_on == "") {
    data$Sex <- as.factor(data$Sex)
  }
  data$APOE4 <- as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == "VISCODE2")] <- "timepoint"
  data$delta_age <- data$Corrected_Predicted_Age - data$Actual_Age
  data <- data %>% filter(timepoint %in% c("bl", "m12", "m24"))
  return(data)
}

# Generate a cross-sectional scatter plot with regression line.
# delta_age is binarised into accelerated (>0) vs. decelerated (<0) aging.
#
# Parameters:
#   data            - prepared data frame
#   formula         - right-hand-side string for lmer (without outcome ~)
#   target          - column name of the outcome variable
#   thresh_year     - maximum year for data filtering
#   target_label    - y-axis label
#   include_scatter - overlay raw data points (default FALSE)
#   p_value_given   - FDR-corrected p-value for the interaction term
make_plots_crossSectional <- function(data, formula, target, thresh_year,
                                      target_label, include_scatter = FALSE,
                                      p_value_given) {

  data <- data %>% filter(year <= thresh_year)

  # Extract variable names referenced in the formula
  columns_needed <- unique(strsplit(formula, " ")[[1]][
    sapply(strsplit(formula, " ")[[1]], function(str) {
      !grepl("[^A-Za-z0-9_ ]", str)
    })
  ])
  columns_needed <- columns_needed[suppressWarnings(is.na(as.numeric(columns_needed)))]
  data <- data[, c(columns_needed, target)]
  data <- na.omit(data)

  # Fit mixed-effects model and extract population-level predictions
  target_formula <- paste(target, formula, sep = " ~ ")
  og <- lmer(as.formula(target_formula), data, REML = FALSE)

  # Assign quartile labels; then collapse to binary (Q1 = decelerated, Q4 = accelerated)
  Quartiles <- quantile(data[, colnames(data) == "delta_age"],
                        probs = c(0, 0.25, 0.5, 0.75, 1))
  data$Quartiles <- cut(data[, colnames(data) == "delta_age"],
                        breaks = Quartiles, labels = c("Q1", "Q2", "Q3", "Q4"),
                        include.lowest = TRUE)

  pred1 <- predict(og, re.form = NA, se.fit = TRUE, nsim = 500)
  temp_model_data <- data
  temp_model_data$boot_fit <- pred1$fit

  data$Quartiles[data$delta_age > 0] <- "Q4"
  data$Quartiles[data$delta_age < 0] <- "Q1"

  plot_Q1Q4 <- temp_model_data %>% filter(Quartiles == "Q1" | Quartiles == "Q4")

  # Population subtitle
  if (add_on == "") {
    title_sub <- "Whole Population"
  } else if (add_on == "__MALE") {
    title_sub <- "Male Population"
  } else {
    title_sub <- "Female Population"
  }

  # Format p-value; bold if significant
  p_val <- p_value_given
  if (p_val < 0.001) {
    p_val_out <- "< 0.001"
  } else {
    p_val_out <- as.character(round(p_val, 4))
  }

  if (p_val < 0.05) {
    p_val_out <- as.formula(paste0("'P-value '~ bold('", p_val_out, "')"))
  } else {
    p_val_out <- paste0("P-value ", p_val_out)
  }

  if (!include_scatter) {
    w <- ggplot(plot_Q1Q4, aes(x = year, y = boot_fit, col = Quartiles)) +
      geom_smooth(aes(year, boot_fit), linewidth = 1.5, method = lm, se = TRUE) +
      ylab(target_label) + xlab("Year") +
      labs(title = title_sub, subtitle = p_val_out) +
      scale_color_manual(values = c("skyblue", "coral1"), name = "",
                         labels = c("metaboAge < 0", "metaboAge > 0"))
  } else {
    w <- ggplot(plot_Q1Q4, aes(x = year, y = !!sym(target), col = Quartiles)) +
      geom_point() +
      geom_smooth(aes(year, boot_fit), linewidth = 1.5, method = lm, se = TRUE) +
      ylab(target_label) + xlab("Year") +
      labs(title = title_sub, subtitle = p_val_out) +
      scale_color_manual(values = c("skyblue", "coral1"), name = "",
                         labels = c("metaboAge < 0", "metaboAge > 0"))
  }

  return(w)
}


# -----------------------------------------------------------------------------
# Data Paths and Preprocessing Functions
# -----------------------------------------------------------------------------

paths <- list(
  "CrossSectionalDataPrepared/CSF_ALL_ADNI",
  "CrossSectionalDataPrepared/FDG_ALL_ADNI",
  "CrossSectionalDataPrepared/MRI_ALL_ADNI",
  "CrossSectionalDataPrepared/MRI_ALL_ADNI",
  "CrossSectionalDataPrepared/Cognitivie_Scores"
)

# Each function accepts a path and uses the outer-loop variable `add_on`
# to select the sex-stratified file variant.
preprocess <- list(

  # CSF biomarkers: ABETA, TAU, PTAU (log-transformed)
  function(x) {
    csf_cross <- read.csv(paste0(x, add_on, ".csv"))
    csf_cross <- prep_factors(csf_cross, add_on)
    csf_cross$ABETA <- log(csf_cross$ABETA)
    csf_cross$TAU   <- log(csf_cross$TAU)
    csf_cross$PTAU  <- log(csf_cross$PTAU)
    # Keep only subjects with more than one observation
    subs_multi <- names(table(csf_cross$RID)[table(csf_cross$RID) != 1])
    csf_cross <- csf_cross %>% filter(RID %in% subs_multi)
    return(csf_cross)
  },

  # FDG-PET brain glucose metabolism
  function(x) {
    fdg_cross <- read.csv(paste0(x, add_on, ".csv"))
    fdg_cross <- prep_factors(fdg_cross, add_on)
    return(fdg_cross)
  },

  # MRI: entorhinal cortex thickness
  function(x) {
    mri_cross <- read.csv(paste0(x, add_on, ".csv"))
    mri_cross <- prep_factors(mri_cross, add_on)
    return(mri_cross)
  },

  # MRI: hippocampal volume
  function(x) {
    mri_cross <- read.csv(paste0(x, add_on, ".csv"))
    mri_cross <- prep_factors(mri_cross, add_on)
    return(mri_cross)
  },

  # Cognitive composite: memory score
  function(x) {
    cogn_cross <- readRDS(paste0(x, add_on, "_ALL.RDS"))
    cogn_cross_main <- merge(
      prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[, c("RID", "timepoint", "ADNI_MEM")],
      prep_factors(cogn_cross$`_ADNI_EF`, add_on),
      by = c("RID", "timepoint")
    )
    cogn_cross_main <- na.omit(cogn_cross_main)
    return(cogn_cross_main)
  }
)


# -----------------------------------------------------------------------------
# Load FDR-Corrected P-Values
# -----------------------------------------------------------------------------

cross_sectional_p_values <- read.csv("FDR_Corrected_All_pvalues.csv")
colnames(cross_sectional_p_values) <- c("Phenotype", "p_value_whole",
                                        "p_value_male", "p_value_female")
cross_sectional_p_values <- cross_sectional_p_values %>%
  dplyr::filter(Phenotype %in% c("CSF PTAU", "PET Angular Mean",
                                 "MRI Entorhinal Thickness",
                                 "MRI Hippocampal Volume",
                                 "Cognition Composite Memory"))

targets      <- c("CSF PTAU", "PET Angular Mean", "MRI Entorhinal Thickness",
                  "MRI Hippocampal Volume", "Cognition Composite Memory")
targets_name <- c("PTAU", "Ang_Mean", "entorhinal_thickness", "HippVol", "ADNI_MEM")
ylabel       <- c("CSF pTau181", "Brain glucose metabolism",
                  "Entorhinal thickness", "Hippocampal volume", "Memory Score")

rownames(cross_sectional_p_values) <- cross_sectional_p_values$Phenotype
fdr_corrected_p_value_whole  <- cross_sectional_p_values[targets, 2]
fdr_corrected_p_value_male   <- cross_sectional_p_values[targets, 3]
fdr_corrected_p_value_female <- cross_sectional_p_values[targets, 4]


# -----------------------------------------------------------------------------
# Generate Cross-Sectional Plots (Whole, Male, Female Populations)
# -----------------------------------------------------------------------------

cross_plots <- list()

for (i in 1:length(targets)) {
  for (add_on in c("", "__MALE", "__FEMALE")) {

    data <- lapply(paths[i], preprocess[[i]])[[1]]
    data <- data %>% filter(delta_age < 30)

    # Select FDR-corrected p-value for the current population stratum
    if (add_on == "") {
      p_val <- fdr_corrected_p_value_whole[i]
    } else if (add_on == "__MALE") {
      p_val <- fdr_corrected_p_value_male[i]
    } else {
      p_val <- fdr_corrected_p_value_female[i]
    }

    # Fit linear model to extract regression coefficient (beta)
    temp_y <- targets_name[i]
    colnames(data)[which(colnames(data) == temp_y)] <- "Temp_Y"
    fit_model     <- lm(Temp_Y ~ delta_age, data = data)
    stats_summary <- summary(fit_model)$coefficients
    beta_val      <- round(stats_summary["delta_age", "Estimate"], 3)
    p_label       <- if (p_val < 0.001) " < 0.001" else paste0(" = ", round(p_val, 4))

    # Build subtitle: bold if significant
    full_text <- paste0("P", p_label, " | \u03b2 = ", beta_val)
    if (p_val < 0.05) {
      stat_lab <- bquote(bold(.(full_text)))
    } else {
      stat_lab <- full_text
    }

    cross_plots[[paste0(targets[i], add_on)]] <-
      ggplot(data, aes(x = delta_age, y = Temp_Y)) +
      geom_point(alpha = 0.3, color = "grey40", size = 1.2) +
      geom_smooth(method = lm, color = "coral1", fill = "skyblue",
                  se = TRUE, linewidth = 1.2) +
      labs(title = if (add_on == "") "Whole Population"
                   else if (add_on == "__MALE") "Male Population"
                   else "Female Population",
           subtitle = stat_lab) +
      xlab("Delta Age") +
      ylab(ylabel[[i]]) +
      theme_ipsum() +
      theme(
        axis.title.y     = element_text(size = 14),
        axis.title.x     = element_text(size = 14),
        plot.title       = element_text(size = 18, face = "bold"),
        plot.subtitle    = element_text(size = 13, color = "black"),
        panel.grid.minor = element_blank()
      )
  }
}


# -----------------------------------------------------------------------------
# Assemble and Save Final Figure (5 phenotypes x 3 sex strata)
# -----------------------------------------------------------------------------

final_fig <- ggarrange(
  plotlist      = cross_plots,
  labels        = c("A", "", "", "B", "", "", "C", "", "", "D", "", "", "E", "", ""),
  ncol          = 3,
  nrow          = 5,
  common.legend = TRUE,
  legend        = "bottom"
)

ggsave("CrossSectionalPlots_Final.png", final_fig,
       width = 15, height = 25, units = "in", dpi = 300)

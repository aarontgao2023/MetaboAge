# =============================================================================
# Longitudinal Analysis: MetaboAge Delta and Progression of Brain/Cognitive Phenotypes
# =============================================================================
#
# Description: Fits linear mixed-effects models to examine whether metabolic
#              age delta (metaboAge - actual age) predicts longitudinal change
#              in brain imaging biomarkers and cognitive composite scores in the
#              ADNI cohort. Analyses are stratified by sex (whole, male, female).
#
# Input files:
#   - LongitudianlDataPrepared/[modality]_Longitudinal_ADNI[_MALE|_FEMALE].csv
#   - LongitudianlDataPrepared/Cognitivie_Scores[_MALE|_FEMALE].RDS
#   - FDR_Corrected_Longitudinal_pvalues.csv
#
# Output:
#   - longitudinal_effect_plots.png
# =============================================================================

source("fun.R")
library(rempsyc)
library(dplyr)
library(ggpubr)
library(lme4)
library(ggtext)


# -----------------------------------------------------------------------------
# Helper Functions
# -----------------------------------------------------------------------------

# Prepare common factor variables (Sex, APOE4) and rename VISCODE2 -> timepoint.
prep_factors <- function(data, add_on) {
  if (add_on == "") {
    data$Sex <- as.factor(data$Sex)
  }
  data$APOE4 <- as.factor(data$APOE4)
  colnames(data)[which(colnames(data) == "VISCODE2")] <- "timepoint"
  return(data)
}

# Fit a longitudinal mixed-effects model and generate trajectory plots
# stratified by metaboAge delta quartile (Q1 vs Q4).
#
# The interaction term delta_age:year is the primary term of interest,
# testing whether metabolic age acceleration modifies the rate of change
# in the outcome over time.
#
# Parameters:
#   data            - prepared data frame
#   formula         - right-hand-side lmer formula string (without outcome ~)
#   target          - column name of the outcome variable
#   thresh_year     - maximum follow-up year to include
#   target_label    - y-axis label
#   include_scatter - overlay observed data points (default FALSE)
#   p_value_given   - FDR-corrected p-value for the delta_age x year interaction
make_plots_longitudinal <- function(data, formula, target, thresh_year,
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

  target_formula <- paste(target, formula, sep = " ~ ")

  # Fit mixed-effects model (ML for comparison, not REML)
  og <- lmer(as.formula(target_formula), data, REML = FALSE)

  # Extract delta_age x year interaction coefficient and standard error
  fixed_effects <- summary(og)$coefficients
  inter_row     <- grep("delta_age:year|year:delta_age",
                        rownames(fixed_effects), value = TRUE)
  beta_val <- round(fixed_effects[inter_row, "Estimate"], 4)
  se_val   <- round(fixed_effects[inter_row, "Std. Error"], 4)

  # Format subtitle with FDR-corrected p-value and interaction beta (HTML for ggtext)
  p_val   <- p_value_given
  p_label <- if (p_val < 0.001) "P < 0.001" else paste0("P = ", round(p_val, 3))
  content <- paste0(p_label, "<br>Interaction \u03b2 = ", beta_val)

  if (p_val < 0.05) {
    stat_lab <- paste0("<b>", content, "</b>")
  } else {
    stat_lab <- content
  }

  # Assign quartile groups; model-predicted trajectories at population level
  Quartiles    <- quantile(data$delta_age, probs = c(0, 0.25, 0.5, 0.75, 1))
  data$Quartiles <- cut(data$delta_age, breaks = Quartiles,
                        labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  data$boot_fit <- predict(og, re.form = NA)

  # Keep only extreme quartiles for visualisation
  plot_Q1Q4 <- data %>%
    mutate(Group = ifelse(Quartiles %in% c("Q1", "Q2"), "sub-median", "super-median")) %>%
    filter(Quartiles %in% c("Q1", "Q4"))

  if (add_on == "") {
    title_sub <- "Whole Population"
  } else if (add_on == "__MALE") {
    title_sub <- "Male Population"
  } else {
    title_sub <- "Female Population"
  }

  p <- ggplot(plot_Q1Q4, aes(x = year, y = boot_fit, col = Group)) +
    geom_smooth(method = lm, linewidth = 1.5, se = TRUE) +
    ylab(target_label) + xlab("Year") +
    labs(title = title_sub, subtitle = stat_lab) +
    scale_color_manual(values = c("skyblue", "coral1"), name = "") +
    theme_minimal() +
    theme(
      plot.title    = element_text(size = 18, face = "bold"),
      plot.subtitle = element_markdown(size = 12, color = "black"),
      axis.title    = element_text(size = 14),
      axis.title.y  = element_text(size = 11, face = "plain", color = "grey20",
                                   family = "sans", margin = margin(r = 10))
    )

  if (include_scatter) {
    p <- p + geom_point(aes(y = !!sym(target)), alpha = 0.2)
  }

  return(p)
}


# -----------------------------------------------------------------------------
# Data Paths and Preprocessing Functions
# -----------------------------------------------------------------------------

longitudinal_p_values <- read.csv("FDR_Corrected_Longitudinal_pvalues.csv")

paths <- list(
  "LongitudianlDataPrepared/Amyloid_Longitudinal_ADNI",
  "LongitudianlDataPrepared/CSF_Longitudinal_ADNI",
  "LongitudianlDataPrepared/MRI_Longitudinal_ADNI",
  "LongitudianlDataPrepared/Cognitivie_Scores",
  "LongitudianlDataPrepared/Cognitivie_Scores"
)

# Each function accepts a path and uses the outer-loop variable `add_on`
# to select the sex-stratified file variant.
preprocess <- list(

  # Amyloid PET (SUVR composite); exclude half-year and irregular timepoints
  function(x) {
    amyloid_cross <- read.csv(paste0(x, add_on, ".csv"))
    amyloid_cross <- amyloid_cross[which(amyloid_cross$year %% 1 != 0.5), ]
    amyloid_cross <- amyloid_cross[which(!amyloid_cross$year %in% c(1, 3, 10)), ]
    amyloid_cross <- prep_factors(amyloid_cross, add_on)
    # Require > 2 observations per subject
    subs_multi <- names(table(amyloid_cross$RID)[table(amyloid_cross$RID) > 2])
    subs_multi <- subs_multi[!subs_multi %in% c(668, 677)]
    amyloid_cross <- amyloid_cross %>% filter(RID %in% subs_multi)
    return(amyloid_cross)
  },

  # CSF biomarkers: ABETA, TAU, PTAU (log-transformed)
  function(x) {
    csf_cross <- read.csv(paste0(x, add_on, ".csv"))
    if (add_on != "") {
      csf_cross <- csf_cross %>% filter(BATCH == "MEDIAN")
    }
    csf_cross <- prep_factors(csf_cross, add_on)
    csf_cross$ABETA <- log(csf_cross$ABETA)
    csf_cross$TAU   <- log(csf_cross$TAU)
    csf_cross$PTAU  <- log(csf_cross$PTAU)
    # Remove subjects with only a single observation
    subs_multi <- names(table(csf_cross$RID)[table(csf_cross$RID) != 1])
    csf_cross <- csf_cross %>% filter(RID %in% subs_multi)
    return(csf_cross)
  },

  # MRI: entorhinal cortex thickness
  function(x) {
    mri_cross <- read.csv(paste0(x, add_on, ".csv"))
    mri_cross <- prep_factors(mri_cross, add_on)
    return(mri_cross)
  },

  # Cognitive composite: executive function
  function(x) {
    cogn_cross <- readRDS(paste0(x, add_on, ".RDS"))
    cogn_cross_main <- merge(
      prep_factors(cogn_cross$`_ADNI_MEM`, add_on)[, c("RID", "timepoint", "ADNI_MEM")],
      prep_factors(cogn_cross$`_ADNI_EF`, add_on),
      by = c("RID", "timepoint")
    )
    cogn_cross_main <- na.omit(cogn_cross_main)
    return(cogn_cross_main)
  },

  # Cognitive composite: memory
  function(x) {
    cogn_cross <- readRDS(paste0(x, add_on, ".RDS"))
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
# Mixed-Effects Model Formulas
# -----------------------------------------------------------------------------

# Whole-population formulas include Sex as a covariate
formula_list_whole <- list(
  "1 + APOE4 + delta_age * year + Sex + Actual_Age + ( 1 + year | RID )",
  "1 + Actual_Age + APOE4 + delta_age * year + Sex + ( 1 + year | RID )",
  "1 + Actual_Age + Sex + APOE4 + delta_age * year + scale( ICV ) + Education + ( 1 + year | RID )",
  "1 + Sex + Actual_Age + Education + delta_age * year + ( 1 + year | RID )",
  "1 + Sex + Actual_Age + Education + delta_age * year + ( 1 + year | RID )"
)

# Sex-stratified formulas omit Sex as a covariate
formula_list_male_female <- list(
  "1 + APOE4 + delta_age * year + Actual_Age + ( 1 + year | RID )",
  "1 + Actual_Age + APOE4 + delta_age * year + ( 1 + year | RID )",
  "1 + Actual_Age + APOE4 + delta_age * year + scale( ICV ) + Education + ( 1 + year | RID )",
  "1 + Actual_Age + Education + delta_age * year + ( 1 + year | RID )",
  "1 + Actual_Age + Education + delta_age * year + ( 1 + year | RID )"
)

years_threshold <- c(4, 3, 8, 8, 8)
targets      <- c("SUMMARYSUVR_COMPOSITE_REFNORM", "PTAU", "entorhinal_thickness",
                  "ADNI_EF", "ADNI_MEM")
target_names <- c("Amyloid SUVR", "CSF PTAU", "Entorhinal Cortical Thickness",
                  "Composite Scores for Executive Function",
                  "Composite Scores for Memory")

# Load FDR-corrected p-values
rownames(longitudinal_p_values) <- longitudinal_p_values$Phenotype
phenotype_keys <- c("Amyloid SUVR Composite Score", "CSF PTAU",
                    "MRI Entorhinal Thickness",
                    "Cognition Composite Executive Function",
                    "Cognition Composite Memory")
fdr_corrected_p_value_whole  <- longitudinal_p_values[phenotype_keys, 2]
fdr_corrected_p_value_male   <- longitudinal_p_values[phenotype_keys, 3]
fdr_corrected_p_value_female <- longitudinal_p_values[phenotype_keys, 4]


# -----------------------------------------------------------------------------
# Generate Longitudinal Trajectory Plots
# -----------------------------------------------------------------------------

long_effect_plots           <- list()
long_effect_plots_w_scatter <- list()

for (i in 1:length(targets)) {
  for (add_on in c("", "__MALE", "__FEMALE")) {

    data    <- lapply(paths[i], preprocess[[i]])[[1]]
    formula <- if (add_on == "") formula_list_whole[[i]] else formula_list_male_female[[i]]

    p_value_corr <- if (add_on == "") {
      fdr_corrected_p_value_whole[i]
    } else if (add_on == "__MALE") {
      fdr_corrected_p_value_male[i]
    } else {
      fdr_corrected_p_value_female[i]
    }

    key <- paste0(target_names[i], add_on)

    long_effect_plots[[key]] <- make_plots_longitudinal(
      data          = data,
      formula       = formula,
      target        = targets[[i]],
      thresh_year   = years_threshold[[i]],
      target_label  = target_names[i],
      p_value_given = p_value_corr
    )

    long_effect_plots_w_scatter[[key]] <- make_plots_longitudinal(
      data            = data,
      formula         = formula,
      target          = targets[[i]],
      thresh_year     = years_threshold[[i]],
      target_label    = target_names[i],
      include_scatter = TRUE,
      p_value_given   = p_value_corr
    )
  }
}


# -----------------------------------------------------------------------------
# Assemble and Save Figures
# -----------------------------------------------------------------------------

# Main figure: MRI entorhinal thickness, executive function, memory (rows 3–5)
annotate_figure(ggarrange(
  long_effect_plots[[7]],  long_effect_plots[[8]],  long_effect_plots[[9]],
  long_effect_plots[[10]], long_effect_plots[[11]], long_effect_plots[[12]],
  long_effect_plots[[13]], long_effect_plots[[14]], long_effect_plots[[15]],
  labels = c("A", "", "", "B", "", "", "C", "", ""),
  common.legend = TRUE, ncol = 3, nrow = 3
))
ggsave("longitudinal_effect_plots.png", width = 10, height = 10, units = "in")


# ============================================================
# STEP 2: MONGOLIA CTmax MODEL SELECTION AND ANALYSIS
# ============================================================
#
# INPUT:
#   CTmax_Merged_Data.csv
#
# MODEL VARIABLES:
#   ctmax                   CTmax (degrees Celsius)
#   weight_g                Body mass (g)
#   parasite_burden_percent Parasite mass / body mass * 100
#   start_temp              Starting temperature (degrees C)
#   species_group           Grayling / Lenok
#   habitat                 Lake / River
#
# NO log transformations.
# NO z-score standardization.
#
# All models use the SAME complete-case dataset.
#
# OUTPUTS:
#   Model rankings and AICc comparisons
#   Competitive models
#   Coefficient estimates and confidence intervals
#   Parasite diagnostics
#   Model diagnostics
#   Figures
#   Journal-style model comparison table
#   Saved model objects
#
# ============================================================


# 1. PACKAGES ------------------------------------------------

required_packages <- c("tidyverse", "gt")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(gt)


# 2. FILE PATHS ----------------------------------------------

input_file <- "CTmax_Merged_Data.csv"
output_dir <- "CTmax_Model_Results"

if (!file.exists(input_file)) {
  stop(
    "Cannot find CTmax_Merged_Data.csv. ",
    "Run Step 1 first and check your working directory."
  )
}

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# 3. IMPORT DATA ---------------------------------------------

raw_data <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

required_columns <- c(
  "trial",
  "fish_id",
  "species",
  "site_type",
  "ctmax",
  "weight_g",
  "start_temp"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!"parasite_burden_percent" %in% names(raw_data)) {

  if (!"parasite_burden" %in% names(raw_data)) {
    stop(
      "Neither parasite_burden_percent nor ",
      "parasite_burden was found."
    )
  }

  raw_data <- raw_data %>%
    mutate(
      parasite_burden_percent =
        as.numeric(parasite_burden) * 100
    )
}


# 4. PREPARE ANALYSIS DATA -----------------------------------

analysis_data <- raw_data %>%
  mutate(
    trial = suppressWarnings(as.numeric(trial)),
    ctmax = suppressWarnings(as.numeric(ctmax)),
    weight_g = suppressWarnings(as.numeric(weight_g)),
    start_temp = suppressWarnings(as.numeric(start_temp)),
    parasite_burden_percent =
      suppressWarnings(as.numeric(parasite_burden_percent)),

    fish_id = trimws(as.character(fish_id)),

    habitat = factor(
      tolower(trimws(as.character(site_type))),
      levels = c("lake", "river")
    ),

    species_group = case_when(
      grepl(
        "grayling",
        as.character(species),
        ignore.case = TRUE
      ) ~ "Grayling",

      grepl(
        "lenok",
        as.character(species),
        ignore.case = TRUE
      ) ~ "Lenok",

      TRUE ~ NA_character_
    ),

    species_group = factor(
      species_group,
      levels = c("Grayling", "Lenok")
    )
  ) %>%
  filter(
    !(trial %in% c(1, 26, 67, 84, 108)),
    is.finite(ctmax),
    !is.na(species_group),
    !is.na(habitat),
    is.finite(weight_g),
    weight_g > 0,
    is.finite(parasite_burden_percent),
    parasite_burden_percent >= 0,
    is.finite(start_temp),
    !is.na(fish_id),
    nzchar(fish_id)
  ) %>%
  droplevels()


# 5. VALIDATE ANALYSIS DATA ----------------------------------

if (nrow(analysis_data) < 30) {
  stop(
    "Fewer than 30 complete observations remain. ",
    "Check the merge and missing-data summary."
  )
}

group_counts <- analysis_data %>%
  count(
    species_group,
    habitat,
    name = "n"
  )

cat("\nSPECIES BY HABITAT SAMPLE SIZES:\n")

print(
  group_counts,
  n = Inf,
  width = Inf
)

if (
  nrow(group_counts) != 4 ||
  any(group_counts$n < 5)
) {
  stop(
    "All four species-by-habitat groups must ",
    "contain at least five fish."
  )
}

duplicate_fish <- analysis_data %>%
  count(habitat, fish_id, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_fish) > 0) {

  print(duplicate_fish, n = Inf)

  stop(
    "Duplicate fish IDs detected within habitat. ",
    "Check the merged dataset."
  )
}

cat(
  "\nTotal observations in all candidate models:",
  nrow(analysis_data),
  "\n"
)

write_csv(
  analysis_data,
  file.path(
    output_dir,
    "CTmax_Analysis_Dataset.csv"
  )
)

write_csv(
  group_counts,
  file.path(
    output_dir,
    "Species_Habitat_Sample_Sizes.csv"
  )
)


# 6. PARASITE DISTRIBUTION DIAGNOSTICS ------------------------

parasite_diagnostics <- analysis_data %>%
  group_by(
    species_group,
    habitat
  ) %>%
  summarise(
    n = n(),

    parasite_positive =
      sum(parasite_burden_percent > 0),

    parasite_zero =
      sum(parasite_burden_percent == 0),

    prevalence_percent =
      100 * mean(parasite_burden_percent > 0),

    mean_burden_percent =
      mean(parasite_burden_percent),

    sd_burden_percent =
      sd(parasite_burden_percent),

    min_burden_percent =
      min(parasite_burden_percent),

    max_burden_percent =
      max(parasite_burden_percent),

    .groups = "drop"
  )

cat("\nPARASITE DISTRIBUTION:\n")

print(
  parasite_diagnostics,
  n = Inf,
  width = Inf
)

write_csv(
  parasite_diagnostics,
  file.path(
    output_dir,
    "Parasite_Distribution_Diagnostics.csv"
  )
)


# 7. DEFINE ALL 19 CANDIDATE MODELS --------------------------
#
# Models preserve the previous candidate-model structure.
# Only original biological units are used.

candidate_formulas <- list(

  M01 = ctmax ~ species_group + habitat,

  M02 = ctmax ~ species_group * habitat,

  M03 = ctmax ~ species_group * habitat +
    weight_g,

  M04 = ctmax ~ species_group * habitat +
    weight_g + habitat:weight_g,

  M05 = ctmax ~ species_group * habitat +
    weight_g + species_group:weight_g,

  M06 = ctmax ~ species_group * habitat * weight_g,

  M07 = ctmax ~ species_group * habitat +
    parasite_burden_percent,

  M08 = ctmax ~ species_group * habitat +
    parasite_burden_percent +
    habitat:parasite_burden_percent,

  M09 = ctmax ~ species_group * habitat +
    weight_g + parasite_burden_percent,

  M10 = ctmax ~ species_group * habitat +
    weight_g + habitat:weight_g +
    parasite_burden_percent,

  M11 = ctmax ~ species_group * habitat +
    weight_g + species_group:weight_g +
    parasite_burden_percent,

  M12 = ctmax ~ species_group * habitat +
    weight_g + parasite_burden_percent +
    habitat:parasite_burden_percent,

  M13 = ctmax ~ species_group * habitat +
    start_temp,

  M14 = ctmax ~ species_group * habitat +
    weight_g + start_temp,

  M15 = ctmax ~ species_group * habitat +
    parasite_burden_percent + start_temp,

  M16 = ctmax ~ species_group * habitat +
    weight_g + parasite_burden_percent +
    start_temp,

  M17 = ctmax ~ species_group * habitat +
    weight_g + habitat:weight_g +
    parasite_burden_percent + start_temp,

  M18 = ctmax ~ species_group * habitat +
    weight_g + species_group:weight_g +
    parasite_burden_percent + start_temp,

  M19 = ctmax ~ species_group * habitat * weight_g +
    parasite_burden_percent + start_temp
)


# 8. FIT ALL ESTIMABLE MODELS --------------------------------
#
# Model IDs and fitted model objects are kept separate.
#
# This avoids the previous:
# Error in formula.character(model)

fitted_models <- list()
skipped_list <- list()

for (model_id in names(candidate_formulas)) {

  current_formula <- candidate_formulas[[model_id]]

  fitted_object <- tryCatch(
    stats::lm(
      formula = current_formula,
      data = analysis_data,
      na.action = na.fail
    ),
    error = function(e) e
  )

  if (inherits(fitted_object, "error")) {

    skipped_list[[model_id]] <- data.frame(
      model_id = model_id,
      reason = conditionMessage(fitted_object)
    )

    next
  }

  if (
    fitted_object$rank <
    length(stats::coef(fitted_object))
  ) {

    skipped_list[[model_id]] <- data.frame(
      model_id = model_id,
      reason = "Rank-deficient model"
    )

    next
  }

  fitted_models[[model_id]] <- fitted_object
}

skipped_models <- if (length(skipped_list) > 0) {

  dplyr::bind_rows(skipped_list)

} else {

  tibble(
    model_id = character(),
    reason = character()
  )
}

write_csv(
  skipped_models,
  file.path(
    output_dir,
    "Skipped_Models.csv"
  )
)

cat(
  "\nSuccessfully fitted models:",
  length(fitted_models),
  "of",
  length(candidate_formulas),
  "\n"
)

if (nrow(skipped_models) > 0) {

  cat("\nSKIPPED MODELS:\n")

  print(
    skipped_models,
    n = Inf,
    width = Inf
  )
}

if (length(fitted_models) == 0) {
  stop("No candidate models fitted successfully.")
}


# 9. AICc CALCULATION ----------------------------------------

calculate_aicc <- function(fitted_object) {

  sample_n <- stats::nobs(fitted_object)

  parameter_k <- attr(
    stats::logLik(fitted_object),
    "df"
  )

  if (sample_n <= parameter_k + 1) {
    return(Inf)
  }

  stats::AIC(fitted_object) +
    (
      2 * parameter_k * (parameter_k + 1)
    ) /
    (
      sample_n - parameter_k - 1
    )
}


# 10. BUILD MODEL RANKING TABLE ------------------------------
#
# IMPORTANT:
# All values are calculated BEFORE data.frame().
#
# This prevents model IDs from being interpreted as formulas.

ranking_list <- list()

for (model_id in names(fitted_models)) {

  fitted_object <- fitted_models[[model_id]]

  formula_text <- paste(
    deparse(
      candidate_formulas[[model_id]],
      width.cutoff = 500
    ),
    collapse = " "
  )

  sample_n <- stats::nobs(fitted_object)

  parameter_k <- attr(
    stats::logLik(fitted_object),
    "df"
  )

  model_aic <- stats::AIC(fitted_object)

  model_aicc <- calculate_aicc(fitted_object)

  model_summary <- summary(fitted_object)

  model_r2 <- model_summary$r.squared

  model_adjusted_r2 <- model_summary$adj.r.squared

  ranking_list[[model_id]] <- data.frame(
    model_id = model_id,
    model_structure = formula_text,
    n = sample_n,
    K = parameter_k,
    AIC = model_aic,
    AICc = model_aicc,
    R2 = model_r2,
    adjusted_R2 = model_adjusted_r2,
    stringsAsFactors = FALSE
  )
}

model_rankings <- dplyr::bind_rows(
  ranking_list
) %>%
  filter(
    is.finite(AICc)
  ) %>%
  arrange(
    AICc
  )

if (nrow(model_rankings) == 0) {
  stop(
    "No fitted candidate models have finite AICc."
  )
}

model_rankings <- model_rankings %>%
  mutate(
    delta_AICc = AICc - min(AICc),

    relative_likelihood =
      exp(-0.5 * delta_AICc),

    akaike_weight =
      relative_likelihood /
      sum(relative_likelihood),

    cumulative_weight =
      cumsum(akaike_weight),

    rank = row_number()
  ) %>%
  select(
    rank,
    model_id,
    model_structure,
    n,
    K,
    AIC,
    AICc,
    delta_AICc,
    akaike_weight,
    cumulative_weight,
    R2,
    adjusted_R2
  )

cat("\nFULL MODEL RANKINGS:\n")

print(
  model_rankings,
  n = Inf,
  width = Inf
)

write_csv(
  model_rankings,
  file.path(
    output_dir,
    "CTmax_Model_AICc_Rankings.csv"
  )
)


# 11. COMPETITIVE MODELS -------------------------------------

competitive_models <- model_rankings %>%
  filter(
    delta_AICc <= 2
  )

cat("\nCOMPETITIVE MODELS (Delta AICc <= 2):\n")

print(
  competitive_models,
  n = Inf,
  width = Inf
)

write_csv(
  competitive_models,
  file.path(
    output_dir,
    "Competitive_Models_Delta_AICc_2.csv"
  )
)


# 12. COEFFICIENTS FOR ALL MODELS -----------------------------

coefficient_list <- list()

for (model_id in names(fitted_models)) {

  fitted_object <- fitted_models[[model_id]]

  coefficient_matrix <-
    summary(fitted_object)$coefficients

  coefficient_list[[model_id]] <- data.frame(
    model_id = model_id,

    term = rownames(coefficient_matrix),

    estimate = coefficient_matrix[, 1],

    std_error = coefficient_matrix[, 2],

    t_value = coefficient_matrix[, 3],

    p_value = coefficient_matrix[, 4],

    row.names = NULL
  )
}

all_coefficients <- dplyr::bind_rows(
  coefficient_list
)

write_csv(
  all_coefficients,
  file.path(
    output_dir,
    "All_Model_Coefficients.csv"
  )
)


# 13. BEST-SUPPORTED MODEL -----------------------------------

best_model_id <- model_rankings$model_id[1]

best_model <- fitted_models[[best_model_id]]

cat("\n========================================\n")
cat("BEST-SUPPORTED MODEL:", best_model_id, "\n")
cat("========================================\n")

print(
  summary(best_model)
)

best_matrix <- summary(best_model)$coefficients

best_ci <- stats::confint(best_model)

best_coefficients <- data.frame(
  term = rownames(best_matrix),

  estimate = best_matrix[, 1],

  std_error = best_matrix[, 2],

  t_value = best_matrix[, 3],

  p_value = best_matrix[, 4],

  CI_lower = best_ci[, 1],

  CI_upper = best_ci[, 2],

  row.names = NULL
)

write_csv(
  best_coefficients,
  file.path(
    output_dir,
    "Best_Model_Coefficients.csv"
  )
)


# 14. DESCRIPTIVE MODEL AVERAGING -----------------------------
#
# Full model averaging across the fitted candidate set.
#
# Absent coefficients contribute zero.
#
# Interpret these descriptively, especially when
# interactions are included.

all_terms <- unique(
  all_coefficients$term
)

averaging_list <- list()

for (term_name in all_terms) {

  weighted_estimate <- 0

  inclusion_weight <- 0

  for (model_id in names(fitted_models)) {

    current_weight <- model_rankings$akaike_weight[
      match(
        model_id,
        model_rankings$model_id
      )
    ]

    current_coefs <- stats::coef(
      fitted_models[[model_id]]
    )

    if (term_name %in% names(current_coefs)) {

      weighted_estimate <-
        weighted_estimate +
        current_weight *
        unname(current_coefs[term_name])

      inclusion_weight <-
        inclusion_weight +
        current_weight
    }
  }

  averaging_list[[term_name]] <- data.frame(
    term = term_name,

    full_averaged_estimate =
      weighted_estimate,

    summed_model_weight =
      inclusion_weight
  )
}

model_averaged_coefficients <- dplyr::bind_rows(
  averaging_list
)

write_csv(
  model_averaged_coefficients,
  file.path(
    output_dir,
    "Model_Averaged_Coefficients.csv"
  )
)


# 15. BEST MODEL DIAGNOSTICS ---------------------------------

diagnostic_data <- analysis_data %>%
  mutate(
    fitted_ctmax =
      as.numeric(stats::fitted(best_model)),

    residual =
      as.numeric(stats::residuals(best_model)),

    standardized_residual =
      as.numeric(stats::rstandard(best_model)),

    leverage =
      as.numeric(stats::hatvalues(best_model)),

    cooks_distance =
      as.numeric(stats::cooks.distance(best_model))
  )

write_csv(
  diagnostic_data,
  file.path(
    output_dir,
    "Best_Model_Diagnostics.csv"
  )
)

influential_observations <- diagnostic_data %>%
  filter(
    cooks_distance >
      4 / nrow(diagnostic_data)
  )

write_csv(
  influential_observations,
  file.path(
    output_dir,
    "Potentially_Influential_Observations.csv"
  )
)


# 16. FIGURE: CTmax VS BODY MASS -----------------------------

theme_set(
  theme_classic(base_size = 13)
)

plot_mass <- ggplot(
  analysis_data,
  aes(
    x = weight_g,
    y = ctmax,
    color = habitat
  )
) +
  geom_point(
    alpha = 0.75,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE
  ) +
  facet_wrap(
    ~ species_group,
    scales = "free_x"
  ) +
  labs(
    x = "Body mass (g)",
    y = "CTmax (°C)",
    color = "Habitat"
  )

ggsave(
  file.path(
    output_dir,
    "CTmax_vs_Body_Mass.png"
  ),
  plot_mass,
  width = 9,
  height = 5,
  dpi = 300
)


# 17. FIGURE: CTmax VS PARASITE BURDEN ------------------------

plot_parasites <- ggplot(
  analysis_data,
  aes(
    x = parasite_burden_percent,
    y = ctmax,
    color = habitat
  )
) +
  geom_point(
    alpha = 0.75,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE
  ) +
  facet_wrap(
    ~ species_group,
    scales = "free_x"
  ) +
  labs(
    x = "Parasite burden (% of body mass)",
    y = "CTmax (°C)",
    color = "Habitat"
  )

ggsave(
  file.path(
    output_dir,
    "CTmax_vs_Parasite_Burden.png"
  ),
  plot_parasites,
  width = 9,
  height = 5,
  dpi = 300
)


# 18. FIGURE: CTmax BY SPECIES AND HABITAT --------------------

plot_habitat <- ggplot(
  analysis_data,
  aes(
    x = habitat,
    y = ctmax,
    fill = species_group
  )
) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.5,
    size = 1.5
  ) +
  labs(
    x = "Habitat",
    y = "CTmax (°C)",
    fill = "Species"
  )

ggsave(
  file.path(
    output_dir,
    "CTmax_by_Species_and_Habitat.png"
  ),
  plot_habitat,
  width = 7,
  height = 5,
  dpi = 300
)


# 19. FIGURE: RESIDUALS VS FITTED -----------------------------

plot_residuals <- ggplot(
  diagnostic_data,
  aes(
    x = fitted_ctmax,
    y = residual
  )
) +
  geom_point(
    alpha = 0.7
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "Fitted CTmax (°C)",
    y = "Residual (°C)"
  )

ggsave(
  file.path(
    output_dir,
    "Best_Model_Residuals.png"
  ),
  plot_residuals,
  width = 7,
  height = 5,
  dpi = 300
)


# 20. FIGURE: NORMAL Q-Q PLOT --------------------------------

plot_qq <- ggplot(
  diagnostic_data,
  aes(
    sample = residual
  )
) +
  stat_qq() +
  stat_qq_line() +
  labs(
    x = "Theoretical quantiles",
    y = "Residual quantiles"
  )

ggsave(
  file.path(
    output_dir,
    "Best_Model_QQ_Plot.png"
  ),
  plot_qq,
  width = 7,
  height = 5,
  dpi = 300
)


# 21. JOURNAL-STYLE MODEL SELECTION TABLE --------------------

journal_table <- model_rankings %>%
  select(
    model_id,
    model_structure,
    K,
    AICc,
    delta_AICc,
    akaike_weight
  ) %>%
  gt::gt() %>%
  gt::tab_header(
    title = gt::md(
      "**Table 1. Candidate models explaining variation in CTmax**"
    ),
    subtitle = paste0(
      "AICc model comparison; n = ",
      nrow(analysis_data)
    )
  ) %>%
  gt::cols_label(
    model_id = "Model",
    model_structure = "Model structure",
    K = "K",
    AICc = "AICc",
    delta_AICc = gt::html("&Delta;AICc"),
    akaike_weight = "Akaike weight"
  ) %>%
  gt::fmt_number(
    columns = c(
      AICc,
      delta_AICc
    ),
    decimals = 2
  ) %>%
  gt::fmt_number(
    columns = akaike_weight,
    decimals = 3
  ) %>%
  gt::tab_style(
    style = gt::cell_text(
      weight = "bold"
    ),
    locations = gt::cells_body(
      rows = 1
    )
  ) %>%
  gt::tab_options(
    table.font.names = "Arial",
    table.font.size = gt::px(12),
    heading.title.font.size = gt::px(15),
    heading.subtitle.font.size = gt::px(11),
    column_labels.font.weight = "bold",
    table.border.top.width = gt::px(2),
    table.border.bottom.width = gt::px(2),
    table_body.hlines.style = "none",
    data_row.padding = gt::px(7)
  ) %>%
  gt::tab_source_note(
    gt::md(
      paste0(
        "All continuous predictors are expressed in ",
        "original units: body mass (g), parasite ",
        "burden (% of body mass), and starting ",
        "temperature (°C). K includes residual ",
        "variance. The lowest-AICc model is bold."
      )
    )
  )

# Display in RStudio Viewer
print(journal_table)

# Save as HTML
gt::gtsave(
  journal_table,
  file.path(
    output_dir,
    "CTmax_Journal_Model_Table.html"
  )
)


# 22. SAVE FITTED MODEL OBJECTS -------------------------------

saveRDS(
  fitted_models,
  file.path(
    output_dir,
    "All_Fitted_Models.rds"
  )
)

saveRDS(
  best_model,
  file.path(
    output_dir,
    "Best_Model.rds"
  )
)


# 23. FINAL SUMMARY ------------------------------------------

cat("\n========================================\n")
cat("CTmax MODEL ANALYSIS COMPLETE\n")
cat("========================================\n")

cat(
  "\nObservations:",
  nrow(analysis_data)
)

cat(
  "\nCandidate models:",
  length(candidate_formulas)
)

cat(
  "\nSuccessfully fitted:",
  length(fitted_models)
)

cat(
  "\nSkipped models:",
  nrow(skipped_models)
)

cat(
  "\nBest-supported model:",
  best_model_id
)

cat(
  "\nBest AICc:",
  round(
    model_rankings$AICc[1],
    2
  )
)

cat(
  "\nBest Akaike weight:",
  round(
    model_rankings$akaike_weight[1],
    3
  )
)

cat(
  "\nCompetitive models (Delta AICc <= 2):",
  nrow(competitive_models)
)

cat(
  "\n\nAll results saved to:",
  output_dir,
  "\n"
)

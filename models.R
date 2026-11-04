
# ============================================================
# STEP 2: MONGOLIA CTmax MODEL COMPARISON
# 20 CANDIDATE MODELS
# NO MASS + FULTON'S K IN THE SAME MODEL
# BASE R ONLY
# ============================================================


# 1. SETTINGS -----------------------------------------------

options(
  na.print = "NA",
  max.print = 10000,
  width = 120
)

input_file <- "CTmax_Merged_Data.csv"
output_dir <- "CTmax_Model_Results_No_Collinearity"

excluded_trials <- c(1, 26, 67, 84, 108)

vif_threshold <- 5

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

output_path <- function(filename) {
  file.path(output_dir, filename)
}

save_csv <- function(x, filename) {
  write.csv(
    x,
    output_path(filename),
    row.names = FALSE,
    na = ""
  )
}

as_number <- function(x) {
  suppressWarnings(
    as.numeric(as.character(x))
  )
}


# 2. IMPORT DATA --------------------------------------------

if (!file.exists(input_file)) {
  stop(
    "Cannot find CTmax_Merged_Data.csv. ",
    "Run Step 1 first."
  )
}

raw_data <- read.csv(
  input_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
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
    "Missing columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

length_candidates <- c(
  "tl_mm",
  "total_length_mm",
  "total_length",
  "length_mm"
)

length_column <- length_candidates[
  length_candidates %in% names(raw_data)
]

if (length(length_column) == 0) {
  stop(
    "No total length column found. ",
    "Total length in mm is required for Fulton's K."
  )
}

length_column <- length_column[1]

if (!"parasite_burden_percent" %in% names(raw_data)) {

  if (!"parasite_burden" %in% names(raw_data)) {
    stop("No parasite burden column found.")
  }

  raw_data$parasite_burden_percent <-
    as_number(raw_data$parasite_burden) * 100
}


# 3. CLEAN VARIABLES ----------------------------------------

raw_data$trial <- as_number(raw_data$trial)
raw_data$ctmax <- as_number(raw_data$ctmax)
raw_data$weight_g <- as_number(raw_data$weight_g)
raw_data$start_temp <- as_number(raw_data$start_temp)

raw_data$tl_mm_analysis <- as_number(
  raw_data[[length_column]]
)

raw_data$parasite_burden_percent <- as_number(
  raw_data$parasite_burden_percent
)

raw_data$fish_id <- trimws(
  as.character(raw_data$fish_id)
)

raw_data$habitat <- tolower(
  trimws(as.character(raw_data$site_type))
)

species_text <- tolower(
  trimws(as.character(raw_data$species))
)

raw_data$species_group <- NA_character_

raw_data$species_group[
  grepl("grayling", species_text)
] <- "Grayling"

raw_data$species_group[
  grepl("lenok", species_text)
] <- "Lenok"


# 4. CALCULATE FULTON'S K -----------------------------------

# K = 100 * weight(g) / length(cm)^3

raw_data$fultons_k <- (
  100 * raw_data$weight_g /
    (raw_data$tl_mm_analysis / 10)^3
)


# 5. CREATE COMMON ANALYSIS DATASET -------------------------

keep_rows <- (
  !(raw_data$trial %in% excluded_trials) &
    is.finite(raw_data$ctmax) &
    is.finite(raw_data$weight_g) &
    raw_data$weight_g > 0 &
    is.finite(raw_data$tl_mm_analysis) &
    raw_data$tl_mm_analysis > 0 &
    is.finite(raw_data$fultons_k) &
    is.finite(raw_data$start_temp) &
    is.finite(raw_data$parasite_burden_percent) &
    raw_data$parasite_burden_percent >= 0 &
    !is.na(raw_data$species_group) &
    raw_data$habitat %in% c("lake", "river") &
    !is.na(raw_data$fish_id) &
    nzchar(raw_data$fish_id)
)

keep_rows[is.na(keep_rows)] <- FALSE

analysis_data <- raw_data[
  keep_rows,
  ,
  drop = FALSE
]

analysis_data$species_group <- factor(
  analysis_data$species_group,
  levels = c("Grayling", "Lenok")
)

analysis_data$habitat <- factor(
  analysis_data$habitat,
  levels = c("lake", "river")
)

rownames(analysis_data) <- NULL

sample_n <- nrow(analysis_data)

if (sample_n < 30) {
  stop(
    "Fewer than 30 complete observations remain."
  )
}

group_counts <- as.data.frame(
  table(
    species_group = analysis_data$species_group,
    habitat = analysis_data$habitat
  )
)

names(group_counts)[3] <- "n"

if (any(group_counts$n < 5)) {
  stop(
    "At least one species-habitat group ",
    "contains fewer than five fish."
  )
}

save_csv(
  analysis_data,
  "CTmax_Analysis_Dataset.csv"
)

save_csv(
  group_counts,
  "Species_Habitat_Sample_Sizes.csv"
)

cat(
  "\nAnalysis sample size:",
  sample_n,
  "\n"
)


# 6. DESCRIPTIVE SUMMARIES ----------------------------------

summary_rows <- list()
counter <- 0

for (sp in levels(analysis_data$species_group)) {

  for (hab in levels(analysis_data$habitat)) {

    counter <- counter + 1

    d <- analysis_data[
      analysis_data$species_group == sp &
        analysis_data$habitat == hab,
      ,
      drop = FALSE
    ]

    summary_rows[[counter]] <- data.frame(
      species = sp,
      habitat = hab,
      n = nrow(d),
      mean_ctmax = mean(d$ctmax),
      sd_ctmax = sd(d$ctmax),
      mean_weight_g = mean(d$weight_g),
      mean_length_mm = mean(d$tl_mm_analysis),
      mean_fultons_k = mean(d$fultons_k),
      mean_parasite_percent =
        mean(d$parasite_burden_percent),
      mean_start_temp = mean(d$start_temp)
    )
  }
}

group_summary <- do.call(
  rbind,
  summary_rows
)

save_csv(
  group_summary,
  "Species_Habitat_Descriptive_Summary.csv"
)


# 7. EXAMINE PREDICTOR CORRELATIONS --------------------------

continuous_predictors <- analysis_data[
  ,
  c(
    "weight_g",
    "tl_mm_analysis",
    "fultons_k",
    "parasite_burden_percent",
    "start_temp"
  )
]

correlation_matrix <- cor(
  continuous_predictors,
  use = "complete.obs"
)

write.csv(
  correlation_matrix,
  output_path("Predictor_Correlation_Matrix.csv")
)

cat("\nPredictor correlations:\n")
print(round(correlation_matrix, 3))


# 8. DEFINE 20 CANDIDATE MODELS -----------------------------
#
# MASS MODELS: M01-M10
# CONDITION MODELS: M11-M20
#
# No model includes weight_g and fultons_k together.
#
# All models use the same complete-case dataset.
#
# Species-habitat interaction retained in all models.

candidate_formulas <- list(

  # BASELINE AND BODY MASS MODELS

  M01 = ctmax ~ species_group * habitat,

  M02 = ctmax ~ species_group * habitat +
    weight_g,

  M03 = ctmax ~ species_group * habitat +
    weight_g + habitat:weight_g,

  M04 = ctmax ~ species_group * habitat +
    weight_g + species_group:weight_g,

  M05 = ctmax ~ species_group * habitat *
    weight_g,

  M06 = ctmax ~ species_group * habitat +
    parasite_burden_percent,

  M07 = ctmax ~ species_group * habitat +
    weight_g + parasite_burden_percent,

  M08 = ctmax ~ species_group * habitat +
    weight_g + habitat:weight_g +
    parasite_burden_percent,

  M09 = ctmax ~ species_group * habitat +
    weight_g + parasite_burden_percent +
    start_temp,

  M10 = ctmax ~ species_group * habitat +
    weight_g + habitat:weight_g +
    parasite_burden_percent + start_temp,


  # FULTON'S CONDITION FACTOR MODELS

  M11 = ctmax ~ species_group * habitat +
    fultons_k,

  M12 = ctmax ~ species_group * habitat +
    fultons_k + habitat:fultons_k,

  M13 = ctmax ~ species_group * habitat +
    fultons_k + species_group:fultons_k,

  M14 = ctmax ~ species_group * habitat *
    fultons_k,

  M15 = ctmax ~ species_group * habitat +
    fultons_k + parasite_burden_percent,

  M16 = ctmax ~ species_group * habitat +
    fultons_k + habitat:fultons_k +
    parasite_burden_percent,

  M17 = ctmax ~ species_group * habitat +
    fultons_k + species_group:fultons_k +
    parasite_burden_percent,

  M18 = ctmax ~ species_group * habitat +
    fultons_k + parasite_burden_percent +
    start_temp,

  M19 = ctmax ~ species_group * habitat +
    fultons_k + habitat:fultons_k +
    parasite_burden_percent + start_temp,

  M20 = ctmax ~ species_group * habitat +
    fultons_k + species_group:fultons_k +
    parasite_burden_percent + start_temp
)

stopifnot(
  length(candidate_formulas) == 20
)


# 9. CHECK FOR FORBIDDEN PREDICTOR COMBINATIONS -------------

for (model_id in names(candidate_formulas)) {

  variables <- all.vars(
    candidate_formulas[[model_id]]
  )

  if (
    "weight_g" %in% variables &&
    "fultons_k" %in% variables
  ) {
    stop(
      "Forbidden mass-K combination in ",
      model_id
    )
  }

  if (
    "tl_mm_analysis" %in% variables &&
    "fultons_k" %in% variables
  ) {
    stop(
      "Forbidden length-K combination in ",
      model_id
    )
  }
}

cat(
  "\nAll 20 models passed the ",
  "mass-condition separation check.\n",
  sep = ""
)


# 10. VIF FUNCTION ------------------------------------------
#
# Calculate VIFs from the model matrix.
#
# For interactions, ordinary term-level VIFs can be
# difficult to interpret. This function calculates
# VIFs for the individual non-intercept columns.
#
# A high VIF may reflect interaction coding as well
# as biological predictor overlap.
#
# Models exceeding the threshold are recorded and
# excluded from the primary ranking.

calculate_vif <- function(fitted_model) {

  X <- model.matrix(fitted_model)

  if ("(Intercept)" %in% colnames(X)) {
    X <- X[
      ,
      colnames(X) != "(Intercept)",
      drop = FALSE
    ]
  }

  if (ncol(X) == 0) {
    return(numeric(0))
  }

  if (ncol(X) == 1) {
    result <- 1
    names(result) <- colnames(X)
    return(result)
  }

  vif_values <- numeric(ncol(X))
  names(vif_values) <- colnames(X)

  for (j in seq_len(ncol(X))) {

    y <- X[, j]
    others <- X[, -j, drop = FALSE]

    if (sd(y) == 0) {
      vif_values[j] <- Inf
      next
    }

    auxiliary_fit <- tryCatch(
      lm.fit(
        x = cbind(
          "(Intercept)" = 1,
          others
        ),
        y = y
      ),
      error = function(e) NULL
    )

    if (is.null(auxiliary_fit)) {
      vif_values[j] <- Inf
      next
    }

    residual_ss <- sum(
      auxiliary_fit$residuals^2
    )

    total_ss <- sum(
      (y - mean(y))^2
    )

    if (total_ss <= 0) {
      vif_values[j] <- Inf
      next
    }

    r_squared <- 1 - residual_ss / total_ss

    if (r_squared >= 1) {
      vif_values[j] <- Inf
    } else {
      vif_values[j] <- 1 / (1 - r_squared)
    }
  }

  vif_values
}


# 11. AICc FUNCTION -----------------------------------------

calculate_aicc <- function(fitted_model) {

  n <- nobs(fitted_model)

  k <- attr(
    logLik(fitted_model),
    "df"
  )

  if (n <= k + 1) {
    return(Inf)
  }

  AIC(fitted_model) +
    2 * k * (k + 1) / (n - k - 1)
}


# 12. FIT AND SCREEN CANDIDATE MODELS -----------------------

fitted_models <- list()
ranking_rows <- list()
skipped_rows <- list()
vif_rows <- list()

for (model_id in names(candidate_formulas)) {

  current_formula <- candidate_formulas[[model_id]]

  fitted_model <- tryCatch(
    lm(
      current_formula,
      data = analysis_data,
      na.action = na.fail
    ),
    error = function(e) e
  )

  if (inherits(fitted_model, "error")) {

    skipped_rows[[model_id]] <- data.frame(
      model_id = model_id,
      reason = conditionMessage(fitted_model)
    )

    next
  }

  # Reject rank-deficient models.

  if (
    fitted_model$rank <
    length(coef(fitted_model))
  ) {

    skipped_rows[[model_id]] <- data.frame(
      model_id = model_id,
      reason = "Rank-deficient model"
    )

    next
  }

  # Calculate VIFs.

  current_vif <- calculate_vif(
    fitted_model
  )

  for (term_name in names(current_vif)) {

    vif_rows[[length(vif_rows) + 1]] <- data.frame(
      model_id = model_id,
      term = term_name,
      VIF = unname(current_vif[term_name])
    )
  }

  max_vif <- if (length(current_vif) > 0) {
    max(current_vif)
  } else {
    1
  }

  # Interaction-column VIF screening.
  #
  # High VIFs involving interaction terms may be
  # caused by coding rather than biological
  # redundancy. Such models are excluded from the
  # primary table but preserved in the screening
  # report for manual review.

  if (
    !is.finite(max_vif) ||
    max_vif > vif_threshold
  ) {

    skipped_rows[[model_id]] <- data.frame(
      model_id = model_id,
      reason = paste0(
        "VIF screening: maximum VIF = ",
        round(max_vif, 3)
      )
    )

    next
  }

  fitted_models[[model_id]] <- fitted_model

  model_summary <- summary(fitted_model)

  formula_text <- paste(
    deparse(
      current_formula,
      width.cutoff = 500
    ),
    collapse = " "
  )

  ranking_rows[[model_id]] <- data.frame(
    model_id = model_id,
    model_structure = formula_text,
    n = nobs(fitted_model),
    K = attr(logLik(fitted_model), "df"),
    AIC = AIC(fitted_model),
    AICc = calculate_aicc(fitted_model),
    max_VIF = max_vif,
    R2 = model_summary$r.squared,
    adjusted_R2 = model_summary$adj.r.squared,
    stringsAsFactors = FALSE
  )
}

skipped_models <- if (length(skipped_rows) > 0) {
  do.call(rbind, skipped_rows)
} else {
  data.frame(
    model_id = character(),
    reason = character()
  )
}

rownames(skipped_models) <- NULL

save_csv(
  skipped_models,
  "Skipped_Models.csv"
)

vif_results <- if (length(vif_rows) > 0) {
  do.call(rbind, vif_rows)
} else {
  data.frame(
    model_id = character(),
    term = character(),
    VIF = numeric()
  )
}

rownames(vif_results) <- NULL

save_csv(
  vif_results,
  "All_Model_VIF_Results.csv"
)

if (length(fitted_models) == 0) {
  stop(
    "All models failed screening. ",
    "Review All_Model_VIF_Results.csv. ",
    "Interaction coding may be inflating VIFs."
  )
}


# 13. AICc MODEL COMPARISON ---------------------------------

model_rankings <- do.call(
  rbind,
  ranking_rows
)

model_rankings <- model_rankings[
  is.finite(model_rankings$AICc),
  ,
  drop = FALSE
]

if (nrow(model_rankings) == 0) {
  stop("No models have finite AICc.")
}

model_rankings <- model_rankings[
  order(model_rankings$AICc),
  ,
  drop = FALSE
]

rownames(model_rankings) <- NULL

model_rankings$delta_AICc <-
  model_rankings$AICc -
  min(model_rankings$AICc)

relative_likelihood <- exp(
  -0.5 * model_rankings$delta_AICc
)

model_rankings$akaike_weight <-
  relative_likelihood /
  sum(relative_likelihood)

model_rankings$cumulative_weight <- cumsum(
  model_rankings$akaike_weight
)

model_rankings$rank <- seq_len(
  nrow(model_rankings)
)

model_rankings <- model_rankings[
  ,
  c(
    "rank",
    "model_id",
    "model_structure",
    "n",
    "K",
    "AIC",
    "AICc",
    "delta_AICc",
    "akaike_weight",
    "cumulative_weight",
    "max_VIF",
    "R2",
    "adjusted_R2"
  )
]

save_csv(
  model_rankings,
  "CTmax_Model_AICc_Rankings.csv"
)

competitive_models <- model_rankings[
  model_rankings$delta_AICc <= 2,
  ,
  drop = FALSE
]

save_csv(
  competitive_models,
  "Competitive_Models_Delta_AICc_2.csv"
)


# 14. COEFFICIENTS FOR ALL MODELS ---------------------------

coefficient_rows <- list()

for (model_id in names(fitted_models)) {

  fitted_model <- fitted_models[[model_id]]

  coefficient_matrix <- summary(
    fitted_model
  )$coefficients

  confidence_intervals <- confint(
    fitted_model,
    level = 0.95
  )

  coefficient_rows[[model_id]] <- data.frame(
    model_id = model_id,
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[, 1],
    std_error = coefficient_matrix[, 2],
    t_value = coefficient_matrix[, 3],
    p_value = coefficient_matrix[, 4],
    CI_lower = confidence_intervals[, 1],
    CI_upper = confidence_intervals[, 2],
    row.names = NULL
  )
}

all_coefficients <- do.call(
  rbind,
  coefficient_rows
)

rownames(all_coefficients) <- NULL

save_csv(
  all_coefficients,
  "All_Model_Coefficients.csv"
)


# 15. BEST-SUPPORTED MODEL ---------------------------------

best_model_id <- model_rankings$model_id[1]

best_model <- fitted_models[[best_model_id]]

best_coefficients <- all_coefficients[
  all_coefficients$model_id == best_model_id,
  ,
  drop = FALSE
]

save_csv(
  best_coefficients,
  "Best_Model_Coefficients.csv"
)

summary_text <- capture.output(
  print(summary(best_model))
)

writeLines(
  summary_text,
  output_path("Best_Model_Summary.txt")
)


# 16. MODEL AVERAGING ---------------------------------------
#
# Full model averaging across retained models.
#
# Absent coefficients contribute zero.
# These are descriptive estimates, not averaged
# standard errors or confidence intervals.

all_terms <- unique(
  all_coefficients$term
)

averaging_rows <- list()

for (term_name in all_terms) {

  weighted_estimate <- 0
  inclusion_weight <- 0

  for (model_id in model_rankings$model_id) {

    current_weight <- model_rankings$akaike_weight[
      match(
        model_id,
        model_rankings$model_id
      )
    ]

    current_coefficients <- coef(
      fitted_models[[model_id]]
    )

    if (
      term_name %in%
      names(current_coefficients)
    ) {

      weighted_estimate <- weighted_estimate +
        current_weight *
        unname(
          current_coefficients[term_name]
        )

      inclusion_weight <- inclusion_weight +
        current_weight
    }
  }

  averaging_rows[[term_name]] <- data.frame(
    term = term_name,
    full_averaged_estimate = weighted_estimate,
    summed_model_weight = inclusion_weight
  )
}

model_averaged_coefficients <- do.call(
  rbind,
  averaging_rows
)

rownames(model_averaged_coefficients) <- NULL

save_csv(
  model_averaged_coefficients,
  "Model_Averaged_Coefficients.csv"
)


# 17. MODEL DIAGNOSTICS -------------------------------------

diagnostic_data <- analysis_data

diagnostic_data$fitted_ctmax <- as.numeric(
  fitted(best_model)
)

diagnostic_data$residual <- as.numeric(
  residuals(best_model)
)

diagnostic_data$standardized_residual <- as.numeric(
  rstandard(best_model)
)

diagnostic_data$leverage <- as.numeric(
  hatvalues(best_model)
)

diagnostic_data$cooks_distance <- as.numeric(
  cooks.distance(best_model)
)

save_csv(
  diagnostic_data,
  "Best_Model_Diagnostics.csv"
)

influential_rows <- (
  diagnostic_data$cooks_distance >
    4 / nrow(diagnostic_data)
)

influential_rows[
  is.na(influential_rows)
] <- FALSE

save_csv(
  diagnostic_data[
    influential_rows,
    ,
    drop = FALSE
  ],
  "Potentially_Influential_Observations.csv"
)


# 18. PLOTTING SETTINGS ------------------------------------

habitat_colors <- c(
  lake = "#2166AC",
  river = "#B2182B"
)

plot_predictor <- function(
    predictor,
    x_label,
    file_name
) {

  png(
    output_path(file_name),
    width = 2400,
    height = 1200,
    res = 300
  )

  old_par <- par(no.readonly = TRUE)

  par(
    mfrow = c(1, 2),
    mar = c(5, 5, 3, 1)
  )

  for (sp in c("Grayling", "Lenok")) {

    d <- analysis_data[
      analysis_data$species_group == sp,
      ,
      drop = FALSE
    ]

    plot(
      d[[predictor]],
      d$ctmax,
      pch = 16,
      col = habitat_colors[
        as.character(d$habitat)
      ],
      xlab = x_label,
      ylab = "CTmax (degrees C)",
      main = sp
    )

    for (hab in c("lake", "river")) {

      group_data <- d[
        d$habitat == hab,
        ,
        drop = FALSE
      ]

      x <- group_data[[predictor]]

      if (
        nrow(group_data) >= 3 &&
        length(unique(x)) >= 2
      ) {

        fit_line <- lm(
          group_data$ctmax ~ x
        )

        abline(
          fit_line,
          col = habitat_colors[hab],
          lwd = 2
        )
      }
    }

    legend(
      "topright",
      legend = c("Lake", "River"),
      col = habitat_colors,
      pch = 16,
      bty = "n",
      cex = 0.8
    )
  }

  par(old_par)
  dev.off()
}


# 19. CTmax VS BODY MASS -----------------------------------

plot_predictor(
  "weight_g",
  "Body mass (g)",
  "CTmax_vs_Body_Mass.png"
)


# 20. CTmax VS PARASITE BURDEN ------------------------------

plot_predictor(
  "parasite_burden_percent",
  "Parasite burden (% body mass)",
  "CTmax_vs_Parasite_Burden.png"
)


# 21. CTmax VS FULTON'S K ----------------------------------

plot_predictor(
  "fultons_k",
  "Fulton's condition factor (K)",
  "CTmax_vs_Fultons_K.png"
)


# 22. CTmax VS STARTING TEMPERATURE ------------------------

plot_predictor(
  "start_temp",
  "Starting temperature (degrees C)",
  "CTmax_vs_Starting_Temperature.png"
)


# 23. FULTON'S K VS LENGTH ---------------------------------

png(
  output_path("Fultons_K_vs_Total_Length.png"),
  width = 2400,
  height = 1200,
  res = 300
)

par(
  mfrow = c(1, 2),
  mar = c(5, 5, 3, 1)
)

for (sp in c("Grayling", "Lenok")) {

  d <- analysis_data[
    analysis_data$species_group == sp,
    ,
    drop = FALSE
  ]

  plot(
    d$tl_mm_analysis,
    d$fultons_k,
    pch = 16,
    col = habitat_colors[
      as.character(d$habitat)
    ],
    xlab = "Total length (mm)",
    ylab = "Fulton's condition factor (K)",
    main = sp
  )

  legend(
    "topright",
    legend = c("Lake", "River"),
    col = habitat_colors,
    pch = 16,
    bty = "n",
    cex = 0.8
  )
}

dev.off()


# 24. CTmax BY SPECIES AND HABITAT -------------------------

png(
  output_path("CTmax_by_Species_and_Habitat.png"),
  width = 2100,
  height = 1500,
  res = 300
)

group_factor <- interaction(
  analysis_data$species_group,
  analysis_data$habitat,
  sep = " - "
)

boxplot(
  analysis_data$ctmax ~ group_factor,
  col = c(
    "#92B6D5",
    "#92B6D5",
    "#D9A4A4",
    "#D9A4A4"
  ),
  ylab = "CTmax (degrees C)",
  xlab = "Species and habitat",
  las = 2,
  main = "CTmax by species and habitat"
)

stripchart(
  analysis_data$ctmax ~ group_factor,
  vertical = TRUE,
  method = "jitter",
  jitter = 0.12,
  pch = 16,
  col = adjustcolor(
    "black",
    alpha.f = 0.45
  ),
  add = TRUE
)

dev.off()


# 25. RESIDUAL DIAGNOSTICS ---------------------------------

png(
  output_path("Best_Model_Residuals.png"),
  width = 1800,
  height = 1400,
  res = 300
)

plot(
  diagnostic_data$fitted_ctmax,
  diagnostic_data$residual,
  pch = 16,
  col = "gray35",
  xlab = "Fitted CTmax (degrees C)",
  ylab = "Residual (degrees C)",
  main = paste(
    "Residual diagnostics:",
    best_model_id
  )
)

abline(
  h = 0,
  lty = 2,
  lwd = 2
)

dev.off()

png(
  output_path("Best_Model_QQ_Plot.png"),
  width = 1800,
  height = 1400,
  res = 300
)

qqnorm(
  diagnostic_data$residual,
  pch = 16,
  col = "gray35",
  main = paste(
    "Normal Q-Q:",
    best_model_id
  )
)

qqline(
  diagnostic_data$residual,
  lwd = 2
)

dev.off()


# 26. FORMATTED HTML MODEL TABLE ---------------------------

html_escape <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "&", "&amp;", x,
    fixed = TRUE
  )

  x <- gsub(
    "<", "&lt;", x,
    fixed = TRUE
  )

  x <- gsub(
    ">", "&gt;", x,
    fixed = TRUE
  )

  x <- gsub(
    '"', "&quot;", x,
    fixed = TRUE
  )

  x
}

html_rows <- character()

for (i in seq_len(nrow(model_rankings))) {

  current_row <- model_rankings[
    i,
    ,
    drop = FALSE
  ]

  row_style <- if (i == 1) {
    ' style="font-weight:bold;background:#f0f0f0;"'
  } else {
    ""
  }

  html_rows <- c(
    html_rows,
    paste0(
      "<tr", row_style, ">",
      "<td>",
      html_escape(current_row$model_id),
      "</td>",
      "<td>",
      html_escape(current_row$model_structure),
      "</td>",
      "<td>",
      current_row$K,
      "</td>",
      "<td>",
      sprintf("%.2f", current_row$AICc),
      "</td>",
      "<td>",
      sprintf("%.2f", current_row$delta_AICc),
      "</td>",
      "<td>",
      sprintf("%.3f", current_row$akaike_weight),
      "</td>",
      "<td>",
      sprintf("%.2f", current_row$max_VIF),
      "</td>",
      "</tr>"
    )
  )
}

html_document <- c(
  "<!DOCTYPE html>",
  "<html>",
  "<head>",
  '<meta charset="UTF-8">',
  "<title>CTmax Model Comparison</title>",
  "<style>",
  "body {font-family: Arial, sans-serif; margin: 35px;}",
  "h2 {font-size: 19px;}",
  "p {font-size: 12px; color: #444;}",
  "table {border-collapse: collapse; width: 100%;}",
  "th, td {padding: 9px; border-bottom: 1px solid #ddd;}",
  "th {text-align: left; border-top: 2px solid black;",
  "border-bottom: 2px solid black;}",
  "td {font-size: 12px;}",
  "tr:last-child td {border-bottom: 2px solid black;}",
  "</style>",
  "</head>",
  "<body>",
  "<h2>Table 1. Candidate models explaining CTmax</h2>",
  paste0(
    "<p>AICc model comparison; n = ",
    sample_n,
    " fish. Models exceeding the VIF threshold",
    " were excluded from this ranking.</p>"
  ),
  "<table>",
  "<thead><tr>",
  "<th>Model</th>",
  "<th>Model structure</th>",
  "<th>K</th>",
  "<th>AICc</th>",
  "<th>&Delta;AICc</th>",
  "<th>Akaike weight</th>",
  "<th>Max VIF</th>",
  "</tr></thead>",
  "<tbody>",
  html_rows,
  "</tbody></table>",
  "<p>Continuous predictors use original biological units.",
  "K is the number of estimated parameters, including",
  "residual variance. Fulton's K and body mass were",
  "not included together in any candidate model.",
  "The best-supported retained model is bold.</p>",
  "</body>",
  "</html>"
)

writeLines(
  html_document,
  output_path("CTmax_Journal_Model_Table.html"),
  useBytes = TRUE
)


# 27. SAVE FITTED MODELS -----------------------------------

saveRDS(
  fitted_models,
  output_path("All_Fitted_Models.rds")
)

saveRDS(
  best_model,
  output_path("Best_Model.rds")
)


# 28. FINAL SUMMARY ----------------------------------------

summary_lines <- c(
  "MONGOLIA CTmax MODEL ANALYSIS",
  "============================",
  paste("Observations:", sample_n),
  paste(
    "Candidate models:",
    length(candidate_formulas)
  ),
  paste(
    "Retained models:",
    length(fitted_models)
  ),
  paste(
    "Skipped models:",
    nrow(skipped_models)
  ),
  paste(
    "VIF threshold:",
    vif_threshold
  ),
  paste(
    "Best model:",
    best_model_id
  ),
  paste(
    "Best AICc:",
    round(model_rankings$AICc[1], 3)
  ),
  paste(
    "Best Akaike weight:",
    round(model_rankings$akaike_weight[1], 4)
  ),
  paste(
    "Competitive models (delta AICc <= 2):",
    nrow(competitive_models)
  )
)

writeLines(
  summary_lines,
  output_path("Analysis_Summary.txt")
)

cat(
  "\nANALYSIS COMPLETE\n",
  "Best model: ", best_model_id, "\n",
  "Observations: ", sample_n, "\n",
  "Retained models: ", length(fitted_models), "\n",
  "Skipped models: ", nrow(skipped_models), "\n",
  "Results folder: ", output_dir, "\n",
  sep = ""
)


# 29. DISPLAY TABLE IN RSTUDIO VIEWER ----------------------

table_file <- normalizePath(
  output_path("CTmax_Journal_Model_Table.html"),
  winslash = "/",
  mustWork = TRUE
)

if (interactive()) {

  viewer <- getOption("viewer")

  if (is.function(viewer)) {
    viewer(table_file)
  } else {
    browseURL(table_file)
  }

}

# ============================================================
# END OF STEP 2
# ============================================================

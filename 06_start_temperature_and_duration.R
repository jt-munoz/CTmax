# 15. TRIAL DURATION AND STARTING TEMPERATURE VS. CTMAX ----
# Trial duration is largely determined by starting temperature,
# ramp rate, and the CTmax endpoint; regressions are descriptive only.

# PREPARE DATA ----
regression_data <- data_final %>%
  mutate(
    species = factor(
      str_to_title(as.character(species)),
      levels = c("Grayling", "Lenok")
    ),
    site_type = factor(
      tolower(as.character(site_type)),
      levels = c("lake", "river")
    ),
    group4 = factor(
      paste(str_to_title(as.character(site_type)), species),
      levels = c(
        "Lake Grayling", "River Grayling",
        "Lake Lenok", "River Lenok"
      )
    ),
    trial_length = suppressWarnings(as.numeric(as.character(trial_length))),
    start_temp = suppressWarnings(as.numeric(as.character(start_temp)))
  )

group_colors <- c(
  "Lake Grayling" = "#3278A3",
  "River Grayling" = "#C45A4A",
  "Lake Lenok" = "#3278A3",
  "River Lenok" = "#C45A4A"
)

# REGRESSION FUNCTION ----
run_ctmax_regression <- function(data, predictor) {

  d <- data %>%
    filter(
      is.finite(ctmax),
      is.finite(.data[[predictor]]),
      !is.na(group4)
    )

  results <- d %>%
    group_by(group4, .drop = TRUE) %>%
    group_modify(~ {
      model <- lm(reformulate(predictor, response = "ctmax"), data = .x)
      coefs <- summary(model)$coefficients
      ci <- confint(model)

      tibble(
        n = nrow(.x),
        intercept = unname(coef(model)[1]),
        slope = unname(coef(model)[2]),
        slope_se = coefs[2, "Std. Error"],
        ci_lower = ci[2, 1],
        ci_upper = ci[2, 2],
        r_squared = summary(model)$r.squared,
        p = coefs[2, "Pr(>|t|)"]
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_holm = p.adjust(p, method = "holm"),
      label = sprintf(
        "n = %d | Slope = %+.3f\nR² = %.3f | p = %s",
        n, slope, r_squared,
        format.pval(p, digits = 2, eps = 0.001)
      )
    )

  list(data = d, results = results)
}

# TRIAL DURATION REGRESSIONS ----
trial_results <- run_ctmax_regression(
  regression_data, "trial_length"
)

cat("\n=== TRIAL DURATION VS. CTMAX ===\n")
print(trial_results$results, n = Inf, width = Inf)

p_trial_ctmax <- ggplot(
  trial_results$data,
  aes(trial_length, ctmax, color = group4)
) +
  geom_point(alpha = 0.75, size = 2.2) +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, linewidth = 1,
    alpha = 0.15
  ) +
  geom_text(
    data = trial_results$results,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.3,
    size = 3.2, color = "black"
  ) +
  facet_wrap(~ group4, ncol = 2, scales = "free_x") +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Trial Duration vs. CTmax",
    subtitle = "Linear regressions with 95% confidence intervals",
    x = "Trial duration (minutes)",
    y = "CTmax (°C)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F5F5F5"),
    panel.spacing = grid::unit(0.8, "lines")
  )

print(p_trial_ctmax)

# STARTING TEMPERATURE REGRESSIONS ----
start_results <- run_ctmax_regression(
  regression_data, "start_temp"
)

cat("\n=== STARTING TEMPERATURE VS. CTMAX ===\n")
print(start_results$results, n = Inf, width = Inf)

p_start_ctmax <- ggplot(
  start_results$data,
  aes(start_temp, ctmax, color = group4)
) +
  geom_point(alpha = 0.75, size = 2.2) +
  geom_smooth(
    method = "lm", formula = y ~ x,
    se = TRUE, linewidth = 1,
    alpha = 0.15
  ) +
  geom_text(
    data = start_results$results,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.3,
    size = 3.2, color = "black"
  ) +
  facet_wrap(~ group4, ncol = 2, scales = "free_x") +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Starting Temperature vs. CTmax",
    subtitle = "Linear regressions with 95% confidence intervals",
    x = "Starting temperature (°C)",
    y = "CTmax (°C)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F5F5F5"),
    panel.spacing = grid::unit(0.8, "lines")
  )

print(p_start_ctmax)

# POOLED MODELS: DIFFERENCES IN SLOPES ----
# Group-specific intercepts are included in all models.
# Compare common slopes against group-specific slopes.

compare_slopes <- function(data, predictor) {

  d <- data %>%
    filter(
      is.finite(ctmax),
      is.finite(.data[[predictor]]),
      !is.na(group4)
    )

  common_model <- lm(
    reformulate(c("group4", predictor), response = "ctmax"),
    data = d
  )

  interaction_model <- lm(
    as.formula(paste("ctmax ~ group4 *", predictor)),
    data = d
  )

  cat("\n=== ", toupper(predictor), ": COMMON SLOPE MODEL ===\n", sep = "")
  print(summary(common_model))

  cat("\n=== ", toupper(predictor), ": GROUP-SPECIFIC SLOPES ===\n", sep = "")
  print(summary(interaction_model))

  cat("\n=== ", toupper(predictor), ": SLOPE COMPARISON ===\n", sep = "")
  print(anova(common_model, interaction_model))

  invisible(list(
    common = common_model,
    interaction = interaction_model
  ))
}

trial_models <- compare_slopes(
  regression_data, "trial_length"
)

start_models <- compare_slopes(
  regression_data, "start_temp"
)

# STARTING TEMPERATURE: ADJUST FOR BODY MASS ----
# Evaluate starting temperature while controlling for
# species-habitat group and log-transformed body mass.

start_adjusted_data <- regression_data %>%
  filter(
    is.finite(ctmax),
    is.finite(start_temp),
    is.finite(weight_g),
    weight_g > 0,
    !is.na(group4)
  ) %>%
  mutate(
    log_mass_c = log(weight_g) - mean(log(weight_g)),
    start_temp_c = start_temp - mean(start_temp)
  )

start_adjusted_lm <- lm(
  ctmax ~ group4 + log_mass_c + start_temp_c,
  data = start_adjusted_data
)

cat("\n=== STARTING TEMPERATURE ADJUSTED FOR GROUP AND MASS ===\n")
print(summary(start_adjusted_lm))
print(confint(start_adjusted_lm))

# TEST WHETHER STARTING TEMPERATURE IMPROVES MODEL FIT ----
start_without_lm <- lm(
  ctmax ~ group4 + log_mass_c,
  data = start_adjusted_data
)

cat("\n=== ADDED VALUE OF STARTING TEMPERATURE ===\n")
print(anova(start_without_lm, start_adjusted_lm))

# DESCRIPTIVE SUMMARY OF EXPERIMENTAL VARIABLES ----
experimental_summary <- regression_data %>%
  group_by(group4, .drop = TRUE) %>%
  summarise(
    n_ctmax = sum(is.finite(ctmax)),
    n_duration = sum(is.finite(trial_length)),
    mean_duration = mean(trial_length, na.rm = TRUE),
    sd_duration = sd(trial_length, na.rm = TRUE),
    n_start = sum(is.finite(start_temp)),
    mean_start = mean(start_temp, na.rm = TRUE),
    sd_start = sd(start_temp, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== EXPERIMENTAL VARIABLE SUMMARY ===\n")
print(experimental_summary, n = Inf, width = Inf)

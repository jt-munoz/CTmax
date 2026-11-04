# 13. CTMAX, EXPERIMENTAL CONDITIONS, AND DIAGNOSTICS ----

# PREPARE DATA ----
analysis_data <- data_final %>%
  mutate(
    species = factor(as.character(species), levels = c("Grayling", "Lenok")),
    site_type = factor(tolower(as.character(site_type)),
                       levels = c("lake", "river")),
    group4 = factor(
      paste(str_to_title(as.character(site_type)), species),
      levels = c("Lake Grayling", "River Grayling",
                 "Lake Lenok", "River Lenok")
    )
  )

analysis_colors <- c(
  Lake = "#3278A3", River = "#C45A4A",
  Grayling = "#3278A3", Lenok = "#C45A4A"
)

# FOUR-PANEL COMPARISON DATA ----
make_comparisons <- function(data, variable) {
  base <- data %>%
    transmute(
      species, site_type,
      value = suppressWarnings(as.numeric(as.character(.data[[variable]])))
    ) %>%
    filter(is.finite(value))

  bind_rows(
    base %>%
      mutate(
        panel = paste("Lake vs. River —", species),
        group = str_to_title(as.character(site_type))
      ),
    base %>%
      mutate(
        panel = paste("Grayling vs. Lenok —",
                      str_to_title(as.character(site_type))),
        group = as.character(species)
      )
  ) %>%
    mutate(
      panel = factor(panel, levels = c(
        "Lake vs. River — Grayling",
        "Lake vs. River — Lenok",
        "Grayling vs. Lenok — Lake",
        "Grayling vs. Lenok — River"
      )),
      group = factor(group, levels = c(
        "Lake", "River", "Grayling", "Lenok"
      ))
    )
}

# DISTRIBUTION PLOTTING FUNCTION ----
plot_distributions <- function(d, title, x_label, type = "ecdf") {
  s <- d %>%
    group_by(panel, group, .drop = TRUE) %>%
    summarise(
      n = n(),
      mean = mean(value),
      sd = if (n() > 1) sd(value) else NA_real_,
      min = min(value),
      max = max(value),
      .groups = "drop"
    ) %>%
    mutate(
      range = max - min,
      cv = if_else(mean == 0, NA_real_, 100 * sd / mean),
      label = sprintf(
        "%s (n=%d)\nMean=%.2f | Range=%.2f | CV=%s",
        group, n, mean, range,
        if_else(is.na(cv), "NA", sprintf("%.1f%%", cv))
      ),
      position = if_else(
        group %in% c("Lake", "Grayling"), 1.1, 3.5
      )
    )

  p <- ggplot(d, aes(value, color = group))

  if (type == "density") {
    valid <- s %>%
      filter(n > 1, range > 0) %>%
      select(panel, group)

    p <- p +
      geom_density(
        data = semi_join(d, valid, by = c("panel", "group")),
        aes(fill = group),
        alpha = 0.25, linewidth = 1.1
      ) +
      geom_vline(
        data = s,
        aes(xintercept = mean, color = group),
        linetype = "dashed", linewidth = 0.7
      ) +
      scale_fill_manual(values = analysis_colors) +
      scale_y_continuous(
        expand = expansion(mult = c(0.03, 0.40))
      ) +
      labs(y = "Density",
           subtitle = "Dashed lines indicate group means")
  } else {
    p <- p +
      stat_ecdf(geom = "step", linewidth = 1.1) +
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0, 1, 0.25),
        expand = expansion(mult = c(0.02, 0.32))
      ) +
      labs(y = "Cumulative proportion")
  }

  p +
    geom_text(
      data = s,
      aes(x = Inf, y = Inf, label = label,
          color = group, vjust = position),
      inherit.aes = FALSE,
      hjust = 1.02, size = 2.9,
      fontface = "bold", lineheight = 0.95
    ) +
    facet_wrap(~ panel, ncol = 2, scales = "free") +
    scale_color_manual(values = analysis_colors) +
    labs(title = title, x = x_label) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "#F5F5F5"),
      panel.spacing = grid::unit(0.8, "lines"),
      plot.title = element_text(face = "bold")
    )
}

# CTMAX DISTRIBUTIONS ----
ctmax_data <- make_comparisons(analysis_data, "ctmax")

ctmax_summary <- ctmax_data %>%
  group_by(panel, group, .drop = TRUE) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    min = min(value),
    max = max(value),
    range = max - min,
    cv = 100 * sd / mean,
    .groups = "drop"
  )

cat("\n=== CTMAX SUMMARY ===\n")
print(ctmax_summary, n = Inf, width = Inf)

p_ctmax_density <- plot_distributions(
  ctmax_data, "CTmax Distributions", "CTmax (°C)", "density"
)
p_ctmax_ecdf <- plot_distributions(
  ctmax_data, "Cumulative CTmax Distributions",
  "CTmax (°C)", "ecdf"
)

print(p_ctmax_density)
print(p_ctmax_ecdf)

# CTMAX INDIVIDUAL OBSERVATIONS ----
p_ctmax_points <- ggplot(
  analysis_data %>% filter(is.finite(ctmax)),
  aes(group4, ctmax, color = group4)
) +
  geom_boxplot(
    width = 0.45, outlier.shape = NA,
    color = "grey35", fill = NA
  ) +
  geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 2) +
  scale_color_manual(values = c(
    "#3278A3", "#C45A4A", "#3278A3", "#C45A4A"
  )) +
  labs(
    title = "Individual CTmax Observations",
    x = NULL, y = "CTmax (°C)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 20, hjust = 1))

print(p_ctmax_points)

# EXPERIMENTAL CONDITIONS ----
condition_vars <- c(
  start_temp = "Starting Temperature",
  trial_length = "Trial Duration",
  `O2 start (%)` = "Starting Oxygen (%)",
  `O2 end (%)` = "Ending Oxygen (%)"
)

condition_vars <- condition_vars[
  names(condition_vars) %in% names(analysis_data)
]

experimental_data <- analysis_data %>%
  select(group4, all_of(names(condition_vars))) %>%
  pivot_longer(
    -group4, names_to = "variable", values_to = "raw_value",
    values_transform = list(raw_value = as.character)
  ) %>%
  mutate(
    value = suppressWarnings(parse_number(raw_value)),
    variable = factor(
      variable,
      levels = names(condition_vars),
      labels = unname(condition_vars)
    )
  )

experimental_summary <- experimental_data %>%
  group_by(variable, group4, .drop = TRUE) %>%
  summarise(
    n = sum(is.finite(value)),
    missing = sum(!is.finite(value)),
    mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
    sd = if (sum(is.finite(value)) > 1) sd(value, na.rm = TRUE) else NA_real_,
    min = if (all(is.na(value))) NA_real_ else min(value, na.rm = TRUE),
    max = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== EXPERIMENTAL CONDITIONS AND MISSINGNESS ===\n")
print(experimental_summary, n = Inf, width = Inf)

p_experimental_box <- ggplot(
  experimental_data %>% filter(is.finite(value)),
  aes(group4, value, color = group4)
) +
  geom_boxplot(
    width = 0.45, outlier.shape = NA,
    color = "grey35", fill = NA
  ) +
  geom_jitter(width = 0.12, alpha = 0.65, size = 1.5) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(
    "#3278A3", "#C45A4A", "#3278A3", "#C45A4A"
  )) +
  labs(
    title = "Experimental Conditions by Fish Group",
    x = NULL, y = "Measured value"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8)
  )

print(p_experimental_box)

p_experimental_missing <- ggplot(
  experimental_summary,
  aes(group4, missing, fill = group4)
) +
  geom_col(width = 0.65) +
  geom_text(aes(label = missing), vjust = -0.25, size = 3) +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c(
    "#3278A3", "#C45A4A", "#3278A3", "#C45A4A"
  )) +
  scale_y_continuous(
    breaks = function(x) pretty(x, n = 5),
    expand = expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Missing Experimental Measurements",
    x = NULL, y = "Number missing"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8)
  )

print(p_experimental_missing)

# PARASITE COUNTS AND PREVALENCE ----
if ("gill_parasite_count" %in% names(analysis_data)) {

  parasite_counts <- analysis_data %>%
    mutate(
      gill_count = suppressWarnings(
        as.numeric(as.character(gill_parasite_count))
      )
    ) %>%
    filter(is.finite(gill_count))

  cat("\n=== RECORDED GILL PARASITE COUNTS ===\n")
  print(
    parasite_counts %>%
      group_by(group4) %>%
      summarise(
        n_measured = n(),
        mean_count = mean(gill_count),
        max_count = max(gill_count),
        zero_count = sum(gill_count == 0),
        .groups = "drop"
      ), n = Inf
  )

  if (nrow(parasite_counts) > 0) {
    p_gill_counts <- ggplot(
      parasite_counts, aes(gill_count)
    ) +
      geom_bar(fill = "#3278A3", width = 0.8) +
      facet_wrap(~ group4, scales = "free_y", ncol = 2) +
      scale_x_continuous(breaks = function(x) pretty(x, n = 8)) +
      labs(
        title = "Gill Parasite Count Distributions",
        x = "Number of gill parasites",
        y = "Number of fish"
      ) +
      theme_bw(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))

    print(p_gill_counts)
  }
}

if ("parasite_present" %in% names(analysis_data)) {

  parasite_status <- analysis_data %>%
    mutate(
      status_raw = tolower(trimws(as.character(parasite_present))),
      positive = case_when(
        status_raw %in% c("1", "true", "yes", "present") ~ 1,
        status_raw %in% c("0", "false", "no", "absent") ~ 0,
        TRUE ~ NA_real_
      )
    )

  parasite_prevalence <- parasite_status %>%
    group_by(group4) %>%
    summarise(
      n_assessed = sum(!is.na(positive)),
      n_missing = sum(is.na(positive)),
      n_positive = sum(positive == 1, na.rm = TRUE),
      prevalence = if_else(
        n_assessed > 0, n_positive / n_assessed, NA_real_
      ),
      .groups = "drop"
    )

  cat("\n=== PARASITE PREVALENCE ===\n")
  print(parasite_prevalence, n = Inf)

  p_parasite_prevalence <- ggplot(
    parasite_prevalence %>% filter(n_assessed > 0),
    aes(group4, prevalence, fill = group4)
  ) +
    geom_col(width = 0.65) +
    geom_text(
      aes(label = sprintf("%d/%d", n_positive, n_assessed)),
      vjust = -0.35, size = 3.5
    ) +
    scale_fill_manual(values = c(
      "#3278A3", "#C45A4A", "#3278A3", "#C45A4A"
    )) +
    scale_y_continuous(
      labels = scales::percent, limits = c(0, 1.12),
      breaks = seq(0, 1, 0.25)
    ) +
    labs(
      title = "Recorded Parasite Prevalence",
      subtitle = "Labels show positive fish / assessed fish",
      x = NULL, y = "Proportion with recorded parasites"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 20, hjust = 1)
    )

  print(p_parasite_prevalence)
}

# CTMAX BY EXPERIMENTAL DATE AND BLOCK ----
if ("date" %in% names(analysis_data)) {

  date_data <- analysis_data %>%
    mutate(
      trial_date = as.Date(as.character(date), format = "%m/%d/%Y")
    ) %>%
    filter(!is.na(trial_date), is.finite(ctmax))

  cat("\n=== CTMAX BY DATE AND GROUP ===\n")
  print(
    date_data %>%
      group_by(trial_date, group4) %>%
      summarise(
        n = n(),
        mean_ctmax = mean(ctmax),
        sd_ctmax = if (n() > 1) sd(ctmax) else NA_real_,
        .groups = "drop"
      ), n = Inf
  )

  if (nrow(date_data) > 0) {
    p_ctmax_date <- ggplot(
      date_data, aes(trial_date, ctmax, color = site_type)
    ) +
      geom_point(
        position = position_jitter(width = 0.18, height = 0),
        alpha = 0.75, size = 2
      ) +
      facet_wrap(~ species, ncol = 1) +
      scale_color_manual(values = c(
        lake = "#3278A3", river = "#C45A4A"
      )) +
      scale_x_date(date_labels = "%b %d") +
      labs(
        title = "CTmax by Experimental Date",
        x = "Experimental date",
        y = "CTmax (°C)", color = "Habitat"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "top",
        strip.text = element_text(face = "bold")
      )

    print(p_ctmax_date)
  }
}

if ("block" %in% names(analysis_data)) {

  block_data <- analysis_data %>%
    mutate(block = factor(block)) %>%
    filter(!is.na(block), is.finite(ctmax))

  if (nrow(block_data) > 0) {
    p_ctmax_block <- ggplot(
      block_data, aes(block, ctmax, color = site_type)
    ) +
      geom_boxplot(
        aes(group = interaction(block, site_type)),
        outlier.shape = NA, color = "grey35", fill = NA
      ) +
      geom_jitter(width = 0.15, alpha = 0.7, size = 1.7) +
      facet_wrap(~ species, ncol = 1) +
      scale_color_manual(values = c(
        lake = "#3278A3", river = "#C45A4A"
      )) +
      labs(
        title = "CTmax by Experimental Block",
        x = "Experimental block",
        y = "CTmax (°C)", color = "Habitat"
      ) +
      theme_bw(base_size = 11) +
      theme(
        legend.position = "top",
        strip.text = element_text(face = "bold")
      )

    print(p_ctmax_block)
  }
}

# PRELIMINARY CTMAX MODEL ----
# Diagnostic model only; not the final candidate model set.
# Uses mass as the body-size predictor and allows species
# and habitat groups to have different mean CTmax values.

model_data <- analysis_data %>%
  filter(
    is.finite(ctmax),
    is.finite(weight_g),
    weight_g > 0,
    !is.na(species),
    !is.na(site_type)
  ) %>%
  mutate(
    log_mass_c = log(weight_g) - mean(log(weight_g))
  )

ctmax_diagnostic_lm <- lm(
  ctmax ~ species * site_type + log_mass_c,
  data = model_data
)

cat("\n=== PRELIMINARY CTMAX MODEL ===\n")
print(summary(ctmax_diagnostic_lm))
cat("\nObservations used:", nobs(ctmax_diagnostic_lm), "\n")
print(confint(ctmax_diagnostic_lm))

# COLLINEARITY ----
cat("\n=== MULTICOLLINEARITY ===\n")

if (requireNamespace("car", quietly = TRUE)) {
  print(car::vif(ctmax_diagnostic_lm))
} else {
  cat("Install the car package to calculate VIFs: install.packages('car')\n")
}

# MODEL DIAGNOSTICS ----
diagnostic_data <- model_data %>%
  mutate(
    fitted_ctmax = fitted(ctmax_diagnostic_lm),
    residual = resid(ctmax_diagnostic_lm),
    std_residual = rstandard(ctmax_diagnostic_lm),
    cooks_d = cooks.distance(ctmax_diagnostic_lm)
  )

p_residual_fitted <- ggplot(
  diagnostic_data, aes(fitted_ctmax, residual, color = site_type)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "black") +
  facet_wrap(~ species) +
  scale_color_manual(values = c(
    lake = "#3278A3", river = "#C45A4A"
  )) +
  labs(
    title = "Residuals vs. Fitted CTmax",
    x = "Fitted CTmax (°C)",
    y = "Residual (°C)", color = "Habitat"
  ) +
  theme_bw(base_size = 11)

p_residual_qq <- ggplot(
  diagnostic_data, aes(sample = std_residual)
) +
  stat_qq(alpha = 0.7) +
  stat_qq_line(color = "#C45A4A") +
  labs(
    title = "Normal Q–Q Plot",
    x = "Theoretical quantiles",
    y = "Standardized residuals"
  ) +
  theme_bw(base_size = 11)

p_residual_groups <- ggplot(
  diagnostic_data, aes(group4, residual, color = group4)
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(width = 0.45, color = "grey35", fill = NA,
               outlier.shape = NA) +
  geom_jitter(width = 0.12, alpha = 0.7, size = 1.7) +
  scale_color_manual(values = c(
    "#3278A3", "#C45A4A", "#3278A3", "#C45A4A"
  )) +
  labs(
    title = "Residuals by Fish Group",
    x = NULL, y = "Residual (°C)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

p_influence <- ggplot(
  diagnostic_data,
  aes(seq_along(cooks_d), cooks_d)
) +
  geom_segment(
    aes(xend = seq_along(cooks_d), yend = 0),
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = 4 / nrow(diagnostic_data),
    linetype = "dashed", color = "#C45A4A"
  ) +
  labs(
    title = "Cook's Distance",
    subtitle = "Dashed line = 4/n screening threshold",
    x = "Model observation number",
    y = "Cook's distance"
  ) +
  theme_bw(base_size = 11)

print(p_residual_fitted)
print(p_residual_qq)
print(p_residual_groups)
print(p_influence)

# RESIDUALS BY BLOCK ----
if ("block" %in% names(diagnostic_data)) {
  p_residual_block <- ggplot(
    diagnostic_data %>% filter(!is.na(block)),
    aes(factor(block), residual, color = site_type)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_jitter(width = 0.15, size = 1.8, alpha = 0.75) +
    facet_wrap(~ species) +
    scale_color_manual(values = c(
      lake = "#3278A3", river = "#C45A4A"
    )) +
    labs(
      title = "Model Residuals by Experimental Block",
      x = "Block", y = "Residual (°C)", color = "Habitat"
    ) +
    theme_bw(base_size = 11)

  print(p_residual_block)
}

# INFLUENTIAL OBSERVATIONS ----
cat("\n=== LARGEST COOK'S DISTANCES ===\n")

print(
  diagnostic_data %>%
    mutate(model_row = row_number()) %>%
    arrange(desc(cooks_d)) %>%
    select(
      model_row, species, site_type,
      ctmax, weight_g, residual, cooks_d
    ) %>%
    slice_head(n = 10),
  n = 10
)


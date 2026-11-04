# 8. LENGTH-WEIGHT RELATIONSHIPS ----

# PREPARE DATA ----
lw_base <- data_final %>%
  filter(!is.na(length_mm), !is.na(weight_g),
         length_mm > 0, weight_g > 0) %>%
  mutate(log_weight = log(weight_g))

lw_data <- bind_rows(
  lw_base %>%
    mutate(
      panel = paste("Lake vs. River —", species),
      comparison = str_to_title(as.character(site_type))
    ),
  lw_base %>%
    mutate(
      panel = paste("Grayling vs. Lenok —",
                    str_to_title(as.character(site_type))),
      comparison = as.character(species)
    )
) %>%
  mutate(
    panel = factor(panel, levels = c(
      "Lake vs. River — Grayling",
      "Lake vs. River — Lenok",
      "Grayling vs. Lenok — Lake",
      "Grayling vs. Lenok — River"
    )),
    comparison = factor(
      comparison, levels = c("Lake", "River", "Grayling", "Lenok")
    )
  )

lw_colors <- c(
  Lake = "#3278A3",
  River = "#C45A4A",
  Grayling = "#3278A3",
  Lenok = "#C45A4A"
)

# FIT LOG-LOG ANCOVA MODELS ----
# Center log-length at the geometric mean length within each panel.

lw_fits <- lw_data %>%
  group_nest(panel) %>%
  mutate(
    ref_length = map_dbl(data, ~ exp(mean(log(.x$length_mm)))),
    model = map2(data, ref_length, \(d, ref) {
      d <- droplevels(d)
      d$log_length_c <- log(d$length_mm / ref)
      lm(log_weight ~ log_length_c * comparison, data = d)
    })
  )

# TEST SLOPE AND WEIGHT DIFFERENCES ----
# Slope: ANCOVA interaction test.
# Weight: ANCOVA group coefficient at the reference length.
# Holm corrections are applied separately across four panels per test.

lw_tests <- lw_fits %>%
  mutate(results = pmap(
    list(data, model, ref_length),
    function(d, model, ref) {
      d <- droplevels(d)
      d$log_length_c <- log(d$length_mm / ref)

      groups <- levels(d$comparison)
      common <- lm(log_weight ~ log_length_c + comparison, data = d)
      ci <- confint(model)

      group_term <- paste0("comparison", groups[2])
      slope_term <- paste0("log_length_c:comparison", groups[2])

      weight_effect <- unname(coef(model)[group_term])
      slope_effect <- unname(coef(model)[slope_term])

      tibble(
        group_1 = groups[1],
        group_2 = groups[2],
        n_1 = sum(d$comparison == groups[1]),
        n_2 = sum(d$comparison == groups[2]),
        ref_length_mm = ref,
        slope_1 = unname(coef(model)["log_length_c"]),
        slope_2 = unname(coef(model)["log_length_c"]) + slope_effect,
        slope_difference = slope_effect,
        slope_lower = ci[slope_term, 1],
        slope_upper = ci[slope_term, 2],
        p_slope = anova(common, model)[["Pr(>F)"]][2],
        weight_difference_pct = 100 * (exp(weight_effect) - 1),
        weight_lower_pct = 100 * (exp(ci[group_term, 1]) - 1),
        weight_upper_pct = 100 * (exp(ci[group_term, 2]) - 1),
        p_weight = summary(model)$coefficients[group_term, "Pr(>|t|)"]
      )
    }
  )) %>%
  select(panel, results) %>%
  unnest(results) %>%
  mutate(
    p_slope_holm = p.adjust(p_slope, method = "holm"),
    p_weight_holm = p.adjust(p_weight, method = "holm"),
    significant = p_slope_holm < 0.05 | p_weight_holm < 0.05,
    slope_ci = sprintf("%+.3f [%.3f, %.3f]",
                       slope_difference, slope_lower, slope_upper),
    weight_ci = sprintf("%+.1f%% [%.1f%%, %.1f%%]",
                        weight_difference_pct,
                        weight_lower_pct, weight_upper_pct),
    annotation = paste0(
      "Slope: ", p_label(p_slope_holm),
      "  |  Weight: ", p_label(p_weight_holm)
    )
  )

# PREDICT POWER-LAW CURVES ----
# Models use log-transformed values; figures use millimeters and grams.

lw_curves <- lw_fits %>%
  mutate(curves = pmap(
    list(data, model, ref_length),
    function(d, model, ref) {
      grid <- d %>%
        group_by(comparison, .drop = TRUE) %>%
        summarise(
          min_length = min(length_mm),
          max_length = max(length_mm),
          .groups = "drop"
        ) %>%
        mutate(
          length_mm = map2(
            min_length, max_length,
            ~ seq(.x, .y, length.out = 100)
          )
        ) %>%
        select(comparison, length_mm) %>%
        unnest(length_mm) %>%
        mutate(log_length_c = log(length_mm / ref))

      grid$predicted_weight <- exp(predict(model, newdata = grid))
      grid
    }
  )) %>%
  select(panel, curves) %>%
  unnest(curves)

# SAMPLE SIZES AND GROUP LABELS ----
lw_counts <- lw_data %>%
  count(panel, comparison, name = "n") %>%
  left_join(
    lw_data %>%
      group_by(panel) %>%
      summarise(
        min_length = min(length_mm),
        max_length = max(length_mm),
        .groups = "drop"
      ),
    by = "panel"
  ) %>%
  mutate(
    is_left = as.character(comparison) %in% c("Lake", "Grayling"),
    x = if_else(is_left, min_length, max_length),
    hjust = if_else(is_left, 0, 1),
    label = paste0(as.character(comparison), "  n=", n)
  )

# VISUALIZE ALL FOUR COMPARISONS ----
p_length_weight <- ggplot(
  lw_data, aes(length_mm, weight_g, color = comparison)
) +
  geom_rect(
    data = filter(lw_tests, significant),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#E0EEE7"
  ) +
  geom_point(size = 1.9, alpha = 0.65) +
  geom_line(
    data = lw_curves,
    aes(x = length_mm, y = predicted_weight,
        color = comparison, group = comparison),
    inherit.aes = FALSE, linewidth = 1.2
  ) +
  geom_text(
    data = lw_tests,
    aes(x = Inf, y = Inf, label = annotation),
    inherit.aes = FALSE,
    hjust = 1.05, vjust = 1.3, size = 3.3
  ) +
  geom_text(
    data = lw_counts,
    aes(x = x, y = -Inf, label = label,
        color = comparison, hjust = hjust),
    inherit.aes = FALSE,
    fontface = "bold", size = 4.3, vjust = -0.4
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "free") +
  scale_color_manual(values = lw_colors) +
  scale_x_continuous(expand = expansion(mult = 0.07)) +
  scale_y_continuous(expand = expansion(mult = c(0.13, 0.15))) +
  labs(
    title = "Length–Weight Relationships",
    subtitle = "Log–log ANCOVA | Holm-adjusted p-values | Shaded: p < 0.05",
    x = "Total length (mm)",
    y = "Body mass (g)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "#F4F4F4"),
    panel.spacing = unit(0.7, "lines"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10)
  )

print(p_length_weight)

# DISPLAY STATISTICAL RESULTS ----
cat("\n=== LENGTH-WEIGHT ANCOVA RESULTS ===\n")
cat("Effects = second group minus first group.\n")
cat("Weight differences evaluated at each panel's reference length.\n")
cat("Holm-adjusted p-values used for significance at alpha = 0.05.\n\n")

lw_tests %>%
  select(
    panel, group_1, group_2, n_1, n_2, ref_length_mm,
    slope_1, slope_2, slope_ci, p_slope, p_slope_holm,
    weight_ci, p_weight, p_weight_holm, significant
  ) %>%
  print(n = Inf, width = Inf)

# INDIVIDUAL MODEL SUMMARIES ----
for (i in seq_len(nrow(lw_fits))) {
  cat("\n========================================\n")
  cat("COMPARISON:", as.character(lw_fits$panel[i]), "\n")
  cat("========================================\n")

  cat("\nMODEL SUMMARY:\n")
  print(summary(lw_fits$model[[i]]))

  cat("\n95% CONFIDENCE INTERVALS:\n")
  print(confint(lw_fits$model[[i]]))

  cat("\nADJUSTED P-VALUES:\n")
  print(
    lw_tests %>%
      filter(panel == lw_fits$panel[i]) %>%
      select(p_slope_holm, p_weight_holm)
  )
}

# 9. SIZE DISTRIBUTIONS BY GROUP ----

# PREPARE DATA ----
size_hist_data <- data_final %>%
  mutate(
    group = factor(
      paste(species, str_to_title(site_type)),
      levels = c("Grayling Lake", "Grayling River",
                 "Lenok Lake", "Lenok River")
    )
  ) %>%
  select(group, length_mm, weight_g) %>%
  pivot_longer(
    c(length_mm, weight_g),
    names_to = "variable", values_to = "value"
  ) %>%
  mutate(
    variable = factor(
      variable,
      levels = c("length_mm", "weight_g"),
      labels = c("Total length (mm)", "Body mass (g)")
    )
  ) %>%
  filter(!is.na(value))

group_colors <- c(
  "Grayling Lake" = "#3278A3",
  "Grayling River" = "#C45A4A",
  "Lenok Lake" = "#3278A3",
  "Lenok River" = "#C45A4A"
)

# SUMMARY STATISTICS ----
size_hist_summary <- size_hist_data %>%
  group_by(variable, group) %>%
  summarise(
    n = n(),
    median = median(value),
    mean = mean(value),
    sd = sd(value),
    min = min(value),
    max = max(value),
    .groups = "drop"
  )

print(size_hist_summary, n = Inf)

# HISTOGRAMS ----
p_size_distributions <- ggplot(size_hist_data, aes(value, fill = group)) +
  geom_histogram(
    data = filter(size_hist_data, variable == "Total length (mm)"),
    binwidth = 20, boundary = 0,
    color = "white", linewidth = 0.3, alpha = 0.8
  ) +
  geom_histogram(
    data = filter(size_hist_data, variable == "Body mass (g)"),
    binwidth = 50, boundary = 0,
    color = "white", linewidth = 0.3, alpha = 0.8
  ) +
  geom_vline(
    data = size_hist_summary,
    aes(xintercept = median),
    linetype = "dashed", linewidth = 0.8,
    color = "#333333"
  ) +
  geom_text(
    data = size_hist_summary,
    aes(x = Inf, y = Inf,
        label = paste0("n = ", n, "\nMedian = ", round(median, 1))),
    inherit.aes = FALSE,
    hjust = 1.1, vjust = 1.3, size = 3.1
  ) +
  facet_grid(variable ~ group, scales = "free_x") +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "Length and Weight Distributions by Sampling Group",
    subtitle = "Dashed lines indicate group medians",
    x = NULL,
    y = "Number of fish"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    strip.text.y = element_text(angle = 0),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(face = "bold", size = 13)
  )

print(p_size_distributions)

# 10. SIZE DISTRIBUTION OVERLAP ----

# PREPARE DATA ----
size_base <- data_final %>%
  select(species, site_type, length_mm, weight_g) %>%
  pivot_longer(c(length_mm, weight_g),
               names_to = "variable", values_to = "value") %>%
  filter(!is.na(value))

size_data <- bind_rows(
  size_base %>%
    mutate(
      panel = paste("Lake vs. River —", species),
      group = str_to_title(as.character(site_type))
    ),
  size_base %>%
    mutate(
      panel = paste("Grayling vs. Lenok —",
                    str_to_title(as.character(site_type))),
      group = as.character(species)
    )
) %>%
  mutate(
    variable = factor(
      variable, levels = c("length_mm", "weight_g"),
      labels = c("Length", "Weight")
    ),
    panel = factor(panel, levels = c(
      "Lake vs. River — Grayling",
      "Lake vs. River — Lenok",
      "Grayling vs. Lenok — Lake",
      "Grayling vs. Lenok — River"
    )),
    group = factor(group, levels = c("Lake", "River", "Grayling", "Lenok"))
  )

size_colors <- c(
  Lake = "#3278A3", River = "#C45A4A",
  Grayling = "#3278A3", Lenok = "#C45A4A"
)

# SUMMARY STATISTICS ----
size_summary <- size_data %>%
  group_by(variable, panel, group, .drop = TRUE) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    min = min(value),
    max = max(value),
    .groups = "drop"
  ) %>%
  mutate(
    range = max - min,
    cv = 100 * sd / mean,
    label = sprintf(
      "%s: n=%d | Mean=%.1f | Range=%.1f | CV=%.1f%%",
      group, n, mean, range, cv
    ),
    label_position = if_else(group %in% c("Lake", "Grayling"), 1.2, 2.8)
  )

cat("\n=== SIZE DISTRIBUTION SUMMARY ===\n")
print(size_summary, n = Inf, width = Inf)

# OBSERVED RANGE OVERLAP ----
size_overlap <- size_summary %>%
  group_by(variable, panel, .drop = TRUE) %>%
  summarise(
    overlap_min = max(min),
    overlap_max = min(max),
    overlap_width = pmax(0, overlap_max - overlap_min),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_min = if_else(overlap_width > 0, overlap_min, NA_real_),
    overlap_max = if_else(overlap_width > 0, overlap_max, NA_real_)
  )

cat("\n=== OBSERVED SIZE OVERLAP ===\n")
print(size_overlap, n = Inf, width = Inf)

# PLOTTING FUNCTION ----
plot_size <- function(type = "density", variable_name = NULL) {

  d <- size_data
  s <- size_summary

  if (!is.null(variable_name)) {
    d <- filter(d, variable == variable_name)
    s <- filter(s, variable == variable_name)
  }

  # Distinct facet IDs keep length and weight separate in overviews.
  facet_levels <- unlist(lapply(levels(size_data$variable),
                                \(v) paste(v, levels(size_data$panel), sep = " | ")))
  d <- d %>% mutate(facet = factor(paste(variable, panel, sep = " | "), levels = facet_levels))
  s <- s %>% mutate(facet = factor(paste(variable, panel, sep = " | "), levels = facet_levels))
  p <- ggplot(d, aes(value, color = group))

  if (type == "density") {
    p <- p +
      geom_density(
        aes(fill = group),
        alpha = 0.25, linewidth = 1.1, adjust = 1
      ) +
      geom_vline(
        data = s,
        aes(xintercept = mean, color = group),
        linetype = "dashed", linewidth = 0.7
      ) +
      scale_fill_manual(values = size_colors) +
      labs(
        y = "Density",
        subtitle = "Dashed lines indicate group means"
      )
  } else {
    p <- p +
      stat_ecdf(geom = "step", linewidth = 1.1) +
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0, 1, 0.25),
        expand = expansion(mult = c(0.02, 0.18))
      ) +
      labs(y = "Cumulative proportion")
  }

  p +
    geom_text(
      data = s,
      aes(x = Inf, y = Inf, label = label,
          color = group, vjust = label_position),
      inherit.aes = FALSE,
      hjust = 1.02, size = 3, fontface = "bold"
    ) +
    facet_wrap(
      ~ facet, ncol = 2, scales = "free",
      labeller = labeller(facet = function(x) {
        if (is.null(variable_name)) x else sub("^.* \\| ", "", x)
      })
    ) +
    scale_color_manual(values = size_colors) +
    labs(
      title = if (is.null(variable_name)) {
        if (type == "density") {
          "Length and Weight Distributions"
        } else {
          "Cumulative Length and Weight Distributions"
        }
      } else {
        paste(
          if (type == "density") "Size Distributions" else "Cumulative Distributions",
          "—", variable_name
        )
      },
      x = if (is.null(variable_name)) {
        "Total length (mm) or body mass (g)"
      } else if (variable_name == "Length") {
        "Total length (mm)"
      } else {
        "Body mass (g)"
      }
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "#F5F5F5"),
      panel.spacing = unit(0.6, "lines"),
      plot.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 9)
    )
}

# EIGHT-PANEL OVERVIEW FIGURES ----
p_size_density <- plot_size("density")
p_size_ecdf <- plot_size("ecdf")

print(p_size_density)
print(p_size_ecdf)

# FOUR-PANEL LENGTH FIGURES ----
p_length_density <- plot_size("density", "Length")
p_length_ecdf <- plot_size("ecdf", "Length")

print(p_length_density)
print(p_length_ecdf)

# FOUR-PANEL WEIGHT FIGURES ----
p_weight_density <- plot_size("density", "Weight")
p_weight_ecdf <- plot_size("ecdf", "Weight")

print(p_weight_density)
print(p_weight_ecdf)


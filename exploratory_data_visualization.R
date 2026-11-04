
# ============================================================
# STEP 2: EXPLORATORY DATA ANALYSIS
# Author: Joseph Munoz
# Descriptive comparisons -> morphology/condition -> experimental checks -> diagnostics.
# Run top to bottom. All figures and test tables remain available as objects.
# ============================================================

# 1. PREPARE DATA ----

library(tidyverse)
set.seed(123)

data_final <- read_csv("CTmax_Merged_Data.csv", show_col_types = FALSE) %>%
  mutate(
    species = factor(species, levels = c("Grayling", "Lenok")),
    site_type = factor(site_type, levels = c("lake", "river"))
  )

habitat_colors <- c(lake = "#3278A3", river = "#C45A4A")
species_colors <- c(Grayling = "#648B62", Lenok = "#9770A8")

variables <- c(
  ctmax = "CTmax (°C)",
  weight_g = "Body mass (g)",
  length_mm = "Total length (mm)",
  fultons_k = "Fulton's K",
  parasite_burden = "Parasite burden (g/g)"
)

p_label <- function(p) {
  ifelse(is.na(p), "Not tested",
         ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))
}

bootstrap_diff <- function(x, y, fun, B = 2000) {
  estimates <- replicate(B,
                         fun(sample(y, length(y), replace = TRUE)) -
                           fun(sample(x, length(x), replace = TRUE)))
  unname(quantile(estimates, c(0.025, 0.975)))
}

# 2. SUMMARY STATISTICS ----

summary_stats <- data_final %>%
  group_by(species, site_type) %>%
  summarise(
    n = n(),
    ctmax_mean = mean(ctmax, na.rm = TRUE),
    ctmax_sd = sd(ctmax, na.rm = TRUE),
    ctmax_min = min(ctmax, na.rm = TRUE),
    ctmax_max = max(ctmax, na.rm = TRUE),
    mass_mean = mean(weight_g, na.rm = TRUE),
    mass_sd = sd(weight_g, na.rm = TRUE),
    length_mean = mean(length_mm, na.rm = TRUE),
    length_sd = sd(length_mm, na.rm = TRUE),
    k_mean = mean(fultons_k, na.rm = TRUE),
    k_sd = sd(fultons_k, na.rm = TRUE),
    parasite_median = median(parasite_burden, na.rm = TRUE),
    parasite_q25 = quantile(parasite_burden, 0.25, na.rm = TRUE),
    parasite_q75 = quantile(parasite_burden, 0.75, na.rm = TRUE),
    parasite_prevalence = 100 * mean(parasite_present == 1, na.rm = TRUE),
    .groups = "drop"
  )

table1 <- summary_stats %>%
  transmute(
    species, site_type,
    `Sample size` = as.character(n),
    `CTmax mean ± SD (°C)` = sprintf("%.2f ± %.2f", ctmax_mean, ctmax_sd),
    `CTmax range (°C)` = sprintf("%.2f–%.2f", ctmax_min, ctmax_max),
    `Body mass mean ± SD (g)` = sprintf("%.1f ± %.1f", mass_mean, mass_sd),
    `Total length mean ± SD (mm)` = sprintf("%.1f ± %.1f", length_mean, length_sd),
    `Fulton's K mean ± SD` = sprintf("%.3f ± %.3f", k_mean, k_sd),
    `Parasite burden median (IQR)` = sprintf(
      "%.4f (%.4f–%.4f)", parasite_median, parasite_q25, parasite_q75),
    `Parasite prevalence (%)` = sprintf("%.1f", parasite_prevalence)
  ) %>%
  pivot_longer(-c(species, site_type), names_to = "Variable", values_to = "Value") %>%
  mutate(group = factor(
    paste(species, str_to_title(site_type)),
    levels = c("Grayling Lake", "Grayling River", "Lenok Lake", "Lenok River")
  )) %>%
  select(Variable, group, Value) %>%
  pivot_wider(names_from = group, values_from = Value)

print(table1, n = Inf, width = Inf)
View(table1)

# 3. STATISTICAL COMPARISON FUNCTIONS ----

# Effect = second group minus first group.
# Continuous: Welch mean difference and 95% CI.
# Parasite burden: median difference and bootstrap 95% CI.
# Prevalence: percentage-point difference and bootstrap 95% CI.

compare_groups <- function(data, within, between) {
  data %>%
    group_by(across(all_of(within))) %>%
    group_modify(~ {
      d <- .x

      tests <- map_dfr(names(variables), function(v) {
        groups <- split(d[[v]], d[[between]], drop = TRUE)
        groups <- lapply(groups, \(x) x[!is.na(x)])

        if (length(groups) != 2 || any(lengths(groups) < 2)) {
          return(tibble(variable = v, test = "Not tested",
                        effect = NA_real_, lower = NA_real_,
                        upper = NA_real_, p = NA_real_))
        }

        x <- groups[[1]]
        y <- groups[[2]]

        if (v == "parasite_burden") {
          ci <- bootstrap_diff(x, y, median)
          p <- suppressWarnings(
            wilcox.test(x, y, exact = FALSE)$p.value)
          effect <- median(y) - median(x)
          test <- "Wilcoxon"
        } else {
          t <- tryCatch(t.test(y, x), error = function(e) NULL)
          ci <- if (is.null(t)) c(NA_real_, NA_real_) else t$conf.int[1:2]
          p <- if (is.null(t)) NA_real_ else t$p.value
          effect <- mean(y) - mean(x)
          test <- "Welch t-test"
        }

        tibble(variable = v, test = test, effect = effect,
               lower = ci[1], upper = ci[2], p = p)
      })

      prev <- d %>%
        filter(!is.na(parasite_present)) %>%
        mutate(present = as.integer(parasite_present == 1))

      groups <- split(prev$present, prev[[between]], drop = TRUE)

      if (length(groups) == 2 && all(lengths(groups) > 0)) {
        x <- groups[[1]]
        y <- groups[[2]]
        ci <- 100 * bootstrap_diff(x, y, mean)
        effect <- 100 * (mean(y) - mean(x))
        p <- fisher.test(table(
          prev[[between]], factor(prev$present, levels = 0:1)
        ))$p.value
      } else {
        ci <- c(NA_real_, NA_real_)
        effect <- p <- NA_real_
      }

      bind_rows(
        tests,
        tibble(variable = "parasite_prevalence",
               test = "Fisher exact", effect = effect,
               lower = ci[1], upper = ci[2], p = p)
      ) %>%
        mutate(p_holm = p.adjust(p, method = "holm"))
    }) %>%
    ungroup() %>%
    mutate(
      variable_label = recode(
        variable, !!!variables,
        parasite_prevalence = "Parasite prevalence"
      ),
      significant = !is.na(p_holm) & p_holm < 0.05,
      p_display = p_label(p_holm)
    )
}

# 4. COMPACT BOXPLOT FUNCTION ----

# CTmax always appears first, followed by the remaining variables.
# Significant panels receive a subtle background tint.

plot_groups <- function(data, within, between, colors, tests, title) {

  panel_levels <- unlist(lapply(
    unname(variables),
    \(v) paste(v, levels(data[[within]]), sep = " | ")
  ))

  d <- data %>%
    select(all_of(c(within, between, names(variables)))) %>%
    pivot_longer(all_of(names(variables)),
                 names_to = "variable", values_to = "value") %>%
    mutate(
      panel = factor(
        paste(recode(variable, !!!variables),
              .data[[within]], sep = " | "),
        levels = panel_levels
      )
    )

  counts <- d %>%
    filter(!is.na(value)) %>%
    count(panel, .data[[between]], name = "n")

  labels <- tests %>%
    filter(variable != "parasite_prevalence") %>%
    mutate(
      panel = factor(
        paste(variable_label, .data[[within]], sep = " | "),
        levels = panel_levels
      ),
      label = paste(test, p_display, sep = " | ")
    )

  ggplot(d, aes(.data[[between]], value, fill = .data[[between]])) +
    geom_rect(
      data = filter(labels, significant),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = "#E0EEE7"
    ) +
    geom_boxplot(width = 0.50, alpha = 0.45,
                 outlier.shape = NA, na.rm = TRUE) +
    geom_jitter(aes(color = .data[[between]]),
                width = 0.10, size = 1.5, alpha = 0.65, na.rm = TRUE) +
    geom_text(
      data = labels,
      aes(x = 1.5, y = Inf, label = label),
      inherit.aes = FALSE, vjust = 1.2, size = 3
    ) +
    geom_text(
      data = counts,
      aes(x = .data[[between]], y = -Inf,
          label = paste0("n=", n)),
      inherit.aes = FALSE, vjust = -0.4, size = 3
    ) +
    facet_wrap(~ panel, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 9),
      strip.background = element_rect(fill = "#F5F5F5"),
      panel.spacing = unit(0.5, "lines"),
      plot.title = element_text(face = "bold", size = 12)
    )
}

# 5. HABITAT COMPARISONS WITHIN SPECIES ----

habitat_tests <- compare_groups(data_final, "species", "site_type")

p_habitat <- plot_groups(
  data_final, "species", "site_type", habitat_colors,
  habitat_tests, "Lake vs. River Within Species"
)

print(habitat_tests, n = Inf, width = Inf)
print(p_habitat)

# 6. SPECIES COMPARISONS WITHIN HABITAT ----

species_tests <- compare_groups(data_final, "site_type", "species")

p_species <- plot_groups(
  data_final, "site_type", "species", species_colors,
  species_tests, "Grayling vs. Lenok Within Habitat"
)

print(species_tests, n = Inf, width = Inf)
print(p_species)

# 7. CTMAX COMPARISONS AMONG SITES ----
# Verify raw site labels against site_name in the merged file before reporting results.

# PREPARE SITES ----
site_data <- data_final %>%
  mutate(
    site = factor(site, levels = c(
      "lake_campsite2", "stomach_bay", "refuge_camp",
      "eg uur camp", "uur spring", "last chance"
    )),
    panel = factor(
      paste(str_to_title(site_type), species, sep = " | "),
      levels = c(
        "Lake | Grayling", "Lake | Lenok",
        "River | Grayling", "River | Lenok"
      )
    )
  )

site_labels <- c(
  lake_campsite2 = "Campsite 2",
  stomach_bay = "Stomach Bay",
  refuge_camp = "Refuge Camp",
  "eg uur camp" = "Eg-Uur Camp",
  "uur spring" = "Uur Spring",
  "last chance" = "Last Chance"
)

# SITE SUMMARY ----
site_summary <- site_data %>%
  group_by(site_type, species, site, .drop = TRUE) %>%
  summarise(
    n = sum(!is.na(ctmax)),
    mean_ctmax = mean(ctmax, na.rm = TRUE),
    sd_ctmax = sd(ctmax, na.rm = TRUE),
    min_ctmax = min(ctmax, na.rm = TRUE),
    max_ctmax = max(ctmax, na.rm = TRUE),
    .groups = "drop"
  )

print(site_summary, n = Inf)

# OVERALL SITE TESTS ----
site_tests <- site_data %>%
  group_by(site_type, species, panel) %>%
  group_modify(~ {
    d <- droplevels(filter(.x, !is.na(ctmax), !is.na(site)))
    sizes <- table(d$site)

    t <- if (length(sizes) >= 2 && all(sizes >= 2)) {
      tryCatch(oneway.test(ctmax ~ site, data = d),
               error = function(e) NULL)
    } else NULL

    tibble(
      n_sites = length(sizes),
      test = "Welch ANOVA",
      p = if (is.null(t)) NA_real_ else t$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    significant = !is.na(p) & p < 0.05,
    label = paste(test, p_label(p), sep = " | ")
  )

print(site_tests, n = Inf)

# PAIRWISE SITE DIFFERENCES AND 95% CI ----
site_pairwise <- site_data %>%
  group_by(site_type, species) %>%
  group_modify(~ {
    d <- droplevels(filter(.x, !is.na(ctmax), !is.na(site)))
    sites <- levels(d$site)

    if (length(sites) < 2) {
      return(tibble(
        site_1 = character(), site_2 = character(),
        effect = numeric(), lower = numeric(),
        upper = numeric(), p = numeric()
      ))
    }

    map_dfr(combn(sites, 2, simplify = FALSE), function(s) {
      x <- d$ctmax[d$site == s[1]]
      y <- d$ctmax[d$site == s[2]]

      t <- if (length(x) >= 2 && length(y) >= 2) {
        tryCatch(t.test(y, x), error = function(e) NULL)
      } else NULL

      tibble(
        site_1 = s[1],
        site_2 = s[2],
        effect = mean(y) - mean(x),
        lower = if (is.null(t)) NA_real_ else t$conf.int[1],
        upper = if (is.null(t)) NA_real_ else t$conf.int[2],
        p = if (is.null(t)) NA_real_ else t$p.value
      )
    })
  }) %>%
  group_by(site_type, species) %>%
  mutate(p_holm = p.adjust(p, method = "holm")) %>%
  ungroup() %>%
  mutate(
    test = "Pairwise Welch t-test",
    effect_ci = sprintf("%.2f [%.2f, %.2f]", effect, lower, upper),
    p_display = p_label(p_holm)
  )

print(site_pairwise, n = Inf, width = Inf)

# SITE BOXPLOTS ----
site_counts <- site_data %>%
  filter(!is.na(ctmax)) %>%
  count(panel, site, name = "n")

p_sites <- ggplot(site_data, aes(site, ctmax, fill = site_type)) +
  geom_rect(
    data = filter(site_tests, significant),
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    inherit.aes = FALSE, fill = "#E0EEE7"
  ) +
  geom_boxplot(width = 0.5, alpha = 0.45,
               outlier.shape = NA, na.rm = TRUE) +
  geom_jitter(aes(color = site_type),
              width = 0.10, size = 1.6, alpha = 0.65, na.rm = TRUE) +
  geom_text(
    data = site_tests,
    aes(x = Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = 1.05, vjust = 1.2, size = 3
  ) +
  geom_text(
    data = site_counts,
    aes(x = site, y = -Inf, label = paste0("n=", n)),
    inherit.aes = FALSE, vjust = -0.4, size = 3
  ) +
  facet_wrap(~ panel, ncol = 2, scales = "free") +
  scale_x_discrete(labels = site_labels, drop = TRUE) +
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  labs(title = "CTmax Across Sampling Sites",
       x = NULL, y = "CTmax (°C)") +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "#F5F5F5"),
    axis.text.x = element_text(size = 8),
    panel.spacing = unit(0.5, "lines"),
    plot.title = element_text(face = "bold", size = 12)
  )

print(p_sites)

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

# 11. PARASITE BURDEN AND FULTON'S K DISTRIBUTIONS ----

# PREPARE DATA ----
condition_base <- data_final %>%
  select(species, site_type, parasite_burden, fultons_k) %>%
  pivot_longer(
    c(parasite_burden, fultons_k),
    names_to = "variable", values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  mutate(
    species = str_to_title(as.character(species)),
    site_type = str_to_title(as.character(site_type)),
    variable = factor(
      variable,
      levels = c("parasite_burden", "fultons_k"),
      labels = c("Parasite Burden", "Fulton's K")
    )
  )

condition_data <- bind_rows(
  condition_base %>%
    mutate(
      panel = paste("Lake vs. River —", species),
      group = site_type
    ),
  condition_base %>%
    mutate(
      panel = paste("Grayling vs. Lenok —", site_type),
      group = species
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

condition_colors <- c(
  Lake = "#3278A3", River = "#C45A4A",
  Grayling = "#3278A3", Lenok = "#C45A4A"
)

# SUMMARY STATISTICS ----
condition_summary <- condition_data %>%
  group_by(variable, panel, group, .drop = TRUE) %>%
  summarise(
    n = n(),
    mean = mean(value),
    sd = if (n() > 1) sd(value) else NA_real_,
    median = median(value),
    min = min(value),
    max = max(value),
    zero_n = sum(value == 0),
    .groups = "drop"
  ) %>%
  mutate(
    range = max - min,
    cv = if_else(mean == 0, NA_real_, 100 * sd / mean),
    zero_pct = 100 * zero_n / n,
    constant_group = n < 2 | range == 0,
    mean_text = if_else(
      variable == "Fulton's K",
      sprintf("%.3f", mean),
      sprintf("%.3g", mean)
    ),
    range_text = if_else(
      variable == "Fulton's K",
      sprintf("%.3f", range),
      sprintf("%.3g", range)
    ),
    cv_text = if_else(is.na(cv), "NA", sprintf("%.1f%%", cv)),
    label = sprintf(
      "%s (n=%d)\nMean=%s | Range=%s | CV=%s",
      group, n, mean_text, range_text, cv_text
    ),
    label_position = if_else(
      group %in% c("Lake", "Grayling"), 1.1, 3.5
    )
  )

cat("\n=== PARASITE BURDEN AND FULTON'S K SUMMARY ===\n")
print(condition_summary, n = Inf, width = Inf)

# OBSERVED RANGE OVERLAP ----
condition_overlap <- condition_summary %>%
  group_by(variable, panel, .drop = TRUE) %>%
  summarise(
    overlap_min = max(min),
    overlap_max = min(max),
    overlap_width = pmax(0, overlap_max - overlap_min),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_min = if_else(
      overlap_max >= overlap_min, overlap_min, NA_real_
    ),
    overlap_max = if_else(
      !is.na(overlap_min), overlap_max, NA_real_
    )
  )

cat("\n=== OBSERVED RANGE OVERLAP ===\n")
print(condition_overlap, n = Inf, width = Inf)

# PLOTTING FUNCTION ----
plot_condition <- function(type = "density", variable_name = NULL) {

  d <- condition_data
  s <- condition_summary

  if (!is.null(variable_name)) {
    d <- d %>% filter(variable == variable_name)
    s <- s %>% filter(variable == variable_name)
  }

  d <- d %>%
    mutate(facet = paste(variable, panel, sep = " | "))

  s <- s %>%
    mutate(facet = paste(variable, panel, sep = " | "))

  facet_levels <- unlist(lapply(
    levels(condition_data$variable),
    function(v) paste(v, levels(condition_data$panel), sep = " | ")
  ))

  d$facet <- factor(d$facet, levels = facet_levels)
  s$facet <- factor(s$facet, levels = facet_levels)

  # Groups with no variation cannot have a density curve.
  density_groups <- s %>%
    filter(!constant_group) %>%
    select(variable, panel, group)

  d_density <- d %>%
    semi_join(density_groups, by = c("variable", "panel", "group"))

  s_constant <- s %>% filter(constant_group)
  s_variable <- s %>% filter(!constant_group)

  p <- ggplot(d, aes(value, color = group))

  if (type == "density") {
    p <- p +
      geom_density(
        data = d_density,
        aes(fill = group),
        alpha = 0.25,
        linewidth = 1.1,
        adjust = 1
      ) +
      geom_vline(
        data = s_variable,
        aes(xintercept = mean, color = group),
        linetype = "dashed",
        linewidth = 0.7
      ) +
      geom_vline(
        data = s_constant,
        aes(xintercept = mean, color = group),
        linetype = "dotdash",
        linewidth = 1.1
      ) +
      scale_fill_manual(values = condition_colors) +
      scale_y_continuous(
        expand = expansion(mult = c(0.03, 0.45))
      ) +
      labs(
        y = "Density",
        subtitle = paste(
          "Dashed lines = group means;",
          "dot-dashed lines = groups with no variation"
        )
      )
  } else {
    p <- p +
      stat_ecdf(geom = "step", linewidth = 1.1) +
      scale_y_continuous(
        labels = scales::percent,
        breaks = seq(0, 1, 0.25),
        expand = expansion(mult = c(0.02, 0.35))
      ) +
      labs(y = "Cumulative proportion")
  }

  plot_title <- if (is.null(variable_name)) {
    if (type == "density") {
      "Parasite Burden and Fulton's K Distributions"
    } else {
      "Cumulative Parasite Burden and Fulton's K Distributions"
    }
  } else {
    paste(
      if (type == "density") {
        "Distribution"
      } else {
        "Cumulative Distribution"
      },
      "—", variable_name
    )
  }

  x_label <- if (is.null(variable_name)) {
    "Parasite burden (g/g) or Fulton's K"
  } else if (variable_name == "Parasite Burden") {
    "Parasite burden (g/g)"
  } else {
    "Fulton's condition factor (K)"
  }

  p +
    geom_text(
      data = s,
      aes(
        x = Inf, y = Inf,
        label = label, color = group,
        vjust = label_position
      ),
      inherit.aes = FALSE,
      hjust = 1.02,
      size = if (is.null(variable_name)) 2.5 else 2.9,
      fontface = "bold",
      lineheight = 0.95
    ) +
    facet_wrap(
      ~ facet, ncol = 2, scales = "free",
      labeller = labeller(facet = function(x) {
        if (is.null(variable_name)) {
          x
        } else {
          sub("^.* \\| ", "", x)
        }
      })
    ) +
    scale_color_manual(values = condition_colors) +
    labs(title = plot_title, x = x_label) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "#F5F5F5"),
      panel.spacing = unit(0.7, "lines"),
      plot.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 9)
    )
}

# EIGHT-PANEL OVERVIEW FIGURES ----
p_condition_density <- plot_condition("density")
p_condition_ecdf <- plot_condition("ecdf")

print(p_condition_density)
print(p_condition_ecdf)

# FOUR-PANEL PARASITE BURDEN FIGURES ----
p_parasite_density <- plot_condition("density", "Parasite Burden")
p_parasite_ecdf <- plot_condition("ecdf", "Parasite Burden")

print(p_parasite_density)
print(p_parasite_ecdf)

# FOUR-PANEL FULTON'S K FIGURES ----
p_fulton_density <- plot_condition("density", "Fulton's K")
p_fulton_ecdf <- plot_condition("ecdf", "Fulton's K")

print(p_fulton_density)
print(p_fulton_ecdf)

# 12. GILL RAKER COUNT VS CTMAX — LAKE GRAYLING ----

# PREPARE DATA ----
raker_data <- data_final %>%
  filter(species == "Grayling", site_type == "lake") %>%
  transmute(
    ctmax = ctmax,
    raker_count = gillraker_count
  ) %>%
  filter(!is.na(ctmax), !is.na(raker_count))

# LINEAR MODEL ----
raker_lm <- lm(ctmax ~ raker_count, data = raker_data)
raker_sum <- summary(raker_lm)

p_val <- coef(raker_sum)[2, 4]
r2 <- raker_sum$r.squared
n <- nrow(raker_data)

p_text <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p_val)))
label_text <- paste0("Linear model | ", p_text, " | R² = ", sprintf("%.3f", r2), " | n = ", n)

# FIGURE ----
p_raker <- ggplot(raker_data, aes(x = raker_count, y = ctmax)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label_text,
    hjust = 1.05, vjust = 1.3,
    size = 4, fontface = "bold"
  ) +
  labs(
    title = "Lake Grayling: Gill Raker Count vs. CTmax",
    x = "Gill raker count",
    y = "CTmax (°C)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold")
  )

print(p_raker)

# MODEL RESULTS ----
cat("\n=== LAKE GRAYLING: GILL RAKER COUNT VS CTMAX ===\n")
print(summary(raker_lm))
cat("\n95% CONFIDENCE INTERVALS:\n")
print(confint(raker_lm))

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

# 14. COOLER COLOR AND CTMAX ----

# PREPARE DATA ----
cooler_data <- data_final %>%
  mutate(
    cooler_raw = tolower(trimws(as.character(testing_chamber))),
    cooler = case_when(
      grepl("^white", cooler_raw) ~ "White",
      cooler_raw == "black" ~ "Black",
      TRUE ~ NA_character_
    ),
    cooler = factor(cooler, levels = c("White", "Black")),
    species = factor(
      as.character(species),
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
    )
  ) %>%
  filter(!is.na(cooler), is.finite(ctmax))

cooler_colors <- c(
  White = "#3278A3",
  Black = "#444444"
)

# SUMMARY STATISTICS ----
cooler_summary <- cooler_data %>%
  group_by(site_type, species, group4, cooler, .drop = TRUE) %>%
  summarise(
    n = n(),
    mean = mean(ctmax),
    sd = if (n() > 1) sd(ctmax) else NA_real_,
    min = min(ctmax),
    max = max(ctmax),
    .groups = "drop"
  ) %>%
  mutate(range = max - min)

cat("\n=== CTMAX BY COOLER COLOR AND FISH GROUP ===\n")
print(cooler_summary, n = Inf, width = Inf)

cat("\n=== COOLER USE BY EXPERIMENTAL BLOCK ===\n")
print(
  cooler_data %>%
    count(site_type, species, block, cooler, .drop = TRUE),
  n = Inf
)

# FOUR-PANEL FIGURE: ALL GROUPS ----
p_cooler_all <- ggplot(
  cooler_data, aes(cooler, ctmax, color = cooler)
) +
  geom_boxplot(
    width = 0.45, outlier.shape = NA,
    fill = NA, color = "grey35"
  ) +
  geom_jitter(
    width = 0.12, height = 0,
    alpha = 0.8, size = 2
  ) +
  geom_text(
    data = cooler_summary,
    aes(x = cooler, y = Inf, label = paste0("n = ", n)),
    inherit.aes = FALSE, vjust = 1.5, size = 3.5
  ) +
  facet_wrap(~ group4, ncol = 2) +
  scale_color_manual(values = cooler_colors) +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = "CTmax by Cooler Color",
    subtitle = "White includes white coolers with mesh",
    x = "Cooler color",
    y = "CTmax (°C)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "#F5F5F5"),
    panel.spacing = grid::unit(0.8, "lines")
  )

print(p_cooler_all)

# TWO-PANEL FIGURE: LAKE FISH ----
lake_cooler <- cooler_data %>%
  filter(site_type == "lake")

p_cooler_lake <- ggplot(
  lake_cooler, aes(cooler, ctmax, color = cooler)
) +
  geom_boxplot(
    width = 0.45, outlier.shape = NA,
    fill = NA, color = "grey35"
  ) +
  geom_jitter(
    width = 0.12, height = 0,
    alpha = 0.8, size = 2.2
  ) +
  geom_text(
    data = cooler_summary %>% filter(site_type == "lake"),
    aes(
      x = cooler, y = Inf,
      label = sprintf("n=%d | Mean=%.2f", n, mean)
    ),
    inherit.aes = FALSE,
    vjust = 1.5, size = 3.2
  ) +
  facet_wrap(~ species, ncol = 2) +
  scale_color_manual(values = cooler_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15))
  ) +
  labs(
    title = "Cooler Color and CTmax — Lake Fish",
    x = "Cooler color",
    y = "CTmax (°C)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "#F5F5F5")
  )

print(p_cooler_lake)

# SPECIES-SPECIFIC WELCH TESTS ----
cooler_tests <- lake_cooler %>%
  group_by(species) %>%
  group_modify(~ {
    d <- droplevels(.x)

    if (n_distinct(d$cooler) < 2 ||
        any(table(d$cooler) < 2)) {
      return(tibble(
        difference_black_minus_white = NA_real_,
        p = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      ))
    }

    test <- t.test(ctmax ~ cooler, data = d)

    tibble(
      difference_black_minus_white =
        mean(d$ctmax[d$cooler == "Black"]) -
        mean(d$ctmax[d$cooler == "White"]),
      p = test$p.value,
      ci_lower = -test$conf.int[2],
      ci_upper = -test$conf.int[1]
    )
  }) %>%
  ungroup() %>%
  mutate(p_holm = p.adjust(p, method = "holm"))

cat("\n=== WHITE VS. BLACK: LAKE FISH ===\n")
print(cooler_tests, n = Inf, width = Inf)

# ADJUSTED MODEL: SPECIES AND MASS ----
lake_model_data <- lake_cooler %>%
  filter(is.finite(weight_g), weight_g > 0) %>%
  mutate(
    log_mass_c = log(weight_g) - mean(log(weight_g))
  )

cooler_lm <- lm(
  ctmax ~ species + log_mass_c + cooler,
  data = lake_model_data
)

cat("\n=== COOLER EFFECT ADJUSTED FOR SPECIES AND MASS ===\n")
print(summary(cooler_lm))
print(confint(cooler_lm))

# ADJUSTED MODEL: ALSO ACCOUNT FOR BLOCK ----
lake_model_data <- lake_model_data %>%
  mutate(block = factor(block))

cooler_block_lm <- lm(
  ctmax ~ species + log_mass_c + cooler + block,
  data = lake_model_data
)

cat("\n=== COOLER EFFECT ADJUSTED FOR SPECIES, MASS, AND BLOCK ===\n")
print(summary(cooler_block_lm))
print(confint(cooler_block_lm))

# COMPARE COOLER COEFFICIENTS ----
cooler_effects <- bind_rows(
  tibble(
    model = "Species + Mass",
    estimate = coef(cooler_lm)["coolerBlack"],
    p = summary(cooler_lm)$coefficients["coolerBlack", "Pr(>|t|)"]
  ),
  tibble(
    model = "Species + Mass + Block",
    estimate = coef(cooler_block_lm)["coolerBlack"],
    p = summary(cooler_block_lm)$coefficients[
      "coolerBlack", "Pr(>|t|)"
    ]
  )
)

cat("\n=== ADJUSTED COOLER EFFECT COMPARISON ===\n")
print(cooler_effects, n = Inf)

cat("\n=== COOLER USE BY BLOCK: LAKE FISH ===\n")
print(with(lake_model_data, table(block, cooler)))

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

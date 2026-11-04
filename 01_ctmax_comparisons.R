
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


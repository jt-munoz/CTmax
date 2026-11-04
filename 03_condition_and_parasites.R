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


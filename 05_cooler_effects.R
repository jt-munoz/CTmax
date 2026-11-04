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


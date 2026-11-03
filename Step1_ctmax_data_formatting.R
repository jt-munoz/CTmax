# ============================================================
# STEP 1: MERGE MONGOLIA CTMAX AND OTOLITH DATA
# Author: Joseph Munoz
# Updated: November 3, 2026
# ============================================================

# 1. PACKAGES AND FILES ----

library(tidyverse)

ctmax_file <- "CTMax_Data_Eg_Uur2"
lake_file <- "Otolith Data Sheet_2026_Lake_Sophia.xlsx - Otoliths_Data_Sheet_2026_Lake.csv"
river_file <- "Otoliths_Data_Sheet_2026_Rivers - Otoliths_Data_Sheet_2026_Rivers.csv"

ctmax_raw <- read_csv(ctmax_file, show_col_types = FALSE)
lake_raw <- read_csv(lake_file, show_col_types = FALSE)
river_raw <- read_csv(river_file, show_col_types = FALSE)

# 2. HELPER FUNCTIONS ----

to_number <- function(x) suppressWarnings(as.numeric(trimws(as.character(x))))

clean_id <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_
  x
}

check_columns <- function(data, required, name) {
  missing <- setdiff(required, names(data))
  
  if (length(missing) > 0) {
    stop(name, " missing columns: ", paste(missing, collapse = ", "))
  } else {
    message(name, ": All required columns found.")
  }
}

# 3. PREPARE CTMAX DATA ----

# Check required columns
check_columns(
  ctmax_raw,
  c("trial", "complete_id", "species", "site", "site_type",
    "ctmax", "start_temp"),
  "CTmax"
)

# Exclude invalid trials
excluded_trials <- c(1, 26, 67, 84, 108)

# Clean and prepare CTmax data
ctmax_data <- ctmax_raw %>%
  mutate(
    trial = to_number(trial),
    fish_id = clean_id(complete_id),
    site_type = str_to_lower(str_trim(as.character(site_type))),
    ctmax = to_number(ctmax),
    start_temp = to_number(start_temp),
    species_group = case_when(
      str_detect(str_to_lower(as.character(species)), "grayling") ~ "Grayling",
      str_detect(str_to_lower(as.character(species)), "lenok") ~ "Lenok",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!trial %in% excluded_trials)

# Check for missing or invalid essential values
if (anyNA(ctmax_data$trial)) stop("Missing trial numbers.")
if (anyNA(ctmax_data$fish_id)) stop("Missing CTmax fish IDs.")
if (anyNA(ctmax_data$ctmax)) stop("Missing or invalid CTmax values.")
if (anyNA(ctmax_data$start_temp)) stop("Missing or invalid starting temperatures.")
#missing start temps for trials 105-107
if (anyNA(ctmax_data$species_group)) stop("Unrecognized species.")
if (anyNA(ctmax_data$site_type) ||
    any(!ctmax_data$site_type %in% c("lake", "river"))) {
  stop("Unexpected habitat values.")
}

# Check for duplicate fish IDs
duplicate_ctmax <- ctmax_data %>%
  count(site_type, fish_id) %>%
  filter(n > 1)

if (nrow(duplicate_ctmax) > 0) {
  print(duplicate_ctmax)
  stop("Duplicate CTmax fish IDs within habitat.")
} else {
  message("No duplicate CTmax fish IDs found.")
}

# 4. PREPARE OTOLITH DATA ----

otolith_required <- c(
  "sample_id", "Q (g)", "tl_mm",
  "Parasite load/in grams (cysts for Lenok)",
  "Nematode parasite in g for Lenok"
)

check_columns(lake_raw, otolith_required, "Lake otolith")
check_columns(river_raw, c(otolith_required, "sample_id_changed"), "River otolith")

prepare_otolith <- function(data, habitat) {
  data %>%
    mutate(across(everything(), as.character)) %>%
    mutate(
      site_type = habitat,
      otolith_row = row_number(),
      sample_id_original = clean_id(sample_id),
      sample_id_changed_clean = if (habitat == "river") clean_id(sample_id_changed) else NA_character_,
      weight_g = to_number(.data[["Q (g)"]]),
      length_mm = to_number(tl_mm),
      parasite_mass_recorded_g = to_number(.data[["Parasite load/in grams (cysts for Lenok)"]]),
      nematode_mass_g = to_number(.data[["Nematode parasite in g for Lenok"]])
    ) %>%
    mutate(
      parasite_mass_g = if_else(
        is.na(parasite_mass_recorded_g) & is.na(nematode_mass_g),
        NA_real_,
        coalesce(parasite_mass_recorded_g, 0) + coalesce(nematode_mass_g, 0)
      ),
      gill_parasite_count = case_when(
        site_type == "lake" & sample_id_original == "BL022" ~ 4,
        site_type == "lake" & sample_id_original == "BL023" ~ 8,
        site_type == "lake" & sample_id_original == "BL024" ~ 3,
        site_type == "river" & sample_id_changed_clean == "CTM_BG008" ~ 10,
        site_type == "river" & sample_id_changed_clean == "CTM_BG014" ~ 3,
        TRUE ~ 0
      )
    )
}

otolith_data <- bind_rows(
  prepare_otolith(lake_raw, "lake"),
  prepare_otolith(river_raw, "river")
) %>%
  mutate(otolith_key = paste(site_type, otolith_row, sep = "_"))

# CHECK FOR MISSING FISH IDS ----

if (anyNA(otolith_data$sample_id_original)) {
  stop("Missing original otolith fish IDs.")
} else {
  message("All otolith fish have original IDs.")
}

# 5. MATCH AND MERGE RECORDS ----

# CREATE MATCHING KEYS ----

changed_keys <- otolith_data %>%
  filter(!is.na(sample_id_changed_clean)) %>%
  transmute(site_type, fish_id = sample_id_changed_clean, otolith_key,
            match_priority = 1L, match_source = "sample_id_changed")

original_keys <- otolith_data %>%
  transmute(site_type, fish_id = sample_id_original, otolith_key,
            match_priority = 2L, match_source = "sample_id")

# MATCH CTMAX FISH TO OTOLITH RECORDS ----

matches <- bind_rows(changed_keys, original_keys) %>%
  semi_join(ctmax_data, by = c("site_type", "fish_id")) %>%
  arrange(match_priority) %>%
  distinct(site_type, fish_id, otolith_key, .keep_all = TRUE) %>%
  group_by(site_type, fish_id) %>%
  filter(match_priority == min(match_priority)) %>%
  ungroup()

# CHECK FOR AMBIGUOUS OR REUSED MATCHES ----

if (anyDuplicated(matches[c("site_type", "fish_id")])) {
  stop("Multiple otolith records match the same CTmax fish.")
}

if (anyDuplicated(matches$otolith_key)) {
  stop("An otolith record matches multiple CTmax fish.")
}

# MERGE DATASETS ----

data_merged <- ctmax_data %>%
  left_join(matches, by = c("site_type", "fish_id")) %>%
  left_join(otolith_data %>% select(-site_type), by = "otolith_key") %>%
  mutate(otolith_matched = !is.na(otolith_key))

# VERIFY MERGE ----

if (nrow(data_merged) != nrow(ctmax_data)) {
  stop("Merge changed the number of CTmax records.")
}

if (any(!data_merged$otolith_matched)) {
  print(data_merged %>%
          filter(!otolith_matched) %>%
          select(trial, fish_id, site_type))
  stop("Some CTmax fish have no matching otolith record.")
}

message("Successfully matched all ", nrow(data_merged), " CTmax fish.")

# VERIFY AND STANDARDIZE SPECIES ----

standardize_species <- function(x) {
  case_when(
    grepl("grayling", x, ignore.case = TRUE) ~ "Grayling",
    grepl("lenok", x, ignore.case = TRUE) ~ "Lenok",
    TRUE ~ NA_character_
  )
}

data_merged <- data_merged %>%
  mutate(
    species = standardize_species(species.y),
    species_ctmax = standardize_species(species.x)
  )
if (anyNA(data_merged$species) ||
    anyNA(data_merged$species_ctmax) ||
    any(data_merged$species != data_merged$species_ctmax)) {
  stop("Species mismatch between CTmax and otolith records.")
}

data_merged <- data_merged %>%
  select(-species_ctmax, -species.x, -species.y)

message("Species verified for all ", nrow(data_merged), " fish.")

# 6. CALCULATE AND VALIDATE VARIABLES ----

# Parasite load = recorded parasite mass + nematode mass (g).
# Parasite burden = parasite load / fish mass (g).
# Fulton K = 100000 * mass (g) / total length (mm)^3.

known_positive <- c(
  "CTM_BL003", "CTM_G019", "CTM_BL028",
  "CTM_BG007", "CTM_BG008", "CTM_BL031", "CTM_BG014"
)

data_merged <- data_merged %>%
  mutate(
    parasite_unknown = fish_id %in% known_positive,
    parasite_mass_recorded_g = if_else(
      otolith_matched & !parasite_unknown,
      coalesce(parasite_mass_recorded_g, 0), NA_real_
    ),
    nematode_mass_g = if_else(
      otolith_matched & !parasite_unknown,
      coalesce(nematode_mass_g, 0), NA_real_
    ),
    parasite_load = parasite_mass_recorded_g + nematode_mass_g,
    parasite_burden = if_else(
      weight_g > 0, parasite_load / weight_g, NA_real_
    ),
    parasite_present = case_when(
      parasite_unknown ~ 1L,
      !is.na(parasite_load) ~ as.integer(parasite_load > 0),
      TRUE ~ NA_integer_
    ),
    fultons_k = if_else(
      weight_g > 0 & length_mm > 0,
      100000 * weight_g / length_mm^3, NA_real_
    )
  )

if (any(data_merged$parasite_load < 0, na.rm = TRUE)) {
  stop("Negative parasite mass detected.")
}

if (any(data_merged$weight_g <= 0 | data_merged$length_mm <= 0,
        na.rm = TRUE)) {
  stop("Nonpositive fish mass or length detected.")
}

if (any(!data_merged$otolith_matched & !is.na(data_merged$parasite_load))) {
  stop("Unmatched fish incorrectly assigned parasite values.")
}

if (any(data_merged$parasite_unknown &
        !is.na(data_merged$parasite_burden), na.rm = TRUE)) {
  stop("Unknown parasite mass incorrectly assigned burden.")
}

message("Variables calculated for ", nrow(data_merged), " fish.")

# 7. EXPORT DATA AND SUMMARIES ----

# CREATE MERGE AUDIT ----

merge_audit <- data_merged %>%
  select(
    trial, fish_id, species, site_type, otolith_key,
    match_source, sample_id_original, sample_id_changed_clean,
    weight_g, length_mm, parasite_unknown,
    parasite_mass_recorded_g, nematode_mass_g,
    parasite_load, parasite_burden, parasite_present
  )

# SUMMARIZE BY SPECIES AND HABITAT ----

group_summary <- data_merged %>%
  group_by(species, site_type) %>%
  summarise(
    n_fish = n(),
    n_parasite_burden = sum(!is.na(parasite_burden)),
    n_parasite_missing = sum(is.na(parasite_burden)),
    n_parasite_present = sum(parasite_present == 1, na.rm = TRUE),
    mean_ctmax = mean(ctmax, na.rm = TRUE),
    .groups = "drop"
  )

# PREPARE FINAL DATASET ----

data_final <- data_merged %>%
  select(where(~ !all(is.na(.)))) %>%
  select(-any_of(c(
    "otolith_row", "otolith_key", "match_priority", "match_source",
    "sample_id_original", "sample_id_changed_clean",
    "species.x", "species.y"
  )))

# EXPORT FILES ----

write_csv(data_final, "CTmax_Merged_Data.csv", na = "")
write_csv(merge_audit, "CTmax_Merge_Audit.csv", na = "")
write_csv(group_summary, "CTmax_Data_Summary.csv", na = "")

cat("\nSTEP 1 COMPLETE:", nrow(data_final), "fish\n")
cat("Final variables:", ncol(data_final), "\n")
print(group_summary)

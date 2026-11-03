
# ============================================================
# STEP 1: MERGE MONGOLIA CTmax AND OTOLITH DATA
#
# Corrected parasite calculations
# Original biological units
#
# OUTPUT:
#   CTmax_Merged_Data.csv
#   CTmax_Merge_Audit.csv
#   CTmax_Merge_Summary.csv
#   CTmax_Missing_Data_Summary.csv
#   CTmax_Otolith_Match_Diagnostics.csv
# ============================================================


# 1. PACKAGES ------------------------------------------------

if (!requireNamespace("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

library(tidyverse)


# 2. FILE NAMES ----------------------------------------------

ctmax_file <- "CTMax_Data_Eg_Uur2.csv"

lake_file <- paste0(
  "Otolith Data Sheet_2026_Lake_Sophia.xlsx",
  " - Otoliths_Data_Sheet_2026_Lake.csv"
)

river_file <- paste0(
  "Otoliths_Data_Sheet_2026_Rivers",
  " - Otoliths_Data_Sheet_2026_Rivers.csv"
)

# Accept the extensionless CTmax file if necessary.
if (!file.exists(ctmax_file)) {
  if (file.exists("CTMax_Data_Eg_Uur2")) {
    ctmax_file <- "CTMax_Data_Eg_Uur2"
  } else {
    stop("Cannot find the CTmax input file.")
  }
}

for (file_path in c(lake_file, river_file)) {
  if (!file.exists(file_path)) {
    stop("Cannot find input file: ", file_path)
  }
}


# 3. READ DATA -----------------------------------------------

ctmax_raw <- readr::read_csv(
  ctmax_file,
  show_col_types = FALSE,
  name_repair = "unique"
)

lake_raw <- readr::read_csv(
  lake_file,
  show_col_types = FALSE,
  name_repair = "unique"
)

river_raw <- readr::read_csv(
  river_file,
  show_col_types = FALSE,
  name_repair = "unique"
)


# 4. HELPER FUNCTIONS ----------------------------------------

# Convert numeric-looking text to numeric.
# Blank cells become NA.
to_number <- function(x) {

  x <- trimws(as.character(x))

  x[x == ""] <- NA_character_

  suppressWarnings(as.numeric(x))
}


# Standardize fish IDs for matching.
clean_id <- function(x) {

  x <- toupper(trimws(as.character(x)))

  x[x %in% c("", "NA", "N/A", "NULL")] <- NA_character_

  x
}


# Require expected columns.
check_columns <- function(data, required, dataset_name) {

  missing <- setdiff(required, names(data))

  if (length(missing) > 0) {
    stop(
      dataset_name,
      " is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}


# Return a column or NA if the column is absent.
optional_column <- function(data, column_name) {

  if (column_name %in% names(data)) {
    return(data[[column_name]])
  }

  rep(NA_character_, nrow(data))
}


# 5. VALIDATE CTmax DATA -------------------------------------

check_columns(
  ctmax_raw,
  c(
    "trial",
    "complete_id",
    "species",
    "site_type",
    "ctmax",
    "start_temp"
  ),
  "CTmax dataset"
)


# 6. PREPARE CTmax DATA --------------------------------------

ctmax_data <- ctmax_raw %>%
  mutate(
    ctmax_row_id = row_number(),

    trial = to_number(trial),

    fish_id = clean_id(complete_id),

    site_type = tolower(
      trimws(as.character(site_type))
    ),

    ctmax = to_number(ctmax),

    start_temp = to_number(start_temp)
  ) %>%
  filter(
    !(trial %in% c(1, 26, 67, 84, 108))
  )

if (any(is.na(ctmax_data$fish_id))) {
  stop(
    "Some CTmax fish have missing IDs. ",
    "Inspect complete_id before merging."
  )
}

if (any(!ctmax_data$site_type %in% c("lake", "river"))) {
  stop(
    "Unexpected site_type values in CTmax data. ",
    "Expected lake or river."
  )
}

# Verify CTmax fish IDs are unique within habitat.
duplicate_ctmax <- ctmax_data %>%
  count(site_type, fish_id) %>%
  filter(n > 1)

if (nrow(duplicate_ctmax) > 0) {

  print(duplicate_ctmax, n = Inf)

  stop(
    "Duplicate CTmax fish IDs detected within habitat. ",
    "Resolve before merging."
  )
}


# 7. VALIDATE OTOLITH DATA -----------------------------------

otolith_required <- c(
  "sample_id",
  "Q (g)",
  "tl_mm",
  "Parasite load/in grams (cysts for Lenok)",
  "Nematode parasite in g for Lenok"
)

check_columns(
  lake_raw,
  otolith_required,
  "Lake otolith dataset"
)

check_columns(
  river_raw,
  otolith_required,
  "River otolith dataset"
)


# 8. PREPARE OTOLITH DATA ------------------------------------

prepare_otolith <- function(data, habitat_name) {

  data %>%
    mutate(
      otolith_row = row_number(),

      site_type = habitat_name,

      sample_id_original = clean_id(sample_id),

      sample_id_changed_clean = clean_id(
        optional_column(data, "sample_id_changed")
      ),

      weight_g = to_number(.data[["Q (g)"]]),

      length_mm = to_number(.data[["tl_mm"]]),

      cyst_mass_g = to_number(
        .data[["Parasite load/in grams (cysts for Lenok)"]]
      ),

      nematode_mass_g = to_number(
        .data[["Nematode parasite in g for Lenok"]]
      )
    ) %>%
    select(
      site_type,
      otolith_row,
      sample_id_original,
      sample_id_changed_clean,
      weight_g,
      length_mm,
      cyst_mass_g,
      nematode_mass_g
    )
}

lake_otolith <- prepare_otolith(
  lake_raw,
  "lake"
)

river_otolith <- prepare_otolith(
  river_raw,
  "river"
)

otolith_data <- bind_rows(
  lake_otolith,
  river_otolith
) %>%
  mutate(
    otolith_key = paste(
      site_type,
      otolith_row,
      sep = "_"
    )
  )


# 9. SPECIAL LAKE RECORDS ------------------------------------
#
# These six lake otolith IDs previously contained
# nonnumeric parasite measurements.
#
# Do not automatically discard them if they correspond
# to CTmax fish.
#
# Their original text values have been converted to NA.
# If matched, blank parasite fields will be zero.

special_lake_ids <- c(
  "G001", "G002", "G003",
  "BL001", "BL002", "BL003"
)

special_lake_records <- lake_otolith %>%
  filter(
    sample_id_original %in% special_lake_ids |
      sample_id_changed_clean %in% special_lake_ids
  )

cat("\nSpecial lake otolith records:\n")

print(
  special_lake_records,
  n = Inf,
  width = Inf
)


# 10. BUILD OTOLITH MATCH KEYS -------------------------------
#
# Priority:
#   1. sample_id_changed
#   2. sample_id
#
# Restrict candidate matches to CTmax fish IDs.
# Avoid treating unrelated duplicate otolith IDs as matches.

ctmax_keys <- ctmax_data %>%
  distinct(site_type, fish_id)

changed_keys <- otolith_data %>%
  filter(!is.na(sample_id_changed_clean)) %>%
  transmute(
    site_type,
    fish_id = sample_id_changed_clean,
    otolith_key,
    match_priority = 1L,
    match_source = "sample_id_changed"
  )

original_keys <- otolith_data %>%
  filter(!is.na(sample_id_original)) %>%
  transmute(
    site_type,
    fish_id = sample_id_original,
    otolith_key,
    match_priority = 2L,
    match_source = "sample_id"
  )

candidate_matches <- bind_rows(
  changed_keys,
  original_keys
) %>%
  semi_join(
    ctmax_keys,
    by = c("site_type", "fish_id")
  ) %>%
  distinct(
    site_type,
    fish_id,
    otolith_key,
    .keep_all = TRUE
  )


# 11. CHECK AMBIGUOUS MATCHES --------------------------------
#
# A CTmax fish cannot have multiple equally preferred
# otolith matches.
#
# An otolith record cannot be assigned to multiple fish.

best_priority <- candidate_matches %>%
  group_by(site_type, fish_id) %>%
  filter(match_priority == min(match_priority)) %>%
  ungroup()

ambiguous_fish <- best_priority %>%
  count(site_type, fish_id, name = "n_matches") %>%
  filter(n_matches > 1)

if (nrow(ambiguous_fish) > 0) {

  print(ambiguous_fish, n = Inf)

  stop(
    "Ambiguous otolith matches detected. ",
    "Resolve the fish IDs before continuing."
  )
}

reused_otolith <- best_priority %>%
  count(otolith_key, name = "n_fish") %>%
  filter(n_fish > 1)

if (nrow(reused_otolith) > 0) {

  print(reused_otolith, n = Inf)

  stop(
    "An otolith record matches multiple CTmax fish. ",
    "Resolve before continuing."
  )
}


# 12. CREATE MATCH TABLE -------------------------------------

match_table <- best_priority %>%
  select(
    site_type,
    fish_id,
    otolith_key,
    match_priority,
    match_source
  )

otolith_measurements <- otolith_data %>%
  select(
    otolith_key,
    otolith_row,
    sample_id_original,
    sample_id_changed_clean,
    weight_g,
    length_mm,
    cyst_mass_g,
    nematode_mass_g
  )


# 13. MERGE CTmax AND OTOLITH DATA ----------------------------

data_merged <- ctmax_data %>%
  left_join(
    match_table,
    by = c("site_type", "fish_id")
  ) %>%
  left_join(
    otolith_measurements,
    by = "otolith_key"
  ) %>%
  mutate(
    otolith_matched = !is.na(otolith_key)
  )

if (nrow(data_merged) != nrow(ctmax_data)) {
  stop(
    "The merge changed the number of CTmax records. ",
    "Check match keys."
  )
}


# 14. CORRECT PARASITE CALCULATIONS ---------------------------
#
# CRITICAL RULE:
#
# Matched otolith fish:
#   blank cyst mass     -> 0
#   blank nematode mass -> 0
#
# Unmatched CTmax fish:
#   parasite measurements remain NA
#
# Never convert unknown parasite status to zero.

data_merged <- data_merged %>%
  mutate(

    cyst_mass_g = if_else(
      otolith_matched,
      coalesce(cyst_mass_g, 0),
      NA_real_
    ),

    nematode_mass_g = if_else(
      otolith_matched,
      coalesce(nematode_mass_g, 0),
      NA_real_
    ),

    parasite_load = if_else(
      otolith_matched,
      cyst_mass_g + nematode_mass_g,
      NA_real_
    ),

    parasite_burden = if_else(
      otolith_matched &
        !is.na(weight_g) &
        weight_g > 0,
      parasite_load / weight_g,
      NA_real_
    ),

    parasite_burden_percent =
      parasite_burden * 100,

    parasite_present = if_else(
      otolith_matched,
      as.integer(parasite_load > 0),
      NA_integer_
    ),

    log_parasite_burden = if_else(
      !is.na(parasite_burden) &
        parasite_burden >= 0,
      log1p(parasite_burden),
      NA_real_
    )
  )


# 15. FULTON'S CONDITION FACTOR ------------------------------
#
# Fulton K using mass in grams and length in mm.
#
# No predicted_weight_g or relative_condition columns.

data_merged <- data_merged %>%
  mutate(
    fultons_k = if_else(
      !is.na(weight_g) &
        weight_g > 0 &
        !is.na(length_mm) &
        length_mm > 0,
      100000 * weight_g / length_mm^3,
      NA_real_
    )
  )


# 16. SPECIES-SPECIFIC CONDITION RESIDUAL --------------------
#
# Residual from:
#   log(mass) ~ log(length)
#
# This is a separate derived condition index.
# It does not transform the weight_g variable
# used by your Step 2 candidate models.

data_merged <- data_merged %>%
  mutate(
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

    condition_residual = NA_real_
  )

for (current_species in c("Grayling", "Lenok")) {

  species_rows <- which(
    data_merged$species_group == current_species &
      is.finite(data_merged$weight_g) &
      data_merged$weight_g > 0 &
      is.finite(data_merged$length_mm) &
      data_merged$length_mm > 0
  )

  if (length(species_rows) < 3) {
    next
  }

  species_data <- data_merged[species_rows, ]

  condition_fit <- lm(
    log(weight_g) ~ log(length_mm),
    data = species_data
  )

  data_merged$condition_residual[species_rows] <-
    as.numeric(residuals(condition_fit))
}


# 17. VALIDATE PARASITE VALUES -------------------------------

invalid_parasites <- data_merged %>%
  filter(
    (!is.na(cyst_mass_g) & cyst_mass_g < 0) |
      (!is.na(nematode_mass_g) & nematode_mass_g < 0) |
      (!is.na(parasite_load) & parasite_load < 0)
  )

if (nrow(invalid_parasites) > 0) {

  print(invalid_parasites, n = Inf)

  stop(
    "Negative parasite measurements detected. ",
    "Inspect original otolith data."
  )
}

# Check that unmatched fish retain unknown parasite status.

unmatched_check <- data_merged %>%
  filter(
    !otolith_matched &
      (
        !is.na(parasite_load) |
          !is.na(parasite_burden) |
          !is.na(parasite_burden_percent)
      )
  )

if (nrow(unmatched_check) > 0) {
  stop(
    "Unmatched fish were incorrectly assigned ",
    "parasite measurements."
  )
}

# Check that matched fish have calculated parasite load.

matched_missing <- data_merged %>%
  filter(
    otolith_matched &
      is.na(parasite_load)
  )

if (nrow(matched_missing) > 0) {
  stop(
    "Matched otolith fish have missing parasite load."
  )
}


# 18. MERGE AUDIT --------------------------------------------

merge_audit <- data_merged %>%
  select(
    ctmax_row_id,
    trial,
    fish_id,
    species,
    site_type,
    otolith_matched,
    match_source,
    otolith_key,
    sample_id_original,
    sample_id_changed_clean,
    weight_g,
    length_mm,
    cyst_mass_g,
    nematode_mass_g,
    parasite_load,
    parasite_burden,
    parasite_burden_percent
  )

write_csv(
  merge_audit,
  "CTmax_Merge_Audit.csv"
)


# 19. MERGE SUMMARY ------------------------------------------

merge_summary <- data_merged %>%
  group_by(species_group, site_type) %>%
  summarise(
    ctmax_fish = n(),

    matched_otolith = sum(otolith_matched),

    unmatched_otolith = sum(!otolith_matched),

    with_weight = sum(!is.na(weight_g)),

    with_length = sum(!is.na(length_mm)),

    with_parasite_burden =
      sum(!is.na(parasite_burden_percent)),

    parasite_positive = sum(
      parasite_load > 0,
      na.rm = TRUE
    ),

    parasite_zero = sum(
      parasite_load == 0,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

print(
  merge_summary,
  n = Inf,
  width = Inf
)

write_csv(
  merge_summary,
  "CTmax_Merge_Summary.csv"
)


# 20. MISSING DATA SUMMARY -----------------------------------

important_variables <- c(
  "ctmax",
  "start_temp",
  "weight_g",
  "length_mm",
  "cyst_mass_g",
  "nematode_mass_g",
  "parasite_load",
  "parasite_burden",
  "parasite_burden_percent",
  "parasite_present",
  "fultons_k",
  "condition_residual"
)

important_variables <- intersect(
  important_variables,
  names(data_merged)
)

missing_summary <- tibble(
  variable = important_variables,

  missing_n = vapply(
    important_variables,
    function(variable_name) {
      sum(is.na(data_merged[[variable_name]]))
    },
    numeric(1)
  )
) %>%
  mutate(
    total_n = nrow(data_merged),
    missing_percent =
      100 * missing_n / total_n
  )

print(
  missing_summary,
  n = Inf
)

write_csv(
  missing_summary,
  "CTmax_Missing_Data_Summary.csv"
)


# 21. MATCH DIAGNOSTICS --------------------------------------

match_diagnostics <- data_merged %>%
  count(
    site_type,
    species_group,
    otolith_matched,
    match_source,
    name = "n"
  )

write_csv(
  match_diagnostics,
  "CTmax_Otolith_Match_Diagnostics.csv"
)


# 22. SAVE FINAL MERGED DATA ---------------------------------

# Remove temporary matching identifiers from the
# main analytical dataset only.
# Retain them in the audit file.

data_final <- data_merged %>%
  select(
    -any_of(c(
      "ctmax_row_id",
      "otolith_row",
      "otolith_key",
      "match_priority",
      "match_source",
      "sample_id_original",
      "sample_id_changed_clean"
    ))
  )

write_csv(
  data_final,
  "CTmax_Merged_Data.csv",
  na = ""
)


# 23. FINAL VERIFICATION -------------------------------------

cat("\n========================================\n")
cat("STEP 1 COMPLETE\n")
cat("========================================\n")

cat(
  "\nCTmax fish retained:",
  nrow(data_final)
)

cat(
  "\nMatched otolith fish:",
  sum(data_final$otolith_matched)
)

cat(
  "\nUnmatched CTmax fish:",
  sum(!data_final$otolith_matched)
)

cat(
  "\nFish with calculated parasite burden:",
  sum(!is.na(data_final$parasite_burden_percent))
)

cat(
  "\nFish with positive parasite mass:",
  sum(data_final$parasite_load > 0, na.rm = TRUE)
)

cat(
  "\nFish with zero parasite mass:",
  sum(data_final$parasite_load == 0, na.rm = TRUE)
)

cat(
  "\n\nMain output: CTmax_Merged_Data.csv\n"
)

cat("\nReady for Step 2 model selection.\n")

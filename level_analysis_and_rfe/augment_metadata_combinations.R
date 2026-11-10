# ==============================================================================
# Metadata-combination data augmentation for archaeological ML workflows
#
# Strategy:
#   1. Define named buckets of metadata columns
#   2. Generate all valid combos (1-3 buckets, exactly 1 column per bucket)
#   3. Split data by target_col (class), augment each class independently,
#      then merge -- prevents mean() blending features across classes
#   4. For each combo within a class: group_by those columns, summarise
#        - numeric features          -> mean() of the group
#        - integer cols (e.g. mni)   -> mean(), then rounded to integer
#        - non-numeric, non-grouping -> NA
#        - grouping columns          -> kept (constant within group)
#        - absent cols (e.g. mni)    -> padded as NA, then imputed to 1L
#   5. Drop any synthetic row that is an exact copy of an original row
#   6. bind_rows(original, all_synthetic) -- same column width throughout
# ==============================================================================

library(dplyr)
library(purrr)
library(stringr)
library(tibble)


# -- 0. Bucket definitions -----------------------------------------------------

buckets <- list(
  temporal     = c("geological_period", "cultural_period", "culture_period_2",
                   "start_date_cal_bp", "end_date_cal_bp", "mis_stage"),
  spatial      = c("region", "country"),
  site_context = c("karstic_system", "open_air_site"),
  biological   = c("age_category", "gender"),
  taxonomy     = c("hominin_species")
)


# -- 1. Generate all valid combinations ----------------------------------------
#
# Valid combo = choose 1-3 buckets (never repeat a bucket),
# then pick exactly one column from each chosen bucket.
#
# Approximate counts:
#   1-bucket combos : 6+2+2+2+1 = 13
#   2-bucket combos : all pairs  = 78
#   3-bucket combos : all triples= 174
#   Total                       ~ 265

generate_combos <- function(buckets,
                             min_buckets = 1L,
                             max_buckets = 3L) {
  bucket_names <- names(buckets)

  bucket_subsets <- unlist(
    lapply(seq(min_buckets, max_buckets), function(k)
      combn(bucket_names, k, simplify = FALSE)
    ),
    recursive = FALSE
  )

  all_combos <- lapply(bucket_subsets, function(bnames) {
    col_lists <- buckets[bnames]
    grid      <- do.call(expand.grid, c(col_lists, stringsAsFactors = FALSE))
    lapply(seq_len(nrow(grid)), function(i)
      as.list(grid[i, , drop = FALSE])
    )
  })

  unlist(all_combos, recursive = FALSE)
}


# -- 2. Augment for one specific combination -----------------------------------
#
# @param df             Tibble for ONE class only
# @param combo          Named list, one column name per bucket
# @param aug_id         Character tag for traceability
# @param feature_cols   Columns used as ML features (subset of names(df))
# @param integer_cols   Cols to round after mean() e.g. c("mni")
# @param required_cols  Cols imputed to 1L when NA after padding (e.g. c("mni"))
# @param min_group_size Groups smaller than this are skipped (default 2)
#
# @return Tibble with same columns as df + .aug_id, or empty tibble

augment_with_combination <- function(df,
                                     combo,
                                     aug_id,
                                     feature_cols,
                                     integer_cols   = character(0),
                                     required_cols  = character(0),
                                     min_group_size = 2L) {

  combo_cols <- unlist(combo, use.names = FALSE)

  # Skip silently if any combo column is absent from df
  if (!all(combo_cols %in% names(df))) {
    return(tibble::tibble())
  }

  # Drop rows where any grouping column is NA
  df_valid <- df |>
    dplyr::filter(dplyr::if_all(dplyr::all_of(combo_cols), ~ !is.na(.x)))

  if (nrow(df_valid) == 0) return(tibble::tibble())

  # Classify feature columns (only those present in df)
  feature_cols_present     <- intersect(feature_cols, names(df))
  numeric_features         <- feature_cols_present[
    vapply(df[feature_cols_present], is.numeric, logical(1))
  ]
  non_numeric_features     <- setdiff(feature_cols_present, numeric_features)
  non_numeric_non_grouping <- setdiff(non_numeric_features, combo_cols)

  # Core: group -> summarise
  #   grouping cols            : kept by summarise automatically
  #   numeric feature cols     : mean()
  #   non-numeric feature cols : NA
  new_rows <- df_valid |>
    dplyr::group_by(dplyr::across(dplyr::all_of(combo_cols))) |>
    dplyr::filter(dplyr::n() >= min_group_size) |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(numeric_features),         ~ mean(.x, na.rm = TRUE)),
      dplyr::across(dplyr::all_of(non_numeric_non_grouping), ~ NA),
      .groups = "drop"
    )

  if (nrow(new_rows) == 0) return(tibble::tibble())

  # Round integer cols that are already in new_rows (i.e. were in feature_cols)
  integer_in_features <- intersect(integer_cols, names(new_rows))
  if (length(integer_in_features) > 0) {
    new_rows <- new_rows |>
      dplyr::mutate(
        dplyr::across(dplyr::all_of(integer_in_features), ~ as.integer(round(.x)))
      )
  }

  new_rows <- new_rows |> dplyr::mutate(.aug_id = aug_id)

  # Pad ALL columns present in df but absent from new_rows -> NA
  # (e.g. mni when it is not in feature_cols arrives here as NA)
  absent_cols <- setdiff(names(df), names(new_rows))
  for (col in absent_cols) {
    new_rows[[col]] <- NA
  }

  # Reorder columns to exactly match df
  new_rows <- new_rows |>
    dplyr::select(dplyr::all_of(names(df)), .aug_id) |>
    tibble::as_tibble()

  # Impute required cols: NA/NaN -> 1L
  # Runs AFTER padding so cols like mni (not in feature_cols) now exist
  required_present <- intersect(required_cols, names(new_rows))
  if (length(required_present) > 0) {
    new_rows <- new_rows |>
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(required_present),
          ~ dplyr::if_else(is.na(.x) | is.nan(.x), 1L, as.integer(.x))
        )
      )
  }

  # Drop exact copies of original rows
  new_rows |>
    dplyr::anti_join(
      dplyr::select(df, dplyr::all_of(feature_cols_present)),
      by = feature_cols_present
    )
}


# -- 3. Apply across all combinations & bind with original ---------------------

build_augmented_dataset <- function(df,
                                    combos,
                                    feature_cols,
                                    integer_cols   = character(0),
                                    required_cols  = character(0),
                                    min_group_size = 2L) {

  original_tagged <- df |> dplyr::mutate(.aug_id = "original")

  # Pre-filter: only run combos whose columns all exist in df
  valid_combos <- Filter(function(combo) {
    all(unlist(combo, use.names = FALSE) %in% names(df))
  }, combos)

  n_dropped <- length(combos) - length(valid_combos)
  if (n_dropped > 0) {
    message(
      n_dropped, " combo(s) skipped -- columns not found in df.\n",
      "  Bucket cols present : ",
      paste(intersect(unlist(buckets), names(df)), collapse = ", "), "\n",
      "  Bucket cols missing : ",
      paste(setdiff(unlist(buckets),   names(df)), collapse = ", ")
    )
  }

  if (length(valid_combos) == 0) {
    message("No valid combos -- returning original data unchanged.")
    return(original_tagged |> tibble::as_tibble())
  }

  synthetic_rows <- purrr::imap_dfr(valid_combos, function(combo, idx) {
    aug_id <- paste0(
      "aug_", stringr::str_pad(idx, width = 4, pad = "0"),
      "_", paste(unlist(combo), collapse = "+")
    )
    augment_with_combination(
      df             = df,
      combo          = combo,
      aug_id         = aug_id,
      feature_cols   = feature_cols,
      integer_cols   = integer_cols,
      required_cols  = required_cols,
      min_group_size = min_group_size
    )
  })

  dplyr::bind_rows(original_tagged, synthetic_rows) |>
    tibble::as_tibble()
}


# -- 4. Factory: returns the augment_fun closure for step_meta_smote ----------
#
# KEY: augmentation is applied per class (target_col) separately, then merged.
# This mirrors SMOTE behaviour -- means are computed within a single class so
# synthetic rows never blend feature distributions across class boundaries.
#
# @param target_col    The classification target column (e.g. "taphonomic_context")
#
# Returned function contract:
#   input  -- full dataframe including bucket cols, mni, target_col
#   output -- bind_rows(original, synthetic), .aug_id dropped, same columns

make_augment_fun <- function(combos,
                              feature_cols,
                              target_col,
                              integer_cols   = character(0),
                              required_cols  = character(0),
                              min_group_size = 2L) {
  function(x) {

    class_levels <- unique(x[[target_col]])

    per_class <- purrr::map(class_levels, function(lvl) {
      df_class <- dplyr::filter(x, .data[[target_col]] == lvl)

      build_augmented_dataset(
        df             = df_class,
        combos         = combos,
        feature_cols   = feature_cols,
        integer_cols   = integer_cols,
        required_cols  = required_cols,
        min_group_size = min_group_size
      )
    })

    dplyr::bind_rows(per_class) |>
      dplyr::select(-.aug_id)
  }
}


# ==============================================================================
# USAGE
# ==============================================================================

# combos <- generate_combos(buckets)
#
# my_augment_fun <- make_augment_fun(
#   combos         = combos,
#   feature_cols   = feature_cols,   # does NOT need to include mni
#   target_col     = target_col,     # split by class before augmenting
#   integer_cols   = c("mni"),       # rounded after mean()
#   required_cols  = c("mni"),       # imputed to 1L after padding if NA
#   min_group_size = 2L
# )
#
# recipe_obj <- build_recipe_for_task(df = model_df, ...) |>
#   step_meta_smote(
#     raw_df_id    = raw_df_id,
#     raw_df_sum   = raw_df_sum,
#     site_id_col  = site_id_col,
#     target_col   = target_variable,
#     mni_col      = mni_col,
#     feature_cols = feature_cols,
#     augment_fun  = my_augment_fun,
#     neighbors    = 1,
#     over_ratio   = 1,
#     skip         = TRUE
#   )

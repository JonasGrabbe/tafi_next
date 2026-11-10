# Metadata-stratified centroids for the training-only interpolation donor pool.
# Original aggregate rows remain the ONLY real training rows. Metadata centroids
# are NOT appended directly to the fitting matrix; they are donors for interpolation.
# Missing metadata: skip only combinations/rows needing that metadata.
# Anatomical/index zeros: valid, even when every dimension is zero.
# Empty metadata groups are skipped; only non-finite GENERATED means discarded.

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
})

buckets <- list(
  temporal = c("geological_period", "cultural_period", "culture_period_2",
               "start_date_cal_bp", "end_date_cal_bp", "mis_stage"),
  spatial = c("region", "country"),
  site_context = c("karstic_system", "open_air_site"),
  biological = c("age_category", "gender"),
  taxonomy = c("hominin_species")
)

generate_combos <- function(buckets, min_buckets = 1L, max_buckets = 3L) {
  min_buckets <- as.integer(min_buckets)
  max_buckets <- min(as.integer(max_buckets), length(buckets))
  if (min_buckets < 1L || max_buckets < min_buckets) stop("Invalid bucket limits")
  subsets <- unlist(lapply(seq.int(min_buckets, max_buckets), function(k) {
    combn(names(buckets), k, simplify = FALSE)
  }), recursive = FALSE)
  unlist(lapply(subsets, function(bnames) {
    grid <- do.call(expand.grid, c(buckets[bnames], stringsAsFactors = FALSE))
    lapply(seq_len(nrow(grid)), function(i) as.list(grid[i, , drop = FALSE]))
  }), recursive = FALSE)
}

# Both real NA and text placeholders indicate unavailable metadata.
# Zero is NOT treated as missing.
is_known_metadata <- function(x) {
  z <- tolower(trimws(as.character(x)))
  !is.na(x) & !is.na(z) & nzchar(z) &
    !(z %in% c("na", "n/a", "nan", "null", "unknown", "unkown",
               "not recorded", "not known", "missing", "-"))
}

# A single metadata combination, applied to ONE class and training IDs only.
augment_with_combination <- function(df, combo, aug_id, feature_cols,
                                     integer_cols = character(0),
                                     required_cols = character(0),
                                     min_group_size = 2L,
                                     min_unique_sites = 2L,
                                     site_id_col = "site_id") {
  combo_cols <- as.character(unlist(combo, use.names = FALSE))
  if (!all(combo_cols %in% names(df))) return(tibble::tibble())
  if (!all(feature_cols %in% names(df))) stop("Metadata donor feature missing")
  if (!all(vapply(df[feature_cols], is.numeric, logical(1)))) {
    stop("Metadata donor predictors must all be numeric")
  }
  if (length(required_cols)) {
    stop("Required-column imputation is disabled: do not invent MNI or anatomy")
  }

  eligible <- Reduce(`&`, lapply(df[combo_cols], is_known_metadata))
  # Keep anatomical zeros; reject only non-finite feature vectors.
  X <- as.matrix(df[, feature_cols, drop = FALSE])
  eligible <- eligible & rowSums(!is.finite(X)) == 0L
  valid <- df[eligible, , drop = FALSE]
  if (nrow(valid) < min_group_size) return(tibble::tibble())

  has_site_id <- site_id_col %in% names(valid)
  if (has_site_id) {
    valid[[site_id_col]] <- as.character(valid[[site_id_col]])
  }
  grouped <- valid |>
    dplyr::group_by(dplyr::across(dplyr::all_of(combo_cols))) |>
    dplyr::filter(dplyr::n() >= min_group_size)

  if (has_site_id) {
    # Avoid synthetic group means made from multiple individual rows of only
    # one archaeological site (that aggregate already exists in db_sum).
    grouped <- grouped |>
      dplyr::filter(dplyr::n_distinct(.data[[site_id_col]]) >= min_unique_sites)
  }
  if (nrow(grouped) == 0L) return(tibble::tibble())

  generated <- grouped |>
    dplyr::summarise(
      dplyr::across(dplyr::all_of(feature_cols), ~ mean(.x)),
      .meta_n_parent_rows = dplyr::n(),
      .meta_source_sites = if (has_site_id) {
        paste(sort(unique(.data[[site_id_col]])), collapse = ";")
      } else NA_character_,
      .groups = "drop"
    )
  if (nrow(generated) == 0L) return(tibble::tibble())

  for (nm in intersect(integer_cols, feature_cols)) {
    generated[[nm]] <- as.numeric(round(generated[[nm]]))
  }

  G <- as.matrix(generated[, feature_cols, drop = FALSE])
  anatomy_cols <- setdiff(feature_cols, integer_cols)
  A <- as.matrix(generated[, anatomy_cols, drop = FALSE])
  # All-zero means may be legitimate (e.g. zero Heip evenness for sparse
  # skeletons). They are never treated as evidence of an empty metadata group:
  # empty groups have already been excluded by n() / n_distinct(site_id).
  nonfinite <- rowSums(!is.finite(G)) > 0L
  zero_vector <- !nonfinite & rowSums(abs(A)) == 0
  n_nonfinite_removed <- sum(nonfinite)
  n_zero_valid_candidates <- sum(zero_vector)
  keep <- !nonfinite
  generated <- generated[keep, , drop = FALSE]
  if (nrow(generated) == 0L) return(tibble::tibble())

  # Preserve original column types for metadata and ordinary unused fields.
  # Missing values in auxiliary metadata columns are allowed: they are NOT ML
  # predictors and NOT reasons to delete generated anatomical centroids.
  out <- df[rep(NA_integer_, nrow(generated)), , drop = FALSE]
  for (nm in intersect(names(generated), names(df))) out[[nm]] <- generated[[nm]]
  out$.aug_id <- rep(aug_id, nrow(generated))
  out$.meta_source_sites <- generated$.meta_source_sites
  out$.meta_n_parent_rows <- generated$.meta_n_parent_rows

  # Synthetic centroid identical to a source row adds no new donor geometry.
  out_unique <- dplyr::anti_join(out, df[, feature_cols, drop = FALSE],
                                 by = feature_cols)
  attr(out_unique, "metadata_stats") <- list(
    nonfinite_discarded = n_nonfinite_removed,
    all_zero_valid_candidates = n_zero_valid_candidates,
    exact_original_duplicate_discarded = nrow(out) - nrow(out_unique)
  )
  out_unique
}

build_augmented_dataset <- function(df, combos, feature_cols,
                                    integer_cols = character(0),
                                    required_cols = character(0),
                                    min_group_size = 2L,
                                    min_unique_sites = 2L,
                                    site_id_col = "site_id") {
  original <- df
  original$.aug_id <- "original"
  original$.meta_source_sites <- NA_character_
  original$.meta_n_parent_rows <- NA_integer_
  valid_combos <- Filter(function(combo) {
    all(unlist(combo, use.names = FALSE) %in% names(df))
  }, combos)
  generated_results <- purrr::imap(valid_combos, function(combo, i) {
    augment_with_combination(
      df = df, combo = combo,
      aug_id = paste0("meta_", i, "_", paste(unlist(combo), collapse = "+")),
      feature_cols = feature_cols,
      integer_cols = integer_cols,
      required_cols = required_cols,
      min_group_size = min_group_size,
      min_unique_sites = min_unique_sites,
      site_id_col = site_id_col
    )
  })
  removed_stat <- function(name) sum(vapply(generated_results, function(z) {
    stat <- attr(z, "metadata_stats")
    if (is.null(stat)) 0 else as.numeric(stat[[name]])
  }, numeric(1)))
  generated <- dplyr::bind_rows(generated_results)
  n_pre_dedup <- nrow(generated)
  if (nrow(generated)) {
    generated <- dplyr::distinct(
      generated, dplyr::across(dplyr::all_of(feature_cols)), .keep_all = TRUE
    )
  }
  out <- dplyr::bind_rows(original, generated)
  attr(out, "metadata_counts") <- list(
    n_combos_total = length(combos),
    n_combos_with_columns = length(valid_combos),
    n_candidates_retained = nrow(generated),
    n_nonfinite_discarded = removed_stat("nonfinite_discarded"),
    n_zero_valid_candidates = removed_stat("all_zero_valid_candidates"),
    n_source_duplicates_discarded = removed_stat("exact_original_duplicate_discarded"),
    n_between_combo_duplicates_discarded = n_pre_dedup - nrow(generated)
  )
  out
}

make_augment_fun <- function(combos, feature_cols, target_col,
                             integer_cols = character(0),
                             required_cols = character(0),
                             min_group_size = 2L,
                             min_unique_sites = 2L,
                             site_id_col = "site_id") {
  force(combos); force(feature_cols); force(target_col)
  function(x) {
    if (!(target_col %in% names(x)) || anyNA(x[[target_col]])) {
      stop("Missing training class labels in donor pool")
    }
    classes <- unique(as.character(x[[target_col]]))
    pieces <- lapply(classes, function(cls) {
      one <- x[as.character(x[[target_col]]) == cls, , drop = FALSE]
      augmented <- build_augmented_dataset(
        df = one, combos = combos, feature_cols = feature_cols,
        integer_cols = integer_cols, required_cols = required_cols,
        min_group_size = min_group_size,
        min_unique_sites = min_unique_sites, site_id_col = site_id_col
      )
      counts <- attr(augmented, "metadata_counts")
      # Explicitly restore class: grouped summaries have no class column.
      augmented[[target_col]] <- factor(
        rep(cls, nrow(augmented)), levels = levels(x[[target_col]])
      )
      list(rows = augmented, audit = tibble::tibble(
        class = cls,
        combinations_total = counts$n_combos_total,
        combinations_available = counts$n_combos_with_columns,
        metadata_centroids_retained = counts$n_candidates_retained,
        nonfinite_centroids_discarded = counts$n_nonfinite_discarded,
        all_zero_candidates_valid = counts$n_zero_valid_candidates,
        exact_source_duplicates_discarded = counts$n_source_duplicates_discarded,
        overlapping_combo_duplicates_discarded = counts$n_between_combo_duplicates_discarded
      ))
    })
    out <- dplyr::bind_rows(lapply(pieces, `[[`, "rows"))
    audit <- dplyr::bind_rows(lapply(pieces, `[[`, "audit"))
    # Retain centroid source-site provenance separately from classifier predictors.
    provenance <- out |>
      dplyr::filter(.data$.aug_id != "original") |>
      dplyr::select(dplyr::all_of(target_col), .aug_id,
                    .meta_source_sites, .meta_n_parent_rows)
    out <- out |>
      dplyr::select(-.aug_id, -.meta_source_sites, -.meta_n_parent_rows)
    attr(out, "metadata_audit") <- audit
    attr(out, "metadata_provenance") <- provenance
    out
  }
}

# Supervised SPR classification: shared model definitions and CV methods.
# Source through run_cv_spr.R or run_projection_spr.R.
# Site-level folds are independent; interpolation donors are restricted to
# training-site identifiers. Fixed hyperparameters; no nested tuning.

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(recipes)
  library(workflows)
  library(parsnip)
  library(withr)
})
tidymodels::tidymodels_prefer()

identity_aug <- function(x) x

default_smote_fun <- function(df, target_col, neighbors, over_ratio) {
  themis::smote(
    df,
    var = target_col,
    k = neighbors,
    over_ratio = over_ratio
  )
}

# -----------------------------------------------------------------------------
# SMOTE-style interpolation that returns ONLY newly generated rows.
# Generates an exact number of new within-class nearest-neighbour interpolations.
# Training folds keep their ORIGINAL aggregated rows separately.
# -----------------------------------------------------------------------------
interpolate_synthetic_only <- function(pool, target_col, feature_cols,
                                       target_new, neighbors = 3L,
                                       distance_cols = feature_cols,
                                       on_degenerate = NULL) {
  lev <- levels(pool[[target_col]])
  pieces <- list()
  for (cls in lev) {
    n_new <- as.integer(target_new[[cls]])
    if (is.na(n_new) || n_new < 0L) stop("Invalid synthetic count for ", cls)
    if (n_new == 0L) next

    donors <- pool[as.character(pool[[target_col]]) == cls, , drop = FALSE]
    donors <- dplyr::distinct(
      donors, dplyr::across(dplyr::all_of(distance_cols)), .keep_all = TRUE
    )
    n <- nrow(donors)
    if (n < 1L) {
      stop("Cannot augment class ", cls, ": no valid donor profile")
    }
    if (n == 1L) {
      # A constant feature within a class cannot yield distinct interpolants.
      # Keep the observed value: deterministic exact-profile oversampling.
      # This is recorded separately from genuine interpolation.
      if (is.function(on_degenerate)) on_degenerate(cls, n_new)
      repeated <- donors[rep(1L, n_new), feature_cols, drop = FALSE]
      repeated[[target_col]] <- factor(rep(cls, n_new), levels = lev)
      pieces[[cls]] <- repeated[, c(target_col, feature_cols), drop = FALSE]
      next
    }
    k <- min(as.integer(neighbors), n - 1L)
    if (k < neighbors) warning("Class ", cls, ": reducing k to ", k)

    mat <- as.matrix(donors[, feature_cols, drop = FALSE])
    distmat <- as.matrix(stats::dist(donors[, distance_cols, drop = FALSE]))
    diag(distmat) <- Inf
    starts <- sample.int(n, size = n_new, replace = TRUE)
    ends <- vapply(starts, function(i) {
      eligible <- order(distmat[i, ])[seq_len(k)]
      eligible[sample.int(length(eligible), size = 1L)]
    }, integer(1))
    alpha <- stats::runif(n_new)
    synth <- mat[starts, , drop = FALSE] +
      (mat[ends, , drop = FALSE] - mat[starts, , drop = FALSE]) * alpha

    # An entirely zero vector can be a valid anatomical/index observation,
    # especially in one-variable abundance or evenness spaces. Reject only
    # genuinely non-finite interpolations, never valid numerical zeros.
    bad_rows <- function(x) {
      rowSums(!is.finite(x)) > 0L
    }
    bad <- bad_rows(synth)
    attempts <- 0L
    while (any(bad) && attempts < 30L) {
      attempts <- attempts + 1L
      n_bad <- sum(bad)
      re_starts <- sample.int(n, size = n_bad, replace = TRUE)
      re_ends <- vapply(re_starts, function(i) {
        eligible <- order(distmat[i, ])[seq_len(k)]
        eligible[sample.int(length(eligible), size = 1L)]
      }, integer(1))
      re_alpha <- stats::runif(n_bad)
      synth[bad, ] <- mat[re_starts, , drop = FALSE] +
        (mat[re_ends, , drop = FALSE] - mat[re_starts, , drop = FALSE]) * re_alpha
      bad <- bad_rows(synth)
    }
    if (any(bad)) stop("Cannot generate requested valid synthetic profiles for ", cls)

    piece <- tibble::as_tibble(as.data.frame(synth, check.names = FALSE))
    piece[[target_col]] <- factor(rep(cls, n_new), levels = lev)
    pieces[[cls]] <- piece[, c(target_col, feature_cols), drop = FALSE]
  }
  if (!length(pieces)) return(pool[FALSE, c(target_col, feature_cols), drop = FALSE])
  dplyr::bind_rows(pieces)
}

make_smote_train_sum_id <- function(model_df_train,
                                    raw_df_id,
                                    raw_df_sum,
                                    site_id_col,
                                    target_col,
                                    mni_col,
                                    feature_cols,
                                    augment_fun = identity_aug,
                                    smote_fun = default_smote_fun, # not used in this interpolation implementation
                                    neighbors = 3L,
                                    over_ratio = 1,  # retained for prior-compatible calling conventions
                                    smote_use_mni = FALSE,
                                    synthetic_strategy = "multiplier",
                                    synthetic_multiplier = 3,
                                    synthetic_target_per_class = 150L,
                                    seed = 42L,
                                    audit_file = NULL,
                                    audit_context = list()) {
  synthetic_strategy <- match.arg(
    synthetic_strategy, c("multiplier", "fixed_total")
  )
  if (!length(feature_cols)) stop("Empty feature set")
  if (mni_col %in% feature_cols) stop("MNI must not be in feature_cols")
  required <- c(site_id_col, target_col)
  for (z in list(model_df_train, raw_df_id, raw_df_sum)) {
    if (!all(required %in% names(z))) stop("Missing ID/outcome columns")
  }
  if (anyDuplicated(as.character(model_df_train[[site_id_col]]))) {
    stop("Main training fold has duplicate site IDs")
  }
  train_ids <- unique(as.character(model_df_train[[site_id_col]]))
  if (anyNA(train_ids)) stop("NA site ID in training fold")
  original <- model_df_train[, c(target_col, feature_cols), drop = FALSE]
  target_levels <- levels(model_df_train[[target_col]])
  if (is.null(target_levels)) stop("Outcome must be a factor")
  original[[target_col]] <- factor(original[[target_col]], levels = target_levels)
  # By dataset definition zeros are valid, and non-finite anatomical values
  # should not occur. Fail conspicuously rather than silently discard sites.
  if (!all(vapply(original[feature_cols], is.numeric, logical(1)))) {
    stop("Original anatomical features must be numeric")
  }
  if (anyNA(original[[target_col]])) stop("Missing original training label")
  X_original <- as.matrix(original[, feature_cols, drop = FALSE])
  if (any(!is.finite(X_original))) {
    stop("Unexpected NA/NaN/Inf in ORIGINAL training anatomical features")
  }

  take_train <- function(x) {
    x[[site_id_col]] <- as.character(x[[site_id_col]])
    dplyr::filter(x, .data[[site_id_col]] %in% train_ids)
  }
  id_train <- take_train(raw_df_id)
  sum_train <- take_train(raw_df_sum)
  if (!nrow(id_train) || !nrow(sum_train)) stop("No auxiliary donors for training sites")
  # These are donors ONLY. They do not go into the fitting data as original rows.
  donors_before <- dplyr::bind_rows(id_train, sum_train)
  donor_pool <- augment_fun(donors_before)
  if (nrow(donor_pool) < nrow(donors_before)) {
    stop("Metadata augmentation cannot remove the original donor profiles")
  }
  donor_pool[[target_col]] <- factor(as.character(donor_pool[[target_col]]),
                                     levels = target_levels)
  if (anyNA(donor_pool[[target_col]])) stop("Missing or invalid donor outcome")

  distance_cols <- c(feature_cols, if (smote_use_mni) mni_col)
  if (!all(distance_cols %in% names(donor_pool))) stop("Missing donor features")
  if (!all(vapply(donor_pool[distance_cols], is.numeric, logical(1)))) {
    stop("Donor features must be numeric")
  }
  donor_pool <- donor_pool[, c(target_col, distance_cols), drop = FALSE]
  valid <- stats::complete.cases(donor_pool) &
    apply(as.matrix(donor_pool[, distance_cols, drop = FALSE]), 1L,
          function(x) all(is.finite(x)))
  n_removed <- sum(!valid)
  if (n_removed > 0L) warning("Dropping ", n_removed, " invalid donor rows")
  donor_pool <- donor_pool[valid, , drop = FALSE]

  original_counts <- table(factor(original[[target_col]], levels = target_levels))
  desired_total <- if (synthetic_strategy == "fixed_total") {
    setNames(rep(as.integer(synthetic_target_per_class), length(target_levels)),
             target_levels)
  } else {
    # Equalize every outcome class to multiplier * the LARGEST ORIGINAL
    # TRAINING-FOLD class. E.g. 32 F / 34 Rest, multiplier 3 -> 102 / 102.
    # Never compute this target using all reference sites or assessment sites.
    shared_target <- as.integer(ceiling(max(as.integer(original_counts)) *
                                        synthetic_multiplier))
    setNames(rep(shared_target, length(target_levels)), target_levels)
  }
  if (synthetic_strategy == "fixed_total" &&
      any(as.integer(original_counts) > desired_total)) {
    stop("fixed_total is smaller than an original training class count; ",
         "cannot balance without deleting original site records")
  }
  to_generate <- pmax(0L, desired_total - as.integer(original_counts))
  names(to_generate) <- target_levels

  # A CSV record documents any exact-profile oversampling required when
  # the class has only one distinct donor vector in the selected feature space.
  write_degenerate_audit <- function(class_name, n_repeated) {
    if (is.null(audit_file) || !nzchar(audit_file)) return(invisible(NULL))
    context_field <- function(key) {
      value <- audit_context[[key]]
      if (is.null(value) || length(value) != 1L) NA_character_
      else as.character(value)
    }
    record <- data.frame(
      level = context_field("level"),
      task = context_field("task"),
      space = context_field("space"),
      model = context_field("model"),
      condition = context_field("condition"),
      fold = context_field("fold"),
      class = class_name,
      distinct_donor_profiles = 1L,
      repeated_training_rows = as.integer(n_repeated),
      method = "exact_profile_oversampling",
      stringsAsFactors = FALSE
    )
    dir.create(dirname(audit_file), recursive = TRUE, showWarnings = FALSE)
    append_existing <- file.exists(audit_file)
    utils::write.table(record, file = audit_file, sep = ",", quote = TRUE,
                       row.names = FALSE, col.names = !append_existing,
                       append = append_existing, na = "")
    invisible(NULL)
  }
  # Deterministic synthetic draws, isolated from the random-forest RNG.
  # Same training sites and same donor pool -> same rows on rerun.
  synth <- withr::with_seed(as.integer(seed), {
    interpolate_synthetic_only(
      donor_pool, target_col, feature_cols,
      target_new = to_generate, neighbors = neighbors,
      distance_cols = distance_cols,
      on_degenerate = write_degenerate_audit
    )
  })
  # The final fitting matrix contains only the ORIGINAL training sites + new
  # synthetic interpolations. No original individual or metadata donor is added.
  result <- dplyr::bind_rows(
    tibble::as_tibble(original), tibble::as_tibble(synth)
  )
  stopifnot(nrow(result) == nrow(original) + nrow(synth))
  result
}

# Recipe construction. Site identifiers are IDs, never predictors.
build_recipe_for_task <- function(df,
                                  target_variable = "taphonomic_context",
                                  site_name_col = "site_name",
                                  raw_df_id = NULL,
                                  raw_df_sum = NULL,
                                  site_id_col = "site_id",
                                  mni_col = "mni",
                                  feature_cols = NULL,
                                  use_smote = FALSE,
                                  use_downsample = FALSE,
                                  use_metadata_aug = FALSE,
                                  smote_use_mni = FALSE,
                                  smote_neighbors = 3L,
                                  smote_over_ratio = 1,
                                  synthetic_strategy = "multiplier",
                                  synthetic_multiplier = 3,
                                  synthetic_target_per_class = 150L,
                                  metadata_min_unique_sites = 2L,
                                  seed = 42L,
                                  normalize = FALSE,
                                  audit_file = NULL,
                                  audit_context = list()) {
  rec <- recipes::recipe(
    as.formula(paste(target_variable, "~ .")), data = df
  ) %>%
    recipes::update_role(!!rlang::sym(site_name_col), new_role = "ID")
  # ReMNAU and MNI-normalized proportional representation remain on their native scales.

  if (use_metadata_aug && !use_smote) {
    stop("Metadata augmentation is currently implemented within the custom SMOTE pool; enable use_smote")
  }
  if (use_smote) {
    if (is.null(raw_df_id) || is.null(raw_df_sum)) stop("Missing auxiliary tables")
    augmentation <- identity_aug
    if (use_metadata_aug) {
      if (!exists("make_augment_fun", mode = "function") ||
          !exists("generate_combos", mode = "function") ||
          !exists("buckets")) {
        stop("Source spr_metadata_augmentation.R before enabling use_metadata_aug")
      }
      augmentation <- make_augment_fun(
        combos = generate_combos(buckets),
        feature_cols = c(feature_cols, if (smote_use_mni) mni_col),
        target_col = target_variable,
        integer_cols = if (smote_use_mni) mni_col else character(),
        required_cols = character(), # NEVER invent MNI = 1
        min_group_size = 2L,
        min_unique_sites = metadata_min_unique_sites,
        site_id_col = site_id_col
      )
    }
    rec <- step_meta_smote(
      recipe = rec,
      raw_df_id = raw_df_id,
      raw_df_sum = raw_df_sum,
      site_id_col = site_id_col,
      target_col = target_variable,
      mni_col = mni_col,
      feature_cols = feature_cols,
      augment_fun = augmentation,
      neighbors = smote_neighbors,
      over_ratio = smote_over_ratio,
      smote_use_mni = smote_use_mni,
      synthetic_strategy = synthetic_strategy,
      synthetic_multiplier = synthetic_multiplier,
      synthetic_target_per_class = synthetic_target_per_class,
      seed = seed,
      audit_file = audit_file,
      audit_context = audit_context,
      skip = TRUE
    )
  }
  if (isTRUE(normalize)) {
    rec <- rec %>% recipes::step_zv(recipes::all_numeric_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric_predictors())
  }
  rec
}

# -----------------------------------------------------------
# --. paasssgt  Model spec + tunable parameters for each classifier
# -----------------------------------------------------------

step_meta_smote <- function(recipe,
                            raw_df_id,
                            raw_df_sum,
                            site_id_col = "site_id",
                            target_col,
                            mni_col = "mni",
                            feature_cols,
                            augment_fun = identity_aug,
                            smote_fun = default_smote_fun,
                            neighbors = 3L,
                            over_ratio = 1,
                            smote_use_mni = FALSE,
                            synthetic_strategy = "multiplier",
                            synthetic_multiplier = 3,
                            synthetic_target_per_class = 150L,
                            seed = 42L,
                            audit_file = NULL,
                            audit_context = list(),
                            skip = TRUE,
                            id = recipes::rand_id("meta_smote")) {
  recipes::add_step(recipe, structure(
    list(raw_df_id = raw_df_id,
         raw_df_sum = raw_df_sum,
         site_id_col = site_id_col,
         target_col = target_col,
         mni_col = mni_col,
         feature_cols = feature_cols,
         augment_fun = augment_fun,
         smote_fun = smote_fun,
         neighbors = neighbors,
         over_ratio = over_ratio,
         smote_use_mni = smote_use_mni,
         synthetic_strategy = synthetic_strategy,
         synthetic_multiplier = synthetic_multiplier,
         synthetic_target_per_class = synthetic_target_per_class,
         seed = as.integer(seed),
         audit_file = audit_file,
         audit_context = audit_context,
         role = NA, trained = FALSE, skip = skip, id = id),
    class = c("step_meta_smote", "step")
  ))
}

prep.step_meta_smote <- function(x, training, info = NULL, ...) {
  x$trained <- TRUE
  x
}

bake.step_meta_smote <- function(object, new_data, ...) {
  # Return original training-fold aggregates PLUS synthetic interpolations.
  # db_id and metadata centroids are used ONLY as interpolation donors.
  make_smote_train_sum_id(
    model_df_train = new_data,
    raw_df_id = object$raw_df_id,
    raw_df_sum = object$raw_df_sum,
    site_id_col = object$site_id_col,
    target_col = object$target_col,
    mni_col = object$mni_col,
    feature_cols = object$feature_cols,
    augment_fun = object$augment_fun,
    smote_fun = object$smote_fun,
    neighbors = object$neighbors,
    over_ratio = object$over_ratio,
    smote_use_mni = object$smote_use_mni,
    synthetic_strategy = object$synthetic_strategy,
    synthetic_multiplier = object$synthetic_multiplier,
    synthetic_target_per_class = object$synthetic_target_per_class,
    seed = object$seed,
    audit_file = object$audit_file,
    audit_context = object$audit_context
  )
}

print.step_meta_smote <- function(x, ...) {
  cat("Meta SMOTE augmentation [", if (x$trained) "trained" else "untrained", "]\n")
  invisible(x)
}


tapho.loader.wfl7 <- function(
  data,
  sum               = TRUE,
  unique_level      = TRUE,
  level             = NULL,
  calc_level        = c("min", "max", "mode"),
  recalc_method     = c("min", "filter", "direct"),
  indet             = TRUE,
  indet_method      = c("up", "dist", "delet"),
  group_by          = NULL,
  select            = NULL,
  fragments_per_bone = 3,
  exclude_sections  = c("teeth"),
  data_variant        = "ProportionalRep_MNIbased"   # <-- new parameter
) {

  # ── Setup ──────────────────────────────────────────────────────────────────
  meta <- meta_set(indet)

  on.exit({
    if ("meta_env" %in% search()) {
      detach("meta_env")
    }
  }, add = TRUE)
  
  calc_level    <- match.arg(calc_level)
  indet_method  <- match.arg(indet_method)
  recalc_method <- match.arg(recalc_method)

  # ── Stage 1: Load and parse all Excel sheets ───────────────────────────────
  db <- read_excel_sheets(
    file_path          = data,
    level              = level,
    calc_level         = calc_level,
    indet              = indet,
    indet_method       = indet_method,
    fragments_per_bone = fragments_per_bone,
    meta               = meta
  )

  # ── Stage 2: Define group_by variants ─────────────────────────────────────
  if (!"site_id" %in% names(db)) {
    stop("The Excel loader did not return the required site_id column")
  }
  raw_id <- trimws(as.character(db$site_id))
  if (anyNA(raw_id) || any(!nzchar(raw_id))) {
    stop("One or more Excel source records have a missing/blank site_id")
  }
  # Each Excel source record is linked to its stable site_id.
  message("Input records: ", nrow(db), " across ",
          length(unique(raw_id)), " distinct site_id values")


  # One aggregate row must represent ONE archaeological site_id.
  # Grouping on site_name can split one site_id if its descriptive name varies
  # across the original records. Never repair this by distinct(site_id) after
  # MNE calculation, because that would discard rather than aggregate counts.
  group_by_sum <- if (is.null(group_by)) "site_id" else union("site_id", group_by)
  group_by_id  <- if ("site_name" %in% group_by) setdiff(group_by, "site_name") else group_by

  # ── Stage 3: Helper — aggregate + MNE pipeline for a given grouping ────────
  process_db <- function(grp) {

    d <- select_and_group_by(
      data         = db,
      group_by     = grp,
      select       = select,
      method       = recalc_method,
      level        = level,
      unique_level = unique_level
    )

    element_lookup <- elements_data %>%
      dplyr::transmute(
        element_key        = gsub(" ", ".", ElementsOne),
        One_Skeleton_count = Bone_Count,
        Skeleton_Section   = .data$Skeleton_Section
      )

    results_list <- lapply(seq_len(nrow(d)), function(i) {

      input_df <- data.frame(
        element = d$element[[i]],
        count   = d$count[[i]]
      )

      bone_data <- input_df %>%
        dplyr::left_join(
          elements_data %>%
            dplyr::select(ElementsOne, Bone_Count, Skeleton_Section),
          by = c("element" = "ElementsOne")
        ) %>%
        dplyr::filter(!Skeleton_Section %in% exclude_sections)

      if (nrow(bone_data) == 0) {
        return(NULL)
      }

      skeleton <- data.frame(
        element = bone_data$element,
        count   = bone_data$Bone_Count
      )

      MNI <- d$mni[[i]]

      var_db <- convertMNEtoVariables_v3(
        MNE         = bone_data$count,
        skeleton    = skeleton,
        individuals = MNI
      )

      out_elements <- var_db$element_2

      if (is.list(out_elements)) {
        out_elements <- out_elements[[1]]
      }

      out_keys <- gsub(" ", ".", out_elements)

      out_one_skeleton_count <- element_lookup$One_Skeleton_count[
        match(out_keys, element_lookup$element_key)
      ]

      out_sections <- element_lookup$Skeleton_Section[
        match(out_keys, element_lookup$element_key)
      ]

      var_db %>%
        dplyr::mutate(
          One_Skeleton_count = as.numeric(out_one_skeleton_count),
          Skeleton_Section   = as.character(out_sections)
        )
    })

    if (any(vapply(results_list, is.null, logical(1)))) {
      stop("Empty anatomical profile after filtering; check input site(s)")
    }
    keys <- as.character(results_list[[1]]$element_2)
    key_clean <- gsub("[^[:alnum:]_]+", "_", keys)
    key_clean <- gsub("^_+|_+$", "", key_clean)
    if (anyDuplicated(key_clean)) stop("Anatomical element names collide after cleaning")
    for (r in results_list) {
      if (!identical(as.character(r$element_2), keys)) {
        stop("Inconsistent anatomical element definitions/order between sites")
      }
    }
    wide_numeric <- function(field, prefix) {
      mat <- do.call(rbind, lapply(results_list, function(x) as.numeric(x[[field]])))
      colnames(mat) <- paste0(prefix, key_clean)
      tibble::as_tibble(mat, .name_repair = "minimal")
    }
    prop <- wide_numeric("ProportionalRep_MNIbased", "ProportionalRep_MNIbased_")
    rem <- wide_numeric("ReMNAU", "ReMNAU_")
    out <- dplyr::select(d, -dplyr::all_of(c("count", "element")))
    out$element <- lapply(results_list, `[[`, "element_2")
    out$One_Skeleton_count <- lapply(results_list, `[[`, "One_Skeleton_count")
    out$Skeleton_Section <- lapply(results_list, `[[`, "Skeleton_Section")
    out$MNE <- lapply(results_list, `[[`, "MNE")
    out$MNAU <- lapply(results_list, `[[`, "MNAU")
    out$PercentMNAU <- lapply(results_list, `[[`, "PercentMNAU")
    out$MNI_used_assemblage <- lapply(results_list, `[[`, "MNI_used_assemblage")
    out$Mean_part_abundance <- vapply(results_list,
      function(x) as.numeric(x$Mean_part_abundance[[1]]), numeric(1))

    derived <- purrr::pmap_dfr(list(
      element = out$element,
      One_Skeleton_count = out$One_Skeleton_count,
      Skeleton_Section = out$Skeleton_Section,
      MNE = out$MNE,
      MNI_used = out$MNI_used_assemblage
    ), make_tapho_derived_features)
    density <- purrr::pmap_dfr(list(
      element = out$element,
      One_Skeleton_count = out$One_Skeleton_count,
      Skeleton_Section = out$Skeleton_Section,
      MNE = out$MNE,
      MNI_used = out$MNI_used_assemblage,
      PercentMNAU = out$PercentMNAU,
      MNAU = out$MNAU
    ), make_tapho_density_features)
    dplyr::bind_cols(out, prop, rem, derived, density,
                     .name_repair = "check_unique")
  }

  # ── Stage 4: Run pipeline for both groupings ───────────────────────────────
  db_sum <- process_db(group_by_sum)
  db_id  <- process_db(group_by_id)

  # ── Stage 5: Return both as a named list ──────────────────────────────────
  list(db_id = db_id, db_sum = db_sum)
}




# ---------------------------------------------------------------------------
# Explicit spaces: anatomical vectors depend on the selected anatomical level.
# ---------------------------------------------------------------------------
ML7_ENGINEERED <- c(
  "mnau_proximal_appendicular_mean",
  "mnau_lower_appendicular_mean",
  "mnau_distal_appendicular_mean",
  "mnau_appendicular_mean",
  "mnau_postcranial_total_mean",
  "mnau_top2_share_mcnaughton",
  "mnau_upper_appendicular_mean",
  "low_density_survival_index",
  "mnau_vertebrae_to_thorax_ratio",
  "all_sections_dominance",
  "mnau_lower_upper_appendicular_diff",
  "mnau_heip_evenness",
  "mnau_girdle_appendicular_logratio",
  "postcranial_continuity_hm",
  "mnau_girdle_to_appendicular_ratio",
  "mnau_gini",
  "density_representation_r2",
  "mnau_cv",
  "mnau_girdle_appendicular_diff",
  "mnau_full_axial_mean",
  "mnau_arms",
  "dominance_gap_all",
  "mnau_girdle_share_of_appendicular",
  "section_completeness_025_all",
  "whole_body_continuity_hm",
  "density_representation_slope",
  "mnau_limb_shaft_proxy_mean",
  "density_residual_sd"
)


# -----------------------------------------------------------------------------
# Resolve structurally undefined DERIVED indices, without modifying observed MNE
# or any valid anatomical proportions. Some mathematical indices are undefined
# for assemblages with zero total representation, one measured skeletal section,
# or insufficient variation for a density regression. An undefined derived index
# is NOT a missing MNE measurement; to retain all valid sites, its operational
# ML encoding is 0, with per-site provenance written to CSV.
#
# For Heip evenness, recompute from the cached original MNE and skeletal-section
# information using the same section pooling as make_tapho_derived_features().
# A zero total or <2 eligible sections has defined ML convention E=0.
# For other engineered indices that are unestimable, 0 is a modelling convention,
# NOT a scientific claim that the physical statistic itself equals zero.
# -----------------------------------------------------------------------------
ml7_heip_from_source <- function(mne, one_skeleton_count,
                                  section, mni) {
  sec_names <- c("cranium", "teeth", "hyoid", "vertebrae", "thorax",
                 "shoulder", "arms", "hands", "pelvic", "legs", "foot")
  mne <- as.numeric(unlist(mne, use.names = FALSE))
  counts <- as.numeric(unlist(one_skeleton_count, use.names = FALSE))
  sec <- as.character(unlist(section, use.names = FALSE))
  mni_values <- as.numeric(unlist(mni, use.names = FALSE))
  # The converter repeats assemblage-level MNI once per skeletal element.
  # All repeated values must agree; they are NOT additional individuals.
  mni_unique <- unique(mni_values)
  if (length(mni_unique) != 1L || any(!is.finite(mni_values))) {
    stop("Cannot repair Heip index: inconsistent/non-finite repeated MNI")
  }
  mni <- mni_unique[[1L]]
  if (length(mne) != length(counts) || length(mne) != length(sec)) {
    stop("Cannot repair Heip index: cached MNE/section vectors have different lengths")
  }
  if (length(mni) != 1L || !is.finite(mni) || mni <= 0) {
    stop("Cannot repair Heip index: missing/nonpositive MNI")
  }
  if (any(!is.finite(mne)) || any(mne < 0, na.rm = TRUE)) {
    stop("Cannot repair Heip index: invalid original MNE counts")
  }
  ok <- is.finite(counts) & counts > 0 & !is.na(sec) & sec %in% sec_names
  parts <- split(which(ok), factor(sec[ok], levels = sec_names), drop = TRUE)
  section_values <- vapply(parts, function(ix) {
    sum(mne[ix]) / sum(counts[ix]) / mni
  }, numeric(1))
  k <- length(section_values)
  total <- sum(section_values)
  if (!is.finite(total)) stop("Cannot repair Heip index: non-finite section total")
  if (k < 2L || total <= 0) {
    reason <- if (k < 2L) "fewer_than_two_eligible_sections" else "zero_total_representation"
    if (k == 0L && any(mne > 0)) {
      stop("Cannot repair Heip: positive MNE with no valid section/count mapping")
    }
    return(list(value = 0, reason = reason, n_sections = k,
                section_total = total))
  }
  p <- section_values / total
  p <- p[p > 0]
  H <- -sum(p * log(p))
  val <- (exp(H) - 1) / (k - 1)
  if (!is.finite(val) || val < -1e-8 || val > 1+1e-8) {
    stop("Cannot repair Heip index: recovered value is out of [0, 1]")
  }
  list(value = max(0, min(1, val)), reason = "recomputed_from_original_MNE",
       n_sections = k, section_total = total)
}

ml7_resolve_derived_indices <- function(tables, audit_file) {
  if (!all(c("db_sum", "db_id") %in% names(tables))) {
    stop("Expected db_sum and db_id before derived-index audit")
  }
  audited <- list()
  for (tblname in c("db_sum", "db_id")) {
    d <- tables[[tblname]]
    required <- c("site_id", "MNE", "One_Skeleton_count",
                  "Skeleton_Section", "MNI_used_assemblage",
                  "Mean_part_abundance", ML7_ENGINEERED)
    if (length(setdiff(required, names(d)))) {
      stop("Missing needed columns for derived-index audit in ", tblname,
           ": ", paste(setdiff(required, names(d)), collapse = ", "))
    }
    core_family <- grep("^(ReMNAU_|ProportionalRep_MNIbased_)",
                        names(d), value = TRUE)
    base <- unique(c("Mean_part_abundance", core_family))
    if (!all(vapply(d[base], is.numeric, logical(1)))) {
      stop("Non-numeric observed anatomical features in ", tblname)
    }
    observed <- as.matrix(d[, base, drop = FALSE])
    if (any(!is.finite(observed))) {
      bad <- which(rowSums(!is.finite(observed)) > 0L)
      stop("Non-finite ORIGINAL skeletal representation in ", tblname,
           "; site_id=", paste(head(d$site_id[bad], 20L), collapse = ", "),
           ". Only derived indices receive operational zero coding.")
    }
    for (feature in ML7_ENGINEERED) {
      x <- d[[feature]]
      if (!is.numeric(x)) stop("Non-numeric derived index: ", feature)
      ix <- which(!is.finite(x))
      if (!length(ix)) next
      for (i in ix) {
        outcome <- if (feature == "mnau_heip_evenness") {
          ml7_heip_from_source(d$MNE[[i]], d$One_Skeleton_count[[i]],
                                d$Skeleton_Section[[i]],
                                d$MNI_used_assemblage[[i]])
        } else {
          list(value = 0, reason = "unestimable_derived_index_encoded_zero",
               n_sections = NA_integer_, section_total = NA_real_)
        }
        audited[[length(audited) + 1L]] <- tibble::tibble(
          source_table = tblname,
          site_id = as.character(d$site_id[[i]]),
          site_name = if ("site_name" %in% names(d))
            paste(as.character(unlist(d$site_name[[i]])), collapse = ";") else NA_character_,
          feature = feature, original_value = as.character(x[[i]]),
          encoded_value = as.numeric(outcome$value),
          reason = outcome$reason,
          n_eligible_sections = as.integer(outcome$n_sections),
          section_representation_total = as.numeric(outcome$section_total)
        )
        d[[feature]][i] <- outcome$value
      }
    }
    remain <- as.matrix(d[, ML7_ENGINEERED, drop = FALSE])
    if (any(!is.finite(remain))) stop("Unresolved derived indices in ", tblname)
    tables[[tblname]] <- d
  }
  audit <- if (length(audited)) dplyr::bind_rows(audited) else tibble::tibble(
    source_table = character(), site_id = character(),
    site_name = character(), feature = character(),
    original_value = character(), encoded_value = double(),
    reason = character(), n_eligible_sections = integer(),
    section_representation_total = double()
  )
  dir.create(dirname(audit_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(audit, audit_file)
  # Preserve source values in the prepared cache and export the applied index conventions.
  # Do not overwrite the original loader cache with harmonised values; retain
  # the untouched cache for provenance and reapply this deterministic policy.
  tables
}

ml7_space_manifest <- function(df) {
  p <- grep("^ProportionalRep_MNIbased_", names(df), value = TRUE)
  r <- grep("^ReMNAU_", names(df), value = TRUE)
  # Match identical underlying elements in the two families.
  if (!length(p) || !length(r)) stop("Both anatomical families must be present")
  if (!setequal(sub("^ProportionalRep_MNIbased_", "", p),
                sub("^ReMNAU_", "", r))) {
    stop("ReMNAU vs ProportionalRep anatomical elements do not match")
  }
  basic <- c("Mean_part_abundance", "mnau_heip_evenness")
  missing <- setdiff(c(basic, ML7_ENGINEERED), names(df))
  if (length(missing)) {
    stop("Missing required engineered/basic feature(s): ", paste(missing, collapse=", "))
  }
  out <- list(
    abundance = basic[1],
    evenness = basic[2],
    abundance_evenness = basic,
    remnau = r,
    pmnau = p,
    engineered = ML7_ENGINEERED,
    combined = unique(c(basic, r, p, ML7_ENGINEERED))
  )
  for (nm in names(out)) {
    feats <- out[[nm]]
    if (anyDuplicated(feats) || anyNA(feats)) stop("Invalid feature manifest ", nm)
    if (!all(vapply(df[feats], is.numeric, logical(1)))) {
      stop("Non-numeric variable in feature manifest: ", nm)
    }
  }
  out
}

ml7_recode <- function(df, task, target = "taphonomic_context") {
  # Harmonise ALL known carnivore subcontexts consistently in db_sum/db_id.
  # The original archaeological label remains in the unrecoded source table.
  original_labels <- trimws(as.character(df[[target]]))
  carnivore_labels <- c("Carnivore", "Carnivore-related",
                        "Carnivore Scavenge (in situ)",
                        "Carnivore den(transport)",
                        "Carnivore den (transport)")
  original_labels[tolower(original_labels) %in% tolower(carnivore_labels)] <-
    "Carnivore"
  df[[target]] <- original_labels
  ref <- c("Funerary Context", "Cannibalism", "Carnivore", "Geological transport")
  if (task == "fc_vs_rest") {
    df <- df[df[[target]] %in% ref, , drop=FALSE]
    y <- ifelse(df[[target]] == "Funerary Context", "Funerary Context", "rest")
    levels <- c("Funerary Context", "rest")
  } else if (task == "cannibalism_vs_fc") {
    df <- df[df[[target]] %in% c("Funerary Context", "Cannibalism"), , drop=FALSE]
    y <- as.character(df[[target]])
    levels <- c("Funerary Context", "Cannibalism")
  } else if (task == "4class") {
    df <- df[df[[target]] %in% ref, , drop=FALSE]
    y <- as.character(df[[target]])
    levels <- ref
  } else stop("Unknown task: ", task)
  df[[target]] <- factor(y, levels=levels)
  if (anyNA(df[[target]])) stop("Invalid supervised label")
  df
}

ml7_validate <- function(df, cols, context, target = "taphonomic_context") {
  if (!all(c("site_id", target, cols) %in% names(df))) {
    stop("Missing site ID, target or features in ", context)
  }
  xx <- as.matrix(df[, cols, drop=FALSE])
  if (!all(is.finite(xx))) {
    ix <- which(rowSums(!is.finite(xx)) > 0L)
    stop("Unexpected NA/NaN/Inf anatomical predictor(s) in ",context,
         ": site_id=",paste(head(df$site_id[ix], 15L),collapse=", "))
  }
  if (anyNA(df[[target]])) stop("Missing labelled outcome in ", context)
  invisible(TRUE)
}

ml7_fixed_model <- function(model, n_features, n_classes, fixed, seed) {
  if (model == "LogisticRegression") {
    # glmnet explicitly requires >= 2 predictor columns. Two spaces in this
    # experiment (abundance / evenness) each have exactly one predictor.
    # Use the SAME nnet multinomial-logistic ridge engine for *all* spaces and
    # both 2-class and 4-class outcomes: no added dummy predictor, no silent
    # change from penalized to unpenalized logistic for one-feature spaces.
    # For 2-class outcomes nnet::multinom fits the binary logistic special case.
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Install the R package 'nnet' for LogisticRegression")
    }
    parsnip::multinom_reg(
      mode = "classification",
      penalty = fixed$logistic_penalty
    ) %>%
      parsnip::set_engine("nnet", MaxNWts = 10000L, maxit = 500L)
  } else if (model == "SVM_RBF") {
    # The SVM is fit by ml7_svm_e1071_fit() using the identical prepared
    # recipe, rather than the kernlab probability-calibration path that
    # Keep the class-probability matrix in its original orientation.
    stop("SVM_RBF must use ml7_svm_e1071_fit(), not ml7_fixed_model()")
  } else if (model == "KNN") {
    parsnip::nearest_neighbor(mode="classification", neighbors=fixed$knn_neighbors,
      weight_func="rectangular", dist_power=2) %>% parsnip::set_engine("kknn")
  } else if (model == "RandomForest") {
    parsnip::rand_forest(mode="classification",
      mtry=min(as.integer(n_features), as.integer(fixed$rf_mtry)),
      min_n=fixed$rf_min_n, trees=fixed$rf_trees) %>%
      parsnip::set_engine("ranger", seed=as.integer(seed), num.threads=1L,
        importance="permutation")
  } else stop("Unknown model: ",model)
}

# -----------------------------------------------------------------------
# RBF-SVM via LIBSVM (e1071) -- consistent engine in *all* folds and fits.
# Training uses the SAME recipe as the other classifiers.  bake(..., NULL)
# retains the training-only synthetic rows; bake(..., test) skips the custom
# augmentation step, preventing assessment/unknown rows from acting as donors.
# This avoids the kernlab 4-class Platt-probability failure seen in v3.6.
# -----------------------------------------------------------------------
ml7_svm_e1071_fit <- function(recipe_obj, train, feature_cols, outcome,
                              fixed, seed) {
  if (!requireNamespace("e1071", quietly = TRUE)) {
    stop("SVM_RBF needs e1071: install.packages('e1071')")
  }
  prepared <- recipes::prep(recipe_obj, training = train, retain = TRUE)
  model_train <- recipes::bake(prepared, new_data = NULL)
  active <- intersect(feature_cols, names(model_train))
  if (length(active) < 1L) {
    stop("SVM_RBF has no predictors after fold-specific zero-variance filtering")
  }
  xx <- as.matrix(model_train[, active, drop = FALSE])
  storage.mode(xx) <- "double"
  yy <- factor(as.character(model_train[[outcome]]),
               levels = levels(train[[outcome]]))
  if (any(!is.finite(xx)) || anyNA(yy) || nlevels(yy) < 2L) {
    stop("SVM_RBF: invalid prepared training predictors or class labels")
  }
  if (any(table(yy) == 0L)) stop("SVM_RBF: empty training class")
  # RBF gamma = 1 / number of retained predictors, with unit-scale input.
  # e1071::svm uses LIBSVM one-against-one multi-class classification and
  # internally estimates its probability model with probability=TRUE.
  set.seed(as.integer(seed))
  fitted <- e1071::svm(
    x = xx, y = yy, type = "C-classification", kernel = "radial",
    cost = as.numeric(fixed$svm_cost),
    gamma = 1 / as.numeric(ncol(xx)),
    scale = FALSE, probability = TRUE,
    tolerance = 0.001, cachesize = 128
  )
  structure(list(
    fitted = fitted, prepared = prepared, features = active,
    class_levels = levels(yy), engine = "e1071::svm (LIBSVM)",
    sigma_or_gamma = 1 / as.numeric(ncol(xx))
  ), class = "ml7_svm_e1071")
}

ml7_svm_e1071_predict <- function(object, new_data, labels) {
  if (!inherits(object, "ml7_svm_e1071")) stop("Not an e1071 SVM fit")
  if (!identical(as.character(labels), as.character(object$class_levels))) {
    stop("SVM_RBF class-order mismatch between fit and prediction")
  }
  if (!nrow(new_data)) {
    mat <- matrix(numeric(), nrow = 0L, ncol = length(labels))
    colnames(mat) <- labels
    return(mat)
  }
  processed <- recipes::bake(object$prepared, new_data = new_data)
  if (nrow(processed) != nrow(new_data)) {
    stop("SVM_RBF assessment rows changed during recipe baking")
  }
  x <- as.matrix(processed[, object$features, drop = FALSE])
  storage.mode(x) <- "double"
  if (any(!is.finite(x))) stop("SVM_RBF non-finite assessment predictor")
  pred <- stats::predict(object$fitted, newdata = x, probability = TRUE)
  prob <- attr(pred, "probabilities")
  if (is.null(prob) || length(dim(prob)) != 2L ||
      nrow(prob) != nrow(new_data) ||
      !all(labels %in% colnames(prob))) {
    stop("LIBSVM did not return a valid probability matrix: ",
         "expected ", nrow(new_data), " rows and classes [",
         paste(labels, collapse=", "), "]")
  }
  prob <- as.matrix(prob[, labels, drop = FALSE])
  storage.mode(prob) <- "double"
  if (any(!is.finite(prob)) || any(prob < -1e-8) ||
      any(prob > 1 + 1e-8) ||
      any(abs(rowSums(prob) - 1) > 1e-4)) {
    stop("LIBSVM returned invalid class probabilities; ",
         "do not treat them as valid scores or replace with 0")
  }
  # Only correct machine-rounding excursions (never substantial invalidity).
  # IMPORTANT: pmin(1, pmax(0, prob)) can drop dim(prob) in base R
  # because its first argument is scalar. In-place indexing retains the
  # assessment-row x class matrix, including 1-row and 1-feature cases.
  prob[prob < 0] <- 0
  prob[prob > 1] <- 1
  prob <- sweep(prob, 1L, rowSums(prob), `/`)
  if (!identical(dim(prob), c(nrow(new_data), length(labels)))) {
    stop("LIBSVM probability dimension changed after normalisation")
  }
  colnames(prob) <- labels
  prob
}

ml7_auc <- function(truth_binary, prob) {
  pos <- which(truth_binary)
  neg <- which(!truth_binary)
  if (!length(pos) || !length(neg)) return(NA_real_)
  ranks <- rank(prob, ties.method="average")
  (sum(ranks[pos]) - length(pos)*(length(pos)+1)/2) / (length(pos)*length(neg))
}

ml7_metrics <- function(truth, pred, prob) {
  labs <- levels(truth)
  nc <- length(labs)
  tr <- factor(truth, levels=labs)
  pr <- factor(pred, levels=labs)
  tp <- vapply(labs, function(k) sum(tr==k & pr==k), numeric(1))
  fp <- vapply(labs, function(k) sum(tr!=k & pr==k), numeric(1))
  fn <- vapply(labs, function(k) sum(tr==k & pr!=k), numeric(1))
  recall <- ifelse(tp+fn > 0, tp/(tp+fn), NA_real_)
  precision <- ifelse(tp+fp > 0, tp/(tp+fp), 0)
  f1 <- ifelse(precision+recall > 0,2*precision*recall/(precision+recall),0)
  aucs <- vapply(seq_along(labs), function(i) ml7_auc(tr==labs[i], prob[,i]), numeric(1))
  yy <- sapply(labs, function(k) as.numeric(tr == k))
  if (nc == 2L) {
    # Binary target class = Funerary Context, in first factor position.
    brier <- mean((prob[,1] - as.numeric(tr==labs[1]))^2)
    # vapply(labs, ...) names the per-class vectors by class. Strip the
    # inherited names before c(precision=...) below; otherwise R creates
    # e.g. `precision.Funerary Context` and pivot_longer(precision) fails.
    pval <- unname(precision[[1L]])
    rval <- unname(recall[[1L]])
    fval <- unname(f1[[1L]])
    auc <- unname(aucs[[1L]])
  } else {
    # Multi-class Brier = mean sum_k (predicted_k - observed_k)^2.
    brier <- mean(rowSums((prob-yy)^2))
    pval <- mean(precision); rval <- mean(recall); fval <- mean(f1)
    auc <- mean(aucs,na.rm=TRUE) # macro one-vs-rest
  }
  # Force every metric to a plain scalar, regardless of named class vectors.
  # The output schema must be identical for binary and four-class tasks.
  output <- c(
    balanced_accuracy = as.numeric(mean(recall)),
    auroc = as.numeric(auc),
    brier = as.numeric(brier),
    precision = as.numeric(pval),
    recall = as.numeric(rval),
    f1 = as.numeric(fval),
    accuracy = as.numeric(mean(pr == tr))
  )
  expected <- c("balanced_accuracy", "auroc", "brier", "precision",
                "recall", "f1", "accuracy")
  if (!identical(names(output), expected) ||
      length(output) != length(expected)) {
    stop("Internal metric naming error: ", paste(names(output), collapse=", "))
  }
  output
}

# nnet::multinom() yields a numeric VECTOR for binary class probabilities,
# representing the SECOND factor level. Some parsnip versions reshape that
# N-length vector as a 1 x N tibble with columns .pred_1 ... .pred_N.
# Only correct this exactly diagnosed shape, never silently assume class order.
ml7_prob_matrix <- function(x, labels, model = NULL, n_assess = NULL,
                            engine_levels = NULL) {
  stopifnot(length(labels) >= 2L)
  if (!is.data.frame(x)) stop("Predicted probabilities are not a data frame")
  expected <- paste0(".pred_", labels)
  if (!all(expected %in% names(x))) {
    alternate <- paste0(".pred_", make.names(labels))
    if (all(alternate %in% names(x))) {
      expected <- alternate
    } else {
      # For 2-class nnet::multinom, the N predicted probabilities of class 2
      # can be presented as columns instead of rows by parsnip. No other
      # model / shape / unknown orientation is eligible for this correction.
      binary_nnet_vector <-
        identical(model, "LogisticRegression") &&
        length(labels) == 2L &&
        !is.null(n_assess) &&
        nrow(x) == 1L &&
        ncol(x) == n_assess &&
        identical(names(x), paste0(".pred_", seq_len(n_assess)))

      if (!binary_nnet_vector) {
        stop("Cannot find probability columns for: ",
             paste(labels, collapse = ", "),
             ". Output dimensions: ", nrow(x), "x", ncol(x),
             "; available columns: ", paste(head(names(x), 40L), collapse = ", "))
      }
      if (is.null(engine_levels) ||
          !identical(as.character(engine_levels), as.character(labels))) {
        stop("Binary nnet model factor levels do not match CV outcome levels. ",
             "Refusing to infer probability-column orientation.")
      }
      prob_second <- as.numeric(unlist(x[1L, , drop = FALSE],
                                       use.names = FALSE))
      if (length(prob_second) != n_assess ||
          any(!is.finite(prob_second)) ||
          any(prob_second < -1e-6) || any(prob_second > 1+1e-6)) {
        stop("Invalid binary nnet predicted probabilities")
      }
      mat <- cbind(1-prob_second, prob_second)
      colnames(mat) <- labels
      return(mat)
    }
  }
  mat <- as.matrix(x[, expected, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- labels
  if (!is.null(n_assess) && nrow(mat) != n_assess) {
    stop("Probability predictions have ", nrow(mat), " rows; expected ", n_assess)
  }
  if (any(!is.finite(mat)) || any(mat < -1e-6) || any(mat > 1+1e-6)) {
    stop("Invalid predicted probabilities")
  }
  if (any(abs(rowSums(mat)-1) > 1e-4)) {
    stop("Predicted probabilities do not sum to one")
  }
  mat
}

ml7_make_splits <- function(df, task, repeats, seed, target="taphonomic_context") {
  # Same folds for all 7 spaces and all 4 algorithms.
  set.seed(as.integer(seed))
  splits <- rsample::vfold_cv(df, v=3L, repeats=as.integer(repeats),
    strata=!!rlang::sym(target))
  for (s in splits$splits) {
    train_ids <- as.character(rsample::analysis(s)$site_id)
    test_ids <- as.character(rsample::assessment(s)$site_id)
    if (length(intersect(train_ids,test_ids))) stop("Fold leakage: overlapping site IDs")
  }
  splits
}

ml7_one_fit <- function(level_tables, task, space, model, split, fold_index,
                       seed, fixed, use_smote, use_metadata_aug,
                       smote_use_mni, smote_neighbors,
                       synthetic_strategy, synthetic_multiplier,
                       synthetic_target_per_class, metadata_min_unique_sites,
                       target="taphonomic_context",
                       audit_file = NULL,
                       audit_context = list()) {
  manifests <- ml7_space_manifest(level_tables$db_sum)
  feats <- manifests[[space]]
  main <- ml7_recode(level_tables$db_sum, task)
  ids <- as.character(main$site_id)
  if (anyDuplicated(ids)) stop("Main aggregated dataset repeats a site_id")
  ml7_validate(main, feats, paste(task,space))
  # The split was built for the same ordered site-level rows.
  split_train_ids <- as.character(rsample::analysis(split)$site_id)
  split_test_ids <- as.character(rsample::assessment(split)$site_id)
  train <- main[match(split_train_ids, as.character(main$site_id)),
                c("site_id",target,feats),drop=FALSE]
  test <- main[match(split_test_ids, as.character(main$site_id)),
               c("site_id",target,feats),drop=FALSE]
  if (anyNA(train$site_id) || anyNA(test$site_id)) stop("Fold / data site ID mismatch")
  aux_id <- ml7_recode(level_tables$db_id, task)
  aux_sum <- main
  # IMPORTANT: metadata needed by augmentation! Use full site-level aux_sum,
  # not the predictor-only `main` subset; it still contains all metadata.
  rec <- build_recipe_for_task(
    df=train, target_variable=target, site_name_col="site_id",
    raw_df_id=aux_id, raw_df_sum=aux_sum,
    site_id_col="site_id", mni_col="mni", feature_cols=feats,
    use_smote=use_smote, use_downsample=FALSE,
    use_metadata_aug=use_metadata_aug,
    smote_use_mni=smote_use_mni, smote_neighbors=smote_neighbors,
    synthetic_strategy=synthetic_strategy,
    synthetic_multiplier=synthetic_multiplier,
    synthetic_target_per_class=synthetic_target_per_class,
    metadata_min_unique_sites=metadata_min_unique_sites,
    seed=as.integer(seed+fold_index*7919L),
    normalize = model %in% c("LogisticRegression","KNN","SVM_RBF"),
    audit_file = audit_file,
    audit_context = audit_context
  )
  if (identical(model, "SVM_RBF")) {
    # A single consistent SVM engine, for binary and multiclass tasks alike.
    fitted <- ml7_svm_e1071_fit(
      rec, train, feats, target, fixed,
      seed = as.integer(seed + fold_index * 83L)
    )
    p <- ml7_svm_e1071_predict(
      fitted, test, labels = levels(train[[target]])
    )
  } else {
  spec <- ml7_fixed_model(model, length(feats), nlevels(train[[target]]),
                          fixed=fixed,seed=as.integer(seed+fold_index))
  wf <- workflows::workflow() %>% workflows::add_recipe(rec) %>% workflows::add_model(spec)
  set.seed(as.integer(seed+fold_index*83L))
  # fit() is a generics S3 method; workflows does not export fit().
  fitted <- generics::fit(wf, data=train)
  probs <- predict(fitted,new_data=test,type="prob")
  # Align binary probabilities using the fitted model class order.
  engine_levels <- if (identical(model, "LogisticRegression") &&
                       nlevels(train[[target]]) == 2L) {
    workflows::extract_fit_engine(fitted)$lev
  } else NULL
  p <- ml7_prob_matrix(probs, levels(train[[target]]),
                       model = model, n_assess = nrow(test),
                       engine_levels = engine_levels)
  }  # end non-SVM workflow path
  if (nrow(p) != nrow(test) || ncol(p) != nlevels(train[[target]])) {
    stop("Probability matrix does not match the assessment rows and outcome classes")
  }
  pred_labels <- colnames(p)[max.col(p, ties.method="first")]
  m <- ml7_metrics(test[[target]], pred_labels, p)
  cols <- as.data.frame(p); names(cols) <- paste0("P_", make.names(colnames(p)))
  oof <- dplyr::bind_cols(
    tibble::tibble(site_id=test$site_id,
      truth=as.character(test[[target]]),
      predicted=pred_labels, fold_index=as.integer(fold_index)),
    tibble::as_tibble(cols)
  )
  list(metrics=m, predictions=oof)
}

ml7_analyze <- function(level_tables, level, tasks, spaces, models, repeats,
                        seed, fixed, use_smote, use_metadata_aug,
                        smote_use_mni, smote_neighbors,
                        synthetic_strategy, synthetic_multiplier,
                        synthetic_target_per_class, metadata_min_unique_sites,
                        output_dir) {
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  manifest <- ml7_space_manifest(level_tables$db_sum)
  level_id <- paste(level,collapse="_")
  write.csv(data.frame(space=rep(names(manifest),lengths(manifest)),
    feature=unlist(manifest,use.names=FALSE)),
    file.path(output_dir,"feature_manifest.csv"),row.names=FALSE)
  all_fold <- list(); all_oof <- list(); summaries <- list()
  for (task in tasks) {
    dat <- ml7_recode(level_tables$db_sum, task)
    if (anyDuplicated(dat$site_id)) stop("Duplicate site_id in task ",task)
    counts <- table(dat$taphonomic_context)
    message("  Task ", task, " (", nrow(dat), " references; ",
            paste(paste(names(counts), counts, sep = "="), collapse = ", "), ")")
    splits <- ml7_make_splits(dat[,c("site_id","taphonomic_context")],
                               task,repeats,seed)
    taskdir <- file.path(output_dir,task)
    dir.create(taskdir,recursive=TRUE,showWarnings=FALSE)
    saveRDS(splits, file.path(taskdir,"cv_splits.rds"))
    audit <- purrr::map_dfr(seq_len(nrow(splits)),function(i){
      a <- rsample::analysis(splits$splits[[i]])
      b <- rsample::assessment(splits$splits[[i]])
      tibble::tibble(fold=i,n_train=nrow(a),n_assess=nrow(b),
        train_sites=paste(a$site_id,collapse=";"),
        assessment_sites=paste(b$site_id,collapse=";"))
    })
    readr::write_csv(audit,file.path(taskdir,"fold_membership.csv"))
    # This file is written only if a training class has a single distinct
    # profile and exact-value oversampling is required. Start clean on rerun.
    fallback_audit <- file.path(taskdir, "degenerate_interpolation_fallbacks.csv")
    if (file.exists(fallback_audit)) file.remove(fallback_audit)
    for (space in spaces) for (model in models) {
      message("    ", space, " | ", model)
      fit_rows <- vector("list",nrow(splits))
      for (i in seq_len(nrow(splits))) {
        ans <- tryCatch(
          ml7_one_fit(level_tables,task,space,model,splits$splits[[i]],
            i,seed,fixed,use_smote,use_metadata_aug,smote_use_mni,
            smote_neighbors,synthetic_strategy,synthetic_multiplier,
            synthetic_target_per_class,metadata_min_unique_sites,
            audit_file = fallback_audit,
            audit_context = list(level=level_id, task=task, space=space,
                                 model=model, condition="cross_validation",
                                 fold=i)),
          error=function(e) {
            stop("Model failed: level=",level_id," task=",task,
                 " space=",space," model=",model," fold=",i,
                 "\n",conditionMessage(e),call.=FALSE)
          })
        row <- tibble::as_tibble_row(as.list(ans$metrics)) %>%
          dplyr::mutate(level=level_id,task=task,space=space,
                        model=model,fold_index=i,n_features=length(manifest[[space]]),
                        n_assess=nrow(ans$predictions),.before=1)
        fit_rows[[i]] <- row
        ans$predictions$level <- level_id
        ans$predictions$task <- task
        ans$predictions$space <- space
        ans$predictions$model <- model
        all_oof[[length(all_oof)+1L]] <- ans$predictions
      }
      folds <- dplyr::bind_rows(fit_rows)
      all_fold[[length(all_fold)+1L]] <- folds
      required_metrics <- c("balanced_accuracy", "auroc", "brier",
                            "precision", "recall", "f1", "accuracy")
      missing_metrics <- setdiff(required_metrics, names(folds))
      if (length(missing_metrics)) {
        stop("Metrics are missing after successful model fitting: ",
             paste(missing_metrics, collapse=", "), "; actual columns: ",
             paste(names(folds), collapse=", "))
      }
      if (!all(vapply(folds[required_metrics], is.numeric, logical(1)))) {
        stop("Fold metrics must all be numeric")
      }
      stats <- folds %>%
        tidyr::pivot_longer(dplyr::all_of(required_metrics),
                            names_to="metric",values_to="estimate") %>%
        dplyr::group_by(level,task,space,model,n_features,metric) %>%
        dplyr::summarise(mean=mean(estimate,na.rm=TRUE),
          sd=stats::sd(estimate,na.rm=TRUE),n_folds=dplyr::n(),.groups="drop")
      summaries[[length(summaries)+1L]] <- stats
    }
  }
  list(fold_metrics=dplyr::bind_rows(all_fold),
       predictions=dplyr::bind_rows(all_oof),
       summaries=dplyr::bind_rows(summaries))
}

ml7_publication_table <- function(summaries, output_dir, repeats,
                                  fixed,seed,use_smote,use_metadata_aug) {
  labels <- c(balanced_accuracy="Balanced accuracy",auroc="AUROC",brier="Brier",
              precision="Precision",recall="Recall",f1="F1",accuracy="Accuracy")
  formatted <- summaries %>%
    dplyr::mutate(metric=unname(labels[metric]),
      score=ifelse(metric=="Brier",
        sprintf("%.3f ± %.3f",mean,sd),
        sprintf("%.1f ± %.1f",mean*100,sd*100))) %>%
    dplyr::select(level,task,space,model,n_features,metric,score) %>%
    tidyr::pivot_wider(names_from=metric,values_from=score)
  # stable sorting and columns
  metriccols <- unname(labels)
  formatted <- formatted %>%
    dplyr::select(level,task,space,model,n_features,dplyr::all_of(metriccols)) %>%
    dplyr::arrange(level,task,match(space,c("abundance","evenness","abundance_evenness",
                                  "remnau","pmnau","engineered","combined")),model)
  header <- c("# Supervised SPR: repeated-CV performance", "",
    paste0("Seed: ",seed,"; CV: 3 folds × ",repeats,
           " repeats; one fixed configuration per model; no independent holdout."),
    paste0("Synthetic interpolation: ",use_smote,
           "; metadata centroid donors: ",use_metadata_aug,"."),
    "All entries are mean ± sample SD across held-out CV folds; repeated folds are not independent replicates.",
    "Accuracy-type measures and AUROC are on a 0–100 scale. Brier is in raw squared-probability units.",
    "Binary precision/recall/F1 and AUROC use Funerary Context as the positive class.",
    "For the 4-class task, precision/recall/F1 and AUROC use macro averaging (AUROC: one-vs-rest).",
    "Binary Brier is mean (P(F)−Y)^2; four-class Brier is mean Σ(P_k−Y_k)^2.",
    "These are fixed-setting cross-validation estimates, not nested-CV estimates or external-test scores.",
    "")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  readr::write_csv(formatted,file.path(output_dir,"publication_table_formatted.csv"))
  lines <- c(header,knitr::kable(formatted,format="pipe",align="l"))
  writeLines(lines,file.path(output_dir,"publication_table_all.md"))
  for (task in unique(formatted$task)) for (lvl in unique(formatted$level)) {
    subset <- formatted[formatted$task==task & formatted$level==lvl, , drop=FALSE]
    if (nrow(subset)) writeLines(c(header,knitr::kable(subset,format="pipe")),
      file.path(output_dir,paste0("table_",task,"_L",lvl,".md")))
  }
  invisible(formatted)
}

ml7_plot_heatmaps <- function(summaries, output_dir) {
  # Publication-ready score comparisons; colours/size are consistent within metric.
  gg <- function(metric, lbl, suffix, multiplier) {
    dat <- summaries[summaries$metric==metric, ,drop=FALSE]
    dat$space <- factor(dat$space,levels=c("abundance","evenness",
      "abundance_evenness","remnau","pmnau","engineered","combined"))
    dat$estimate <- dat$mean * multiplier
    for (task in unique(dat$task)) for (lvl in unique(dat$level)) {
      sub <- dat[dat$task==task & dat$level==lvl,,drop=FALSE]
      if (!nrow(sub)) next
      plot <- ggplot2::ggplot(sub,ggplot2::aes(x=space,y=model,fill=estimate))+
        ggplot2::geom_tile(color="white", linewidth=0.3)+
        ggplot2::geom_text(ggplot2::aes(label=ifelse(metric=="brier",
          sprintf("%.3f",estimate),sprintf("%.1f",estimate))),size=3)+
        ggplot2::labs(x="Feature space",y="Classifier",fill=lbl,
          title=paste(lbl,task,paste0("L",lvl),sep=" | "))+
        ggplot2::theme_minimal(base_size=11)+
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=35,hjust=1))
      stem <- file.path(output_dir,paste0("heatmap_",suffix,"_",task,"_L",lvl))
      ggplot2::ggsave(paste0(stem,".png"),plot,width=10,height=4.5,dpi=300)
      ggplot2::ggsave(paste0(stem,".pdf"),plot,width=10,height=4.5)
    }
  }
  gg("balanced_accuracy","Balanced accuracy (%)","BA",100)
  gg("auroc","AUROC (%)","AUROC",100)
  gg("brier","Brier score","Brier",1)
  invisible(NULL)
}

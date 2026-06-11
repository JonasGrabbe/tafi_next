# __________________ tapho_loader_wfl_features ____________________
#
# Combines:
#   * the dual-grouping structure of `tapho.loader.wfl`
#       -> computes an "_id" grouping and a "_sum" grouping and returns
#          both as a named list( db_id = ..., db_sum = ... )
#   * the feature-engineering tail of `tapho.loader`
#       -> for EACH grouping, after the MNE pipeline it appends the
#          engineered taphonomic features and the density-informed features.
#
# Net effect: the full MNE + feature-engineering pipeline is run twice
# (once per grouping) inside `process_db`, and both results are returned.


source("select_data_new.R")
source("data_loader_new.R")
source("archo_preload_2.R")
source("collapse_down_step.R")
source("indet.R")
source("meta/meta_set.R")
source("derived_features.R")
source("density_features_fixed.R")

source("convertMNEtovariables_new.R")

library(dplyr)
library(tibble)
library(readxl)
library(purrr)

tapho.loader.wfl.features <- function(
  data,
  sum                = TRUE,
  unique_level       = TRUE,
  level              = NULL,
  calc_level         = c("min", "max", "mode"),
  recalc_method      = c("min", "filter", "direct"),
  indet              = TRUE,
  indet_method       = c("up", "dist", "delet"),
  group_by           = NULL,
  select             = NULL,
  fragments_per_bone = 3,
  exclude_sections   = c("teeth")
) {

  # ── Setup ──────────────────────────────────────────────────────────────────
  meta_set(indet)
  on.exit(detach("meta_env"), add = TRUE)

  calc_level    <- match.arg(calc_level)
  indet_method  <- match.arg(indet_method)
  recalc_method <- match.arg(recalc_method)   # (the .wfl version dropped this)

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

  # ── Stage 2: Define the two group_by variants ──────────────────────────────
  group_by_sum <- if (is.null(group_by)) "site_name" else union("site_name", group_by)
  group_by_id  <- if ("site_name" %in% group_by) setdiff(group_by, "site_name") else group_by

  # ── Stage 3: Helper — aggregate + MNE + feature-engineering for a grouping ──
  process_db <- function(grp) {

    d <- select_and_group_by(
      data         = db,
      group_by     = grp,
      select       = select,
      method       = recalc_method,
      level        = level,
      unique_level = unique_level
    )

    # Lookup used later to attach section and one-skeleton count
    # to the final output elements.
    element_lookup <- elements_data %>%
      dplyr::transmute(
        element_key        = gsub(" ", ".", ElementsOne),
        One_Skeleton_count = Bone_Count,
        Skeleton_Section   = Skeleton_Section
      )

    # ── Compute MNE-derived variables for each row ───────────────────────────
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

      if (nrow(bone_data) == 0) return(NULL)

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
      if (is.list(out_elements)) out_elements <- out_elements[[1]]

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

    # ── Unpack MNE results into top-level list-columns ───────────────────────
    out <- d %>%
      dplyr::select(-c(count, element)) %>%
      dplyr::mutate(mne_features = results_list) %>%
      dplyr::mutate(
        element                           = lapply(mne_features, `[[`, "element_2"),
        One_Skeleton_count                = lapply(mne_features, `[[`, "One_Skeleton_count"),
        Skeleton_Section                  = lapply(mne_features, `[[`, "Skeleton_Section"),
        MNE                               = lapply(mne_features, `[[`, "MNE"),
        MNAU                              = lapply(mne_features, `[[`, "MNAU"),
        MAU                               = lapply(mne_features, `[[`, "MAU"),
        PercentMAU                        = lapply(mne_features, `[[`, "PercentMAU"),
        PercentMNAU                       = lapply(mne_features, `[[`, "PercentMNAU"),
        ReMNAU                            = lapply(mne_features, `[[`, "ReMNAU"),
        AcReMNAU                          = lapply(mne_features, `[[`, "AcReMNAU"),
        CumulativeReMNAU                  = lapply(mne_features, `[[`, "CumulativeReMNAU"),
        ProportionalRep_MNIbased          = lapply(mne_features, `[[`, "ProportionalRep_MNIbased"),
        MNI_used_assemblage               = lapply(mne_features, `[[`, "MNI_used_assemblage"),
        Mean_MNAU_assemblage              = lapply(mne_features, `[[`, "Mean_MNAU_assemblage"),
        Mean_part_abundance               = lapply(mne_features, `[[`, "Mean_part_abundance"),
        OverallProportionalRep_assemblage = lapply(mne_features, `[[`, "OverallProportionalRep_assemblage")
      ) %>%
      dplyr::select(-mne_features)

    # ── Engineered taphonomic features ───────────────────────────────────────
    derived_features <- purrr::pmap_dfr(
      list(
        element            = out$element,
        One_Skeleton_count = out$One_Skeleton_count,
        Skeleton_Section   = out$Skeleton_Section,
        MNE                = out$MNE,
        MNI_used           = out$MNI_used_assemblage
      ),
      make_tapho_derived_features
    )

    # ── Density-informed engineered features ─────────────────────────────────
    density_features <- purrr::pmap_dfr(
      list(
        element            = out$element,
        One_Skeleton_count = out$One_Skeleton_count,
        Skeleton_Section   = out$Skeleton_Section,
        MNE                = out$MNE,
        MNI_used           = out$MNI_used_assemblage,
        PercentMNAU        = out$PercentMNAU,
        MNAU               = out$MNAU
      ),
      make_tapho_density_features
    )

    # ── Bind engineered features onto the MNE table for this grouping ────────
    dplyr::bind_cols(out, derived_features, density_features)
  }

  # ── Stage 4: Run the full pipeline for both groupings ──────────────────────
  db_sum <- process_db(group_by_sum)
  db_id  <- process_db(group_by_id)

  # ── Stage 5: Return both as a named list ───────────────────────────────────
  list(db_id = db_id, db_sum = db_sum)
}

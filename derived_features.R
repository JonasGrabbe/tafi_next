make_tapho_derived_features <- function(
  element,
  One_Skeleton_count,
  Skeleton_Section,
  MNE,
  MNI_used,
  epsilon = 1e-6
) {

  sections_order <- c(
    "cranium", "teeth", "hyoid", "vertebrae", "thorax", "shoulder",
    "arms", "hands", "pelvic", "legs", "foot"
  )

  groups <- list(
    cranial_complex       = c("cranium", "teeth", "hyoid"),
    axial_main            = c("vertebrae", "thorax"),
    full_axial            = c("cranium", "teeth", "hyoid", "vertebrae", "thorax"),
    postcranial_total     = c("vertebrae", "thorax", "shoulder", "arms", "hands", "pelvic", "legs", "foot"),
    column_trunk          = c("vertebrae", "thorax", "pelvic"),
    appendicular          = c("shoulder", "arms", "hands", "pelic", "legs", "foot"),
    upper_appendicular    = c("shoulder", "arms", "hands"),
    lower_appendicular    = c("pelvic", "legs", "foot"),
    girdle                = c("shoulder", "pelvic"),
    limb_shaft_proxy      = c("arms", "legs"),
    distal_appendicular   = c("hands", "foot"),
    proximal_appendicular = c("shoulder", "arms", "pelvic", "legs")
  )

  # Correct possible typo before use
  groups$appendicular <- c("shoulder", "arms", "hands", "pelvic", "legs", "foot")

  as_first_number <- function(x) {
    x <- suppressWarnings(as.numeric(unlist(x)))
    x <- x[is.finite(x) & !is.na(x)]
    if (length(x) == 0) return(NA_real_)
    x[1]
  }

  sum_ignore_missing <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) return(NA_real_)
    sum(x, na.rm = TRUE)
  }

  mean_ignore_missing <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }

  safe_ratio <- function(num, den) {
    num <- as.numeric(num)
    den <- as.numeric(den)

    if (length(num) == 0 || length(den) == 0) return(NA_real_)
    if (is.na(num) || is.na(den)) return(NA_real_)

    (num + epsilon) / (den + epsilon)
  }

  harmonic_mean <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[!is.na(x)]

    if (length(x) == 0) return(NA_real_)
    if (all(x == 0)) return(0)

    length(x) / sum(1 / (x + epsilon))
  }

  evenness <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[!is.na(x)]

    k <- length(x)
    if (k < 2) return(NA_real_)

    total <- sum(x, na.rm = TRUE)
    if (!is.finite(total) || total <= 0) return(NA_real_)

    p <- x / total
    p <- p[p > 0]

    H <- -sum(p * log(p))
    H / log(k)
  }

  dominance <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x <- x[!is.na(x)]

    total <- sum(x, na.rm = TRUE)
    if (length(x) == 0 || !is.finite(total) || total <= 0) return(NA_real_)

    max(x / total, na.rm = TRUE)
  }

  presence_count <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    sum(x > 0, na.rm = TRUE)
  }

  MNI_used <- as_first_number(MNI_used)

  df <- tibble::tibble(
    element            = as.character(unlist(element)),
    One_Skeleton_count = as.numeric(unlist(One_Skeleton_count)),
    Skeleton_Section   = as.character(unlist(Skeleton_Section)),
    MNE                = as.numeric(unlist(MNE))
  ) %>%
    dplyr::filter(!is.na(Skeleton_Section))

  if (nrow(df) == 0 || is.na(MNI_used) || MNI_used <= 0) {
    return(tibble::tibble())
  }

  section_df <- tibble::tibble(Skeleton_Section = sections_order) %>%
    dplyr::left_join(
      df %>%
        dplyr::group_by(Skeleton_Section) %>%
        dplyr::summarise(
          mne_sum        = sum(MNE, na.rm = TRUE),
          one_count_sum  = sum(One_Skeleton_count, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "Skeleton_Section"
    ) %>%
    dplyr::mutate(
      mne_sum       = dplyr::if_else(is.na(mne_sum), NA_real_, mne_sum),
      one_count_sum = dplyr::if_else(is.na(one_count_sum) | one_count_sum <= 0, NA_real_, one_count_sum),
      mnau_prop     = (mne_sum / one_count_sum) / MNI_used
    )

  mne_section <- stats::setNames(section_df$mne_sum, section_df$Skeleton_Section)
  mnau_section <- stats::setNames(section_df$mnau_prop, section_df$Skeleton_Section)

  section_mne_features <- stats::setNames(
    as.list(mne_section),
    paste0("mne_", names(mne_section))
  )

  section_mnau_features <- stats::setNames(
    as.list(mnau_section),
    paste0("mnau_", names(mnau_section))
  )

  composite_mne <- purrr::imap_dbl(groups, function(secs, nm) {
    sum_ignore_missing(mne_section[secs])
  })

  composite_mnau <- purrr::imap_dbl(groups, function(secs, nm) {
    mean_ignore_missing(mnau_section[secs])
  })

  composite_mne_features <- stats::setNames(
    as.list(composite_mne),
    paste0("mne_", names(composite_mne), "_sum")
  )

  composite_mnau_features <- stats::setNames(
    as.list(composite_mnau),
    paste0("mnau_", names(composite_mnau), "_mean")
  )

  mnau_all_sections_mean <- mean_ignore_missing(mnau_section[sections_order])
  mne_all_sections_sum   <- sum_ignore_missing(mne_section[sections_order])

  marrow_poor_mnau <- mean_ignore_missing(mnau_section[c("hands", "foot", "thorax")])
  marrow_poor_mne  <- sum_ignore_missing(mne_section[c("hands", "foot", "thorax")])

  ratio_features <- list(
    mnau_appendicular_to_axial_main_ratio =
      safe_ratio(composite_mnau["appendicular"], composite_mnau["axial_main"]),

    mnau_column_trunk_to_appendicular_ratio =
      safe_ratio(composite_mnau["column_trunk"], composite_mnau["appendicular"]),

    mnau_upper_to_lower_appendicular_ratio =
      safe_ratio(composite_mnau["upper_appendicular"], composite_mnau["lower_appendicular"]),

    mnau_hand_to_foot_ratio =
      safe_ratio(mnau_section["hands"], mnau_section["foot"]),

    mnau_cranial_to_postcranial_ratio =
      safe_ratio(composite_mnau["cranial_complex"], composite_mnau["postcranial_total"]),

    mnau_cranial_to_axial_main_ratio =
      safe_ratio(composite_mnau["cranial_complex"], composite_mnau["axial_main"]),

    mnau_cranial_to_limb_shaft_proxy_ratio =
      safe_ratio(composite_mnau["cranial_complex"], composite_mnau["limb_shaft_proxy"]),

    mnau_vertebrae_to_thorax_ratio =
      safe_ratio(mnau_section["vertebrae"], mnau_section["thorax"]),

    mnau_arms_to_legs_ratio =
      safe_ratio(mnau_section["arms"], mnau_section["legs"]),

    mnau_distal_to_proximal_appendicular_ratio =
      safe_ratio(composite_mnau["distal_appendicular"], composite_mnau["proximal_appendicular"]),

    mnau_girdle_to_limb_shaft_proxy_ratio =
      safe_ratio(composite_mnau["girdle"], composite_mnau["limb_shaft_proxy"]),

    mnau_distal_to_limb_shaft_proxy_ratio =
      safe_ratio(composite_mnau["distal_appendicular"], composite_mnau["limb_shaft_proxy"]),

    mnau_limb_shaft_proxy_to_axial_main_ratio =
      safe_ratio(composite_mnau["limb_shaft_proxy"], composite_mnau["axial_main"]),

    mnau_limb_shaft_proxy_to_column_trunk_ratio =
      safe_ratio(composite_mnau["limb_shaft_proxy"], composite_mnau["column_trunk"]),

    mnau_marrow_rich_proxy_to_marrow_poor_proxy_ratio =
      safe_ratio(composite_mnau["limb_shaft_proxy"], marrow_poor_mnau),

    mnau_cranial_proportion =
      safe_ratio(composite_mnau["cranial_complex"], mnau_all_sections_mean),

    mnau_limb_shaft_proxy_proportion =
      safe_ratio(composite_mnau["limb_shaft_proxy"], mnau_all_sections_mean),

    mne_appendicular_to_axial_main_ratio =
      safe_ratio(composite_mne["appendicular"], composite_mne["axial_main"]),

    mne_column_trunk_to_appendicular_ratio =
      safe_ratio(composite_mne["column_trunk"], composite_mne["appendicular"]),

    mne_upper_to_lower_appendicular_ratio =
      safe_ratio(composite_mne["upper_appendicular"], composite_mne["lower_appendicular"]),

    mne_hand_to_foot_ratio =
      safe_ratio(mne_section["hands"], mne_section["foot"]),

    mne_cranial_to_postcranial_ratio =
      safe_ratio(composite_mne["cranial_complex"], composite_mne["postcranial_total"]),

    mne_cranial_to_axial_main_ratio =
      safe_ratio(composite_mne["cranial_complex"], composite_mne["axial_main"]),

    mne_cranial_to_limb_shaft_proxy_ratio =
      safe_ratio(composite_mne["cranial_complex"], composite_mne["limb_shaft_proxy"]),

    mne_vertebrae_to_thorax_ratio =
      safe_ratio(mne_section["vertebrae"], mne_section["thorax"]),

    mne_arms_to_legs_ratio =
      safe_ratio(mne_section["arms"], mne_section["legs"]),

    mne_distal_to_proximal_appendicular_ratio =
      safe_ratio(composite_mne["distal_appendicular"], composite_mne["proximal_appendicular"]),

    mne_girdle_to_limb_shaft_proxy_ratio =
      safe_ratio(composite_mne["girdle"], composite_mne["limb_shaft_proxy"]),

    mne_distal_to_limb_shaft_proxy_ratio =
      safe_ratio(composite_mne["distal_appendicular"], composite_mne["limb_shaft_proxy"]),

    mne_limb_shaft_proxy_to_axial_main_ratio =
      safe_ratio(composite_mne["limb_shaft_proxy"], composite_mne["axial_main"]),

    mne_limb_shaft_proxy_to_column_trunk_ratio =
      safe_ratio(composite_mne["limb_shaft_proxy"], composite_mne["column_trunk"]),

    mne_marrow_rich_proxy_to_marrow_poor_proxy_ratio =
      safe_ratio(composite_mne["limb_shaft_proxy"], marrow_poor_mne),

    mne_cranial_proportion =
      safe_ratio(composite_mne["cranial_complex"], mne_all_sections_sum),

    mne_limb_shaft_proxy_proportion =
      safe_ratio(composite_mne["limb_shaft_proxy"], mne_all_sections_sum)
  )

  continuity_features <- list(
    upper_appendicular_continuity_hm =
      harmonic_mean(mnau_section[c("shoulder", "arms", "hands")]),

    lower_appendicular_continuity_hm =
      harmonic_mean(mnau_section[c("pelvic", "legs", "foot")]),

    appendicular_continuity_hm =
      harmonic_mean(mnau_section[c("shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    axial_main_continuity_hm =
      harmonic_mean(mnau_section[c("vertebrae", "thorax")]),

    column_trunk_continuity_hm =
      harmonic_mean(mnau_section[c("vertebrae", "thorax", "pelvic")]),

    full_axial_continuity_hm =
      harmonic_mean(c(composite_mnau["cranial_complex"], mnau_section[c("vertebrae", "thorax")])),

    postcranial_continuity_hm =
      harmonic_mean(mnau_section[c("vertebrae", "thorax", "shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    proximal_appendicular_continuity_hm =
      harmonic_mean(mnau_section[c("shoulder", "arms", "pelvic", "legs")]),

    distal_appendicular_continuity_hm =
      harmonic_mean(mnau_section[c("hands", "foot")]),

    whole_body_continuity_hm =
      harmonic_mean(mnau_section[sections_order])
  )

  evenness_features <- list(
    all_sections_evenness =
      evenness(mnau_section[sections_order]),

    postcranial_evenness =
      evenness(mnau_section[c("vertebrae", "thorax", "shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    appendicular_evenness =
      evenness(mnau_section[c("shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    column_trunk_evenness =
      evenness(mnau_section[c("vertebrae", "thorax", "pelvic")]),

    upper_appendicular_evenness =
      evenness(mnau_section[c("shoulder", "arms", "hands")]),

    lower_appendicular_evenness =
      evenness(mnau_section[c("pelvic", "legs", "foot")]),

    distal_appendicular_evenness =
      evenness(mnau_section[c("hands", "foot")])
  )

  dominance_features <- list(
    all_sections_dominance =
      dominance(mnau_section[sections_order]),

    postcranial_dominance =
      dominance(mnau_section[c("vertebrae", "thorax", "shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    appendicular_dominance =
      dominance(mnau_section[c("shoulder", "arms", "hands", "pelvic", "legs", "foot")])
  )

  share_features <- list(
    mnau_cranial_share =
      safe_ratio(composite_mnau["cranial_complex"], mnau_all_sections_mean),

    mnau_postcranial_share =
      safe_ratio(composite_mnau["postcranial_total"], mnau_all_sections_mean),

    mnau_limb_shaft_proxy_share =
      safe_ratio(composite_mnau["limb_shaft_proxy"], mnau_all_sections_mean),

    mnau_distal_appendicular_share =
      safe_ratio(composite_mnau["distal_appendicular"], composite_mnau["appendicular"]),

    mnau_girdle_share_of_appendicular =
      safe_ratio(composite_mnau["girdle"], composite_mnau["appendicular"]),

    mne_cranial_share =
      safe_ratio(composite_mne["cranial_complex"], mne_all_sections_sum),

    mne_postcranial_share =
      safe_ratio(composite_mne["postcranial_total"], mne_all_sections_sum),

    mne_limb_shaft_proxy_share =
      safe_ratio(composite_mne["limb_shaft_proxy"], mne_all_sections_sum),

    mne_distal_appendicular_share =
      safe_ratio(composite_mne["distal_appendicular"], composite_mne["appendicular"]),

    mne_girdle_share_of_appendicular =
      safe_ratio(composite_mne["girdle"], composite_mne["appendicular"])
  )

  richness_features <- list(
    section_richness_all =
      presence_count(mne_section[sections_order]),

    section_richness_postcranial =
      presence_count(mne_section[c("vertebrae", "thorax", "shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    section_richness_appendicular =
      presence_count(mne_section[c("shoulder", "arms", "hands", "pelvic", "legs", "foot")]),

    cranial_complex_richness =
      presence_count(mne_section[c("cranium", "teeth", "hyoid")])
  )

  tibble::as_tibble(
    c(
      section_mne_features,
      section_mnau_features,
      composite_mne_features,
      composite_mnau_features,
      ratio_features,
      continuity_features,
      evenness_features,
      dominance_features,
      share_features,
      richness_features
    )
  )
}
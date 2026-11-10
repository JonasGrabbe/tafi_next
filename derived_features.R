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

  # ---- NEW helpers (engineered feature expansion) -------------------------

  # Signed difference of two scalars, NA-safe.
  safe_diff <- function(a, b) {
    a <- as.numeric(a); b <- as.numeric(b)
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    if (is.na(a) || is.na(b)) return(NA_real_)
    a - b
  }

  # log( (num+eps)/(den+eps) ), NA-safe. Stable alternative to plain ratios.
  safe_logratio <- function(num, den) {
    num <- as.numeric(num); den <- as.numeric(den)
    if (length(num) == 0 || length(den) == 0) return(NA_real_)
    if (is.na(num) || is.na(den)) return(NA_real_)
    log((num + epsilon) / (den + epsilon))
  }

  max_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    max(x)
  }

  sum_finite <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    sum(x)
  }

  # Proportional shares p_i = x_i / sum(x); NULL if no usable total.
  proportions <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    total <- sum(x)
    if (length(x) == 0 || !is.finite(total) || total <= 0) return(NULL)
    x / total
  }

  # Berger-Parker gap: largest share minus second-largest share.
  dominance_gap <- function(x) {
    p <- proportions(x)
    if (is.null(p) || length(p) < 2) return(NA_real_)
    ps <- sort(p, decreasing = TRUE)
    ps[1] - ps[2]
  }

  # McNaughton two-section dominance: sum of two largest shares.
  top2_share <- function(x) {
    p <- proportions(x)
    if (is.null(p)) return(NA_real_)
    ps <- sort(p, decreasing = TRUE)
    sum(ps[seq_len(min(2, length(ps)))])
  }

  # Simpson dominance lambda = sum(p^2).
  simpson_dominance <- function(x) {
    p <- proportions(x)
    if (is.null(p)) return(NA_real_)
    sum(p^2)
  }

  # Inverse Simpson = 1 / lambda  (effective number of sections, Simpson).
  inverse_simpson <- function(x) {
    d <- simpson_dominance(x)
    if (is.na(d) || d <= 0) return(NA_real_)
    1 / d
  }

  # Hill number of order 1 = exp(Shannon H)  (effective number of sections).
  hill1 <- function(x) {
    p <- proportions(x)
    if (is.null(p)) return(NA_real_)
    p <- p[p > 0]
    if (length(p) < 1) return(NA_real_)
    exp(-sum(p * log(p)))
  }

  # Gini coefficient of inequality (0 = perfectly even, ->1 = concentrated).
  gini <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    n <- length(x)
    if (n < 2 || sum(x) <= 0) return(NA_real_)
    x <- sort(x)
    sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
  }

  # Name of the section holding the maximum value.
  dominant_section <- function(v) {
    v <- v[is.finite(v)]
    if (length(v) == 0) return(NA_character_)
    names(v)[which.max(v)]
  }

  # Smith & Wilson Evar evenness (independent of richness).
  evar_smith_wilson <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x) & x > 0]
    S <- length(x)
    if (S < 2) return(NA_real_)
    lx <- log(x)
    v <- sum((lx - mean(lx))^2) / S
    1 - (2 / pi) * atan(v)
  }

  # Heip evenness = (exp(H) - 1) / (k - 1).
  heip_evenness <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    total <- sum(x)
    k <- length(x)
    if (k < 2 || !is.finite(total) || total <= 0) return(NA_real_)
    p <- x / total; p <- p[p > 0]
    H <- -sum(p * log(p))
    (exp(H) - 1) / (k - 1)
  }

  # Simpson evenness = (1 / lambda) / k.
  simpson_evenness <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    total <- sum(x); k <- length(x)
    if (k < 1 || !is.finite(total) || total <= 0) return(NA_real_)
    p <- x / total
    D <- sum(p^2)
    if (D <= 0) return(NA_real_)
    (1 / D) / k
  }

  # Pielou evenness using only present (>0) sections as k.
  evenness_present_only <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x) & x > 0]
    kp <- length(x)
    if (kp < 2) return(NA_real_)
    p <- x / sum(x)
    H <- -sum(p * log(p))
    H / log(kp)
  }

  # Coefficient of variation.
  cv <- function(x) {
    x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
    if (length(x) < 2) return(NA_real_)
    m <- mean(x)
    if (m == 0) return(NA_real_)
    stats::sd(x) / m
  }

  # Count of sections at or above a representation threshold.
  completeness_count <- function(x, thr) {
    x <- suppressWarnings(as.numeric(x))
    sum(x >= thr, na.rm = TRUE)
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

  # ---- NEW: MNAU contrast features (differences + log-ratios) -------------
  # Built on the regions that survived feature elimination. Differences and
  # log-ratios are used because every plain ratio was eliminated in the RFE run;
  # contrasts on the surviving MNAU scale behave more stably. The three named
  # plain ratios are kept because they were explicitly proposed.

  get_mnau <- function(key) {
    if (key %in% names(composite_mnau)) return(unname(composite_mnau[key]))
    if (key %in% names(mnau_section))   return(unname(mnau_section[key]))
    NA_real_
  }

  contrast_pairs <- list(
    lower_upper_appendicular   = c("lower_appendicular",   "upper_appendicular"),
    proximal_distal_appendicular = c("proximal_appendicular", "distal_appendicular"),
    legs_appendicular          = c("legs",                 "appendicular"),
    legs_postcranial           = c("legs",                 "postcranial_total"),
    girdle_appendicular        = c("girdle",               "appendicular"),
    appendicular_axial_main    = c("appendicular",         "axial_main"),
    postcranial_cranial        = c("postcranial_total",    "cranial_complex"),
    limb_shaft_axial_main      = c("limb_shaft_proxy",     "axial_main"),
    arms_legs                  = c("arms",                 "legs")
  )

  contrast_features <- list()
  for (nm in names(contrast_pairs)) {
    a <- get_mnau(contrast_pairs[[nm]][1])
    b <- get_mnau(contrast_pairs[[nm]][2])
    contrast_features[[paste0("mnau_", nm, "_diff")]]     <- safe_diff(a, b)
    contrast_features[[paste0("mnau_", nm, "_logratio")]] <- safe_logratio(a, b)
  }
  contrast_features$mnau_legs_to_appendicular_ratio <-
    safe_ratio(mnau_section["legs"], composite_mnau["appendicular"])
  contrast_features$mnau_legs_to_postcranial_ratio <-
    safe_ratio(mnau_section["legs"], composite_mnau["postcranial_total"])
  contrast_features$mnau_girdle_to_appendicular_ratio <-
    safe_ratio(composite_mnau["girdle"], composite_mnau["appendicular"])

  # ---- NEW: dominance / concentration features ----------------------------

  dom_argmax <- dominant_section(mnau_section[sections_order])

  dominance_features_new <- list(
    dominance_gap_all =
      dominance_gap(mnau_section[sections_order]),

    mnau_top2_share_mcnaughton =
      top2_share(mnau_section[sections_order]),

    mnau_simpson_dominance =
      simpson_dominance(mnau_section[sections_order]),

    mnau_inverse_simpson =
      inverse_simpson(mnau_section[sections_order]),

    mnau_gini =
      gini(mnau_section[sections_order]),

    dominant_section_is_legs =
      if (is.na(dom_argmax)) NA_real_ else as.numeric(dom_argmax == "legs"),

    dominant_section_is_appendicular =
      if (is.na(dom_argmax)) NA_real_ else as.numeric(dom_argmax %in% groups$appendicular),

    dominant_section_is_axial =
      if (is.na(dom_argmax)) NA_real_ else as.numeric(dom_argmax %in% groups$full_axial),

    appendicular_dominance_refined =
      dominance(mnau_section[groups$appendicular]),

    postcranial_dominance_refined =
      dominance(mnau_section[groups$postcranial_total]),

    lower_appendicular_dominance =
      dominance(mnau_section[groups$lower_appendicular]),

    upper_appendicular_dominance =
      dominance(mnau_section[groups$upper_appendicular]),

    dominant_appendicular_share =
      safe_ratio(max_finite(mnau_section[groups$appendicular]),
                 sum_finite(mnau_section[sections_order])),

    dominant_axial_share =
      safe_ratio(max_finite(mnau_section[groups$full_axial]),
                 sum_finite(mnau_section[sections_order])),

    mnau_dominance_of_region_means =
      dominance(composite_mnau)
  )

  # ---- NEW: diversity / evenness features ---------------------------------

  diversity_evenness_features_new <- list(
    mnau_evar_smith_wilson =
      evar_smith_wilson(mnau_section[sections_order]),

    mnau_heip_evenness =
      heip_evenness(mnau_section[sections_order]),

    mnau_simpson_evenness =
      simpson_evenness(mnau_section[sections_order]),

    mnau_evenness_present_only =
      evenness_present_only(mnau_section[sections_order]),

    appendicular_evenness_refined =
      evenness(mnau_section[groups$appendicular]),

    postcranial_evenness_refined =
      evenness(mnau_section[groups$postcranial_total]),

    lower_appendicular_evenness_refined =
      evenness(mnau_section[groups$lower_appendicular]),

    upper_appendicular_evenness_refined =
      evenness(mnau_section[groups$upper_appendicular]),

    upper_lower_evenness =
      evenness(c(composite_mnau["upper_appendicular"], composite_mnau["lower_appendicular"])),

    proximal_distal_evenness =
      evenness(c(composite_mnau["proximal_appendicular"], composite_mnau["distal_appendicular"])),

    axial_appendicular_evenness =
      evenness(c(composite_mnau["axial_main"], composite_mnau["appendicular"])),

    mnau_region_mean_evenness =
      evenness(composite_mnau),

    effective_sections_hill1 =
      hill1(mnau_section[sections_order]),

    mnau_cv =
      cv(mnau_section[sections_order]),

    mnau_evenness_minus_dominance =
      safe_diff(evenness(mnau_section[sections_order]),
                dominance(mnau_section[sections_order]))
  )

  # ---- NEW: completeness features (graded replacement for richness) -------

  completeness_features <- list(
    section_completeness_025_all =
      completeness_count(mnau_section[sections_order], 0.25),

    section_completeness_050_all =
      completeness_count(mnau_section[sections_order], 0.50),

    section_completeness_075_all =
      completeness_count(mnau_section[sections_order], 0.75),

    postcranial_completeness_050 =
      completeness_count(mnau_section[groups$postcranial_total], 0.50),

    appendicular_completeness_050 =
      completeness_count(mnau_section[groups$appendicular], 0.50)
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
      richness_features,
      contrast_features,
      dominance_features_new,
      diversity_evenness_features_new,
      completeness_features
    )
  )
}
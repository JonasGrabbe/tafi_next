make_tapho_density_features <- function(
  element,
  One_Skeleton_count,
  Skeleton_Section,
  MNE,
  MNI_used,
  PercentMNAU = NULL,
  MNAU = NULL,
  epsilon = 1e-6,
  fragile_threshold = 0,
  teeth_density_placeholder = 0.721100,
  audit_tolerance = 1e-6
) {

  density_feature_names <- c(
    "density_weighted_mean_mnau",
    "density_weighted_mean_mne",
    "density_representation_spearman",
    "density_representation_slope",
    "high_density_to_low_density_ratio",
    "high_density_share",
    "low_density_share",
    "low_density_survival_index",
    "fragile_section_survival_count",
    "density_evenness",
    "density_dominance_index",
    "density_adjusted_cranial_postcranial_contrast",
    "density_adjusted_limb_trunk_contrast",
    "density_adjusted_appendicular_axial_contrast",
    "density_adjusted_distal_proximal_contrast",
    "density_corrected_marrow_proxy_contrast",
    "high_minus_low_density_survival",
    "high_to_low_density_logratio",
    "low_density_to_postcranial_ratio",
    "low_density_deficit_from_appendicular",
    "low_density_axial_survival",
    "fragile_survival_fraction",
    "low_density_presence_050",
    "vertebrae_thorax_low_density_balance",
    "density_representation_r2",
    "density_residual_sd",
    "density_residual_appendicular",
    "density_residual_cranial",
    "density_residual_axial",
    "dominant_low_density_flag",
    "mnau_dominant_section_density",
    "mnau_low_density_dominance"
  )

  empty_density_features <- function() {
    x <- as.list(rep(NA_real_, length(density_feature_names)))
    names(x) <- density_feature_names
    tibble::as_tibble(x)
  }

  as_chr <- function(x) {
    if (is.null(x)) return(character(0))
    as.character(unlist(x, use.names = FALSE))
  }

  as_num <- function(x) {
    if (is.null(x)) return(numeric(0))
    suppressWarnings(as.numeric(unlist(x, use.names = FALSE)))
  }

  first_positive <- function(x) {
    x <- as_num(x)
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0) return(NA_real_)
    x[1]
  }

  recycle_or_na <- function(x, n) {
    x <- as_num(x)
    if (length(x) == n) return(x)
    rep(NA_real_, n)
  }

  sections_order <- c(
    "cranium", "teeth", "hyoid", "vertebrae", "thorax", "shoulder",
    "arms", "hands", "pelvic", "legs", "foot"
  )

  section_density <- c(
    cranium   = 0.721100,
    teeth     = teeth_density_placeholder,
    hyoid     = 0.252500,
    vertebrae = 0.154900,
    thorax    = 0.214000,
    shoulder  = 0.190450,
    arms      = 0.576000,
    hands     = 0.489167,
    pelvic    = 0.326500,
    legs      = 0.451025,
    foot      = 0.479167
  )[sections_order]

  density_bins <- list(
    low  = c("vertebrae", "shoulder", "thorax", "hyoid"),
    mid  = c("pelvic", "legs", "foot"),
    high = c("hands", "arms", "cranium", "teeth")
  )

  groups <- list(
    cranial_complex       = c("cranium", "teeth", "hyoid"),
    axial_main            = c("vertebrae", "thorax"),
    full_axial            = c("cranium", "teeth", "hyoid", "vertebrae", "thorax"),
    postcranial_total     = c("vertebrae", "thorax", "shoulder", "arms", "hands", "pelvic", "legs", "foot"),
    column_trunk          = c("vertebrae", "thorax", "pelvic"),
    appendicular          = c("shoulder", "arms", "hands", "pelvic", "legs", "foot"),
    upper_appendicular    = c("shoulder", "arms", "hands"),
    lower_appendicular    = c("pelvic", "legs", "foot"),
    girdle                = c("shoulder", "pelvic"),
    limb_shaft_proxy      = c("arms", "legs"),
    distal_appendicular   = c("hands", "foot"),
    proximal_appendicular = c("shoulder", "arms", "pelvic", "legs")
  )

  sec <- as_chr(Skeleton_Section)
  mne <- as_num(MNE)
  one_count <- as_num(One_Skeleton_count)
  n <- length(sec)

  if (n == 0 || length(mne) != n) {
    return(empty_density_features())
  }

  if (length(one_count) != n) {
    one_count <- rep(NA_real_, n)
  }

  mni <- first_positive(MNI_used)

  df <- tibble::tibble(
    Skeleton_Section   = sec,
    MNE                = mne,
    One_Skeleton_count = one_count,
    MNAU_supplied      = recycle_or_na(MNAU, n),
    PercentMNAU_supplied = recycle_or_na(PercentMNAU, n)
  ) %>%
    dplyr::filter(!is.na(Skeleton_Section), Skeleton_Section %in% sections_order)

  if (nrow(df) == 0 || !is.finite(mni) || mni <= 0) {
    return(empty_density_features())
  }

  # Element-level audit values. These are not averaged for the density features;
  # they are used only to check that convertMNEtoVariables_v3 is consistent with
  # MNE / One_Skeleton_count and MNI.
  df <- df %>%
    dplyr::mutate(
      element_mnau_from_mne = dplyr::if_else(
        is.finite(One_Skeleton_count) & One_Skeleton_count > 0,
        MNE / One_Skeleton_count,
        NA_real_
      ),
      element_prop_from_mne = element_mnau_from_mne / mni,
      element_percent_from_mne = 100 * element_prop_from_mne,
      mnau_mismatch = is.finite(MNAU_supplied) &
        is.finite(element_mnau_from_mne) &
        abs(MNAU_supplied - element_mnau_from_mne) > audit_tolerance,
      percent_mnau_mismatch = is.finite(PercentMNAU_supplied) &
        is.finite(element_percent_from_mne) &
        abs(PercentMNAU_supplied - element_percent_from_mne) > audit_tolerance
    )

  if (any(df$mnau_mismatch, na.rm = TRUE)) {
    warning("Supplied MNAU differs from MNE / One_Skeleton_count for at least one element.", call. = FALSE)
  }

  if (any(df$percent_mnau_mismatch, na.rm = TRUE)) {
    warning("Supplied PercentMNAU differs from 100 * MNE / One_Skeleton_count / MNI_used for at least one element.", call. = FALSE)
  }

  # Section-level representation: merge section bones first, then calculate.
  # This matches make_tapho_derived_features():
  # section_prop = (sum(MNE) / sum(One_Skeleton_count)) / MNI_used
  section_df <- tibble::tibble(Skeleton_Section = sections_order) %>%
    dplyr::left_join(
      df %>%
        dplyr::group_by(Skeleton_Section) %>%
        dplyr::summarise(
          mne_sum = if (all(is.na(MNE))) NA_real_ else sum(MNE, na.rm = TRUE),
          one_count_sum = if (all(is.na(One_Skeleton_count))) NA_real_ else sum(One_Skeleton_count, na.rm = TRUE),
          .groups = "drop"
        ),
      by = "Skeleton_Section"
    ) %>%
    dplyr::mutate(
      one_count_sum = dplyr::if_else(
        is.na(one_count_sum) | one_count_sum <= 0,
        NA_real_,
        one_count_sum
      ),
      section_prop = (mne_sum / one_count_sum) / mni
    )

  section_rep <- stats::setNames(section_df$section_prop, section_df$Skeleton_Section)
  section_mne <- stats::setNames(section_df$mne_sum, section_df$Skeleton_Section)

  safe_weighted_density <- function(weights) {
    ok <- is.finite(weights) & is.finite(section_density) & weights > 0
    if (!any(ok)) return(NA_real_)
    sum(weights[ok] * section_density[ok], na.rm = TRUE) / sum(weights[ok], na.rm = TRUE)
  }

  safe_share <- function(x, total) {
    if (!is.finite(total) || total <= epsilon) return(NA_real_)
    x / total
  }

  safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) return(NA_real_)
    mean(x)
  }

  d <- as.numeric(section_density)
  r <- as.numeric(section_rep)
  names(r) <- sections_order
  regression_ok <- is.finite(d) & is.finite(r)

  if (sum(regression_ok) < 3 || stats::sd(r[regression_ok]) <= epsilon) {
    density_representation_spearman <- NA_real_
  } else {
    density_representation_spearman <- suppressWarnings(
      stats::cor(d[regression_ok], r[regression_ok], method = "spearman")
    )
  }

  residuals <- stats::setNames(rep(NA_real_, length(sections_order)), sections_order)

  if (sum(regression_ok) < 2) {
    density_representation_slope <- NA_real_
  } else {
    d_ok <- d[regression_ok]
    r_ok <- r[regression_ok]
    d_bar <- mean(d_ok)
    r_bar <- mean(r_ok)
    slope_denominator <- sum((d_ok - d_bar)^2)

    if (!is.finite(slope_denominator) || slope_denominator <= epsilon) {
      density_representation_slope <- NA_real_
    } else {
      density_representation_slope <- sum((d_ok - d_bar) * (r_ok - r_bar)) / slope_denominator
      intercept <- r_bar - density_representation_slope * d_bar
      residuals[regression_ok] <- r_ok - (intercept + density_representation_slope * d_ok)
    }
  }

  low_sections  <- density_bins$low
  mid_sections  <- density_bins$mid
  high_sections <- density_bins$high

  total_rep <- sum(section_rep, na.rm = TRUE)
  if (!is.finite(total_rep) || total_rep <= epsilon) {
    total_rep <- NA_real_
  }

  high_total <- sum(section_rep[high_sections], na.rm = TRUE)
  mid_total  <- sum(section_rep[mid_sections], na.rm = TRUE)
  low_total  <- sum(section_rep[low_sections], na.rm = TRUE)

  density_bin_totals <- c(low = low_total, mid = mid_total, high = high_total)

  if (sum(density_bin_totals, na.rm = TRUE) <= epsilon) {
    density_evenness <- NA_real_
    density_dominance_index <- NA_real_
  } else {
    p <- density_bin_totals / sum(density_bin_totals, na.rm = TRUE)
    p_nonzero <- p[p > 0 & is.finite(p)]
    density_evenness <- -sum(p_nonzero * log(p_nonzero)) / log(3)
    density_dominance_index <- max(p, na.rm = TRUE)
  }

  resid_mean <- function(section_names) {
    safe_mean(residuals[section_names])
  }

  # ---- NEW density features -----------------------------------------------
  # Reuse section_rep (section MNAU proportions), section_density, the
  # regression residuals, and the low/mid/high density bins already computed.

  high_mean_rep        <- safe_mean(section_rep[high_sections])
  low_mean_rep         <- safe_mean(section_rep[low_sections])
  postcranial_mean_rep <- safe_mean(section_rep[groups$postcranial_total])
  appendicular_mean_rep <- safe_mean(section_rep[groups$appendicular])

  high_minus_low_density_survival <-
    if (is.na(high_mean_rep) || is.na(low_mean_rep)) NA_real_ else high_mean_rep - low_mean_rep

  high_to_low_density_logratio <-
    if (is.na(high_mean_rep) || is.na(low_mean_rep)) NA_real_ else
      log((high_mean_rep + epsilon) / (low_mean_rep + epsilon))

  low_density_to_postcranial_ratio <-
    if (is.na(low_mean_rep) || is.na(postcranial_mean_rep)) NA_real_ else
      (low_mean_rep + epsilon) / (postcranial_mean_rep + epsilon)

  low_density_deficit_from_appendicular <-
    if (is.na(appendicular_mean_rep) || is.na(low_mean_rep)) NA_real_ else
      appendicular_mean_rep - low_mean_rep

  low_density_axial_survival <-
    safe_mean(section_rep[c("vertebrae", "thorax", "hyoid")])

  fragile_survival_fraction <- {
    v <- section_rep[low_sections]
    if (all(is.na(v))) NA_real_ else sum(v > fragile_threshold, na.rm = TRUE) / length(low_sections)
  }

  low_density_presence_050 <-
    sum(section_rep[low_sections] >= 0.50, na.rm = TRUE)

  vertebrae_thorax_low_density_balance <- {
    xv <- unname(section_rep["vertebrae"]); xt <- unname(section_rep["thorax"])
    if (is.na(xv) || is.na(xt)) NA_real_ else 1 - abs(xv - xt) / (xv + xt + epsilon)
  }

  density_representation_r2 <- {
    ok <- regression_ok & is.finite(residuals)
    if (sum(ok) < 2) {
      NA_real_
    } else {
      ss_res <- sum(residuals[ok]^2)
      ss_tot <- sum((r[ok] - mean(r[ok]))^2)
      if (!is.finite(ss_tot) || ss_tot <= epsilon) NA_real_ else 1 - ss_res / ss_tot
    }
  }

  density_residual_sd <- {
    rr <- residuals[is.finite(residuals)]
    if (length(rr) < 2) NA_real_ else stats::sd(rr)
  }

  density_residual_appendicular <- resid_mean(groups$appendicular)
  density_residual_cranial      <- resid_mean(groups$cranial_complex)
  density_residual_axial        <- resid_mean(groups$axial_main)

  dom_section_rep <- {
    v <- section_rep[is.finite(section_rep)]
    if (length(v) == 0) NA_character_ else names(v)[which.max(v)]
  }

  dominant_low_density_flag <-
    if (is.na(dom_section_rep)) NA_real_ else as.numeric(dom_section_rep %in% low_sections)

  mnau_dominant_section_density <-
    if (is.na(dom_section_rep)) NA_real_ else unname(section_density[dom_section_rep])

  mnau_low_density_dominance <- {
    v <- section_rep[low_sections]; v <- v[is.finite(v)]
    tot <- sum(v)
    if (length(v) == 0 || !is.finite(tot) || tot <= 0) NA_real_ else max(v / tot)
  }

  tibble::tibble(
    density_weighted_mean_mnau =
      safe_weighted_density(section_rep),

    density_weighted_mean_mne =
      safe_weighted_density(section_mne),

    density_representation_spearman =
      density_representation_spearman,

    density_representation_slope =
      density_representation_slope,

    high_density_to_low_density_ratio =
      (safe_mean(section_rep[high_sections]) + epsilon) /
      (safe_mean(section_rep[low_sections]) + epsilon),

    high_density_share =
      safe_share(high_total, total_rep),

    low_density_share =
      safe_share(low_total, total_rep),

    low_density_survival_index =
      safe_mean(section_rep[low_sections]),

    fragile_section_survival_count =
      sum(section_rep[low_sections] > fragile_threshold, na.rm = TRUE),

    density_evenness =
      density_evenness,

    density_dominance_index =
      density_dominance_index,

    density_adjusted_cranial_postcranial_contrast =
      resid_mean(groups$cranial_complex) -
      resid_mean(groups$postcranial_total),

    density_adjusted_limb_trunk_contrast =
      resid_mean(groups$limb_shaft_proxy) -
      resid_mean(groups$column_trunk),

    density_adjusted_appendicular_axial_contrast =
      resid_mean(groups$appendicular) -
      resid_mean(groups$axial_main),

    density_adjusted_distal_proximal_contrast =
      resid_mean(groups$distal_appendicular) -
      resid_mean(groups$proximal_appendicular),

    density_corrected_marrow_proxy_contrast =
      resid_mean(groups$limb_shaft_proxy) -
      resid_mean(c("hands", "foot", "thorax")),

    high_minus_low_density_survival       = high_minus_low_density_survival,
    high_to_low_density_logratio          = high_to_low_density_logratio,
    low_density_to_postcranial_ratio      = low_density_to_postcranial_ratio,
    low_density_deficit_from_appendicular = low_density_deficit_from_appendicular,
    low_density_axial_survival            = low_density_axial_survival,
    fragile_survival_fraction             = fragile_survival_fraction,
    low_density_presence_050              = low_density_presence_050,
    vertebrae_thorax_low_density_balance  = vertebrae_thorax_low_density_balance,
    density_representation_r2             = density_representation_r2,
    density_residual_sd                   = density_residual_sd,
    density_residual_appendicular         = density_residual_appendicular,
    density_residual_cranial              = density_residual_cranial,
    density_residual_axial                = density_residual_axial,
    dominant_low_density_flag             = dominant_low_density_flag,
    mnau_dominant_section_density         = mnau_dominant_section_density,
    mnau_low_density_dominance            = mnau_low_density_dominance
  )
}

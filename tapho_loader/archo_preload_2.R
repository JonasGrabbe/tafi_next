#' Preload/validate a two-column site table (meta rows 1:31 + skeletal rows)
#'
#' Expected format:
#'   - Column 1: element names
#'   - Column 2: one site/individual column (ind_name)
#'   - Rows 1:31: meta fields
#'   - Rows 32+: skeletal elements (counts)
#'
#' This is hard-coded and does NOT depend on external meta.xlsx at runtime.
archo_preload_site_table <- function(x, strict = FALSE) {

  issues <- character(0)
  add_issue <- function(z) issues <<- c(issues, z)
  warn_or_stop <- function(z) {
    if (isTRUE(strict)) stop(z, call. = FALSE)
    add_issue(z)
  }

  if (!is.data.frame(x)) stop("archo_preload_site_table(): input must be a data.frame.", call. = FALSE)
  if (ncol(x) < 2) stop("archo_preload_site_table(): input must have at least 2 columns.", call. = FALSE)
  if (nrow(x) < 32) stop("archo_preload_site_table(): expected >= 32 rows (31 meta + skeletal).", call. = FALSE)

  # Keep only first 2 columns
  x <- x[, 1:2, drop = FALSE]

  # Preserve 2nd column name = ind_name
  ind_name <- colnames(x)[2]
  if (is.null(ind_name) || !nzchar(ind_name)) ind_name <- "sample_1"

  colnames(x) <- c("element", ind_name)
  x$element <- trimws(as.character(x$element))

  meta <- x[1:31, , drop = FALSE]
  raw  <- x[-c(1:31), , drop = FALSE]

  # --- Required meta keys (hard-coded) ---
  required_meta <- c("site_id", "individual_id", "site_name", "count_type",
                     "geological_period", "cultural_period", "reference")
  missing_meta <- setdiff(required_meta, meta$element)
  if (length(missing_meta) > 0) {
    warn_or_stop(paste0("Missing required meta row(s): ", paste(missing_meta, collapse = ", ")))
  }

  # --- Skeletal count checks ---
  raw_num <- suppressWarnings(as.numeric(raw[[ind_name]]))
  if (any(is.na(raw_num) & !is.na(raw[[ind_name]]))) {
    warn_or_stop("Non-numeric skeletal values found; coerced to NA then set to 0.")
  }

  # NA -> 0 
  if (any(is.na(raw_num))) {
    warn_or_stop("NA skeletal counts replaced with 0 (NA = not identified).")
    raw_num[is.na(raw_num)] <- 0
  }

  # integer coercion
  non_int <- abs(raw_num - round(raw_num)) > 1e-9
  if (any(non_int)) {
    warn_or_stop("Non-integer skeletal counts detected; rounded to nearest integer.")
    raw_num <- round(raw_num)
  }

  # non-negative
  if (any(raw_num < 0, na.rm = TRUE)) {
    warn_or_stop("Negative skeletal counts detected; set to 0.")
    raw_num[raw_num < 0] <- 0
  }

  raw[[ind_name]] <- as.integer(raw_num)

  out <- rbind(meta, raw)

  if (length(issues) > 0 && !isTRUE(strict)) {
    warning(paste0("archo_preload_site_table validation report:\n- ",
                   paste(unique(issues), collapse = "\n- ")), call. = FALSE)
  }

  attr(out, "validation_report") <- unique(issues)
  return(out)
}

# Optional convenience wrapper: if you call archo_preload() it will validate your 2-column site table format
archo_preload <- function(x, strict = FALSE) {
  if (is.data.frame(x) && ncol(x) >= 2) {
    return(archo_preload_site_table(x, strict = strict))
  }
  stop("archo_preload(): unsupported input format (expected a two-column site table).", call. = FALSE)
}

# Helper: parse "1-0" or "1-0-3" individual ids
parse_individual_id <- function(id) {
  id <- trimws(as.character(id))
  parts <- strsplit(id, "-", fixed = TRUE)[[1]]
  parts_num <- suppressWarnings(as.integer(parts))
  list(
    site_id = if (length(parts_num) >= 1) parts_num[1] else NA_integer_,
    sample_no = if (length(parts_num) >= 2) parts_num[2] else NA_integer_,
    individual_no = if (length(parts_num) >= 3) parts_num[3] else NA_integer_
  )


}
# not used
tapho.load_2 <- function(new_site_data, strict = FALSE) {

  #  1) Validate first 
  new_site_data <- archo_preload(new_site_data, strict = strict)

  ind_name <- colnames(new_site_data)[2]

  colnames(new_site_data) <- c("element", "count")
  new_site_meta <- new_site_data[1:31, , drop = FALSE]
  new_site_raw  <- new_site_data[-c(1:31), , drop = FALSE]

  target_levels <- summarize_bone_sections_zero_na(new_site_raw)

  # NA skeletal counts => 0 (already applied by preload, but keep safe)
  new_site_raw$count[is.na(new_site_raw$count)] <- 0

  recalc_new_site_raw <- recal_mne_six(
    as.numeric(new_site_raw$count),
    target_levels,
    skeleton_hierarchy,
    new_site_raw$element,
    sections_order
  )

  one_skeleton <- get_bones_dataframe(target_levels)

  #old func
  data_aligned <- align_loaded_data(
    one_skeleton$element,
    recalc_new_site_raw$element,
    recalc_new_site_raw$count
  )

  tapho_loaded <- loaded_tapho(
    ind_name     = ind_name,
    meta         = new_site_meta,
    level        = target_levels,
    col_element  = data_aligned$element,
    col_count    = data_aligned$count,
    strict       = strict
  )

  return(tapho_loaded)
}

loaded_tapho <- function(ind_name, meta, level, col_element, col_count, strict = FALSE) {

  issues <- character(0)
  add_issue <- function(x) issues <<- c(issues, x)
  warn_or_stop <- function(x) {
    if (isTRUE(strict)) stop(x, call. = FALSE)
    add_issue(x)
  }

  # ---- level validation ----
  if (length(level) != 11) {
    warn_or_stop(sprintf("loaded_tapho(): level must have length 11 (got %d). Padding/truncating.", length(level)))
    if (length(level) < 11) level <- c(level, rep(NA_real_, 11 - length(level)))
    if (length(level) > 11) level <- level[1:11]
  }

  # ---- skeletal vectors ----
  if (length(col_element) != length(col_count)) {
    stop("loaded_tapho(): col_element and col_count must have the same length.", call. = FALSE)
  }
  col_element <- trimws(as.character(col_element))

  col_count_num <- suppressWarnings(as.numeric(col_count))
  if (any(is.na(col_count_num) & !is.na(col_count))) {
    warn_or_stop("loaded_tapho(): Non-numeric skeletal counts found; coerced to NA then set to 0.")
  }
  col_count_num[is.na(col_count_num)] <- 0

  non_int <- abs(col_count_num - round(col_count_num)) > 1e-9
  if (any(non_int)) {
    warn_or_stop("loaded_tapho(): Non-integer skeletal counts detected; rounded to nearest integer.")
    col_count_num <- round(col_count_num)
  }

  if (any(col_count_num < 0, na.rm = TRUE)) {
    warn_or_stop("loaded_tapho(): Negative skeletal counts detected; set to 0.")
    col_count_num[col_count_num < 0] <- 0
  }
  col_count_num <- as.integer(col_count_num)

  # ---- meta parsing (31-row element/count block) ----
  if (!is.data.frame(meta) || !all(c("element", "count") %in% colnames(meta))) {
    stop("loaded_tapho(): meta must be a data.frame with columns element + count.", call. = FALSE)
  }

  meta$element <- trimws(as.character(meta$element))
  meta_named <- stats::setNames(as.list(meta$count), meta$element)

  get_meta <- function(k) {
    if (!is.null(meta_named[[k]])) return(meta_named[[k]])
    NA
  }

  # ---- hard-coded type coercers ----
  parse_int <- function(x, field) {
    out <- suppressWarnings(as.integer(as.numeric(x)))
    if (is.na(out) && !is.na(x) && nzchar(as.character(x))) {
      warn_or_stop(sprintf("Meta field '%s' expected integer, got '%s'.", field, as.character(x)))
    }
    out
  }

  parse_chr <- function(x) {
    if (is.null(x) || is.na(x)) return(NA_character_)
    x <- as.character(x)
    if (!nzchar(trimws(x))) return(NA_character_)
    x
  }

  parse_lgl <- function(x, field) {
    if (is.null(x) || is.na(x)) return(NA)
    if (is.logical(x)) return(x)
    s <- tolower(trimws(as.character(x)))
    if (s %in% c("true","t","1","yes","y")) return(TRUE)
    if (s %in% c("false","f","0","no","n")) return(FALSE)
    warn_or_stop(sprintf("Meta field '%s' expected boolean, got '%s' (kept as NA).", field, as.character(x)))
    NA
  }

  # ---- meta fields (from meta.xlsx schema; cleaned) ----
  site_id            <- parse_int(get_meta("site_id"), "site_id")
  individual_id      <- parse_chr(get_meta("individual_id"))
  site_name          <- parse_chr(get_meta("site_name"))
  layer_name         <- parse_chr(get_meta("layer_name"))
  count_type         <- parse_lgl(get_meta("count_type"), "count_type")

  geological_period  <- parse_chr(get_meta("geological_period"))
  cultural_period    <- parse_chr(get_meta("cultural_period"))
  culture_period_2   <- parse_chr(get_meta("culture_period_2"))

  start_date_cal_bp  <- parse_int(get_meta("start_date_cal_bp"), "start_date_cal_bp")
  end_date_cal_bp    <- parse_int(get_meta("end_date_cal_bp"), "end_date_cal_bp")

  mis_stage          <- parse_chr(get_meta("mis_stage"))

  karstic_system     <- parse_chr(get_meta("karstic_system"))
  karstic_system_2   <- parse_chr(get_meta("karstic_system_2"))
  open_air_site      <- parse_chr(get_meta("open_air_site"))
  open_air_site_2    <- parse_chr(get_meta("open_air_site_2"))

  region             <- parse_chr(get_meta("region"))
  country            <- parse_chr(get_meta("country"))

  mni                <- parse_int(get_meta("mni"), "mni")
  age_category       <- parse_chr(get_meta("age_category"))
  n_infant           <- parse_int(get_meta("n_infant"), "n_infant")
  n_juvenil          <- parse_int(get_meta("n_juvenil"), "n_juvenil")
  n_subadult         <- parse_int(get_meta("n_subadult"), "n_subadult")
  n_adult            <- parse_int(get_meta("n_adult"), "n_adult")

  gender             <- parse_chr(get_meta("gender"))

  taphonomic_context <- parse_chr(get_meta("taphonomic_context"))
  funerary_context   <- parse_chr(get_meta("funerary_context"))
  discovery_context  <- parse_chr(get_meta("discovery_context"))
  degree_of_excavation <- parse_chr(get_meta("degree_of_excavation"))

  hominin_species    <- parse_chr(get_meta("hominin_species"))
  observation        <- parse_chr(get_meta("observation"))
  reference          <- parse_chr(get_meta("reference"))

  # ---- meta restrictions (hard-coded, warning-first) ----
  # start >= end (BP)
  if (!is.na(start_date_cal_bp) && !is.na(end_date_cal_bp) && start_date_cal_bp < end_date_cal_bp) {
    warn_or_stop("Meta restriction: start_date_cal_bp should be >= end_date_cal_bp (BP).")
  }

  # karstic vs open air exclusivity
  if (!is.na(karstic_system) && nzchar(karstic_system) && !is.na(open_air_site) && nzchar(open_air_site)) {
    warn_or_stop("Meta restriction: karstic_system and open_air_site should not both be filled (mutually exclusive).")
  }

  # basic categorical checks (soft)
  allowed_geol <- c("Lower Pleistocene","Middle Pleistocene","Late Pleistocene","Holocene")
  if (!is.na(geological_period) && !(geological_period %in% allowed_geol)) {
    warn_or_stop(sprintf("Meta warning: geological_period '%s' not in allowed set.", geological_period))
  }

  allowed_karst <- c("Cave","Rockshelter")
  if (!is.na(karstic_system) && !(karstic_system %in% allowed_karst)) {
    warn_or_stop(sprintf("Meta warning: karstic_system '%s' unexpected (allowed: Cave, Rockshelter, NA).", karstic_system))
  }

  if (!is.na(open_air_site) && !(open_air_site %in% c("Open Air"))) {
    warn_or_stop(sprintf("Meta warning: open_air_site '%s' unexpected (allowed: Open Air, NA).", open_air_site))
  }

  allowed_gender <- c("Female","Male","Indet/Mixed")
  if (!is.na(gender) && !(gender %in% allowed_gender)) {
    warn_or_stop(sprintf("Meta warning: gender '%s' not in allowed set.", gender))
  }

  allowed_agecat <- c("Infantile","Juvenile","Subadult","Adult","Indet")
  if (!is.na(age_category) && !(age_category %in% allowed_agecat)) {
    warn_or_stop(sprintf("Meta warning: age_category '%s' not in allowed set.", age_category))
  }

  # individual_id format check (string is fine!)
  if (!is.na(individual_id)) {
    ok <- grepl("^\\d+-\\d+(?:-\\d+)?$", individual_id)
    if (!ok) {
      warn_or_stop(sprintf("Meta warning: individual_id '%s' does not match expected pattern like '1-0' or '1-0-3'.", individual_id))
    }
  }

  # derive numeric parts from individual_id if possible
  parsed_id <- if (!is.na(individual_id)) parse_individual_id(individual_id) else list(site_id=NA, sample_no=NA, individual_no=NA)

  if (!is.na(site_id) && !is.na(parsed_id$site_id) && site_id != parsed_id$site_id) {
    warn_or_stop(sprintf("Meta warning: site_id (%d) does not match individual_id prefix (%d).", site_id, parsed_id$site_id))
  }

  # compute sum from counts
  sum_count <- sum(col_count_num, na.rm = TRUE)

  # Emit aggregated warnings
  if (length(issues) > 0 && !isTRUE(strict)) {
    warning(paste0("loaded_tapho validation:\n- ", paste(unique(issues), collapse = "\n- ")), call. = FALSE)
  }

  #  UPDATED OUTPUT: include ALL meta columns (site_id ... reference)
  dataset <- tibble::tibble(
    ind_name = ind_name,

    # Common identifiers
    #site_index = site_id,
    #sample_index = individual_id,

    # Meta columns (1:31) as explicit fields
    site_id = site_id,
    individual_id = individual_id,
    site_name = site_name,
    layer_name = layer_name,
    count_type = count_type,
    geological_period = geological_period,
    cultural_period = cultural_period,
    culture_period_2 = culture_period_2,
    start_date_cal_bp = start_date_cal_bp,
    end_date_cal_bp = end_date_cal_bp,
    mis_stage = mis_stage,
    karstic_system = karstic_system,
    karstic_system_2 = karstic_system_2,
    open_air_site = open_air_site,
    open_air_site_2 = open_air_site_2,
    region = region,
    country = country,
    mni = mni,
    age_category = age_category,
    n_infant = n_infant,
    n_juvenil = n_juvenil,
    n_subadult = n_subadult,
    n_adult = n_adult,
    gender = gender,
    taphonomic_context = taphonomic_context,
    funerary_context = funerary_context,
    discovery_context = discovery_context,
    degree_of_excavation = degree_of_excavation,
    hominin_species = hominin_species,
    observation = observation,
    reference = reference,

    # derived numeric parts for easy later calculations
    #sample_site_id = parsed_id$site_id,
    #sample_no = parsed_id$sample_no,
    #individual_no = parsed_id$individual_no,

    # summary
    sum_count = sum_count,

    # list-cols for analysis
    bone.level = list(level),
    element = list(col_element),
    count = list(col_count_num)
  )

  attr(dataset, "validation_issues") <- unique(issues)
  return(dataset)
}
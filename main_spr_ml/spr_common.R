# Shared data preparation, provenance and site-level description.
# Source after the corrected HSPA loader and spr_ml_core.R.

spr_dependency_paths <- function(data_root, model_root, input_file) {
  relative <- c(
    "tapho_loader_new.R", "try_shared_allocator_proposal.R",
    "select_data_new.R", "data_loader_new.R", "archo_preload_2.R",
    "collapse_down_step.R", "indet.R", "meta/meta_set.R",
    "derived_features.R", "density_features_fixed.R",
    "convertMNEtovariables_new.R"
  )
  files <- c(input_file,
             file.path(data_root, relative),
             file.path(model_root, "spr_ml_core.R"))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("Missing HSPA input or source file(s):\n",
                            paste(missing, collapse = "\n"))
  normalizePath(files, winslash = "/", mustWork = TRUE)
}

spr_source_manifest <- function(data_root, model_root, input_file) {
  paths <- spr_dependency_paths(data_root, model_root, input_file)
  data.frame(
    file = paths,
    md5 = unname(as.character(tools::md5sum(paths))),
    stringsAsFactors = FALSE
  )
}

spr_prepare_level <- function(level, data_root, model_root,
                              input_file, cache_root, source_manifest) {
  if (length(level) != 11L) stop("Each anatomical level must have 11 integers")
  level_id <- paste(level, collapse = "_")
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  cache_path <- file.path(cache_root,
    paste0("prepared_sum_and_ind_by_site_id_L", level_id, ".rds"))
  cache_key <- list(
    schema = "SPR_publication_v1_corrected_loader",
    level = as.integer(level),
    source_md5 = setNames(source_manifest$md5, source_manifest$file)
  )
  stored <- if (file.exists(cache_path))
    tryCatch(readRDS(cache_path), error = function(e) NULL) else NULL

  if (is.list(stored) && identical(stored$key, cache_key)) {
    tables <- stored$tables
    message("  Prepared anatomical representations: cached")
  } else {
    message("  Preparing anatomical representations from HSPA")
    tables <- tapho.loader.wfl7(
      data = input_file, level = level, indet = FALSE,
      indet_method = "delet", exclude_sections = c("teeth")
    )
  }
  if (!is.list(tables) ||
      !all(c("db_sum", "db_id") %in% names(tables))) {
    stop("Loader must return db_sum and db_id")
  }
  sum <- tables$db_sum
  ind <- tables$db_id
  if (!all(c("site_id", "taphonomic_context") %in% names(sum)) ||
      !"site_id" %in% names(ind)) {
    stop("Required site_id or class column missing from prepared data")
  }
  ids <- trimws(as.character(sum$site_id))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("db_sum must have exactly one row per nonmissing site_id")
  }
  donor_ids <- trimws(as.character(ind$site_id))
  if (anyNA(donor_ids) || any(!nzchar(donor_ids))) {
    stop("Missing site_id in individual donor records")
  }
  if (!is.list(stored) || !identical(stored$key, cache_key)) {
    saveRDS(list(key = cache_key, tables = tables), cache_path)
  }
  tables
}

spr_write_site_coordinates <- function(tables, level_dir) {
  x <- tables$db_sum
  need <- c("site_id", "site_name", "taphonomic_context",
            "Mean_part_abundance", "mnau_heip_evenness")
  if (!all(need %in% names(x))) stop("Missing site coordinate column")
  dat <- x[, need, drop = FALSE]
  dat$site_id <- as.character(dat$site_id)
  # The original archaeological context is retained alongside the modelling class.
  ref <- ml7_recode(x, "fc_vs_rest")
  labels <- setNames(as.character(ref$taphonomic_context),
                     as.character(ref$site_id))
  dat$model_reference_label <- unname(labels[dat$site_id])
  dat$sample_group <- ifelse(!is.na(dat$model_reference_label),
                           "Reference",
                           ifelse(tolower(trimws(dat$taphonomic_context)) %in%
                             c("unknown", "unkown"), "Unknown", "Other"))
  readr::write_csv(dat, file.path(level_dir,"site_coordinates_C_Heip.csv"))
  readr::write_csv(dat[dat$sample_group == "Reference",,drop=FALSE],
                   file.path(level_dir,"reference_C_Heip.csv"))
  readr::write_csv(dat[dat$sample_group == "Unknown",,drop=FALSE],
                   file.path(level_dir,"unknown_C_Heip.csv"))
  classes <- table(as.character(ml7_recode(x, "4class")$taphonomic_context))
  readr::write_csv(data.frame(class=names(classes), n_sites=as.integer(classes)),
                   file.path(level_dir,"reference_class_counts.csv"))
  invisible(dat)
}

spr_record_run <- function(out, settings, data_root, model_root, input_file) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  manifest <- spr_source_manifest(data_root, model_root, input_file)
  readr::write_csv(manifest, file.path(out,"input_and_code_checksums.csv"))
  release_scripts <- file.path(model_root, c(
    "run_cv_spr.R", "run_projection_spr.R", "spr_common.R",
    "spr_ml_core.R", "spr_projection_core.R",
    "spr_metadata_augmentation.R"))
  if (any(!file.exists(release_scripts))) {
    stop("The six publication scripts must be installed together")
  }
  readr::write_csv(data.frame(
    file = normalizePath(release_scripts, winslash = "/", mustWork = TRUE),
    md5 = unname(as.character(tools::md5sum(release_scripts)))),
    file.path(out, "analysis_code_checksums.csv"))
  saveRDS(settings, file.path(out,"analysis_settings.rds"))
  dput(settings, file = file.path(out,"analysis_settings.R"))
  writeLines(capture.output(utils::sessionInfo()),
             con = file.path(out,"session_info.txt"))
  manifest
}

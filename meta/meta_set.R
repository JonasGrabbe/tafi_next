
#_________________________________   META ENV   _____________________________________

meta_set <- function(new_indet) {

  if ("meta_env" %in% search()) {
    detach("meta_env")
  }

  meta_file <- paste0(
    "meta_indet_",
    toupper(new_indet),
    ".R"
  )

  source(
    file.path("meta", meta_file),
    local = TRUE
  )

  meta <- list(
    skeleton_hierarchy = skeleton_hierarchy,
    elements_data      = elements_data,
    category_hierarchy = category_hierarchy,
    category_mappings  = category_mappings
  )

  meta_env <- list2env(
    meta,
    parent = emptyenv()
  )

  attach(
    meta_env,
    name = "meta_env",
    warn.conflicts = FALSE
  )

  invisible(meta)
}
#_______________________________________________________________________________________
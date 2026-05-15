

collapse_down_step_special <- function(data, part, current_level, target_level, section_name, MNI) {
  key         <- paste0(part, "_", target_level, "_to_", current_level)
  if (!key %in% names(category_hierarchy)) stop("Missing mapping: ", key)

  mapping      <- category_hierarchy[[key]]
  level_label  <- paste0("Level ", target_level)
  level_maxima <- skeleton_hierarchy[[section_name]][[level_label]]

  # Pre-compute which _indet children are consumed (parent_count > 0)
  consumed_indet_names <- purrr::map(names(mapping), function(parent_elem) {
    parent_clean   <- gsub(" ", ".", parent_elem)
    children       <- gsub(" ", ".", mapping[[parent_elem]])
    indet_children <- children[grepl("_indet$", children)]

    parent_count <- data %>%
      filter(std_element == parent_clean) %>%
      pull(count) %>%
      { if (length(.) == 0) 0L else as.integer(.[1]) }  # ← take first, guard duplicates

    if (parent_count > 0 && length(indet_children) > 0) indet_children else character(0)
  }) %>% unlist() %>% unique()

  new_rows <- purrr::map_dfr(names(mapping), function(parent_elem) {

    parent_clean  <- gsub(" ", ".", parent_elem)
    children      <- gsub(" ", ".", mapping[[parent_elem]])
    real_children <- drop_indet_if_possible(children)

    # Own count at parent level — always scalar
    parent_count <- data %>%
      filter(std_element == parent_clean) %>%
      pull(count) %>%
      { if (length(.) == 0) 0L else as.integer(.[1]) }  # ← take first, guard duplicates

    # _indet children in this mapping
    indet_children <- children[grepl("_indet$", children)]
    indet_counts <- purrr::map_int(indet_children, \(ch) {
      val <- data %>% filter(std_element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val[1])  # ← take first, guard duplicates
    })

    # Only absorb _indet into parent if parent itself has a count
    effective_parent_count <- if (parent_count > 0) {
      if (part %in% c("pelvic", "cranium")) {
        max(c(parent_count, indet_counts))
      } else {
        parent_count + sum(indet_counts)
      }
    } else {
      0L
    }

    indet_consumed <- parent_count > 0 && length(indet_children) > 0

    # Existing counts on real (non-indet) children
    existing_counts <- purrr::map_int(real_children, \(ch) {
      val <- data %>% filter(std_element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val[1])  # ← take first
    })

    max_counts <- purrr::map_int(real_children, \(ch) {
      base <- level_maxima[[ch]]
      if (is.null(base)) 0L else as.integer(base * MNI)
    })

    # Self-loop: single child same as parent
    if (length(real_children) == 1 && real_children[1] == parent_clean) {
      indet_rows <- if (indet_consumed) {
        tibble(std_element = indet_children, count = 0L)
      } else {
        tibble(std_element = character(0), count = integer(0))
      }
      return(bind_rows(
        tibble(std_element = real_children, count = effective_parent_count),
        indet_rows
      ))
    }

    # ── pelvic / cranium logic ─────────────────────────────────────────────
    if (part %in% c("pelvic", "cranium")) {
      if (effective_parent_count > 0 && all(existing_counts == 0)) {
        chosen       <- sample(real_children, 1)
        final_counts <- integer(length(real_children))
        final_counts[real_children == chosen] <- effective_parent_count
      } else {
        final_counts <- existing_counts
      }
    } else {
    # ── original logic ───────────────────────────────────────────────────────
      surplus      <- max(0L, effective_parent_count)
      final_counts <- distribute_surplus(existing_counts, max_counts, surplus)
    }

    indet_rows <- if (indet_consumed) {
      tibble(std_element = indet_children, count = 0L)
    } else {
      tibble(std_element = character(0), count = integer(0))
    }

    bind_rows(
      tibble(std_element = real_children, count = final_counts),
      indet_rows
    )
  })

  parent_names          <- gsub(" ", ".", names(mapping))
  all_child_names       <- unique(gsub(" ", ".", unlist(mapping)))
  non_indet_child_names <- all_child_names[!grepl("_indet$", all_child_names)]

  data %>%
    filter(!std_element %in% c(parent_names, non_indet_child_names, consumed_indet_names)) %>%
    bind_rows(new_rows) %>%
    distinct(std_element, .keep_all = TRUE)  # ← deduplicate as final safety net
}


collapse_down_step_direct_special <- function(data, part, current_level, target_level, section_name, MNI) {
  key <- paste0(part, "_", target_level, "_to_", current_level)
  if (!key %in% names(category_hierarchy)) stop("Missing mapping: ", key)

  mapping      <- category_hierarchy[[key]]
  level_label  <- paste0("Level ", target_level)
  level_maxima <- skeleton_hierarchy[[section_name]][[level_label]]

  new_rows <- purrr::map_dfr(names(mapping), function(parent_elem) {

    parent_clean  <- gsub(" ", ".", parent_elem)
    children      <- gsub(" ", ".", mapping[[parent_elem]])
    real_children <- drop_indet_if_possible(children)
    all_children  <- children  # includes _indet

    parent_count <- data %>%
      filter(std_element == parent_clean) %>%
      pull(count) %>%
      { if (length(.) == 0) 0L else as.integer(.[1]) }

    # _indet children start at 0 — they are new at this level
    indet_children <- children[grepl("_indet$", children)]

    existing_counts <- purrr::map_int(real_children, \(ch) {
      val <- data %>% filter(std_element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val[1])
    })

    max_counts <- purrr::map_int(real_children, \(ch) {
      base <- level_maxima[[ch]]
      if (is.null(base)) 0L else as.integer(base * MNI)
    })

    # Self-loop
    if (length(real_children) == 1 && real_children[1] == parent_clean) {
      return(bind_rows(
        tibble(std_element = real_children, count = as.integer(parent_count)),
        tibble(std_element = indet_children, count = 0L)
      ))
    }

    # ── pelvic / cranium logic ─────────────────────────────────────────────
    if (part %in% c("pelvic", "cranium")) {
      if (parent_count > 0 && all(existing_counts == 0)) {
        chosen       <- sample(real_children, 1)
        final_counts <- integer(length(real_children))
        final_counts[real_children == chosen] <- parent_count
      } else {
        final_counts <- existing_counts
      }
    } else {
    # ── original logic ───────────────────────────────────────────────────────
      surplus      <- max(0L, parent_count)
      final_counts <- distribute_surplus(existing_counts, max_counts, surplus)
    }

    # Return real children with distributed counts + _indet children as 0
    bind_rows(
      tibble(std_element = real_children, count = final_counts),
      tibble(std_element = indet_children, count = 0L)
    )
  })

  parent_names    <- gsub(" ", ".", names(mapping))
  all_child_names <- unique(gsub(" ", ".", unlist(mapping)))  # raw, includes _indet

  data %>%
    filter(!std_element %in% c(parent_names, all_child_names)) %>%
    bind_rows(new_rows) %>%
    distinct(std_element, .keep_all = TRUE)
}


collapse_up_step_special <- function(data, part, current_level, target_level) {
  key <- paste0(part, "_", current_level, "_to_", target_level)
  if (!key %in% names(category_hierarchy)) stop("Missing mapping: ", key)
  mapping <- category_hierarchy[[key]]

  new_rows <- purrr::map_dfr(names(mapping), function(parent_elem) {
    parent_clean  <- gsub(" ", ".", parent_elem)
    children      <- gsub(" ", ".", mapping[[parent_elem]])
    real_children <- children  # includes _indet — they must contribute

    child_counts <- purrr::map_int(real_children, \(ch) {
      val <- data %>% filter(std_element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val[1])
    })

    existing_parent <- data %>% filter(std_element == parent_clean) %>% pull(count)
    existing_parent <- if (length(existing_parent) == 0) 0L else as.integer(existing_parent[1])

    parent_count <- if (part %in% c("pelvic", "cranium")) {
      if (length(child_counts) == 0) existing_parent
      else max(c(child_counts, existing_parent))
    } else {
      # sum children + existing parent, but avoid double-counting self-loops
      self_loop <- real_children == parent_clean
      sum(child_counts[!self_loop]) + existing_parent
    }

    tibble(std_element = parent_clean, count = as.integer(parent_count))
  })

  parent_names    <- gsub(" ", ".", names(mapping))
  all_child_names <- unique(gsub(" ", ".", unlist(mapping)))

  data %>%
    filter(!std_element %in% c(parent_names, all_child_names)) %>%
    bind_rows(new_rows) %>%
    distinct(std_element, .keep_all = TRUE)
}
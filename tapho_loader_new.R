
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

####################
make_tapho_icons <- function(
  icons = c("🦴", "💀", "🔥", "🪨"),
  weights = c(4, 3, 2, 1)
) {

  # --------------------------------------------------
  # Preserve analytical RNG state
  # --------------------------------------------------

  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv
    )
  }

  # Independent temporary seed for cosmetic randomness
  cosmetic_seed <- as.integer(
    (
      as.numeric(Sys.time()) * 1000 +
        Sys.getpid()
    ) %% .Machine$integer.max
  )

  set.seed(cosmetic_seed)


  # --------------------------------------------------
  # Pick 2 DISTINCT icon types, weighted
  # --------------------------------------------------

  selected <- sample(
    icons,
    size = 2,
    replace = FALSE,
    prob = weights
  )


  # --------------------------------------------------
  # Pick the dominant one, also weighted
  # --------------------------------------------------

  dominant <- sample(
    selected,
    size = 1,
    prob = weights[
      match(selected, icons)
    ]
  )

  minor <- setdiff(
    selected,
    dominant
  )


  # --------------------------------------------------
  # Exactly 9 dominant + 1 minor,
  # randomly positioned
  # --------------------------------------------------

  icon_sequence <- c(
    rep(dominant, 9),
    minor
  )

  icon_sequence <- sample(
    icon_sequence
  )


  # --------------------------------------------------
  # Restore analytical RNG exactly as it was
  # --------------------------------------------------

  if (had_seed) {

    assign(
      ".Random.seed",
      old_seed,
      envir = .GlobalEnv
    )

  } else if (
    exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  ) {

    rm(
      ".Random.seed",
      envir = .GlobalEnv
    )
  }


  return(icon_sequence)
}

##### 
make_tapho_status <- function(
  i,
  n,
  site = "",
  icons,
  label = "Tapho-loading"
) {

  progress <- i / n
  percent <- round(progress * 100)


  # --------------------------------------------------
  # 0% = 1 icon
  # 100% = all 10 icons
  # --------------------------------------------------

  icon_steps <- 1 + floor(
    progress * (length(icons) - 1)
  )

  icon_steps <- max(
    1,
    min(
      length(icons),
      icon_steps
    )
  )


  # --------------------------------------------------
  # Visible icons
  # --------------------------------------------------

  icon_progress <- paste0(
    icons[seq_len(icon_steps)],
    collapse = ""
  )


  # --------------------------------------------------
  # Keep icon field width constant
  # --------------------------------------------------

  full_icons <- paste0(
    icons,
    collapse = ""
  )

  full_width <- nchar(
    full_icons,
    type = "width"
  )

  current_width <- nchar(
    icon_progress,
    type = "width"
  )

  icon_field <- paste0(
    icon_progress,
    strrep(
      " ",
      max(
        0,
        full_width - current_width
      )
    )
  )


  # --------------------------------------------------
  # Status
  # --------------------------------------------------

  sprintf(
    "%s (%3d%%):  %s  %d / %d sites · %s",
    label,
    percent,
    icon_field,
    i,
    n,
    site
  )
}
####################


summarize_bone_sections_zero_na <- function(input_df,
  method = "min",
  meta) {
  
    list2env(meta, envir = environment())

    # Ensure method is a single string
    method <- as.character(method)
    if (length(method) != 1) {
    stop("'method' must be a single value: 'min', 'max', or 'mode'. Got length = ", length(method))
    }

    input_df <- input_df %>%
    dplyr::filter(!is.na(count) & count != 0)

    bone_info_df <- add_bone_info(input_df)
    
    sections_order <- c(
    "cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
    "arms", "hands", "pelvic", "legs", "foot"
    )

    stat_mode <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
    }

    stat_fn <- switch(method,
    "min"  = function(x) min(x, na.rm = TRUE),
    "max"  = function(x) max(x, na.rm = TRUE),
    "mode" = stat_mode,
    stop("method must be 'min', 'max', or 'mode'. Got: ", method)
    )

    summary_df <- bone_info_df %>%
    dplyr::group_by(Skeleton_Section) %>%
    dplyr::summarize(min_level = stat_fn(max_Level), .groups = "drop")

    most_detailed_level <- c(5,5,1,3,4,3,3,5,4,3,5)
    result_vector <- stats::setNames(most_detailed_level, sections_order)

    for (section in sections_order) {
    if (section %in% summary_df$Skeleton_Section) {
    level <- summary_df$min_level[summary_df$Skeleton_Section == section]
    result_vector[section] <- if (is.na(level)) 1 else level
    }
    }

    unname(result_vector)
}


distribute_surplus <- function(counts, max_counts, n_indet) {

  counts     <- as.integer(counts)
  max_counts <- as.integer(max_counts)
  n_indet    <- as.integer(n_indet)

  if (length(counts) != length(max_counts)) {
    stop(
      "counts and max_counts must have the same length.",
      call. = FALSE
    )
  }

  if (anyNA(counts) || anyNA(max_counts) || is.na(n_indet)) {
    stop(
      "counts, max_counts and n_indet must not contain NA.",
      call. = FALSE
    )
  }

  if (n_indet < 0L) {
    stop(
      "n_indet cannot be negative.",
      call. = FALSE
    )
  }

  # Existing data must already fit anatomical capacity.
  over <- which(counts > max_counts)

  if (length(over) > 0L) {
    stop(
      "\nExisting count exceeds anatomical capacity.",
      "\ncounts: ", paste(counts, collapse = ", "),
      "\nmax_counts: ", paste(max_counts, collapse = ", "),
      "\nover-capacity positions: ", paste(over, collapse = ", "),
      call. = FALSE
    )
  }

  if (n_indet == 0L) {
    return(counts)
  }

  open_slots <- pmax(
    0L,
    max_counts - counts
  )

  if (sum(open_slots) < n_indet) {
    stop(
      "\nNot enough anatomical capacity to distribute indeterminate specimens.",
      "\ncounts: ", paste(counts, collapse = ", "),
      "\nmax_counts: ", paste(max_counts, collapse = ", "),
      "\nopen_slots: ", paste(open_slots, collapse = ", "),
      "\nn_indet: ", n_indet,
      "\ntotal_open_slots: ", sum(open_slots),
      call. = FALSE
    )
  }

  slot_pool <- rep(
    seq_along(counts),
    times = open_slots
  )

  chosen <- sample(
    slot_pool,
    size = n_indet,
    replace = FALSE
  )

  for (i in chosen) {
    counts[i] <- counts[i] + 1L
  }

  counts
}
collapse_up_step <- function(data, part, current_level, target_level) {
  key <- paste0(part, "_", current_level, "_to_", target_level)
  if (!key %in% names(category_hierarchy)) stop("Missing mapping: ", key)

  element_map <- stack(category_hierarchy[[key]]) %>%
    setNames(c("child", "parent")) %>%
    mutate(across(everything(), \(x) gsub(" ", ".", x)))

  mappable     <- inner_join(data, element_map, by = c("std_element" = "child"))
  non_mappable <- anti_join(data, element_map, by = c("std_element" = "child"))

  aggregated <- mappable %>%
    group_by(std_element = parent) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

 data <-  full_join(non_mappable, aggregated, by = "std_element") %>%
    mutate(count = coalesce(count.x, 0) + coalesce(count.y, 0)) %>%
    select(std_element, count)


  full_join(non_mappable, aggregated, by = "std_element") %>%
    mutate(count = coalesce(count.x, 0) + coalesce(count.y, 0)) %>%
    select(std_element, count)
  }




drop_indet_if_possible <- function(x) {
  x <- gsub(" ", ".", x)
  if (length(x) > 1) {
    tmp <- x[!grepl("_indet$", x)]
    if (length(tmp) > 0) return(tmp)
  }
  x
}




collapse_down_step <- function(
  data,
  part,
  current_level,
  target_level,
  section_name,
  MNI
) {

  key <- paste0(
    part, "_", target_level, "_to_", current_level
  )

  if (!key %in% names(category_hierarchy)) {
    stop("Missing mapping: ", key, call. = FALSE)
  }

  mapping <- category_hierarchy[[key]]

  level_label <- paste0("Level ", target_level)

  level_maxima <- skeleton_hierarchy[[section_name]][[level_label]]

  if (is.null(level_maxima)) {
    stop(
      "Missing anatomical maxima for ",
      section_name, " / ", level_label,
      call. = FALSE
    )
  }

  if (length(MNI) != 1L || is.na(MNI) ||
      !is.finite(MNI) || MNI < 1) {
    stop("Invalid MNI: ", paste(MNI, collapse = ", "),
         call. = FALSE)
  }

  # Read one element count safely.
  # Do not silently ignore duplicated input rows.
  get_count <- function(element_name) {

    val <- data$count[
      data$std_element == element_name
    ]

    if (length(val) == 0L) {
      return(0L)
    }

    if (length(val) != 1L || is.na(val[1])) {
      stop(
        "Missing or duplicated count for element: ",
        element_name,
        "\nkey: ", key,
        call. = FALSE
      )
    }

    val <- as.integer(val[1])

    if (val < 0L) {
      stop(
        "Negative count for element: ",
        element_name,
        call. = FALSE
      )
    }

    val
  }

  new_rows <- purrr::map_dfr(
    names(mapping),
    function(parent_elem) {

      parent_clean <- gsub(" ", ".", parent_elem)

      children <- gsub(
        " ", ".", mapping[[parent_elem]]
      )

      real_children <- drop_indet_if_possible(children)

      # --------------------------------------------------
      # 1. Existing parent and child observations
      # --------------------------------------------------

      parent_count <- get_count(parent_clean)

      existing_counts <- vapply(
        real_children,
        get_count,
        integer(1)
      )

      max_counts <- purrr::map_int(
        real_children,
        function(ch) {

          base <- level_maxima[[ch]]

          if (is.null(base)) {
            stop(
              "Missing anatomical capacity for: ", ch,
              "\nkey: ", key,
              call. = FALSE
            )
          }

          as.integer(base * MNI)
        }
      )

      # --------------------------------------------------
      # 2. Self-loop: element remains unchanged
      # --------------------------------------------------

      if (
        length(real_children) == 1L &&
        real_children[1] == parent_clean
      ) {

        if (parent_count > max_counts[1]) {
          stop(
            "Existing element exceeds anatomical capacity.",
            "\nkey: ", key,
            "\nelement: ", parent_clean,
            "\ncount: ", parent_count,
            "\nmaximum: ", max_counts[1],
            "\nMNI: ", MNI,
            call. = FALSE
          )
        }

        return(
          tibble::tibble(
            std_element = real_children,
            count = as.integer(parent_count)
          )
        )
      }

      # --------------------------------------------------
      # 3. Check already identified children
      # --------------------------------------------------

      over <- which(existing_counts > max_counts)

      if (length(over) > 0L) {
        stop(
          "\nExisting child count exceeds anatomical capacity.",
          "\nkey: ", key,
          "\npart: ", part,
          "\ncurrent_level: ", current_level,
          "\ntarget_level: ", target_level,
          "\nparent_elem: ", parent_clean,
          "\nchildren: ",
          paste(real_children, collapse = ", "),
          "\nexisting_counts: ",
          paste(existing_counts, collapse = ", "),
          "\nmax_counts: ",
          paste(max_counts, collapse = ", "),
          "\nover-capacity children: ",
          paste(real_children[over], collapse = ", "),
          "\nMNI: ", MNI,
          call. = FALSE
        )
      }

      # --------------------------------------------------
      # 4. Calculate genuinely available slots
      # --------------------------------------------------

      open_slots <- pmax(
        0L,
        max_counts - existing_counts
      )

      # IMPORTANT CHANGE:
      # Parent observations are additional bones.
      # Do NOT subtract identified child counts.
      #
      # OLD:
      # surplus <- max(
      #   0L,
      #   parent_count - sum(existing_counts)
      # )
      # OLD 2: 
      #surplus <- max(
      #  0L,
      #  as.integer(parent_count) -
      #    sum(existing_counts, na.rm = TRUE)
      #)

      surplus <- as.integer(parent_count)

      # --------------------------------------------------
      # 5. Validate the additional parent count
      # --------------------------------------------------

      if (sum(open_slots) < surplus) {
        stop(
          "\nCannot distribute parent count.",
          "\nkey: ", key,
          "\npart: ", part,
          "\ncurrent_level: ", current_level,
          "\ntarget_level: ", target_level,
          "\nparent_elem: ", parent_clean,
          "\nchildren: ",
          paste(real_children, collapse = ", "),
          "\nparent_count: ", parent_count,
          "\nexisting_counts: ",
          paste(existing_counts, collapse = ", "),
          "\nmax_counts: ",
          paste(max_counts, collapse = ", "),
          "\nopen_slots: ",
          paste(open_slots, collapse = ", "),
          "\navailable_slots: ", sum(open_slots),
          "\nsurplus: ", surplus,
          "\nMNI: ", MNI,
          call. = FALSE
        )
      }

      # --------------------------------------------------
      # 6. Randomly assign additional parent bones
      #    ONLY into available anatomical slots
      # --------------------------------------------------

      final_counts <- existing_counts

      if (surplus > 0L) {

        slot_pool <- rep(
          seq_along(real_children),
          times = open_slots
        )

        chosen <- slot_pool[
          sample.int(
            n = length(slot_pool),
            size = surplus,
            replace = FALSE
          )
        ]

        additions <- tabulate(
          chosen,
          nbins = length(real_children)
        )

        final_counts <- existing_counts + additions
      }

      # --------------------------------------------------
      # 7. Final anatomical capacity validation
      # --------------------------------------------------

      if (any(final_counts > max_counts)) {
        stop(
          "\nInternal error: distribution exceeded capacity.",
          "\nkey: ", key,
          "\nparent_elem: ", parent_clean,
          "\nfinal_counts: ",
          paste(final_counts, collapse = ", "),
          "\nmax_counts: ",
          paste(max_counts, collapse = ", "),
          call. = FALSE
        )
      }

      tibble::tibble(
        std_element = real_children,
        count = as.integer(final_counts)
      )
    }
  )

  # ----------------------------------------------------
  # 8. Remove replaced parent/child rows and insert
  #    recalculated children
  # ----------------------------------------------------

  parent_names <- gsub(
    " ", ".", names(mapping)
  )

  child_names <- unique(
    gsub(" ", ".", unlist(mapping))
  )

  child_names <- drop_indet_if_possible(
    child_names
  )

  data %>%
    dplyr::filter(
      !std_element %in% c(parent_names, child_names)
    ) %>%
    dplyr::bind_rows(new_rows)
}


is_special <- function(part) part %in% c("pelvic", "cranium")

recal_mne_six <- function(current_mne, bone_names, target_levels, MNI) {

    max_levels <- c(5, 5, 1, 3, 4, 3, 3, 5, 4, 3, 5)
    min_levels <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
    working_base <- tibble(count = current_mne, element = bone_names) %>%
    add_bone_info() %>%
    mutate(std_element = gsub(" ", ".", element))

    purrr::imap(sections_order, \(part, i) {
    data <- working_base %>%
    filter(Skeleton_Section == part) %>%
    select(std_element, count)

    current_max <- max_levels[i]
    current_min <- min_levels[i]
    tgt <- target_levels[i]

    # Pass 1 — collapse from finest detail down to target
    while (current_max > tgt) {
    step_up_level <- current_max - 1

      if (is_special(part)){
        data <- collapse_up_step_special(data, part,
          current_level = current_max,
          target_level  = step_up_level)
      } 
      else{
          data <- collapse_up_step(data, part,
                current_level = current_max,
                target_level  = step_up_level)}
    current_max <- current_max - 1

    }

    # Pass 2 — expand from coarsest detail up to target
    while (current_min < tgt) {

    if (is_special(part)){   
          data <- collapse_down_step_special(data, part,
                        current_level = current_min,
                        target_level  = current_min + 1,
                        section_name  = part,
          MNI           = MNI)}
      else{
        data <- collapse_down_step(data, part,
                current_level = current_min,
                target_level  = current_min + 1,
                section_name  = part,
                MNI           = MNI)}
    current_min <- current_min + 1
    

    }

    data
    }) %>%
    bind_rows() %>%
    rename(element = std_element)
}


align_counts_to_skeleton <- function( one_skeleton_element,
                                      recalc_element,
                                      recalc_count,
                                      fill_value = NA,
                                      aggregate = c("sum", "first")) {
aggregate <- match.arg(aggregate)

# Build a lookup table (handle duplicate elements in recalc data)
recalc_df <- data.frame(element = recalc_element, count = recalc_count)

if (aggregate == "sum") {
recalc_df <- aggregate(count ~ element, data = recalc_df, FUN = sum)
} else if (aggregate == "first") {
recalc_df <- recalc_df[!duplicated(recalc_df$element), ]
}

# Match skeleton elements to recalc counts
aligned_counts <- recalc_df$count[match(one_skeleton_element, recalc_df$element)]

# Fill missing matches if requested
aligned_counts[is.na(aligned_counts)] <- fill_value

# Output aligned dataframe
align_loaded_db <- data.frame(
element = one_skeleton_element,
count = aligned_counts
)

return(align_loaded_db)
}
########################################### tapho-r loader function sum=0 ##############
tapho.load <- function( data, 
                        target_levels        = NULL,
                        calc_level   = c("min", "max", "mode"),
                        indet        = TRUE,
                        indet_method = c("up", "dist", "delet"),
                        fragments_per_bone        = 3,
                        meta,
                        strict = FALSE) {

    
  #Validate & coerce first (warn now, not later)
  new_site_data <- archo_preload_site_table(data, strict = strict)

  ind_name <- colnames(new_site_data)[2]

  colnames(new_site_data) <- c("element", "count")

  new_site_meta <- new_site_data[1:31, , drop = FALSE]
  new_site_raw  <- new_site_data[-c(1:31), , drop = FALSE]

  # NA skeletal values
  new_site_raw$count <- as.numeric(new_site_raw$count)
  new_site_raw$count[is.na(new_site_raw$count)] <- 0
  #MNI
  MNI <- as.numeric(new_site_meta[18,]$count)


  #_indet/_fragments out/in
  new_site_raw <- resolve_indeterminate(new_site_raw         = new_site_raw, 
                                        MNI          = MNI,
                                        indet        = indet,
                                        indet_method = indet_method,
                                        fragments_per_bone = fragments_per_bone)



  # Step 3: Resolve level — compute or use provided
  if (is.null(target_levels)) {
    target_levels <- summarize_bone_sections_zero_na(
      input_df       = new_site_raw,
      method = calc_level,
      meta
    )
  }

  recalc_new_site_raw <- recal_mne_six(
    current_mne = as.numeric(new_site_raw$count),
    bone_names = new_site_raw$element,
    target_levels = target_levels,
    MNI = MNI
  )



  one_skeleton <- get_bones_dataframe(target_levels)

  data_aligned <- align_counts_to_skeleton(
    one_skeleton$element,
    recalc_new_site_raw$element,
    recalc_new_site_raw$count
  )

  #  FIXED CALL SIGNATURE: pass meta DF, not meta$count
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

########################################### tapho-r loader function sum=False, ie db and sheets ##############
tapho.load.db <- function(data,
                          level        = NULL,
                          calc_level   = c("min", "max", "mode"),
                          indet        = TRUE,
                          indet_method = c("up", "dist", "delet"),
                          fragments_per_bone        = 3,
                          meta
) {
  tb <- tapho.load( data         = data[, 1:2],
                    target_levels        = level,
                    calc_level   = calc_level,
                    indet        = indet,
                    indet_method = indet_method,
                    fragments_per_bone        = fragments_per_bone,
                    meta
  
  )

  if (tb$count_type == FALSE) { 
  for (i in 3:ncol(data)) {
    one_column <- data[, c(1,i)]
    #one_column <- one_column %>% rename(count = !!paste0("count", i-1))
    tb_one <- tapho.load( data         = one_column,
                          target_levels        = level,
                          calc_level   = calc_level,
                          indet        = indet,
                          indet_method = indet_method,
                          fragments_per_bone        = fragments_per_bone, 
                          meta
    )
    tb <- tb %>% add_row(tb_one)
  }}

  return(tb)
}



library(dplyr)

# Define the function to read an Excel file with multiple sheets
read_excel_sheets <- function(
  file_path,
  level = NULL,
  calc_level = c("min", "max", "mode"),
  indet = TRUE,
  indet_method = c("up", "dist", "delet"),
  fragments_per_bone = 3,
  meta,
  verbose = TRUE
) {

  # --------------------------------------------------
  # Get Excel sheets
  # --------------------------------------------------

  sheet_names <- readxl::excel_sheets(
    file_path
  )

  n_sheets <- length(
    sheet_names
  )

  total_start <- Sys.time()


  # --------------------------------------------------
  # Create ONE random icon sequence for this run
  # --------------------------------------------------

  run_icons <- make_tapho_icons()


  # --------------------------------------------------
  # Initial progress line
  # --------------------------------------------------

  if (verbose) {

    status_line <- make_tapho_status(
      i = 0,
      n = n_sheets,
      site = "starting...",
      icons = run_icons
    )

    cat(
      "\r",
      status_line,
      sep = ""
    )

    flush.console()
  }


  # --------------------------------------------------
  # Process sheets
  # --------------------------------------------------

  processed_list <- lapply(
    seq_along(sheet_names),
    function(i) {

      sname <- sheet_names[i]


      # ----------------------------------------------
      # Read Excel sheet
      # ----------------------------------------------

      sheet_data <- readxl::read_excel(
        file_path,
        sheet = sname
      )


      # ----------------------------------------------
      # TAPHO processing
      # ----------------------------------------------

      processed_sheet <- tapho.load.db(
        data = sheet_data,
        level = level,
        calc_level = calc_level,
        indet = indet,
        indet_method = indet_method,
        fragments_per_bone =
          fragments_per_bone,
        meta = meta
      )


      # ----------------------------------------------
      # Update SAME console line
      # ----------------------------------------------

      if (verbose) {

        status_line <- make_tapho_status(
          i = i,
          n = n_sheets,
          site = sname,
          icons = run_icons
        )

        # \r returns cursor to beginning
        # No \n here: stay on same line
        cat(
          "\r",
          status_line,
          sep = ""
        )

        flush.console()
      }


      processed_sheet
    }
  )


  # --------------------------------------------------
  # Total elapsed time
  # --------------------------------------------------

  total_time <- as.numeric(
    difftime(
      Sys.time(),
      total_start,
      units = "secs"
    )
  )


  # --------------------------------------------------
  # Final permanent line
  # --------------------------------------------------

  if (verbose) {

    final_icons <- paste0(
      run_icons,
      collapse = ""
    )

    final_line <- sprintf(
      paste0(
        "Tapho-loaded: ",
        "%s  ",
        "%d / %d sites · ",
        "finished in %.1f sec"
      ),
      final_icons,
      n_sheets,
      n_sheets,
      total_time
    )

    # Overwrite current progress line
    # and finally add newline
    cat(
      "\r",
      final_line,
      "\n",
      sep = ""
    )

    flush.console()
  }


  # --------------------------------------------------
  # Combine sheets
  # --------------------------------------------------

  combined_df <- dplyr::bind_rows(
    processed_list
  )

  return(
    combined_df
  )
}

#______ append mne derived var data __________

# ... inside 
#___________________________________________



#tapho.loader <- function(
  # 1. DATA INPUT — always first, it's the primary required argument
#  data,

  # 2. LEVEL CONFIGURATION — core processing logic, used deepest in the call stack
#  level        = NULL,
#  calc_level   = c("min", "max", "mode"),

  # 3. INDETERMINATE HANDLING — another core processing group
#  indet        = TRUE,
#  indet_method = c("up", "dist", "delet"),

  # 4. FILTERING / SELECTION — applied after data is loaded
# group_by     = NULL,
#  select       = NULL,

  # 5. OUTPUT OPTIONS — last, as it affects the final result shape
#  sum          = TRUE
#){
  
  #1 read_excel_sheets(data, calc_level, level, indet, indet_method)

  #2 elect & group_by(data, group_by, select, level, sum, indet)

  #3 convertMNEtovariables(data, indet)

  # return data
#}

tapho.loader.v1 <- function(
  data,
  sum          = TRUE,
  unique_level= TRUE,
  level        = NULL,
  calc_level   = c("min", "max", "mode"),
  recalc_method =  c("min", "filter", "direct"), 
  indet        = TRUE,
  indet_method = c("up", "dist", "delet"),
  group_by     = NULL,
  select       = NULL,
  fragments_per_bone = 3, # all_levels = FALSE /TRUE, 
  exclude_sections = c("teeth") 
) {



#_________________________________   META ENV   _____________________________________
  
meta <- meta_set(indet)

on.exit({
  if ("meta_env" %in% search()) {
    detach("meta_env")
  }
}, add = TRUE)
 
#_______________________________________________________________________________________




  calc_level   <- match.arg(calc_level)
  indet_method <- match.arg(indet_method)

  # Step 1: Load and process all sheets
  db <- read_excel_sheets(
    file_path    = data,
    level        = level,
    calc_level   = calc_level,
    indet        = indet,
    indet_method = indet_method,
    fragments_per_bone = fragments_per_bone,
    meta
  )
  
  # If sum = TRUE and "site_id" not already in group_by, prepend it
if (sum) {
  if (is.null(group_by)) {
    group_by <- "site_id"
  } else if (!"site_id" %in% group_by) {
    group_by <- c("site_id", group_by)
  }
}

  # Step 2: Apply grouping and selection
  db <- select_and_group_by(
    data     = db,
    group_by = group_by,
    select   = select,
    method =  recalc_method,
    level    = level,
    unique_level= unique_level
  )

  # need one_skeleton element at all times 

  # dist down via pesos!!   /solve indet TRUE -> del.(check)

  # regroup / ratio. /...

  # get section info then get one element sum up count and one_element and rename to ... -> feature engineering. 
  # updata $count $elemet $bone.level 
  # select bone section or eliminate eg -c("teeth") take section 
  # eliminate each features corresponding to section named in slote
  # for each site_i -> take element count add skeleton_section via add_bone_info.
  # add to db[sit_i, ] <- skeleton_section 
  # eliminate each section in elimininate_section eg -c("teeth")

  # feature creation. sum "element" (def). sum AND one_skeleton. 
  # -> function(element, count , one_skeleton, skeleton_section)
  # -> 1. add "vertebrae", "thorax", "pelvic", and "shoulder", "arms", "legs", 
  # -> 2. add extrimities / post-cranium. 
  # -> 3. add ratios. list. + same * bone denity * onehot-metas  

  # replace elemtent count, one-skeleton count. 
  # -> add convertMNEtoVar. AND one_skeleton, skeleton_section 
  
  input_df <- data.frame(element = db[1,]$element[[1]] , count = db[1,]$count[[1]])
  data <- input_df %>%
    left_join(
      elements_data %>% select(ElementsOne, Bone_Count, Skeleton_Section),   
      by = c("element" = "ElementsOne")  
    ) 
  
 db <- data %>%
  filter(!Skeleton_Section %in% exclude_sections)

  # take index of teeth and delet in list element and list count and list one_skeleton
  
  # need count, element, OneElement-count, -> vie add_bone_info()$one_element?
  # ...  
  #MNI = 2
  #skeleton <- data.frame(element = data$element , count = data$Bone_Count)
  # Step 3: Convert MNE to variables
  # for each site_i take mne, one_skeleton ->  
  #var_db <- convertMNEtoVariables_v3(
  # MNE = data$count, 
  #skeleton = skeleton, 
  #individuals = MNI
  #)
  
  
   results_list <- lapply(seq_len(nrow(db)), function(i) {
  
    # Build input_df for this row
    input_df <- data.frame(
      element = db[i, ]$element[[1]],
      count   = db[i, ]$count[[1]]
    )

    # Join bone info
    data <- input_df %>%
      left_join(
        elements_data %>% select(ElementsOne, Bone_Count, Skeleton_Section),
        by = c("element" = "ElementsOne")
      ) #%>%
      #filter(!Skeleton_Section %in% exclude_sections)
    
    # Guard: if data is empty after filtering, return NULL
    if (nrow(data) == 0) return(NULL)
     
    MNI      <- db[i, ]$mni  # use real MNI from db instead of hardcoded 2
    skeleton <- data.frame(element = data$element, count = data$Bone_Count)
    
    var_db <- convertMNEtoVariables_v3(
      MNE        = data$count,
      skeleton   = skeleton,
      individuals = MNI
    )

    var_db
   })
  db <- db %>%
  mutate(mne_features = results_list)
  db <- db %>%
    mutate(
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
    select(-mne_features)  # drop the nested tibble column

  return(db)
}

# --- > time calculator sum over 
# plot mne vs desity (need fuente). ->
# plot mne vs %
# plot pca -> mda
# -> ml pipeline -> level 1,2,3, -1,-2 + ratios. -> recursive feature eleimination.

# select feature count , and label ie mne and/or ratios

# smote ind -sum
# delet section,eg -c(teeth)
# how to add ratios and bone desities /  same * bone denity * onehot-metas  
# need one_skeleton$count and $element 
# all level at once ? -> $element $count + one_skeleton$element $count + ratios + bone_desity --- outside: one_hot_metas.
#  -> after  ->
# ....






tapho.loader.v2 <- function(
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
  exclude_sections  = c("teeth")
) {

  # ── Setup ──────────────────────────────────────────────────────────────────
  meta_set(indet)
  on.exit(detach("meta_env"), add = TRUE)

  calc_level   <- match.arg(calc_level)
  indet_method <- match.arg(indet_method)

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

  # ── Stage 2: Group and aggregate by site (and optional extra variables) ────
  if (sum) {
    group_by <- if (is.null(group_by)) "site_id" else union("site_id", group_by)
  }

  db <- select_and_group_by(
    data         = db,
    group_by     = group_by,
    select       = select,
    method       = recalc_method,
    level        = level,
    unique_level = unique_level
  )

  # ── Stage 3: Compute MNE-derived variables for each site row ──────────────
  results_list <- lapply(seq_len(nrow(db)), function(i) {

    input_df <- data.frame(
      element = db$element[[i]],
      count   = db$count[[i]]
    )
  
    # Build bone_data as its own pipeline — no chain into the if
    bone_data <- input_df %>%
      left_join(
        elements_data %>% select(ElementsOne, Bone_Count, Skeleton_Section),
        by = c("element" = "ElementsOne")
      ) %>%
      filter(!Skeleton_Section %in% exclude_sections)
  
    # Guard is now a plain if, not part of the pipe
    if (nrow(bone_data) == 0) return(NULL)
  
    skeleton <- data.frame(element = bone_data$element, count = bone_data$Bone_Count)
    MNI      <- db$mni[[i]]
  
    convertMNEtoVariables_v3(
      MNE         = bone_data$count,
      skeleton    = skeleton,
      individuals = MNI
    )
  })

  # ── Stage 4: Unpack results into top-level columns and return ──────────────
  db %>% select(-c(count,element)) %>%
    mutate(mne_features = results_list) %>%
    mutate(
      # update sum_count /bone.level / 
      element                           = lapply(mne_features, `[[`, "element_2"),
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
    select(-mne_features)
}

# __________________ tapho_loader_ml_wfl ____________________

tapho.loader.wfl <- function(
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
  group_by_sum <- if (is.null(group_by)) "site_name" else union("site_name", group_by)
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

    # Find task_name column index dynamically
    #task_col <- which(names(d) == "task_name")

    d %>%
      select(-c(count, element)) %>%
      mutate(mne_features = results_list) %>%
      mutate(
        element                           = lapply(mne_features, `[[`, "element_2"),
        One_Skeleton_count                = lapply(mne_features, `[[`, "One_Skeleton_count"),   # new
        Skeleton_Section                  = lapply(mne_features, `[[`, "Skeleton_Section"),     # new
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
      select(-mne_features) %>%
      # ── Stage 4: Unnest task_name using dynamic column position ─────────────
    {
      tbl <- .
      task_col_final <- which(names(tbl) == data_variant)
      tbl %>%
        select(1:37,task_col_final) %>%
        mutate(!!data_variant := lapply(seq_along(.data[[data_variant]]), function(i)
          setNames(.data[[data_variant]][[i]], tbl$element[[i]])
        )) %>%
        unnest_wider(!!sym(data_variant))
    }
  }

  # ── Stage 4: Run pipeline for both groupings ───────────────────────────────
  db_sum <- process_db(group_by_sum)
  db_id  <- process_db(group_by_id)

  # ── Stage 5: Return both as a named list ──────────────────────────────────
  list(db_id = db_id, db_sum = db_sum)
}

# ============================================================
# HIERARCHY-AWARE ANATOMICAL ALLOCATION
# ============================================================

# 1. Load physical-slot model and matching algorithms
source("anatomical_allocator_proposed.R")

# 2. Preserve the original indeterminate distributor
#    for special pelvic/cranium handling
distribute_indet_legacy <- distribute_indet

# 3. Install the hierarchy-aware replacements
source("shared_allocator_function_replacements.R")
# ============================================================
# HIERARCHY-AWARE ANATOMICAL ALLOCATION
# ============================================================

tapho.loader.v3 <- function(
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
  verbose            = TRUE
) {

  # ── Setup ──────────────────────────────────────────────────────────────────
  meta_set(indet)
  on.exit(detach("meta_env"), add = TRUE)

  calc_level   <- match.arg(calc_level)
  indet_method <- match.arg(indet_method)

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


  # ── Stage 2: Group and aggregate by site and optional extra variables ──────
  if (sum) {
    group_by <- if (is.null(group_by)) "site_id" else union("site_id", group_by)
  }


  db <- select_and_group_by(
    data         = db,
    group_by     = group_by,
    select       = select,
    method       = recalc_method,
    level        = level,
    unique_level = unique_level
  )

  # Lookup used later to attach section and one-skeleton count
  # to the final output elements.
  element_lookup <- elements_data %>%
    dplyr::transmute(
      element_key          = gsub(" ", ".", ElementsOne),
      One_Skeleton_count   = Bone_Count,
      Skeleton_Section     = Skeleton_Section
    )

  # ── Stage 3: Compute MNE-derived variables for each site row ───────────────
  results_list <- lapply(seq_len(nrow(db)), function(i) {

    input_df <- data.frame(
      element = db$element[[i]],
      count   = db$count[[i]]
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
  
    MNI <- db$mni[[i]]
  
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
  
    var_db <- var_db %>%
      dplyr::mutate(
        One_Skeleton_count = as.numeric(out_one_skeleton_count),
        Skeleton_Section   = as.character(out_sections)
      )
  
    var_db
  })

    # ── Stage 4: Unpack results into top-level columns ──────────────────────────
    out <- db %>%
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
  
    derived_features <- purrr::pmap_dfr(
      list(
        element             = out$element,
        One_Skeleton_count  = out$One_Skeleton_count,
        Skeleton_Section    = out$Skeleton_Section,
        MNE                 = out$MNE,
        MNI_used            = out$MNI_used_assemblage
      ),
      make_tapho_derived_features
    )
  
    dplyr::bind_cols(out, derived_features)
  }



  tapho.loader <- function(
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
    verbose            = TRUE
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
      meta               = meta,
      verbose = verbose
    )
  
###################
library(dplyr)
library(tidyr)
print(colnames(db))
print(db$bone.level)
bone_level_out <- do.call(rbind, db$bone.level)

bone_level_out <- as.data.frame(bone_level_out)

names(bone_level_out) <- paste0("bone.level_", 1:11)

write.csv(
  bone_level_out,
  "bone_levels_2.csv",
  row.names = FALSE
)
###################    
  
    # ── Stage 2: Group and aggregate by site and optional extra variables ──────
    if (sum) {
      group_by <- if (is.null(group_by)) "site_id" else union("site_id", group_by)
    }
  

  
    db <- select_and_group_by(
      data         = db,
      group_by     = group_by,
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
  
    # ── Stage 3: Compute MNE-derived variables for each site row ───────────────
    results_list <- lapply(seq_len(nrow(db)), function(i) {
  
      input_df <- data.frame(
        element = db$element[[i]],
        count   = db$count[[i]]
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
  
      MNI <- db$mni[[i]]
  
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
  
      var_db <- var_db %>%
        dplyr::mutate(
          One_Skeleton_count = as.numeric(out_one_skeleton_count),
          Skeleton_Section   = as.character(out_sections)
        )
  
      var_db
    })
  
    # ── Stage 4: Unpack results into top-level list-columns ────────────────────
    out <- db %>%
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
  
    # ── Stage 5: Existing engineered taphonomic features ───────────────────────
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
  
    # ── Stage 6: New 16 density-informed engineered features ───────────────────
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
  
    # ── Stage 7: Final output ──────────────────────────────────────────────────
    dplyr::bind_cols(
      out,
      derived_features,
      density_features
    )
  }
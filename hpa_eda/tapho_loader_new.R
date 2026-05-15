
source("select_data_new.R")
source("data_loader_new.R")
source("archo_preload_2.R")
source("collapse_down_step.R")
source("indet.R")
source("meta/meta_set.R")

source("convertMNEtovariables_new.R")

library(dplyr)
library(tibble)
library(readxl)
library(purrr)

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
  if (n_indet == 0) return(counts)
  open_slots <- max_counts - counts
  if (sum(open_slots, na.rm=TRUE) < as.numeric(n_indet))
    stop("Not enough open slots to place all indeterminate specimens.")
  slot_pool <- rep(seq_along(counts), times = open_slots)
  chosen    <- sample(slot_pool, size = n_indet, replace = FALSE)
  for (i in chosen) counts[i] <- counts[i] + 1
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



collapse_down_step <- function(data, part, current_level, target_level,section_name, MNI) {
  key         <- paste0(part, "_", target_level, "_to_", current_level)
  if (!key %in% names(category_hierarchy)) stop("Missing mapping: ", key)

  mapping      <- category_hierarchy[[key]]
  level_label  <- paste0("Level ", target_level)
  level_maxima <- skeleton_hierarchy[[section_name]][[level_label]]
#print(level_maxima)
  new_rows <- purrr::map_dfr(names(mapping), function(parent_elem) {

    parent_clean <- gsub(" ", ".", parent_elem)
    children     <- gsub(" ", ".", mapping[[parent_elem]])
    real_children <- drop_indet_if_possible(children)

    parent_count <- data %>%
      filter(std_element == parent_clean) %>%
      pull(count) %>%
      { if (length(.) == 0) 0L else . }

    existing_counts <- purrr::map_int(real_children, \(ch) {
      val <- data %>% filter(std_element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val)
    })

    max_counts <- purrr::map_int(real_children, \(ch) {
      base <- level_maxima[[ch]]
     # print("###############################    WTF      ########################################")
     # print("base :")
     # print(str(base))
     # print(base)
     # print("MNI : ")
     # print(str(MNI))
     # print(MNI)
     # print("#######################################################################")
      if (is.null(base)) 0L else as.integer(base * MNI)
    })


  if (length(real_children) == 1 && real_children[1] == parent_clean) {
    #if (part == "pelvic") {print("Pelvic: length(real_children) == 1 && real_children[1] == parent_clean") 
     # print(real_children)}
            return(tibble(std_element = real_children, count = as.integer(parent_count)))
      }

    surplus      <- max(0L, parent_count) # no double - sum(existing_counts))
    final_counts <- distribute_surplus(existing_counts, max_counts, surplus)
    #if (part == "pelvic") {}

    tibble(std_element = real_children, count = final_counts)
  })

  parent_names <- gsub(" ", ".", names(mapping))
  child_names <- unique(gsub(" ", ".", unlist(mapping)))
  child_names <- drop_indet_if_possible(child_names)
  #if (part == "pelvic") {print("Pelvic: child_names & parent_names ") 
  #print(parent_names)
  #  print(child_names)
  #}
  data %>%
    filter(!std_element %in% c(parent_names, child_names)) %>%
    bind_rows(new_rows)

}

is_special <- function(part) part %in% c("pelvic", "cranium")

recal_mne_six <- function(current_mne, bone_names, target_levels, MNI) {
  print("six")
  print(names(elements_data))
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
    #print("up: from to ")
    #print(current_max)
    #print(step_up_level)
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
    #print(data, n=42)
    }

    # Pass 2 — expand from coarsest detail up to target
    while (current_min < tgt) {
      #print("down: from to ")
       #   print(current_min)
    #print(current_min + 1)
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
    
    #print(data, n=42)
   # if (part == "pelvic") { print(data$)}
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

  #print(new_site_raw)

  # Step 3: Resolve level — compute or use provided
  if (is.null(target_levels)) {
    target_levels <- summarize_bone_sections_zero_na(
      input_df       = new_site_raw,
      method = calc_level,
      meta
    )
  }


  #print("########################################### tapho-r loader function sum=0 ##############")
  #print(new_site_raw$count)
  #print("########################################### tapho-r loader function sum=0 ##############")

  #mni <- new_site_data[18,2]
  #mni <- suppressWarnings(as.numeric(new_site_data$count[18]))

  recalc_new_site_raw <- recal_mne_six(
    current_mne = as.numeric(new_site_raw$count),
    bone_names = new_site_raw$element,
    target_levels = target_levels,
    MNI = MNI
  )

  #print("########################################### tapho-r loader function sum=0 ##############")
  #print(recalc_new_site_raw)
  #print("########################################### tapho-r loader function sum=0 ##############")

  one_skeleton <- get_bones_dataframe(target_levels)
  #(one_skeleton)
  # print(recalc_new_site_raw)
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
read_excel_sheets <- function(file_path, 
                              level        = NULL,
                              calc_level   = c("min", "max", "mode"),
                              indet        = TRUE,
                              indet_method = c("up", "dist", "delet"),
                              fragments_per_bone = 3,
                              meta
) {
  
  # Get all sheet names in the Excel file
  sheet_names <- readxl::excel_sheets(file_path)
  #print(sheet_names)
  # Read each sheet, process it, and store in a list
  processed_list <- lapply(sheet_names, function(sname) {
    
    # Read the current sheet
    sheet_data <- readxl::read_excel(file_path, sheet = sname)
    
    # Process it using your tapho.load.db() function
    processed_sheet <- tapho.load.db(
                                        data         = sheet_data,
                                        level        = level,
                                        calc_level   = calc_level,
                                        indet        = indet,
                                        indet_method = indet_method,
                                        fragments_per_bone = fragments_per_bone,
                                        meta
    )
      
    
    
    # Optionally, you might want to keep track of which sheet it came from
    # processed_sheet <- processed_sheet %>% 
    #   mutate(sheet_name = sname)
    
    processed_sheet
  })
  
  # Combine all processed sheets into one tibble
  combined_df <- bind_rows(processed_list)
  
  return(combined_df)
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
  
meta_set(indet)          # initial load
on.exit(detach("meta_env"), add = TRUE)  # keep this for cleanup
 
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
  #print(db[,c(1:5,33:36)], n=26)
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

  
  input_df <- data.frame(element = db[1,]$element[[1]] , count = db[1,]$count[[1]])
  data <- input_df %>%
    left_join(
      elements_data %>% select(ElementsOne, Bone_Count, Skeleton_Section),   
      by = c("element" = "ElementsOne")  
    ) 
  
 db <- data %>%
  filter(!Skeleton_Section %in% exclude_sections)
  
  
   results_list <- lapply(seq_len(nrow(db)), function(i) {
  
    # Build input_df for this row
    input_df <- data.frame(
      element = db[i, ]$element[[1]],
      count   = db[i, ]$count[[1]]
    )
    #print(length(input_df$element))
    #print(length(input_df$count))
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
    print("test lenght elemen & count : ")
    print(length(var_db$MNE[[1]]))
     print(length(data$element))
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
 #print(results_list)
  return(db)
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
print(db)
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

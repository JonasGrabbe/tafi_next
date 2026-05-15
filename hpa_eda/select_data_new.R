source("recalMNE_new.R")

process_data <- function(data) {

  data_sum0 <- data %>% filter(sum == 0)
  data_sum1 <- data %>% filter(sum == 1)

  # If there is sum=0 data, process it
  if(nrow(data_sum0) > 0) {

  check_mixed_cols <- c( "time_range", "geological_context",
                          "accumulation_cause", "gender", "age", "hominin_species")

  data_sum0 <- data %>% filter(sum == 0)
  data_sum1 <- data %>% filter(sum == 1)

  data_sum0 <- recal_mne_selected_data(data_sum0)

  data_sum0_summarized <- data_sum0 %>%
    group_by(site_index) %>%
      summarise(
        site_index = first(site_index),
        sample_index = first(sample_index),
        site_name = first(site_name),
        sum = 1,
        across(all_of(check_mixed_cols), ~ if (n_distinct(.x) == 1) first(.x) else "mixed"),
        number_of_individuals = sum(number_of_individuals, na.rm = TRUE),
        bone.level = list(first(bone.level)),
        element = list(first(element)),
        count = list(reduce(count, function(x, y) {
          x_vec <- unname(unlist(unname(x)))
          y_vec <- unname(unlist(unname(y)))
          # Get all unique names (indices)
          max_length <- max(length(x_vec), length(y_vec))
          # Pad shorter vector with zeros
          x_vec <- c(x_vec, rep(0, max_length - length(x_vec)))
          y_vec <- c(y_vec, rep(0, max_length - length(y_vec)))
          # Sum the vectors
          unname(x_vec + y_vec)
        })),
        .groups = 'drop'
      )
  data_sum0_summarized$element <- lapply(data_sum0_summarized$element, function(x) {
    sub(".*_", "", x)  # Remove everything before the underscore
  })
  data_combined <- bind_rows(data_sum1, data_sum0_summarized)
  data_combined$count <- lapply(data_combined$count, function(x) {
    unname(x)
  })
  } else {

    data_combined <- data_sum1

  }
  return(data_combined)
  }





select_data <- function(
  data,
  time_range_filter = NULL,
  geological_context_filter = NULL,
  number_of_individuals_filter = NULL,
  age_filter = NULL,
  gender_filter = NULL,
  accumulation_cause_filter = NULL,
  bone_parts_filter = NULL,
  site_name_filter = NULL,
  sum_filter = 1,
  hominin_species_filter = NULL
) {
  library(dplyr)

  data_filtered <- data

  # Filter based on metadata fields
  if (!is.null(time_range_filter)) {
    data_filtered <- data_filtered %>%
      filter(time_range %in% time_range_filter)
  }

  if (!is.null(geological_context_filter)) {
    data_filtered <- data_filtered %>%
      filter(geological_context %in% geological_context_filter)
  }

  if (!is.null(number_of_individuals_filter)) {
    if (length(number_of_individuals_filter) == 1) {
      data_filtered <- data_filtered %>%
        filter(number_of_individuals == number_of_individuals_filter)
    } else {
      data_filtered <- data_filtered %>%
        filter(number_of_individuals >= min(number_of_individuals_filter),
               number_of_individuals <= max(number_of_individuals_filter))
    }
  }

  if (!is.null(age_filter)) {
    data_filtered <- data_filtered %>%
      filter(age %in% age_filter)
  }

  if (!is.null(gender_filter)) {
    data_filtered <- data_filtered %>%
      filter(gender %in% gender_filter)
  }

  if (!is.null(accumulation_cause_filter)) {
    data_filtered <- data_filtered %>%
      filter(accumulation_cause %in% accumulation_cause_filter)
  }

  if (!is.null(site_name_filter)) {
    data_filtered <- data_filtered %>%
      filter(site_name %in% site_name_filter)
  }



  if (!is.null(hominin_species_filter)) {
    if ("hominin_species" %in% colnames(data_filtered)) {
      data_filtered <- data_filtered %>%
        filter(hominin_species %in% hominin_species_filter)
    } else {
      warning("The 'hominin_species' field is not present in the data.")
    }
  }

  # Filter based on 'element' without unnesting 'bone.level'
  if (!is.null(bone_parts_filter)) {
    data_filtered <- data_filtered %>%
      rowwise() %>%
      filter(any(element %in% bone_parts_filter)) %>%
      ungroup()
  }

  if (!is.null(sum_filter)) {
    if (sum_filter == 0) {
    data_filtered <- data_filtered %>%
      filter(sum %in% sum_filter)
  }
  else {
    data_filtered <- process_data(data_filtered)
    data_filtered <- data_filtered  %>%
      mutate(bone.level = map(bone.level, function(x) {
        if (is.list(x) && !is.null(names(x))) {
          return(x)
        } else if (is.numeric(x)) {
          # Create named list with standard bone categories
          bone_categories <- c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
          "arms", "hands", "pelvic", "legs", "foot")
          return(setNames(as.list(x), bone_categories))
        } else {
          return(x)
        }
      }))
    }
  }


  return(data_filtered)
}


############################


get_min <- function(selected_data) {

  # Convert the 'bone.level' list-column to a data frame
  bone_levels_df <- selected_data %>%
    # Select only the 'bone.level' column to simplify
    select(bone.level) %>%
    # Unnest the 'bone.level' list into wider format
    unnest_wider(bone.level, names_sep = "_")

  # Calculate the minimum for each bone section
  min_bone_level <- bone_levels_df %>%
    summarise(across(everything(), ~ min(.x, na.rm = TRUE)))

  # Convert to a named vector for easier use
  min_bone_level_vec <- as.numeric(min_bone_level)
  names(min_bone_level_vec) <- colnames(min_bone_level)

  # Display the result
  return (unname(min_bone_level_vec))
  }


##########################


recalculate_counts <- function(data_filtered) {
  library(dplyr)
  library(tidyr)
  library(purrr)

  #### Integration of recal_mne function ####

  # Calculate the minimum bone level across the filtered data
  min_bone_level <- get_min(data_filtered)

  # Recalculate 'bone.level' to the minimum bone level
  data_recalculated <- data_filtered #%>%
    #mutate(bone.level = list(min_bone_level))

  # Recalculate counts and elements for each row
  recalculated_list <- pmap(
    list(data_recalculated$count, data_recalculated$bone.level, data_recalculated$element),
    function(count_i, bone_level_i, element_i) {
      recal_mne(
        current_mne = unlist(count_i),
        current_levels = bone_level_i,
        target_levels = min_bone_level,
        skeleton_hierarchy = skeleton_hierarchy,
        bone_names = unlist(element_i),
        sections_order = sections_order
      )
    }
  )

  # Extract elements and counts from recalculated results
  data_recalculated <- data_recalculated %>%
    mutate(
      recalculated = recalculated_list,
      bone.level = list(min_bone_level),
      element = map(recalculated, "element"),
      count = map(recalculated, "count")
    ) %>%
    select(-recalculated)

  # Ensure that all samples have counts for the same set of elements
  # Collect all unique elements across the data
  unique_elements <- data_recalculated %>%
    pull(element) %>%
    unlist() %>%
    unique()

  # Sort the unique elements for consistency
  unique_elements <- sort(unique_elements)

  # Align counts across all rows by creating a consistent element list
  data_recalculated <- data_recalculated %>%
    mutate(
      count_vector = map2(element, count, ~ setNames(.y, .x)),
      aligned_counts = map(count_vector, ~ {
        counts <- .x[unique_elements]
        counts[is.na(counts)] <- 0
        counts
      })
    )

  # Create a wide-format data frame with elements as columns
  counts_df <- data_recalculated %>%
    select(aligned_counts) %>%
    mutate(sample_id = row_number()) %>%
    unnest_wider(aligned_counts)

  # Combine the counts with the original data
  data_recalculated <- data_recalculated %>%
    select(-count_vector, -aligned_counts, -count, -element) %>%
    bind_cols(counts_df %>% select(-sample_id))

  # Now, 'data_recalculated' has counts for each element aligned across all samples
  return(data_recalculated)
}


##################################################


group_and_summarize <- function(data, group_cols){


  check_mixed_cols <- c("sum", "time_range", "geological_context",
                        "accumulation_cause", "gender", "age", "hominin_species") %>%
    setdiff(group_cols)


  merged_df <- data %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      across(all_of(check_mixed_cols), ~ if (n_distinct(.x) == 1) first(.x) else "mixed"),
      number_of_individuals = sum(number_of_individuals, na.rm = TRUE),
      bone.level = list(first(bone.level)),
      element = list(first(element)),
      count = list(reduce(count, ~ map2(.x, .y, `+`))),
     .groups = 'drop'
    )

  return(merged_df)
}



#### new ind -> to whole site data

sumby_site <- function(data){   
  # merging all ind to one sum=0 but need to calculcate min bone.level firest then recal all column first
  recalc_data <- data %>%
  group_by(site_id) %>%
  group_modify(~ recal_mne_selected_data(.x)) %>%  # apply function to each site_id subset
  ungroup()
  
  #print("recalc_data[,c(1:5,34:36)]:")
  #print(recalc_data[,c(1:5,34:36)], n=26)
  #print(recalc_data[19,34:36]$element)
  #print(recalc_data[20,34:36]$element)
  # check that all bone.level are the same  otherwise throuh error 
  sumby_data <- recalc_data %>% 
  # Create the composite string (Row-wise) before grouping
  mutate(obs_temp = paste(ind_name, individual_id, observation, sep = " ")) %>%
  
  group_by(site_id) %>%
  summarise(
      ind_name = paste(ind_name, collapse = " "),
      site_id = first(site_id), 
      individual_id = paste(individual_id, collapse = " "), 
      site_name = first(site_name), 
      
      # Layer logic: Unique or list of unique
      layer_name = if (n_distinct(layer_name) == 1) first(layer_name) else paste(unique(layer_name), collapse = ", "),
      
      count_type = TRUE, 
      geological_period = first(geological_period), 
      cultural_period = first(cultural_period),
      culture_period_2 = first(culture_period_2), 
      start_date_cal_bp = first(start_date_cal_bp),
      end_date_cal_bp = first(end_date_cal_bp), 
      mis_stage = first(mis_stage),
      karstic_system = first(karstic_system), 
      karstic_system_2 = first(karstic_system_2), 
      open_air_site = first(open_air_site), 
      open_air_site_2 = first(open_air_site_2), 
      region = first(region), 
      country = first(country), 
      mni = sum(mni, na.rm = TRUE), 
      
      # across(age_category, ~ if (n_distinct(.x) == 1) first(.x) else "Mixed"),
      # [UPDATED] Gender Logic: Uniform -> Indet presence -> Mixed
      across(age_category, ~ case_when(
      n_distinct(.x) == 1 ~ first(.x),       # All same? Keep it.
      "Indet" %in% .x     ~ "Indet",         # Any Indet? Force Indet.
      TRUE                ~ "Mixed"          # Diff values, no Indet? Mixed.
      )),

      n_infant = sum(n_infant),  
      n_juvenil = sum(n_juvenil),  
      n_subadult = sum(n_subadult),   
      n_adult  = sum(n_adult),  
      
      # [UPDATED] Gender Logic: Uniform -> Indet presence -> Mixed
      across(gender, ~ case_when(
      n_distinct(.x) == 1 ~ first(.x),       # All same? Keep it.
      "Indet" %in% .x     ~ "Indet",         # Any Indet? Force Indet.
      TRUE                ~ "Mixed"          # Diff values, no Indet? Mixed.
      )),
      
      taphonomic_context = first(taphonomic_context), 
      funerary_context = first(funerary_context), 
      discovery_context = first(discovery_context), 
      degree_of_excavation = first(degree_of_excavation), 
      hominin_species = first(hominin_species),
      
      # Collapse the pre-calculated row strings
      observation = paste(obs_temp, collapse = "; "),
      
      reference = first(reference), 
      sum_count = sum(sum_count, na.rm = TRUE), 
      bone.level = list(first(bone.level)), 
      element = list(first(element)),
      count = list(Reduce("+", count)) 
  )
  
  return(sumby_data)
  
}



sumby_feature <- function(data, group_by_col, filters = NULL) {
  
  group_sym <- rlang::sym(group_by_col)
  
  # Step 1: Apply filters if provided
  if (!is.null(filters)) {
    for (col in names(filters)) {
      col_sym <- rlang::sym(col)
      vals    <- filters[[col]]
      data    <- data %>% dplyr::filter(!!col_sym %in% vals)
    }
    
    if (nrow(data) == 0) {
      warning("No rows remain after filtering. Check your filter values.")
      return(tibble())
    }
  }
  
  # Step 2: Recalculate MNE within each group
  recalc_data <- data %>%
    group_by(!!group_sym) %>%
    group_modify(~ recal_mne_selected_data(.x)) %>%
    ungroup()
  
  # Step 3: Summarise by group
  sumby_data <- recalc_data %>%
    mutate(obs_temp = paste(ind_name, individual_id, observation, sep = " ")) %>%
    group_by(!!group_sym) %>%
    summarise(
      ind_name             = paste(ind_name, collapse = " "),
      site_id              = if (n_distinct(site_id) == 1) first(site_id) 
                             else paste(unique(site_id), collapse = ", "),
      individual_id        = paste(individual_id, collapse = " "),
      site_name            = if (n_distinct(site_name) == 1) first(site_name)
                             else paste(unique(site_name), collapse = ", "),
      layer_name           = if (n_distinct(layer_name) == 1) first(layer_name)
                             else paste(unique(layer_name), collapse = ", "),
      count_type           = TRUE,
      geological_period    = if (n_distinct(geological_period) == 1) first(geological_period)
                             else paste(unique(geological_period), collapse = ", "),
      cultural_period      = if (n_distinct(cultural_period) == 1) first(cultural_period)
                             else paste(unique(cultural_period), collapse = ", "),
      culture_period_2     = first(culture_period_2),
      start_date_cal_bp    = first(start_date_cal_bp),
      end_date_cal_bp      = first(end_date_cal_bp),
      mis_stage            = if (n_distinct(mis_stage) == 1) first(mis_stage)
                             else paste(unique(mis_stage), collapse = ", "),
      karstic_system       = first(karstic_system),
      karstic_system_2     = first(karstic_system_2),
      open_air_site        = first(open_air_site),
      open_air_site_2      = first(open_air_site_2),
      region               = if (n_distinct(region) == 1) first(region)
                             else paste(unique(region), collapse = ", "),
      country              = if (n_distinct(country) == 1) first(country)
                             else paste(unique(country), collapse = ", "),
      mni                  = sum(mni, na.rm = TRUE),
      
      age_category = case_when(
        n_distinct(age_category) == 1 ~ first(age_category),
        "Indet" %in% age_category     ~ "Indet",
        TRUE                          ~ "Mixed"
      ),
      
      n_infant   = sum(n_infant,   na.rm = TRUE),
      n_juvenil  = sum(n_juvenil,  na.rm = TRUE),
      n_subadult = sum(n_subadult, na.rm = TRUE),
      n_adult    = sum(n_adult,    na.rm = TRUE),
      
      gender = case_when(
        n_distinct(gender) == 1 ~ first(gender),
        "Indet" %in% gender     ~ "Indet",
        TRUE                    ~ "Mixed"
      ),
      
      taphonomic_context   = first(taphonomic_context),
      funerary_context     = first(funerary_context),
      discovery_context    = first(discovery_context),
      degree_of_excavation = first(degree_of_excavation),
      hominin_species      = if (n_distinct(hominin_species) == 1) first(hominin_species)
                             else paste(unique(hominin_species), collapse = ", "),
      observation          = paste(obs_temp, collapse = "; "),
      reference            = first(reference),
      sum_count            = sum(sum_count, na.rm = TRUE),
      bone.level           = list(first(bone.level)),
      element              = list(first(element)),
      count                = list(Reduce("+", count)),
      .groups = "drop"
    )
  
  return(sumby_data)
}



select_and_group_by <- function(data, 
                                group_by = NULL, 
                                select = NULL, 
                                level = NULL, 
                                sum = TRUE, 
                                indet = TRUE) {
  
  # ── Step 1: Apply filters (select) if provided ──────────────────────────────
  if (!is.null(select)) {
    for (col in names(select)) {
      col_sym <- rlang::sym(col)
      vals    <- select[[col]]
      data    <- data %>% dplyr::filter(!!col_sym %in% vals)
    }
    
    if (nrow(data) == 0) {
      warning("No rows remain after filtering. Check your select values.")
      return(tibble::tibble())
    }
  }
  
  # ── Step 2: Group and summarise ─────────────────────────────────────────────
  if (!sum) {
    # No summarisation requested — return filtered data as-is
    return(data)
  }
  
  if (is.null(group_by)) {
    # Default: group by site_id
    db <- sumby_site(data)
    
  } else if (length(group_by) == 1) {
    # Single grouping variable
    db <- sumby_feature(data, group_by_col = group_by)
    
  } else {
    # Multiple grouping variables: create a composite key, group, then drop it
    composite_col <- paste(group_by, collapse = "__x__")
    
    data <- data %>%
      dplyr::mutate(
        !!composite_col := paste(!!!rlang::syms(group_by), sep = " | ")
      )
    
    db <- sumby_feature(data, group_by_col = composite_col)
    
    # Split composite key back into individual columns
    db <- db %>%
      tidyr::separate(
        col  = !!composite_col,
        into = group_by,
        sep  = " \\| ",
        extra = "merge",
        fill  = "right"
      )
  }
  
  return(db)
}


select_and_group_by <- function(data,
  group_by = NULL,
  select   = NULL,
  method = c("min", "filter", "direct"), 
  level    = NULL, 
  unique_level= TRUE) {

# ── Step 1: Filter rows (select) ────────────────────────────────────────────
# Only `select` → filter, no aggregation
# Both NULL     → data passes through unchanged
if (!is.null(select)) {
for (col in names(select)) {
col_sym <- rlang::sym(col)
vals    <- select[[col]]
data    <- data %>% dplyr::filter(!!col_sym %in% vals)
}
if (nrow(data) == 0) {
warning("No rows remain after filtering. Check your select values.")
return(tibble::tibble())
}
}

# ── Step 2: Group and summarise ─────────────────────────────────────────────
# NULL group_by → return filtered (or unchanged) data as-is
if (is.null(group_by)) {

  if (unique_level) {
    return (recal_mne_selected_data(data, method = method, target_level = level))}
  else{
return(data)}
}

# Single or multiple group_by columns → create composite key if needed
if (length(group_by) > 1) {
composite_col <- paste(group_by, collapse = "__x__")
data <- data %>%
dplyr::mutate(
!!composite_col := paste(!!!rlang::syms(group_by), sep = " | ")
)
group_col <- composite_col
} else {
group_col <- group_by
}

group_sym <- rlang::sym(group_col)

# ── Step 3: Recalculate MNE within each group ───────────────────────────────

# if unique_level = TRUE get method level
  
# level and indet are passed down to recal_mne_selected_data()
if (unique_level) {
  # Apply recal_mne_selected_data to the whole dataset at once
  recalc_data <- recal_mne_selected_data(data, method = method, target_level = level)
  
} else {

    recalc_data <- data %>%
      dplyr::group_by(!!group_sym) %>%
      dplyr::group_modify(~ recal_mne_selected_data(.x, method = method, target_level = level)) %>%
      dplyr::ungroup()
  
}

# ── Step 4: Summarise ───────────────────────────────────────────────────────
sumby_data <- recalc_data %>%
dplyr::mutate(obs_temp = paste(ind_name, individual_id, observation, sep = " ")) %>%
dplyr::group_by(!!group_sym) %>%
dplyr::summarise(

ind_name          = paste(ind_name, collapse = " "),
site_id           = if (n_distinct(site_id) == 1) first(site_id)
else paste(unique(site_id), collapse = ", "),
individual_id     = paste(individual_id, collapse = " "),
site_name         = if (n_distinct(site_name) == 1) first(site_name)
else paste(unique(site_name), collapse = ", "),
layer_name        = if (n_distinct(layer_name) == 1) first(layer_name)
else paste(unique(layer_name), collapse = ", "),
count_type        = TRUE,
geological_period = if (n_distinct(geological_period) == 1) first(geological_period)
else paste(unique(geological_period), collapse = ", "),
cultural_period   = if (n_distinct(cultural_period) == 1) first(cultural_period)
else paste(unique(cultural_period), collapse = ", "),
culture_period_2  = first(culture_period_2),
start_date_cal_bp = first(start_date_cal_bp),
end_date_cal_bp   = first(end_date_cal_bp),
mis_stage         = if (n_distinct(mis_stage) == 1) first(mis_stage)
else paste(unique(mis_stage), collapse = ", "),
karstic_system    = first(karstic_system),
karstic_system_2  = first(karstic_system_2),
open_air_site     = first(open_air_site),
open_air_site_2   = first(open_air_site_2),
region            = if (n_distinct(region) == 1) first(region)
else paste(unique(region), collapse = ", "),
country           = if (n_distinct(country) == 1) first(country)
else paste(unique(country), collapse = ", "),
mni               = sum(mni, na.rm = TRUE),

age_category = dplyr::case_when(
n_distinct(age_category) == 1 ~ first(age_category),
"Indet" %in% age_category     ~ "Indet",
TRUE                          ~ "Mixed"
),

n_infant   = sum(n_infant,   na.rm = TRUE),
n_juvenil  = sum(n_juvenil,  na.rm = TRUE),
n_subadult = sum(n_subadult, na.rm = TRUE),
n_adult    = sum(n_adult,    na.rm = TRUE),

gender = dplyr::case_when(
n_distinct(gender) == 1 ~ first(gender),
"Indet" %in% gender     ~ "Indet",
TRUE                    ~ "Mixed"
),

taphonomic_context   = first(taphonomic_context),
funerary_context     = first(funerary_context),
discovery_context    = first(discovery_context),
degree_of_excavation = first(degree_of_excavation),
hominin_species      = if (n_distinct(hominin_species) == 1) first(hominin_species)
else paste(unique(hominin_species), collapse = ", "),
observation          = paste(obs_temp, collapse = "; "),
reference            = first(reference),
sum_count            = sum(Reduce("+", count), na.rm = TRUE) ,
bone.level           = list(first(bone.level)),
element              = list(first(element)),
count                = list(Reduce("+", count)),

.groups = "drop"
)

# ── Step 5: Restore multiple group_by columns from composite key ─────────────
if (length(group_by) > 1) {
sumby_data <- sumby_data %>%
tidyr::separate(
col   = !!composite_col,
into  = group_by,
sep   = " \\| ",
extra = "merge",
fill  = "right"
)
}

return(sumby_data)
}
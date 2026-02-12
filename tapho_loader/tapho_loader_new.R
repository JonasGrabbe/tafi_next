source("meta_new.R")
source("NewSiteDataLoader_new.R")
source("select_data_new.R")
source("data_loader_new.R")
source("archo_preload_2.R")

library(dplyr)
library(tibble)
library(readxl)

summarize_bone_sections_zero_na <- function(input_df, reference_df = elements_data) {
  
  # Filter out rows where input_df$element is 0 or NA

  input_df <- input_df %>%
    filter(!is.na(input_df$count) & input_df$count != 0)


  # Add bone info using the filtered input_df
  bone_info_df <- add_bone_info(input_df, reference_df)

  # Define the ordered list of skeleton sections
    sections_order <- c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
    "arms", "hands", "pelvic", "legs", "foot")

  # Summarize the minimum Level for each section
  summary_df <- bone_info_df %>%
    group_by(Skeleton_Section) %>%
    summarize(min_level = min(max_Level, na.rm = TRUE)) %>%
    ungroup()

  # Create a named vector with all sections, initializing with 1 (instead of NA)
  #result_vector <- setNames(rep(1, length(sections_order)), sections_order)
  most_detailed_level <- c(5,5,1,3,4,3,3,5,4,3,5)
  result_vector <- setNames(most_detailed_level, sections_order)

  # Fill in the result_vector with available data
  for (section in sections_order) {
    if (section %in% summary_df$Skeleton_Section) {
      level <- summary_df$min_level[summary_df$Skeleton_Section == section]
      # If the level is NA, keep it as 1, otherwise use the actual level
      result_vector[section] <- if(is.na(level)) 1 else level
    }else {
      print(section)
    # If the section is not in summary_df, it will remain 1
  }}

  # Return just the numeric vector without names
  return(unname(result_vector))
}

recal_mne_six <- function(current_mne,  target_levels, skeleton_hierarchy, bone_names, sections_order) {
  current_levels = c(5,5,1,3,4,3,3,5,4,3,5)


  # Initialize results list
  new_mne <- list()

  # Create initial tibble with bone info
  current_mne_named <- tibble(count = current_mne, element = bone_names) %>%
    add_bone_info() %>%
    mutate(std_element = gsub(" ", ".", element))

  # Process each skeleton section
  for (part in sections_order) {
    part_idx <- which(sections_order == part)
    current_level <- current_levels[part_idx]
    target_level <- target_levels[part_idx]

    cat("\nProcessing", part, ": Starting at level", current_level, "-> Target", target_level, "\n")

    # Initialize working data
    working_data <- current_mne_named %>%
      filter(Skeleton_Section == part) %>%
      select(std_element, count)

    iteration <- 1
    while(current_level > target_level) {
      intermediate_level <- max(target_level, current_level - 1)
      mapping_name <- paste0(part, "_", current_level, "_to_", intermediate_level)

      if (mapping_name %in% names(category_hierarchy)) {
        # Create element-category mapping
        element_map <- stack(category_hierarchy[[mapping_name]]) %>%
          setNames(c("child", "parent")) %>%
          mutate(
            child = gsub(" ", ".", child),
            parent = gsub(" ", ".", parent)
          )
         

        # Split data into mappable and non-mappable elements
        mappable <- working_data %>%
          inner_join(element_map, by = c("std_element" = "child"))

      
        non_mappable <- working_data %>%
          anti_join(element_map, by = c("std_element" = "child"))


        # Aggregate mappable elements
        aggregated <- mappable %>%
          group_by(parent) %>%
          summarize(count = sum(count, na.rm = TRUE)) %>%
          rename(std_element = parent)


        # Combine with non-mappable elements
        #working_data <- bind_rows(non_mappable, aggregated)
                # SPECIAL SUMMATION LOGIC

                  working_data <- non_mappable %>%
                    full_join(aggregated, by = "std_element") %>%
                    mutate(count = coalesce(count.x, 0) + coalesce(count.y, 0)) %>%
                    select(std_element, count)

        

        current_level <- intermediate_level
        cat("  Step", iteration, "| Level:", current_level, "| Elements:", nrow(working_data), "\n")
        iteration <- iteration + 1
      } else {
        warning(paste("Missing mapping:", mapping_name))
        break
      }
    }
 
    # Store final results
    if(nrow(working_data) > 0) {
      for(i in seq_len(nrow(working_data))) {
        #element_name <- paste(part, working_data$std_element[i], sep = "_")
        element_name <-  working_data$std_element[i]
        new_mne[[element_name]] <- working_data$count[i]
      }
    }
  }

  # Convert to tibble with proper ordering
  new_mne_df <- tibble(
    element = names(new_mne),
    count = unlist(new_mne)
  )

  return(new_mne_df)
}



align_counts_to_skeleton <- function(one_skeleton_element,
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
tapho.load <- function(new_site_data, strict = FALSE) {

  # Validate & coerce first (warn now, not later)
  new_site_data <- archo_preload_site_table(new_site_data, strict = strict)

  ind_name <- colnames(new_site_data)[2]

  colnames(new_site_data) <- c("element", "count")

  new_site_meta <- new_site_data[1:31, , drop = FALSE]
  new_site_raw  <- new_site_data[-c(1:31), , drop = FALSE]

  target_levels <- summarize_bone_sections_zero_na(new_site_raw)

  # NA skeletal values are already set to 0 by validator, but keep this safe:
  new_site_raw$count[is.na(new_site_raw$count)] <- 0



  recalc_new_site_raw <- recal_mne_six(
    as.numeric(new_site_raw$count),
    target_levels,
    skeleton_hierarchy,
    new_site_raw$element,
    sections_order
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
tapho.load.db <- function(data) {
  tb <- tapho.load(data[, 1:2])

  if (tb$count_type == FALSE) { 
  for (i in 3:ncol(data)) {
    one_column <- data[, c(1,i)]
    #one_column <- one_column %>% rename(count = !!paste0("count", i-1))
    tb_one <- tapho.load(one_column)
    tb <- tb %>% add_row(tb_one)
  }}

  return(tb)
}



library(dplyr)

# Define the function to read an Excel file with multiple sheets
read_excel_sheets <- function(file_path) {
  
  # Get all sheet names in the Excel file
  sheet_names <- readxl::excel_sheets(file_path)
  print(sheet_names)
  # Read each sheet, process it, and store in a list
  processed_list <- lapply(sheet_names, function(sname) {
    
    # Read the current sheet
    sheet_data <- readxl::read_excel(file_path, sheet = sname)
    
    # Process it using your tapho.load.db() function
    processed_sheet <- tapho.load.db(sheet_data)
    
    # Optionally, you might want to keep track of which sheet it came from
    # processed_sheet <- processed_sheet %>% 
    #   mutate(sheet_name = sname)
    
    processed_sheet
  })
  
  # Combine all processed sheets into one tibble
  combined_df <- bind_rows(processed_list)
  
  return(combined_df)
}
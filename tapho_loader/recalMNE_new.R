#source("NewSiteDataLoader_new.R")
#source("data_loader_new.R")
source("meta_new.R")
#source("select_data_new.R")

library(dplyr)

sections_order <- c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
"arms", "hands", "pelvic", "legs", "foot")

recal_mne <- function(current_mne, current_levels, target_levels, skeleton_hierarchy, bone_names, sections_order) {
  
  # Initialize an empty list to store new MNE results
  new_mne <- list()
  
  # Create the initial tibble with counts and elements
  current_mne_named <- tibble(count = current_mne, element = bone_names)
  
  # Add bone info to the tibble (assuming add_bone_info function exists)
  current_mne_named <- add_bone_info(current_mne_named)
  
  # Helper function to standardize names
  standardize_name <- function(name) {
    gsub(" ", ".", name)
  }
 
  # Loop over each part/section
  for (part_index in seq_along(sections_order)) {
    part <- sections_order[part_index]
    current_level <- current_levels[[part_index]]
    target_level <- target_levels[[part_index]]
    print(part_index)

    # If current level matches target level, no aggregation is needed; just copy the original counts
    if (current_level == target_level) {
      # Extract counts for this part and add directly to new_mne
      part_data <- current_mne_named %>% filter(Skeleton_Section == part)
      for (i in seq_len(nrow(part_data))) {
        element <- standardize_name(part_data$element[i])
        count <- part_data$count[i]
        new_mne[[paste(part, element, sep = "_")]] <- count
        #new_mne[[category]] <- count
      }
      next
    }
    

    # Create the appropriate mapping name for this part and level transition
    mapping_name <- paste0(part, "_", current_level, "_to_", target_level)
   
    # Check if the mapping exists in the hierarchy
    if (mapping_name %in% names(category_hierarchy)) {
      # Get the mapping
      mapping <- category_hierarchy[[mapping_name]]
      
      # Build a lookup from elements to target categories
      element_to_category <- unlist(lapply(names(mapping), function(category) {
        std_category <- standardize_name(category)
        elements <- mapping[[category]]
        std_elements <- sapply(elements, standardize_name)
        setNames(rep(std_category, length(std_elements)), std_elements)
      }))
  
      # Filter data for this part
      part_data <- current_mne_named %>% filter(Skeleton_Section == part)
      
      # Standardize element names
      part_data$std_element <- standardize_name(part_data$element)
      
      # Map elements to target categories
      part_data$target_category <- element_to_category[part_data$std_element]
      
      # Handle elements not in mapping (optional)
      # For now, exclude them
      part_data <- part_data %>% filter(!is.na(target_category))
    

      # Aggregate counts by target category
      aggregated_counts <- part_data %>%
        group_by(target_category) %>%
        summarize(total_count = sum(count, na.rm = TRUE))

      aggregated_counts$target_category <- factor( aggregated_counts$target_category, levels = names(mapping))
      aggregated_counts <- aggregated_counts %>%
      arrange(target_category)
      
      # Store the counts in new_mne
      for (i in seq_len(nrow(aggregated_counts))) {
        category <- aggregated_counts$target_category[i]
        total_count <- aggregated_counts$total_count[i]
        new_mne[[paste(part, category, sep = "_")]] <- total_count
        #new_mne[[category]] <- total_count
      }
    } else {
      # Mapping not found, handle accordingly
      warning(paste("Mapping not found for", mapping_name))
    }
  }
  
  # Convert new_mne to data frame
  new_mne_df <- tibble(element = names(new_mne), count = unlist(new_mne))
  
  # Return the new MNE results
  return(new_mne_df)
}




recal_mne_selected_data <- function(selected_data) {


  #calc min bone level
  min_bone_level <- get_min(selected_data)
print(min_bone_level)
  # Recalculate counts and elements for each row
  recalculated_list <- pmap(
    list(selected_data$count, selected_data$bone.level, selected_data$element),
    function(count_i, bone_level_i, element_i) {
      recal_mne(
        current_mne = count_i,
        current_levels = bone_level_i,
        target_levels = min_bone_level,
        skeleton_hierarchy = skeleton_hierarchy,
        bone_names = element_i,
        sections_order = sections_order
      )
    }
  )
#print(recalculated_list)
  data_recalculated <- selected_data %>%
    mutate(
      recalculated = recalculated_list,
      bone.level = list(min_bone_level),
      element = map(recalculated, "element"),
      count = map(recalculated, "count")
    ) %>%
    select(-recalculated)



  return(data_recalculated)

}

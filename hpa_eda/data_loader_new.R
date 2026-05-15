library(dplyr)
library(tidyr)
library(stringr)

add_bone_info <- function(input_df, reference_df = elements_data) {
  print(names(elements_data))


  # Step 1: Ensure correct input structure
  if (!all(c("element", "count") %in% names(input_df))) {
    stop("Input dataframe must have 'element' and 'count' columns")
  }
  
  # Step 2: Unnest and clean the element names
  input_df <- input_df %>%
    unnest(element) %>%    # Unnest 'element' if it is a list column
    mutate(element = as.character(element), # Ensure 'element' is a character column
           element = str_replace(element, "^[^.]+\\.", "")) # Clean 'element'
  
  # Step 3: Prepare reference dataframe

  reference_df <- reference_df %>%
    select(ElementsOne, max_Level, Skeleton_Section) %>%
    distinct()
  
  # Step 4: Match bone names and add information
  result_df <- input_df %>%
    left_join(reference_df, by = c("element" = "ElementsOne"))
  
  # Step 5: Handle unmatched bones
  unmatched <- is.na(result_df$max_Level) | is.na(result_df$Skeleton_Section)
  if (any(unmatched)) {
    warning(paste("Unmatched bones:",
                  paste(result_df$element[unmatched], collapse = ", ")))
  }
  
  # Step 6: Ensure consistent column order
  result_df <- result_df %>%
    select(element, count, max_Level, Skeleton_Section)
 
  return(result_df)
}



summarize_bone_sections <- function(input_df, reference_df = elements_data) {
  bone_info_df <- add_bone_info(input_df, reference_df )

  # Define the ordered list of skeleton sections
  sections_order <- c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
  "arms", "hands", "pelvic", "legs", "foot")

  # Summarize the minimum Level for each section
  summary_df <- bone_info_df %>%
    group_by(Skeleton_Section) %>%
    summarize(min_level = min(Level, na.rm = TRUE)) %>%
    ungroup()

  # Create a named vector with all sections, initializing with 1 (instead of NA)
  result_vector <- setNames(rep(1, length(sections_order)), sections_order)

  # Fill in the result_vector with available data
  for (section in sections_order) {
    if (section %in% summary_df$Skeleton_Section) {
      level <- summary_df$min_level[summary_df$Skeleton_Section == section]
      # If the level is NA, keep it as 1, otherwise use the actual level
      result_vector[section] <- if(is.na(level)) 1 else level
    }
    # If the section is not in summary_df, it will remain 1
  }

  # Return just the numeric vector without names
  return(unname(result_vector))
}


# Function to look up bones and return a dataframe with correct element-count mapping
get_bones_dataframe <- function(levels_vector) {

  # Ensure the levels_vector has exactly 10 elements
  if (length(levels_vector) != 11) {
    stop("The input vector must have exactly 11 elements.")
  }

  # Sections in the skeleton hierarchy
  sections <-  c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
  "arms", "hands", "pelvic", "legs", "foot")

  # Initialize a list to store the elements and counts
  elements <- c()
  counts <- c()

  # Loop through each section and level
  for (i in seq_along(sections)) {
    section <- sections[i]
    level <- paste0("Level ", levels_vector[i])

    # Extract the corresponding bones and counts at the specified level
    bones_at_level <- skeleton_hierarchy[[section]][[level]]

    # If the result is a list, unlist it; otherwise, just use it as is
    if (is.list(bones_at_level)) {
      bones_flat <- unlist(bones_at_level)
    } else {
      bones_flat <- bones_at_level
    }

    # Append names (elements) and counts to the vectors
    elements <- c(elements, names(bones_flat))
    counts <- c(counts, as.numeric(bones_flat))
  }

  # Convert the vectors to a dataframe
  df <- data.frame(
    element = elements,
    count = counts,
    stringsAsFactors = FALSE
  )

  return(df)
}


compareBoneLevels <- function(...) {
  # Collect all arguments into a list
  results_list <- list(...)

  # Ensure that at least two results are provided
  if (length(results_list) < 2) {
    stop("At least two results from convertMNEtoVariables must be provided.")
  }

  # Extract and flatten the bone_level vectors from each result
  bone_levels <- lapply(results_list, function(res) {
    if (!is.list(res) || !"bone.level" %in% names(res)) {
      stop("Each input must be a list output from convertMNEtoVariables containing a bone_level vector.")
    }
    # Apply unlist to flatten the bone.level list into a numeric vector
    return(unlist(res$bone.level, use.names = FALSE))
  })

  # Compute the element-wise minimum across all bone_level vectors
  min_bone_level <- Reduce(pmin, bone_levels)

  return(min_bone_level)
}



## recal MNE -> take min (c()) -> get one skeleton + sum up -> recalc
source("meta/meta_set.R")


add_bone_info_indet <- function(data) {
  data %>%
    left_join(
      elements_data %>% select(ElementsOne, Level, Skeleton_Section),   # elements_data_indet
      by = c("element" = "ElementsOne")   # adjust "element" to your col name if different
    ) %>%
    mutate(Level = as.integer(Level))
}


get_indet_siblings <- function(indet_elem, part, current_level, section_name, MNI) {
  # Step 1: find parent
  parent_elem <- get_indet_parent(indet_elem, part, current_level)
  if (length(parent_elem) == 0) stop("No parent found for: ", indet_elem)
  
  # Step 2: get all children of that parent at current level
  key_down <- paste0(part, "_", current_level, "_to_", current_level - 1)
  mapping  <- category_hierarchy[[key_down]]    #category_hierarchy_indet
  
  # find the entry in mapping whose key is the parent
  parent_key <- names(mapping)[gsub(" ", ".", names(mapping)) == parent_elem]
  if (length(parent_key) == 0) stop("Parent not found in mapping: ", parent_elem)
  
  children <- gsub(" ", ".", mapping[[parent_key]])
    if(indet_elem == "lumbar_indet" ){ #print(children)
    }
  # Step 3: drop _indet siblings (only keep real elements)
  drop_indet_if_possible(children)
}

get_indet_parent <- function(indet_elem, part, current_level) {
  key <- paste0(part, "_", current_level, "_to_", current_level - 1) #,"_indet")
  #print(key)
  if (!key %in% names(category_hierarchy)) stop("Missing mapping: ", key) # skeleton_hierarchy_indet
  
  element_map <- stack(category_hierarchy[[key]]) %>%                     # skeleton_hierarchy_indet
    setNames(c("child", "parent")) %>%
    mutate(across(everything(), \(x) gsub(" ", ".", x)))
  #print(element_map)
  element_map %>% filter(child == indet_elem) %>% pull(parent)
}

# ── Function 2: distribute _indet randomly among siblings ────────────────────
distribute_indet <- function(data, part, MNI) {
  data <- add_bone_info_indet(data)
  
  indet_rows <- data %>% filter(grepl("_indet$", element), count > 0)
  if (nrow(indet_rows) == 0) {
    return(data %>% filter(!grepl("_indet$", element)) %>% select(element, count))
  }
  
  result <- data
  
  for (i in seq_len(nrow(indet_rows))) {
    indet_elem    <- indet_rows$element[i]
    indet_count   <- indet_rows$count[i]
    current_level <- indet_rows$Level[i]           # ← straight from add_bone_info_indet
    section       <- indet_rows$Skeleton_Section[i]
    
    siblings <- get_indet_siblings(indet_elem, part, current_level, section, MNI)
    
    level_label  <- paste0("Level ", current_level)
    level_maxima <- skeleton_hierarchy[[section]][[level_label]]  #skeleton_hierarchy_indet
    
    existing_counts <- purrr::map_int(siblings, \(ch) {
      val <- result %>% filter(element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val)
    })
    
    max_counts <- purrr::map_int(siblings, \(ch) {
      base <- level_maxima[[ch]]
      if (is.null(base)) 0L else as.integer(base * MNI)
    })
    #print(siblings)
    #print(existing_counts)
    #print(max_counts)
    #print(indet_count)

    final_counts <- distribute_surplus(existing_counts, max_counts, indet_count)
    #print("final_counts[j]: ")
    for (j in seq_along(siblings)) {
      sib <- siblings[j]
      #cat("sib: ")
      #print(sib)
      #cat("final_counts[j]: ", as.integer(final_counts[j]))
      #cat("count: ", typeof(count))
      result <- result %>%
        mutate(count = if_else(element == sib, as.double(final_counts[j]), count))
    }
  }
  
  result %>%
    filter(!grepl("_indet$", element)) %>%
    select(element, count)
}

# ── Function 1: collapse _indet UP to parent ─────────────────────────────────
collapse_indet_up <- function(data, part) {
  data <- add_bone_info_indet(data)
  
  indet_rows <- data %>% filter(grepl("_indet$", element), count > 0)
  if (nrow(indet_rows) == 0) {
    return(data %>% filter(!grepl("_indet$", element)) %>% select(element, count))
  }
  
  result <- data

  for (i in seq_len(nrow(indet_rows))) {
    indet_elem    <- indet_rows$element[i]
    indet_count   <- indet_rows$count[i]
    current_level <- indet_rows$Level[i]           # ← straight from add_bone_info_indet
    section       <- indet_rows$Skeleton_Section[i]
    print(indet_elem)
    print(part)
    print(current_level)
    parent_elem <- get_indet_parent(indet_elem, part, current_level)
    if (length(parent_elem) == 0) stop("No parent found for: ", indet_elem)
    
    result <- result %>%
      mutate(count = case_when(
        element == parent_elem ~ count + indet_count,
        TRUE ~ count
      ))
  }
  
  result %>%
    filter(!grepl("_indet$", element)) %>%
    select(element, count)
}

# ── Helper: convert rib_fragments count to estimated whole bones ─────────────
convert_fragments <- function(count, fragments_per_bone = 3) {
  as.integer(floor(count / fragments_per_bone))
}

# ── Function 1: collapse rib_fragments UP to parent ──────────────────────────
collapse_fragments_up <- function(data, part, fragments_per_bone = 3) {
  data <- add_bone_info_indet(data)
 
  frag_rows <- data %>% filter(element == "rib_fragments", count > 0)
  if (nrow(frag_rows) == 0) {
    return(data %>% filter(element != "rib_fragments") %>% select(element, count))
  }
  
  result <- data
  #print("frag_rows: ")
   #print(frag_rows)
  for (i in seq_len(nrow(frag_rows))) {
    frag_count    <- frag_rows$count[i]
    current_level <- frag_rows$Level[i]
    section       <- frag_rows$Skeleton_Section[i]
    
    # Convert fragments → estimated whole bones
    whole_count <- convert_fragments(frag_count, fragments_per_bone)
        #print("whole_count: ")
    #print(whole_count)
    # Find parent (same logic as _indet: use category_hierarchy mapping)
    parent_elem <- get_indet_parent("rib_fragments", part, current_level)
        #print("parent_elem : ")
    #print(parent_elem )
    if (length(parent_elem) == 0) stop("No parent found for rib_fragments")
    
    result <- result %>%
      mutate(count = case_when(
        element == parent_elem ~ count + whole_count,
        TRUE ~ count
      ))
  }
  
  result %>%
    filter(element != "rib_fragments") %>%
    select(element, count)
}

# ── Function 2: distribute rib_fragments randomly among siblings ──────────────
distribute_fragments <- function(data, part, MNI, fragments_per_bone = 3) {
  data <- add_bone_info_indet(data)
  
  frag_rows <- data %>% filter(element == "rib_fragments", count > 0)
  if (nrow(frag_rows) == 0) {
    return(data %>% filter(element != "rib_fragments") %>% select(element, count))
  }
  
  result <- data
  
  for (i in seq_len(nrow(frag_rows))) {
    frag_count    <- frag_rows$count[i]
    current_level <- frag_rows$Level[i]
    section       <- frag_rows$Skeleton_Section[i]
    
    # Convert fragments → estimated whole bones
    whole_count <- convert_fragments(frag_count, fragments_per_bone)
    #print("whole_count: ")
    #print(whole_count)
    # Get siblings (excludes rib_fragments itself via drop_indet_if_possible equivalent)
    siblings <- get_indet_siblings("rib_fragments", part, current_level, section, MNI)
    #print("siblings: ")
    #print(siblings)
    level_label  <- paste0("Level ", current_level)
    level_maxima <- skeleton_hierarchy[[section]][[level_label]]  #skeleton_hierarchy_indet
    
    existing_counts <- purrr::map_int(siblings, \(ch) {
      val <- result %>% filter(element == ch) %>% pull(count)
      if (length(val) == 0) 0L else as.integer(val)
    })
    
    max_counts <- purrr::map_int(siblings, \(ch) {
      base <- level_maxima[[ch]]
      if (is.null(base)) 0L else as.integer(base * MNI)
    })
    
    final_counts <- distribute_surplus(existing_counts, max_counts, whole_count)
    #paste("final_counts[j]: ",final_counts)
    for (j in seq_along(siblings)) {
      sib <- siblings[j]
      #paste("sib: ", sib)
      #paste("final_counts[j]: ",final_counts[j])
      #paste("count: ", count)
      result <- result %>%
        mutate(count = if_else(element == sib, as.integer(final_counts[j]), count))
    }
  }
  
  result %>%
    filter(element != "rib_fragments") %>%
    select(element, count)
}

# ── Extended sibling helper that also excludes specific element names ─────────
get_siblings_excluding <- function(elem, part, current_level, section_name, MNI, 
  exclude = NULL) {
parent_elem <- get_indet_parent(elem, part, current_level)
if (length(parent_elem) == 0) stop("No parent found for: ", elem)

key_down <- paste0(part, "_", current_level, "_to_", current_level - 1) #,"_indet")
mapping  <- category_hierarchy[[key_down]]   # skeleton_hierarchy_indet

parent_key <- names(mapping)[gsub(" ", ".", names(mapping)) == parent_elem]
children   <- gsub(" ", ".", mapping[[parent_key]])

# Drop _indet AND any explicitly excluded elements (eg rib_fragments itself)
children <- drop_indet_if_possible(children)
children <- children[!children %in% c(elem, exclude)]
children
}

# ── Master indet resolution function ─────────────────────────────────────────

resolve_indeterminate <- function(
  new_site_raw, MNI,
  indet              = TRUE,
  indet_method       = c("up", "dist", "delet"),
  fragments_per_bone = 3
) {
indet_method <- match.arg(indet_method)

# ── Early exit ──────────────────────────────────────────────────────────────
if (indet) return(new_site_raw)

sections_order <- c(
  "cranium", "teeth", "hyoid", "vertebrae", "thorax", "shoulder",
  "arms", "hands", "pelvic", "legs", "foot"
)

meta_set(!indet)  #  < ─────reset ENV──────────
  
# ── Enrich with bone info to get Skeleton_Section per element ───────────────
data_enriched <- add_bone_info_indet(new_site_raw)

  
  

# ── Final cleanup applied to ALL methods ────────────────────────────────────
clean_indet_fragments <- function(df) {
  df %>%
    filter(
      !grepl("_indet$",     element),
      !grepl("_fragments$", element)
    )
}

# ── Delet: drop immediately, no processing needed ───────────────────────────
if (indet_method == "delet") {
  return(clean_indet_fragments(new_site_raw))
}

# ── For "up" and "dist": process section by section ─────────────────────────
result <- new_site_raw

for (section in sections_order) {
   #print(data_enriched)
  #print("section: ")
  #print(section)
  section_elements <- data_enriched %>%
    filter(Skeleton_Section == section) %>%
    pull(element)
  
  section_data <- result %>% filter(element %in% section_elements)
  
  has_indet     <- any(grepl("_indet$",     section_data$element) & section_data$count > 0)
  has_fragments <- any(grepl("_fragments$", section_data$element) & section_data$count > 0)
  
  if (!has_indet && !has_fragments) next
  
  if (indet_method == "up") {
    
    if (has_indet) {
      section_data <- collapse_indet_up(section_data, part = section)
    }
    if (has_fragments && section == "thorax") {
      section_data <- collapse_fragments_up(section_data, part = section,
                                            fragments_per_bone = fragments_per_bone)
    }
    
  } else if (indet_method == "dist") {
    
    if (has_indet) {
      #print(section)
      #print(MNI)
      section_data <- distribute_indet(section_data, part = section, MNI = MNI)
    }
    if (has_fragments && section == "thorax") {
      section_data <- distribute_fragments(section_data, part = section, MNI = MNI,
                                           fragments_per_bone = fragments_per_bone)
    }
  }
  
  meta_set(!indet)  #  < ─────reset ENV──────────

  result <- result %>%
    filter(!element %in% section_elements) %>%
    bind_rows(section_data) %>%
    arrange(match(element, data_enriched$element))
}

# ── Final safety net: drop any _indet / _fragments rows with count == 0 ─────
# Catches sections skipped because count == 0 but row still present
result %>% clean_indet_fragments()
}
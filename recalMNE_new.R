#source("NewSiteDataLoader_new.R")
source("data_loader_new.R")
#source("meta_new_indet.R")
source("collapse_down_step.R")
#source("select_data_new.R")

library(dplyr)


sections_order <- c("cranium","teeth","hyoid", "vertebrae", "thorax", "shoulder",
"arms", "hands", "pelvic", "legs", "foot")

is_special <- function(part) part %in% c("pelvic", "cranium")

recal_mne <- function(current_mne, current_levels, target_levels, skeleton_hierarchy,
  bone_names, sections_order, MNI) {

working_base <- tibble(count = current_mne, element = bone_names) %>%
add_bone_info() %>%
mutate(std_element = gsub(" ", ".", element))

purrr::imap(sections_order, \(part, i) {
data <- working_base %>%
filter(Skeleton_Section == part) %>%
select(std_element, count)


cur <- current_levels[i]
tgt <- target_levels[i]

if (cur > tgt) {
  if (is_special(part)){
    data <- collapse_up_step_special(data, part,
      current_level = cur,
      target_level  = tgt)}
  else{
data <- collapse_up_step(data, part,
        current_level = cur,
        target_level  = tgt)}


} else if (cur < tgt) {
  if (is_special(part)){
    data <- collapse_down_step_direct_special(data, part,
              current_level = cur,
              target_level  = tgt,
              section_name  = part,
              MNI           = MNI)}
  else{data <- collapse_down_step(data, part,
    current_level = cur,
    target_level  = tgt)}

}
  
 # after collapse_up/down (or even if no collapse happened)
level_label <- paste0("Level ", tgt)
expected <- names(skeleton_hierarchy[[part]][[level_label]])
expected <- gsub(" ", ".", expected) 
if (part == "pelvic") {
  #print("#######################  level_label ##################################")
  #print(part)
  #print(level_label)
  #cat("pelvic cur/tgt:", cur, tgt, "\n")
  #print(setdiff(expected, data$std_element))
  #print("#######################  level_label  ##################################")
}


data
}) %>%
bind_rows() %>%
rename(element = std_element)
}




recal_mne_selected_data <- function(selected_data, 
  method = c("min", "filter", "direct"), 
  target_level = NULL) {

method <- match.arg(method)

# --- Determine target level based on method ---
if (method == "min") {
# Original behaviour: compute minimum across all rows
min_bone_level <- get_min(selected_data)
target_levels_vec <- min_bone_level

} else if (method == "filter") {
# Given a target level, filter out rows where any element-wise 
# bone level is NOT <= target_level, then recalculate on survivors
if (is.null(target_level)) stop("'target_level' must be provided for method = 'filter'")

# Keep only rows where every bone level value <= corresponding target
selected_data <- selected_data %>%
  dplyr::filter(purrr::map_lgl(bone.level, function(bl) {
    bl_vec <- as.numeric(unlist(bl))
    tgt <- target_level[seq_along(bl_vec)]
    all(is.na(bl_vec) | (bl_vec >= tgt))
  }))

if (nrow(selected_data) == 0) {
warning("No rows remain after filtering. Consider a less strict target_level.")
return(selected_data)
}

target_levels_vec <- target_level

} else if (method == "direct") {
# Use the provided target level directly, no filtering, no min calculation
if (is.null(target_level)) stop("'target_level' must be provided for method = 'direct'")
target_levels_vec <- target_level
}

# --- Recalculate MNE for each row using the resolved target level ---
recalculated_list <- pmap(
list(selected_data$count, selected_data$bone.level, selected_data$element, selected_data$mni),
function(count_i, bone_level_i, element_i, mni) {
recal_mne(
current_mne     = count_i,
current_levels  = bone_level_i,
target_levels   = target_levels_vec,
skeleton_hierarchy = skeleton_hierarchy,
bone_names      = element_i,
sections_order  = sections_order,
MNI             = mni
)
}
)

data_recalculated <- selected_data %>%
mutate(
recalculated = recalculated_list,
bone.level   = list(target_levels_vec),
element      = map(recalculated, "element"),
count        = map(recalculated, "count")
) %>%
select(-recalculated)

return(data_recalculated)
}
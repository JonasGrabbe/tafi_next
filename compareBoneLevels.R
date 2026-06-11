#' Compare and Compute Element-wise Minimum of Bone Level Vectors
#'
#' This function takes multiple outputs from the \code{\link{convertMNEtoVariables}} function,
#' extracts the \code{$bone_level} vectors, and computes a new vector containing the element-wise
#' minimum values across all provided bone level vectors.
#'
#' @param ... Outputs from \code{\link{convertMNEtoVariables}} function. At least two outputs
#' should be provided, each containing a \code{$bone_level} vector.
#'
#' @return A numeric vector representing the element-wise minimum values of the provided
#' \code{$bone_level} vectors.
#'
#' @examples
#' # Assuming result1 and result2 are outputs from convertMNEtoVariables function
#' result1 <- convertMNEtoVariables(MNE1, skeleton1)
#' result2 <- convertMNEtoVariables(MNE2, skeleton2)
#'
#' # Compute the element-wise minimum of the bone level vectors
#' min_bone_level <- compareBoneLevels(result1, result2)
#' print(min_bone_level)
#'
#' @export
compareBoneLevels <- function(...) {
  # Collect all arguments into a list
  results_list <- list(...)

  # Ensure that at least two results are provided
  if (length(results_list) < 2) {
    stop("At least two results from convertMNEtoVariables must be provided.")
  }

  # Extract the bone_level vectors from the results
  bone_levels <- lapply(results_list, function(res) {
    if (!is.list(res) || !"bone_level" %in% names(res)) {
      stop("Each input must be a list output from convertMNEtoVariables containing a bone_level vector.")
    }
    return(res$bone_level)
  })

  # Compute the element-wise minimum across all bone_level vectors
  min_bone_level <- Reduce(pmin, bone_levels)

  return(min_bone_level)
}


# Define a modified version of compareBoneLevels to handle dataframe column directly
compareBoneLevelsDF <- function(data) {
  # Extract the bone_level vectors from the dataframe
  bone_levels <- data$bone.level  # This extracts all bone.level vectors
  
  # Ensure there are at least two bone_level vectors
  if (length(bone_levels) < 2) {
    stop("At least two bone.level vectors must be provided.") 
  }   # !! Adjust alert message bone.level -> data? !!
  
  # Compute the element-wise minimum across all bone_level vectors
  min_bone_level <- Reduce(pmin, bone_levels)
  
  return(min_bone_level)
}
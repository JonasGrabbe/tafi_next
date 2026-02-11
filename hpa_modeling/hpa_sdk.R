generate_synthetic_hpa_data <- function(n_samples, features, classes) {
  n_classes <- length(levels(classes))
  x_data <- matrix(sample(0:100, n_samples * length(features), replace = TRUE),
                   nrow = n_samples, ncol = length(features))
  colnames(x_data) <- features
  y_labels <- sample(classes, n_samples, replace = TRUE) 
  synthetic_data <- as_tibble(x_data) %>%
    mutate(!!TARGET_VARIABLE := y_labels, .before = 1) #%>%  
  colnames(synthetic_data) <- make.names(colnames(synthetic_data))
  return(synthetic_data)
}


compute_classification_report_from_confmats <- function(conf_matrices) {
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # Helper function to extract metrics per class from a confusion matrix
  extract_metrics <- function(cm, id, id2) {
    long_cm <- as.data.frame(cm$table) %>%
      rename(Prediction = Prediction, Truth = Truth, count = Freq)
  
    class_labels <- unique(long_cm$Truth)
  
    map_dfr(class_labels, function(cls) {
      tp <- long_cm %>% filter(Prediction == cls, Truth == cls) %>% pull(count)
      fp <- long_cm %>% filter(Prediction == cls, Truth != cls) %>% summarise(fp = sum(count)) %>% pull(fp)
      fn <- long_cm %>% filter(Prediction != cls, Truth == cls) %>% summarise(fn = sum(count)) %>% pull(fn)
      support <- long_cm %>% filter(Truth == cls) %>% summarise(support = sum(count)) %>% pull(support)
  
      precision <- if ((tp + fp) == 0) NA else tp / (tp + fp)
      recall    <- if ((tp + fn) == 0) NA else tp / (tp + fn)
      f1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) NA else {
        2 * precision * recall / (precision + recall)
      }
  
      tibble(
        id = id,
        id2 = id2,
        Class = cls,
        precision = precision,
        recall = recall,
        f1 = f1,
        support = support
      )
    })
  }
  
  
  
  # Apply across all (id, id2) folds
  all_metrics <- pmap_dfr(
    list(conf_matrices$conf_matrix, conf_matrices$id, conf_matrices$id2),
    extract_metrics
  )
  
  # Aggregate mean and sd per class
  per_class_summary <- all_metrics %>%
    group_by(Class) %>%
    summarise(
      precision_mean = mean(precision, na.rm = TRUE),
      precision_sd   = sd(precision, na.rm = TRUE),
      recall_mean    = mean(recall, na.rm = TRUE),
      recall_sd      = sd(recall, na.rm = TRUE),
      f1_mean        = mean(f1, na.rm = TRUE),
      f1_sd          = sd(f1, na.rm = TRUE),
      support        = mean(support, na.rm = TRUE),
      support_sd     = sd(support, na.rm = TRUE),
      .groups = "drop"
    )

  # Macro average (simple mean of classes)
  macro_avg <- per_class_summary %>%
    summarise(
      precision_mean = mean(precision_mean, na.rm = TRUE),
      precision_sd   = mean(precision_sd, na.rm = TRUE),
      recall_mean    = mean(recall_mean, na.rm = TRUE),
      recall_sd      = mean(recall_sd, na.rm = TRUE),
      f1_mean        = mean(f1_mean, na.rm = TRUE),
      f1_sd          = mean(f1_sd, na.rm = TRUE),
      support        = sum(support, na.rm = TRUE),
      support_sd     = NA_real_  # Not applicable
    ) %>%
    mutate(Class = "macro_avg") %>%
    select(Class, everything())

  # Weighted average
  total_support <- sum(per_class_summary$support, na.rm = TRUE)

  weighted_avg <- per_class_summary %>%
    summarise(
      precision_mean = weighted.mean(precision_mean, support, na.rm = TRUE),
      precision_sd   = weighted.mean(precision_sd, support, na.rm = TRUE),
      recall_mean    = weighted.mean(recall_mean, support, na.rm = TRUE),
      recall_sd      = weighted.mean(recall_sd, support, na.rm = TRUE),
      f1_mean        = weighted.mean(f1_mean, support, na.rm = TRUE),
      f1_sd          = weighted.mean(f1_sd, support, na.rm = TRUE),
      support        = total_support,
      support_sd     = NA_real_
    ) %>%
    mutate(Class = "weighted_avg") %>%
    select(Class, everything())

  # Combine full report
  full_report <- bind_rows(per_class_summary, macro_avg, weighted_avg)

  # Pretty formatting
  pretty_report <- full_report %>%
    mutate(
      Precision = sprintf("%.3f (%.3f)", precision_mean, precision_sd),
      Recall    = sprintf("%.3f (%.3f)", recall_mean, recall_sd),
      F1        = sprintf("%.3f (%.3f)", f1_mean, f1_sd),
      Support   = ifelse(
        is.na(support_sd),
        sprintf("%.1f", support),
        sprintf("%.1f (%.1f)", support, support_sd)
      )
    ) %>%
    select(Class, Precision, Recall, F1, Support)

  list(
    raw = full_report,
    pretty = pretty_report
  )
}

#### taphor ####

convertMNEtoVariables <- function(MNE, skeleton, individuals = NULL) {

  if (length(MNE) != nrow(skeleton)) {
    stop("Length of MNE vector must match the number of skeletal elements")
  }

  # Convert MNE to a data frame if it's not already
  if (!is.data.frame(MNE)) {
    MNE <- as.data.frame(MNE)
  }

  Elements <- skeleton$element

  # Calculate MNI and Individuals if not provided
  if (is.null(individuals)) {
    MNI <- ceiling(MNE$MNE / skeleton$count)
    individuals <- max(MNI)
  }

  # Compute Abundance for each part
  abundance <- MNE$MNE / (skeleton$count * individuals)

  # Compute Mean Abundance (MA)
  ma <- sum(MNE$MNE, na.rm = TRUE) / (sum(skeleton$count, na.rm = TRUE) * individuals)

  # Compute Mean Part Abundance (MPA)
  mpa <- sum(abundance, na.rm = TRUE) / length(abundance)

  # Additional metrics calculations
  MAU <- abundance
  MNAU <- abundance
  PMAU <- MAU / max(MNE$MNE / skeleton$count)
  PMNAU <- MNAU / individuals
  ReMNAU <- MNAU / sum(MNAU, na.rm = TRUE)

  # Calculate Accumulative ReMNAU
  AcReMNAU <- cumsum(ReMNAU)

  # Bone level vector calculation using summarize_bone_sections
  bone_info_df <- data.frame(element = skeleton$element, count = MNE$MNE)
  #bone_level <- summarize_bone_sections(bone_info_df)

  # Return a list of all calculated metrics
  return(list(
    MNE = as.vector(MNE$MNE),
    Elements = Elements,
    abundance = abundance,
    MNI = individuals,
    MAU = MAU,
    MNAU = MNAU,
    PercentMAU = PMAU,
    PercentMNAU = PMNAU,
    ReMNAU = ReMNAU,
    AcReMNAU = AcReMNAU,
    MeanAbundance = ma,
    MeanPartAbundance = mpa#,
    #bone_level = bone_level
  ))
}


# --- Revised convertMNEtoVariables_v2 Function ---
convertMNEtoVariables_v15 <- function(MNE_input_vector,
  skeleton_ref_df, # Changed name for clarity
  individuals_site_mni = NULL) { # Changed name

# Minimal argument validation (as requested, but some is still good)
if (length(MNE_input_vector) != nrow(skeleton_ref_df)) {
stop("MNE vector length must match skeleton_ref_df rows.")
}
MNE_values <- as.numeric(MNE_input_vector)
MNE_values[is.na(MNE_values)] <- 0 # Consistent with prior NA filling

skeletal_elements_names <- as.character(skeleton_ref_df$element)
skeletal_element_counts <- as.numeric(skeleton_ref_df$count)
# Ensure counts are positive for division, treat 0 or negative counts as effectively 1 (or a tiny number)
# to avoid Inf/NaN issues that are hard to recover from if the data is truly 0.
# A count of 0 for an element means it shouldn't exist, so MNAU should be 0 or Inf if MNE > 0.
# Let's make safe_counts very small if count is 0, so MNE/safe_count becomes large (Inf-like) if MNE>0
# and then we set non-finite to 0. Or, more simply, if count is 0, MNAU is 0 (unless MNE>0 which is an issue).
# A practical approach: if a skeletal element count is 0, MNAU is undefined or should be 0 if MNE is 0.
# For calculations, use a small positive number to avoid division by zero, then manage Infs.
safe_skeletal_counts <- ifelse(skeletal_element_counts <= 0, 1e-9, skeletal_element_counts)


# MNI for the assemblage
MNI_used <- individuals_site_mni
valid_individuals_provided <- FALSE
if (!is.null(MNI_used)) {
if(is.numeric(MNI_used) && length(MNI_used) == 1 && MNI_used > 0 && !is.na(MNI_used)) {
valid_individuals_provided <- TRUE
} else {
MNI_used <- NULL # Force recalculation
}
}
if (is.null(MNI_used)) {
mni_calc_values <- ceiling(MNE_values / safe_skeletal_counts)
MNI_used <- max(mni_calc_values[is.finite(mni_calc_values)], 0, na.rm = TRUE)
}
MNI_used <- ifelse(MNI_used <= 0, 1, MNI_used) # Ensure MNI is at least 1

# MNAU = MNE / count_in_skeleton
MNAU_calc <- MNE_values / safe_skeletal_counts
MNAU_calc[!is.finite(MNAU_calc)] <- 0 # MNE > 0 and count = 0 would be Inf, set to 0 (or handle as error)

# MAU = MNE / count_in_skeleton (same as MNAU per discussion)
MAU_calc <- MNAU_calc

# %MAU = (MAU_i / max(MAU)) * 100
max_MAU_val <- max(MAU_calc, na.rm = TRUE)
PercentMAU_calc <- if (max_MAU_val > 0) (MAU_calc / max_MAU_val) * 100 else 0
PercentMAU_calc[!is.finite(PercentMAU_calc)] <- 0

# %MNAU (same as %MAU here)
PercentMNAU_calc <- PercentMAU_calc

# ReMNAU (observed) = MNAU_i / sum(MNAU)
sum_MNAU_val <- sum(MNAU_calc, na.rm = TRUE)
ReMNAU_observed_calc <- if (sum_MNAU_val > 0) MNAU_calc / sum_MNAU_val else 0
ReMNAU_observed_calc[!is.finite(ReMNAU_observed_calc)] <- 0

# AcReMNAU = ReMNAU_observed / ReMNAU_expected
# ReMNAU_expected based on skeletal counts proportions
sum_total_skeletal_parts_in_skeleton <- sum(skeletal_element_counts, na.rm = TRUE)
ReMNAU_expected_calc <- if (sum_total_skeletal_parts_in_skeleton > 0) {
skeletal_element_counts / sum_total_skeletal_parts_in_skeleton
} else {
rep(0, length(skeletal_element_counts)) # Vector of zeros
}
ReMNAU_expected_calc[!is.finite(ReMNAU_expected_calc)] <- 0

AcReMNAU_calc <- ReMNAU_observed_calc / ifelse(ReMNAU_expected_calc == 0, 1e-9, ReMNAU_expected_calc) # Avoid div by zero for expected=0
AcReMNAU_calc[!is.finite(AcReMNAU_calc)] <- 0 # if observed was 0, or if expected was 0 and observed >0 -> Inf like

# Cumulative ReMNAU (observed)
CumulativeReMNAU_calc <- cumsum(ReMNAU_observed_calc)

# Proportional Representation (based on MNI) = MNAU / MNI_used
ProportionalRep_MNIbased_calc <- if (MNI_used > 0) MNAU_calc / MNI_used else 0
ProportionalRep_MNIbased_calc[!is.finite(ProportionalRep_MNIbased_calc)] <- 0

# --- Assemblage-wide summary stats ---
Mean_MNAU_assemblage_calc <- mean(MNAU_calc, na.rm = TRUE)
Mean_MNAU_assemblage_calc <- ifelse(is.nan(Mean_MNAU_assemblage_calc), 0, Mean_MNAU_assemblage_calc)

OverallProportionalRep_assemblage_calc <- if (sum_total_skeletal_parts_in_skeleton > 0 && MNI_used > 0) {
sum(MNE_values, na.rm = TRUE) / (sum_total_skeletal_parts_in_skeleton * MNI_used)
} else {0}
OverallProportionalRep_assemblage_calc <- ifelse(is.nan(OverallProportionalRep_assemblage_calc),0,OverallProportionalRep_assemblage_calc)

output <- tibble::tibble(
element = skeletal_elements_names,
MNE = MNE_values,
MNAU = MNAU_calc,
MAU = MAU_calc,
PercentMAU = PercentMAU_calc,
PercentMNAU = PercentMNAU_calc,
ReMNAU = ReMNAU_observed_calc,
AcReMNAU = AcReMNAU_calc,
CumulativeReMNAU = CumulativeReMNAU_calc,
ProportionalRep_MNIbased = ProportionalRep_MNIbased_calc,
# Assemblage-wide values (repeated for each element)
MNI_used_assemblage = MNI_used,
Mean_MNAU_assemblage = Mean_MNAU_assemblage_calc,
OverallProportionalRep_assemblage = OverallProportionalRep_assemblage_calc
)
return(output)
}

convertMNEtoVariables_v2 <- function(MNE,  skeleton, individuals = NULL) {

# Minimal argument validation
if (length(MNE) != nrow(skeleton)) { 
stop("MNE vector length must match skeleton data frame rows.")
}
MNE_values <- as.numeric(MNE) 
MNE_values[is.na(MNE_values)] <- 0

skeletal_elements_names <- as.character(skeleton$element) 
skeletal_element_counts <- as.numeric(skeleton$count)   
safe_skeletal_counts <- ifelse(skeletal_element_counts <= 0, 1e-9, skeletal_element_counts) 

# MNI for the assemblage
MNI_used <- individuals 
valid_individuals_provided <- FALSE
if (!is.null(MNI_used)) {
  # given MNI value should be an non-empty numeric scalar value at least 1!
if(is.numeric(MNI_used) && length(MNI_used) == 1 && MNI_used > 0 && !is.na(MNI_used)) {
valid_individuals_provided <- TRUE
} else {
MNI_used <- NULL
}
}
if (is.null(MNI_used)) {
mni_calc_values <- ceiling(MNE_values / safe_skeletal_counts)
MNI_used <- max(mni_calc_values[is.finite(mni_calc_values)], 0, na.rm = TRUE)
}
MNI_used <- ifelse(MNI_used <= 0, 1, MNI_used)

# MNAU = MNE / count_in_skeleton
MNAU_calc <- MNE_values / safe_skeletal_counts
MNAU_calc[!is.finite(MNAU_calc)] <- 0

MAU_calc <- MNAU_calc # MAU = MNAU

# %MAU
max_MAU_val <- max(MAU_calc, na.rm = TRUE)
PercentMAU_calc <- if (max_MAU_val > 0) {
  (MAU_calc / max_MAU_val) * 100
} else {
  rep(0, length(MAU_calc))
}
PercentMAU_calc[!is.finite(PercentMAU_calc) | is.na(PercentMAU_calc)] <- 0
PercentMNAU_calc <- PercentMAU_calc

# ReMNAU (observed)
sum_MNAU_val <- sum(MNAU_calc, na.rm = TRUE)
ReMNAU_observed_calc <- if (sum_MNAU_val > 0) {
  MNAU_calc / sum_MNAU_val
} else {
  rep(0, length(MNAU_calc))
}
ReMNAU_observed_calc[!is.finite(ReMNAU_observed_calc) | is.na(ReMNAU_observed_calc)] <- 0

# AcReMNAU = ReMNAU_observed / ReMNAU_expected
sum_total_skeletal_parts_in_skeleton <- sum(skeletal_element_counts, na.rm = TRUE)
ReMNAU_expected_calc <- if (sum_total_skeletal_parts_in_skeleton > 0) {
skeletal_element_counts / sum_total_skeletal_parts_in_skeleton
} else { rep(0, length(skeletal_element_counts)) }
ReMNAU_expected_calc[!is.finite(ReMNAU_expected_calc)] <- 0

AcReMNAU_calc <- ReMNAU_observed_calc / ifelse(ReMNAU_expected_calc == 0, 1e-9, ReMNAU_expected_calc)
AcReMNAU_calc[!is.finite(AcReMNAU_calc)] <- 0

CumulativeReMNAU_calc <- cumsum(ReMNAU_observed_calc)

ProportionalRep_MNIbased_calc <- if (MNI_used > 0) MNAU_calc / MNI_used else 0
ProportionalRep_MNIbased_calc[!is.finite(ProportionalRep_MNIbased_calc)] <- 0

Mean_MNAU_assemblage_calc <- mean(MNAU_calc, na.rm = TRUE)
Mean_MNAU_assemblage_calc <- ifelse(is.nan(Mean_MNAU_assemblage_calc), 0, Mean_MNAU_assemblage_calc)

OverallProportionalRep_assemblage_calc <- if (sum_total_skeletal_parts_in_skeleton > 0 && MNI_used > 0) {
sum(MNE_values, na.rm = TRUE) / (sum_total_skeletal_parts_in_skeleton * MNI_used)
} else {0}
OverallProportionalRep_assemblage_calc <- ifelse(is.nan(OverallProportionalRep_assemblage_calc),0,OverallProportionalRep_assemblage_calc)

output <- tibble::tibble(
element = skeletal_elements_names,
MNE = MNE_values, MNAU = MNAU_calc, MAU = MAU_calc,
PercentMAU = PercentMAU_calc, PercentMNAU = PercentMNAU_calc,
ReMNAU = ReMNAU_observed_calc, AcReMNAU = AcReMNAU_calc,
CumulativeReMNAU = CumulativeReMNAU_calc,
ProportionalRep_MNIbased = ProportionalRep_MNIbased_calc,
MNI_used_assemblage = MNI_used,
Mean_MNAU_assemblage = Mean_MNAU_assemblage_calc,
OverallProportionalRep_assemblage = OverallProportionalRep_assemblage_calc
)
return(output)
}

convertMNEtoVariables_v3 <- function(MNE, skeleton, individuals = NULL) {
  stopifnot(length(MNE) == nrow(skeleton))

  # Inputs
  MNE  <- as.numeric(MNE);  MNE[is.na(MNE)] <- 0
  c_i  <- as.numeric(skeleton$count)
  elnm <- as.character(skeleton$element)

  # MNI (use provided if valid; else compute)
  if (!is.null(individuals) && is.finite(individuals) && individuals > 0) {
    MNI_used <- as.numeric(individuals)
  } else {
    MNI_used <- max(ceiling(MNE / c_i), na.rm = TRUE)
    if (!is.finite(MNI_used) || is.na(MNI_used) || MNI_used < 1) MNI_used <- 1
  }

  # Core quantities
  MNAU <- MNE / c_i
  MAU  <- MNAU

  maxMAU <- suppressWarnings(max(MAU, na.rm = TRUE))
  #max-normalization (Binford)
  PercentMAU  <- if (is.finite(maxMAU) && maxMAU > 0) 100 * MAU / maxMAU else rep(0, length(MAU)) #
  # %MNAU MNAU / MNI
  PercentMNAU <- 100 * MNAU / MNI_used   # alternative definition - 


  totalMNAU <- sum(MNAU, na.rm = TRUE)
  ReMNAU    <- if (totalMNAU > 0) MNAU / totalMNAU else rep(0, length(MNAU))

  expected  <- c_i / sum(c_i, na.rm = TRUE)
  AcReMNAU  <- ReMNAU / expected

  CumulativeReMNAU <- cumsum(ReMNAU)  # note: depends on current row order

  ProportionalRep_MNIbased <- MNAU / MNI_used

  Mean_MNAU_assemblage <- mean(MNAU, na.rm = TRUE)

  Mean_part_abundance <- mean(ProportionalRep_MNIbased , na.rm = TRUE)

  OverallProportionalRep_assemblage <-
    sum(MNE, na.rm = TRUE) / (sum(c_i, na.rm = TRUE) * MNI_used)

  tibble::tibble(
    element = elnm,
    MNE = MNE,
    MNAU = MNAU,
    MAU = MAU,
    PercentMAU = PercentMAU,
    PercentMNAU = PercentMNAU,
    ReMNAU = ReMNAU,
    AcReMNAU = AcReMNAU,
    CumulativeReMNAU = CumulativeReMNAU,
    ProportionalRep_MNIbased = ProportionalRep_MNIbased,
    MNI_used_assemblage = MNI_used,
    Mean_MNAU_assemblage = Mean_MNAU_assemblage,
    Mean_part_abundance = Mean_part_abundance, 
    OverallProportionalRep_assemblage = OverallProportionalRep_assemblage
  )
}


plot_pca_3d <- function(input_tibble, title_text) {
  
  # Select only numeric predictor columns (heuristic: all numeric except 'type' if it were numeric, but it's factor)
  # Exclude site_name and type columns
  data_for_pca <- input_tibble %>%
    select(-all_of(c(site_name_col, type_col))) %>%
    # Ensure all selected columns are numeric; PCA will fail otherwise
    # This step might not be necessary if data is already clean numeric
    mutate(across(everything(), as.numeric)) %>% 
    # Remove columns that are all NA or have zero variance after numeric conversion
    select_if(~!all(is.na(.))) %>%
    select_if(~var(., na.rm = TRUE) > 0)

  if (ncol(data_for_pca) < 3) {
    warning(paste("Not enough columns for 3D PCA in", title_text, "(need at least 3). Skipping plot."))
    return(NULL)
  }
  
  # Perform PCA
  # scale. = TRUE is important for variables on different scales
  # center = TRUE is default
  pca_results <- prcomp(data_for_pca, scale. = TRUE, center = TRUE)
  
  # Create a data frame for plotting
  pca_plot_data <- tibble(
    PC1 = pca_results$x[,1],
    PC2 = pca_results$x[,2],
    PC3 = pca_results$x[,3],
    type = input_tibble[[type_col]] # Get the original type column for coloring
  )
  
  # Create 3D scatter plot
  fig <- plot_ly(pca_plot_data, x = ~PC1, y = ~PC2, z = ~PC3, 
                 color = ~type, colors = "viridis", # Or any other color palette
                 marker = list(size = 5)) %>%
    add_markers() %>%
    plotly::layout(title = title_text,
           scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')))
  
  return(fig)
}

plot_pca_2d_and_save <- function(input_tibble, title_text, filename_prefix) {
  
  # Select only numeric predictor columns
  data_for_pca <- input_tibble %>%
    select(-all_of(c(site_name_col, type_col))) %>%
    mutate(across(everything(), as.numeric)) %>%
    select_if(~!all(is.na(.))) %>%
    select_if(~var(., na.rm = TRUE) > 0) # Remove zero-variance columns

  if (ncol(data_for_pca) < 2) { # Need at least 2 columns for 2D PCA
    warning(paste("Not enough columns for 2D PCA in", title_text, "(need at least 2). Skipping plot."))
    return(NULL)
  }
  
  # Perform PCA
  pca_results <- prcomp(data_for_pca, scale. = TRUE, center = TRUE)
  
  # Get variance explained
  pca_summary <- summary(pca_results)
  variance_explained <- pca_summary$importance[2, 1:2] * 100 # PC1 and PC2 variance
  
  # Create a data frame for plotting
  pca_plot_data <- tibble(
    PC1 = pca_results$x[,1],
    PC2 = pca_results$x[,2],
    type = factor(input_tibble[[type_col]]) # Ensure 'type' is a factor for coloring
  )
  
  # Create 2D scatter plot with ggplot2
  pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = type, shape = type)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = title_text,
         x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
         y = paste0("PC2 (", round(variance_explained[2], 1), "%)"),
         color = "Type",
         shape = "Type") +
    theme_bw() +
    scale_color_viridis_d() # Using viridis discrete color scale

  # Print the plot to the R graphics device
  print(pca_plot)
  
  # Save the plot
  output_filename <- paste0(filename_prefix, "_2D_PCA.png")
  ggsave(output_filename, plot = pca_plot, width = 8, height = 6, dpi = 300)
  cat(paste("\nPCA plot saved as:", output_filename, "\n"))
  
  return(pca_plot) # Return the ggplot object
}


cosine_similarity_robust <- function(vec1, vec2) {
  vec1_num <- as.numeric(vec1)
  vec2_num <- as.numeric(vec2)
  
  norm1 <- sqrt(sum(vec1_num^2))
  norm2 <- sqrt(sum(vec2_num^2))
  
  # If either vector has zero norm (e.g., all zeros), similarity is 0
  # (or undefined/NaN - returning 0 assumes no similarity)
  if (norm1 == 0 || norm2 == 0) {
    return(0) 
  }
  sum(vec1_num * vec2_num) / (norm1 * norm2)
}

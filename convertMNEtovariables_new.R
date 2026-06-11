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
    element_2 = elnm,
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
# ==============================================================================
# Supervised ML: 7 feature spaces x 4 algorithms x 3 classification tasks
# Direct tapho_db_sum analysis + bilateral indeterminate redistribution
# Repeated grouped CV + unknown-site projection + cross-conformal summaries
#
# Tasks:
#   1) 4class:
#        Cannibalism / Carnivore / Funerary Context / Geological transport
#   2) cannibalism_vs_fc:
#        Cannibalism / Funerary Context
#   3) fc_vs_rest:
#        Funerary Context / Rest
#        Rest = Cannibalism + Carnivore + Geological transport
#
# IMPORTANT:
# - Input is final/tapho_db_sum.csv. 
# - Unknown/Unkown and Sepulchral/Cuevas Sepulcrales never enter fitting,
#   preprocessing, balancing, CV performance estimation, or conformal calibration.
# - All data-dependent preprocessing is fitted inside each training fold only.
# - Training classes are balanced by random oversampling WITHIN the training fold.
# - No hyperparameter search is performed. Settings are fixed a priori.
# - Repeated grouped CV: 5 repeats x 3 folds = 15 held-out assessments per
#   model-feature-task combination. Each reference record has one OOF prediction
#   per CV repeat (5 total), never from a model that trained on its own site.
# - Unknown-site probabilities are affinities to the reference classes under these
#   balanced classification tasks; they are NOT literal archaeological posterior
#   probabilities.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

INPUT_FILE <- "tapho_db_sum.csv"
OUT_DIR <- "supervised_ml_7spaces_3tasks_R_outputs"

SEED <- 42
CV_FOLDS <- 3L
CV_REPEATS <- 5L

RF_TREES <- 80L
RF_MIN_NODE <- 2L
SVM_COST <- 1
KNN_K <- 7L
LOGISTIC_LAMBDA <- 1

CONFORMAL_ALPHA <- 0.10

REDISTRIBUTE_BILATERAL_INDET <- TRUE
DROP_REDISTRIBUTED_INDET_FROM_FEATURES <- TRUE
PCTMNAU_MAX <- 1

# Current request/example includes Sima de los Huesos PNAS.
# Change to TRUE if you later want to omit it.
EXCLUDE_SIMA_PNAS <- FALSE
SIMA_PNAS_NAMES <- c("Sima de los Huesos PNAS")

# For the compact unknown table that mirrors the table in the prompt.
COMPACT_AFFINITY_SPACES <- c(
  "C + Heip evenness",
  "%MNAU raw",
  "ReMNAU raw",
  "Engineered",
  "Combined"
)

# Descriptive bands for "how many feature spaces are F / middle / NF".
# These are reporting bands, NOT classification thresholds.
AFFINITY_LOW <- 0.40
AFFINITY_HIGH <- 0.60

# Hard binary classification remains P(F) >= 0.5.
CLASS_THRESHOLD <- 0.50

# ------------------------------------------------------------------------------
# 1. Startup and package checks
# ------------------------------------------------------------------------------

cat("\n=== SUPERVISED 7-SPACE / 3-TASK ML PIPELINE ===\n")
cat("Working directory:", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "\n")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(OUT_DIR)) {
  stop("Could not create output directory: ", OUT_DIR)
}

ROOT_OUT <- normalizePath(OUT_DIR, winslash = "/", mustWork = FALSE)
cat("Output directory :", ROOT_OUT, "\n")

writeLines(
  c(
    paste0("Run started: ", Sys.time()),
    paste0("Working directory: ", getwd()),
    paste0("Input: ", INPUT_FILE),
    paste0("Output: ", ROOT_OUT)
  ),
  file.path(ROOT_OUT, "RUN_STARTED.txt")
)

required_packages <- c("tidyverse", "glmnet", "ranger", "e1071")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  writeLines(
    paste("Missing packages:", paste(missing_packages, collapse = ", ")),
    file.path(ROOT_OUT, "RUN_STOPPED_MISSING_PACKAGES.txt")
  )
  stop(
    "Install the following packages and rerun:\n  ",
    paste(missing_packages, collapse = ", ")
  )
}

library(tidyverse)

set.seed(SEED)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file was not found: ", INPUT_FILE,
    "\nCurrent working directory: ", getwd()
  )
}

d_raw <- readr::read_csv(INPUT_FILE, show_col_types = FALSE)

required_columns <- c(
  "site_name", "taphonomic_context", "funerary_context",
  "Mean_part_abundance", "mnau_heip_evenness"
)
missing_required <- setdiff(required_columns, names(d_raw))
if (length(missing_required) > 0) {
  stop(
    "Missing required columns in tapho_db_sum: ",
    paste(missing_required, collapse = ", ")
  )
}

if (!"layer_name" %in% names(d_raw)) d_raw$layer_name <- NA_character_

d <- d_raw %>%
  mutate(record_id = row_number(), .before = 1)

writeLines(
  c(
    paste0("Input rows: ", nrow(d)),
    paste0("Input columns: ", ncol(d)),
    paste0("Unique site_name values: ", n_distinct(d$site_name))
  ),
  file.path(ROOT_OUT, "RUN_INPUT_OK.txt")
)

# ------------------------------------------------------------------------------
# 2. Bilateral indeterminate redistribution
# ------------------------------------------------------------------------------

find_bilateral_triplets <- function(data_names, prefix, feature_type) {
  indet <- data_names[
    stringr::str_starts(data_names, prefix) &
      stringr::str_ends(data_names, "_indet")
  ]

  if (length(indet) == 0) return(tibble())

  tibble(
    feature_type = feature_type,
    indet_col = indet,
    base = stringr::str_remove(indet, "_indet$"),
    left_col = paste0(base, "_left"),
    right_col = paste0(base, "_right")
  ) %>%
    filter(left_col %in% data_names, right_col %in% data_names) %>%
    mutate(
      anatomical_part = stringr::str_remove(base, paste0("^", prefix))
    )
}

# Raw ReMNAU columns only. If adjusted copies exist, do not redistribute them here.
raw_rem_names <- names(d)[
  str_starts(names(d), "ReMNAU_") &
    !str_starts(names(d), "ReMNAU_adjusted_")
]

indet_pairs <- bind_rows(
  find_bilateral_triplets(
    names(d),
    "ProportionalRep_MNIbased_",
    "%MNAU"
  ),
  find_bilateral_triplets(
    raw_rem_names,
    "ReMNAU_",
    "ReMNAU"
  )
)

redistribute_bilateral_pair <- function(
  dat, indet_col, left_col, right_col,
  feature_type, anatomical_part, max_value = Inf
) {
  i_before <- dat[[indet_col]]
  l_before <- dat[[left_col]]
  r_before <- dat[[right_col]]

  i_after <- i_before
  l_after <- l_before
  r_after <- r_before

  audit_rows <- list()
  ai <- 1L

  rows <- which(is.finite(i_before) & i_before > 0)

  for (j in rows) {
    amount <- i_before[j]

    l0 <- if (is.na(l_before[j])) 0 else l_before[j]
    r0 <- if (is.na(r_before[j])) 0 else r_before[j]

    if (is.finite(max_value)) {
      cap_l <- max(0, max_value - l0)
      cap_r <- max(0, max_value - r0)
    } else {
      cap_l <- Inf
      cap_r <- Inf
    }

    # Start with the most even possible 50:50 split.
    add_l <- if (is.finite(cap_l)) min(amount / 2, cap_l) else amount / 2
    add_r <- if (is.finite(cap_r)) min(amount / 2, cap_r) else amount / 2
    remaining <- amount - add_l - add_r

    # Reallocate any remainder to whichever side still has capacity.
    if (remaining > 0) {
      cap_l2 <- if (is.finite(cap_l)) max(0, cap_l - add_l) else Inf
      take_l <- if (is.finite(cap_l2)) min(remaining, cap_l2) else remaining
      add_l <- add_l + take_l
      remaining <- remaining - take_l
    }

    if (remaining > 0) {
      cap_r2 <- if (is.finite(cap_r)) max(0, cap_r - add_r) else Inf
      take_r <- if (is.finite(cap_r2)) min(remaining, cap_r2) else remaining
      add_r <- add_r + take_r
      remaining <- remaining - take_r
    }

    l_after[j] <- l0 + add_l
    r_after[j] <- r0 + add_r
    i_after[j] <- max(0, remaining)

    if (is.finite(max_value)) {
      if (l0 <= max_value + 1e-12) l_after[j] <- min(l_after[j], max_value)
      if (r0 <= max_value + 1e-12) r_after[j] <- min(r_after[j], max_value)
    }

    audit_rows[[ai]] <- tibble(
      record_id = dat$record_id[j],
      site_name = dat$site_name[j],
      taphonomic_context = dat$taphonomic_context[j],
      feature_type = feature_type,
      anatomical_part = anatomical_part,
      indet_col = indet_col,
      left_col = left_col,
      right_col = right_col,
      indet_before = i_before[j],
      left_before = l_before[j],
      right_before = r_before[j],
      assigned_left = add_l,
      assigned_right = add_r,
      residual_indet = i_after[j],
      left_after = l_after[j],
      right_after = r_after[j],
      pct_cap = if (is.finite(max_value)) max_value else NA_real_
    )
    ai <- ai + 1L
  }

  dat[[indet_col]] <- i_after
  dat[[left_col]] <- l_after
  dat[[right_col]] <- r_after

  list(
    data = dat,
    audit = if (length(audit_rows) > 0) bind_rows(audit_rows) else tibble()
  )
}

indet_audit_parts <- list()

if (REDISTRIBUTE_BILATERAL_INDET && nrow(indet_pairs) > 0) {
  for (ii in seq_len(nrow(indet_pairs))) {
    pair <- indet_pairs[ii, ]
    max_val <- if (pair$feature_type == "%MNAU") PCTMNAU_MAX else Inf

    rr <- redistribute_bilateral_pair(
      dat = d,
      indet_col = pair$indet_col,
      left_col = pair$left_col,
      right_col = pair$right_col,
      feature_type = pair$feature_type,
      anatomical_part = pair$anatomical_part,
      max_value = max_val
    )

    d <- rr$data
    indet_audit_parts[[paste(pair$feature_type, pair$anatomical_part, sep = "_")]] <- rr$audit
  }
}

indet_redistribution_audit <- bind_rows(indet_audit_parts)

if (nrow(indet_redistribution_audit) > 0) {
  indet_redistribution_summary <- indet_redistribution_audit %>%
    group_by(feature_type, anatomical_part, indet_col, left_col, right_col) %>%
    summarise(
      n_records_redistributed = n(),
      total_indet_before = sum(indet_before, na.rm = TRUE),
      total_assigned_left = sum(assigned_left, na.rm = TRUE),
      total_assigned_right = sum(assigned_right, na.rm = TRUE),
      total_residual_indet = sum(residual_indet, na.rm = TRUE),
      n_with_residual = sum(residual_indet > 1e-12, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  indet_redistribution_summary <- indet_pairs %>%
    transmute(
      feature_type, anatomical_part, indet_col, left_col, right_col,
      n_records_redistributed = 0L,
      total_indet_before = 0,
      total_assigned_left = 0,
      total_assigned_right = 0,
      total_residual_indet = 0,
      n_with_residual = 0L
    )
}

write_csv(indet_pairs, file.path(ROOT_OUT, "indet_bilateral_pairs_detected.csv"))
write_csv(indet_redistribution_audit, file.path(ROOT_OUT, "indet_redistribution_audit.csv"))
write_csv(indet_redistribution_summary, file.path(ROOT_OUT, "indet_redistribution_summary.csv"))

# ------------------------------------------------------------------------------
# 3. Archaeological labels and target sets
# ------------------------------------------------------------------------------

collapse_class4 <- function(x) {
  case_when(
    x == "Funerary Context" ~ "Funerary Context",
    x == "Cannibalism" ~ "Cannibalism",
    x %in% c(
      "Carnivore",
      "Carnivore-related",
      "Carnivore Scavenge (in situ)",
      "Carnivore den(transport)",
      "Carnivore den (transport)"
    ) ~ "Carnivore",
    x == "Geological transport" ~ "Geological transport",
    TRUE ~ NA_character_
  )
}

d <- d %>%
  mutate(class4 = collapse_class4(taphonomic_context))

unknown <- d %>%
  filter(taphonomic_context %in% c("Unknown", "Unkown"))

if (EXCLUDE_SIMA_PNAS) {
  unknown <- unknown %>% filter(!site_name %in% SIMA_PNAS_NAMES)
}

sepulchral <- d %>%
  filter(taphonomic_context %in% c("Cuevas Sepulcrales", "Sepulchral Cave"))

known_reference <- d %>%
  filter(!is.na(class4))

partition <- tibble(
  dataset = c("Known reference", "Unknown/Unkown projection", "Sepulchral projection"),
  n_rows = c(nrow(known_reference), nrow(unknown), nrow(sepulchral)),
  unique_site_names = c(
    n_distinct(known_reference$site_name),
    n_distinct(unknown$site_name),
    n_distinct(sepulchral$site_name)
  )
)

write_csv(partition, file.path(ROOT_OUT, "sample_partition.csv"))
write_csv(
  known_reference %>% count(class4, name = "n"),
  file.path(ROOT_OUT, "known_reference_class_counts.csv")
)

# ------------------------------------------------------------------------------
# 4. Seven feature spaces
# ------------------------------------------------------------------------------

C_var <- "Mean_part_abundance"
E_var <- "mnau_heip_evenness"

prop_cols <- names(d)[
  str_starts(names(d), "ProportionalRep_MNIbased_")
]

rem_cols <- names(d)[
  str_starts(names(d), "ReMNAU_") &
    !str_starts(names(d), "ReMNAU_adjusted_")
]

prop_indet_exclude <- indet_pairs %>%
  filter(feature_type == "%MNAU") %>%
  pull(indet_col)

rem_indet_exclude <- indet_pairs %>%
  filter(feature_type == "ReMNAU") %>%
  pull(indet_col)

if (DROP_REDISTRIBUTED_INDET_FROM_FEATURES) {
  prop_cols <- setdiff(prop_cols, prop_indet_exclude)
  rem_cols <- setdiff(rem_cols, rem_indet_exclude)
}

eng_cols_requested <- c(
  "mnau_proximal_appendicular_mean",
  "mnau_lower_appendicular_mean",
  "mnau_distal_appendicular_mean",
  "mnau_appendicular_mean",
  "mnau_postcranial_total_mean",
  "mnau_top2_share_mcnaughton",
  "mnau_upper_appendicular_mean",
  "low_density_survival_index",
  "mnau_vertebrae_to_thorax_ratio",
  "all_sections_dominance",
  "mnau_lower_upper_appendicular_diff",
  "mnau_heip_evenness",
  "mnau_girdle_appendicular_logratio",
  "postcranial_continuity_hm",
  "mnau_girdle_to_appendicular_ratio",
  "mnau_gini",
  "density_representation_r2",
  "mnau_cv",
  "mnau_girdle_appendicular_diff",
  "mnau_full_axial_mean",
  "mnau_arms",
  "dominance_gap_all",
  "mnau_girdle_share_of_appendicular",
  "section_completeness_025_all",
  "whole_body_continuity_hm",
  "density_representation_slope",
  "mnau_limb_shaft_proxy_mean",
  "density_residual_sd"
)

eng_cols <- intersect(eng_cols_requested, names(d))

feature_sets <- list(
  `C only` = C_var,
  `Heip evenness only` = E_var,
  `C + Heip evenness` = c(C_var, E_var),
  `%MNAU raw` = prop_cols,
  `ReMNAU raw` = rem_cols,
  `Engineered` = eng_cols,
  `Combined` = unique(c(prop_cols, rem_cols, eng_cols))
)

if (any(lengths(feature_sets) == 0)) {
  stop(
    "At least one feature space contains zero variables: ",
    paste(names(feature_sets)[lengths(feature_sets) == 0], collapse = ", ")
  )
}

feature_manifest <- imap_dfr(
  feature_sets,
  ~ tibble(
    feature_space = .y,
    n_features_requested = length(.x),
    variable = .x
  )
)
write_csv(feature_manifest, file.path(ROOT_OUT, "feature_manifest.csv"))

# ------------------------------------------------------------------------------
# 5. Classification tasks
# ------------------------------------------------------------------------------

task_classes <- list(
  `4class` = c(
    "Cannibalism",
    "Carnivore",
    "Funerary Context",
    "Geological transport"
  ),
  `cannibalism_vs_fc` = c(
    "Cannibalism",
    "Funerary Context"
  ),
  `fc_vs_rest` = c(
    "Rest",
    "Funerary Context"
  )
)

make_task_data <- function(task_name, dat) {
  if (task_name == "4class") {
    return(
      dat %>%
        filter(!is.na(class4)) %>%
        mutate(task_label = class4)
    )
  }

  if (task_name == "cannibalism_vs_fc") {
    return(
      dat %>%
        filter(class4 %in% c("Cannibalism", "Funerary Context")) %>%
        mutate(task_label = class4)
    )
  }

  if (task_name == "fc_vs_rest") {
    return(
      dat %>%
        filter(!is.na(class4)) %>%
        mutate(
          task_label = if_else(
            class4 == "Funerary Context",
            "Funerary Context",
            "Rest"
          )
        )
    )
  }

  stop("Unknown task: ", task_name)
}

# ------------------------------------------------------------------------------
# 6. Grouped stratified fold map
# ------------------------------------------------------------------------------

make_grouped_stratified_folds <- function(
  dat, v = 3L, repeats = 5L, seed = 42L
) {
  if (any(is.na(dat$site_name) | dat$site_name == "")) {
    stop("site_name contains missing/blank values; grouped CV requires site_name.")
  }

  group_labels <- dat %>%
    distinct(site_name, task_label)

  inconsistent <- group_labels %>%
    count(site_name, name = "n_labels") %>%
    filter(n_labels > 1)

  if (nrow(inconsistent) > 0) {
    stop(
      "At least one site_name has multiple task labels within the same task: ",
      paste(head(inconsistent$site_name, 10), collapse = ", ")
    )
  }

  min_groups <- group_labels %>%
    count(task_label, name = "n_groups") %>%
    summarise(min_n = min(n_groups)) %>%
    pull(min_n)

  if (min_groups < v) {
    stop(
      "Cannot make ", v, " grouped folds because the smallest class has only ",
      min_groups, " unique site_name groups."
    )
  }

  maps <- vector("list", repeats)

  for (r in seq_len(repeats)) {
    set.seed(seed + r * 1009L)

    one <- group_labels %>%
      group_by(task_label) %>%
      group_modify(~ {
        g <- sample(.x$site_name, length(.x$site_name), replace = FALSE)
        tibble(
          site_name = g,
          fold = rep(seq_len(v), length.out = length(g))
        )
      }) %>%
      ungroup() %>%
      mutate(repeat_id = r)

    maps[[r]] <- one
  }

  bind_rows(maps)
}

# ------------------------------------------------------------------------------
# 7. Fold-specific preprocessing
# ------------------------------------------------------------------------------

fit_preprocessor <- function(train, vars) {
  med <- map_dbl(vars, function(v) {
    z <- train[[v]]
    m <- suppressWarnings(median(z, na.rm = TRUE))
    if (!is.finite(m)) 0 else m
  })
  names(med) <- vars

  impute <- function(dat) {
    x <- as.data.frame(dat[, vars, drop = FALSE])

    for (v in vars) {
      bad <- !is.finite(x[[v]]) | is.na(x[[v]])
      x[[v]][bad] <- med[[v]]
    }

    as.matrix(x)
  }

  xtr0 <- impute(train)

  mu <- colMeans(xtr0)
  sdv <- apply(xtr0, 2, sd)

  keep <- is.finite(sdv) & sdv > 0

  # If all variables are constant in a fold, stop rather than fabricate a model.
  if (!any(keep)) {
    stop("All variables are zero-variance in one training fold.")
  }

  keep_vars <- vars[keep]
  mu <- mu[keep]
  sdv <- sdv[keep]

  transform_raw <- function(dat) {
    xx <- impute(dat)
    xx[, keep, drop = FALSE]
  }

  transform_scaled <- function(dat) {
    xx <- transform_raw(dat)
    sweep(sweep(xx, 2, mu, "-"), 2, sdv, "/")
  }

  list(
    vars = vars,
    keep_vars = keep_vars,
    medians = med,
    means = mu,
    sds = sdv,
    train_raw = transform_raw(train),
    train_scaled = transform_scaled(train),
    transform_raw = transform_raw,
    transform_scaled = transform_scaled
  )
}

# ------------------------------------------------------------------------------
# 8. Training-fold oversampling
# ------------------------------------------------------------------------------

balance_indices <- function(y, seed) {
  y <- as.character(y)
  classes <- sort(unique(y))
  counts <- table(y)
  target_n <- max(counts)

  set.seed(seed)

  idx <- unlist(
    lapply(classes, function(cl) {
      ii <- which(y == cl)

      if (length(ii) < target_n) {
        c(ii, sample(ii, target_n - length(ii), replace = TRUE))
      } else {
        ii
      }
    }),
    use.names = FALSE
  )

  sample(idx, length(idx), replace = FALSE)
}

# ------------------------------------------------------------------------------
# 9. Four fixed classifiers
# ------------------------------------------------------------------------------

align_probability_matrix <- function(P, classes) {
  P <- as.matrix(P)

  if (is.null(colnames(P))) {
    if (ncol(P) != length(classes)) {
      stop("Probability matrix has no class names and unexpected number of columns.")
    }
    colnames(P) <- classes
  }

  missing <- setdiff(classes, colnames(P))
  if (length(missing) > 0) {
    stop(
      "Probability output is missing class(es): ",
      paste(missing, collapse = ", ")
    )
  }

  P <- P[, classes, drop = FALSE]

  # Numerical guard.
  P[P < 0] <- 0
  rs <- rowSums(P)
  rs[!is.finite(rs) | rs <= 0] <- 1
  P / rs
}

predict_logistic <- function(X_train, y_train, X_new, classes) {
  family_name <- if (length(classes) == 2) "binomial" else "multinomial"

   if (ncol(X_train) == 1L) {
    X_train_glmnet <- cbind(
      X_train,
      .glmnet_zero_dummy = rep(0, nrow(X_train))
    )
    X_new_glmnet <- cbind(
      X_new,
      .glmnet_zero_dummy = rep(0, nrow(X_new))
    )
  } else {
    X_train_glmnet <- X_train
    X_new_glmnet <- X_new
  }

  fit <- glmnet::glmnet(
    x = X_train_glmnet,
    y = factor(y_train, levels = classes),
    family = family_name,
    alpha = 0,
    lambda = LOGISTIC_LAMBDA,
    standardize = FALSE
  )

  pr <- predict(
    fit,
    newx = X_new_glmnet,
    type = "response",
    s = LOGISTIC_LAMBDA
  )

  if (length(classes) == 2) {
    # glmnet binomial response is P(second factor level).
    p2 <- as.numeric(pr)
    P <- cbind(1 - p2, p2)
    colnames(P) <- classes
  } else {
    # n x class x lambda
    P <- pr[, , 1, drop = FALSE]
    P <- P[, , 1]
    if (is.vector(P)) P <- matrix(P, nrow = nrow(X_new))
  }

  align_probability_matrix(P, classes)
}

predict_rf <- function(X_train, y_train, X_new, classes) {
  tr <- as.data.frame(X_train)
  names(tr) <- paste0("X", seq_len(ncol(tr)))
  tr$.y <- factor(y_train, levels = classes)

  te <- as.data.frame(X_new)
  names(te) <- paste0("X", seq_len(ncol(te)))

  mtry_val <- max(1L, floor(sqrt(ncol(X_train))))

  fit <- ranger::ranger(
    .y ~ .,
    data = tr,
    probability = TRUE,
    num.trees = RF_TREES,
    mtry = mtry_val,
    min.node.size = RF_MIN_NODE,
    seed = SEED
  )

  P <- predict(fit, data = te)$predictions
  align_probability_matrix(P, classes)
}

predict_svm <- function(X_train, y_train, X_new, classes) {
  fit <- e1071::svm(
    x = X_train,
    y = factor(y_train, levels = classes),
    kernel = "radial",
    cost = SVM_COST,
    gamma = 1 / ncol(X_train),
    scale = FALSE,
    probability = TRUE
  )

  pred <- predict(
    fit,
    X_new,
    probability = TRUE,
    decision.values = FALSE
  )

  P <- attr(pred, "probabilities")
  align_probability_matrix(P, classes)
}

predict_weighted_knn <- function(X_train, y_train, X_new, classes, k = 7L) {
  k_use <- min(k, nrow(X_train))
  P <- matrix(
    0,
    nrow = nrow(X_new),
    ncol = length(classes),
    dimnames = list(NULL, classes)
  )

  for (i in seq_len(nrow(X_new))) {
    dif <- sweep(X_train, 2, X_new[i, ], "-")
    dst <- sqrt(rowSums(dif^2))
    ord <- order(dst)[seq_len(k_use)]

    d <- dst[ord]
    yy <- y_train[ord]

    # Exact matches dominate without infinite numeric weights.
    if (any(d <= 1e-12)) {
      exact <- yy[d <= 1e-12]
      wtab <- table(factor(exact, levels = classes))
      P[i, ] <- as.numeric(wtab) / sum(wtab)
    } else {
      w <- 1 / d
      score <- vapply(
        classes,
        function(cl) sum(w[yy == cl]),
        numeric(1)
      )
      P[i, ] <- score / sum(score)
    }
  }

  P
}

predict_model <- function(
  model_name,
  X_train_scaled,
  X_train_raw,
  y_train,
  X_new_scaled,
  X_new_raw,
  classes
) {
  if (model_name == "Logistic regression") {
    return(predict_logistic(X_train_scaled, y_train, X_new_scaled, classes))
  }

  if (model_name == "Random forest") {
    return(predict_rf(X_train_raw, y_train, X_new_raw, classes))
  }

  if (model_name == "RBF-SVM") {
    return(predict_svm(X_train_scaled, y_train, X_new_scaled, classes))
  }

  if (model_name == "KNN") {
    return(
      predict_weighted_knn(
        X_train_scaled,
        y_train,
        X_new_scaled,
        classes,
        k = KNN_K
      )
    )
  }

  stop("Unknown model: ", model_name)
}

model_names <- c(
  "Logistic regression",
  "Random forest",
  "RBF-SVM",
  "KNN"
)

# ------------------------------------------------------------------------------
# 10. Metrics
# ------------------------------------------------------------------------------

binary_auc <- function(y01, score) {
  ok <- is.finite(score) & !is.na(y01)
  y01 <- y01[ok]
  score <- score[ok]

  n1 <- sum(y01 == 1)
  n0 <- sum(y01 == 0)

  if (n1 == 0 || n0 == 0) return(NA_real_)

  r <- rank(score, ties.method = "average")
  (sum(r[y01 == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

calc_metrics <- function(truth, P, classes) {
  truth <- factor(truth, levels = classes)
  pred <- factor(classes[max.col(P, ties.method = "first")], levels = classes)

  tab <- table(truth, pred)

  recalls <- vapply(seq_along(classes), function(j) {
    den <- sum(tab[j, ])
    if (den == 0) return(NA_real_)
    tab[j, j] / den
  }, numeric(1))

  precisions <- vapply(seq_along(classes), function(j) {
    den <- sum(tab[, j])
    if (den == 0) return(NA_real_)
    tab[j, j] / den
  }, numeric(1))

  balanced_accuracy <- mean(recalls, na.rm = TRUE)
  macro_precision <- mean(precisions, na.rm = TRUE)
  macro_recall <- mean(recalls, na.rm = TRUE)

  if (length(classes) == 2) {
    positive <- "Funerary Context"
    if (!positive %in% classes) positive <- classes[2]

    y01 <- as.integer(truth == positive)
    roc_auc <- binary_auc(y01, P[, positive])
    brier <- mean((P[, positive] - y01)^2)

    pos_j <- which(classes == positive)
    precision <- precisions[pos_j]
    recall <- recalls[pos_j]
  } else {
    aucs <- vapply(classes, function(cl) {
      binary_auc(
        as.integer(truth == cl),
        P[, cl]
      )
    }, numeric(1))

    roc_auc <- mean(aucs, na.rm = TRUE)

    Y <- sapply(classes, function(cl) as.integer(truth == cl))
    brier <- mean(rowSums((P - Y)^2))

    precision <- macro_precision
    recall <- macro_recall
  }

  tibble(
    n = length(truth),
    balanced_accuracy = balanced_accuracy,
    roc_auc = roc_auc,
    precision = precision,
    recall = recall,
    brier = brier
  )
}

prob_to_long <- function(
  P, meta, classes,
  task_name, feature_space, model_name,
  repeat_id, fold_id, dataset_role
) {
  as_tibble(P, .name_repair = "minimal") %>%
    setNames(classes) %>%
    bind_cols(
      meta %>%
        select(
          record_id, site_name, layer_name,
          taphonomic_context, funerary_context,
          any_of(c("task_label", "class4"))
        )
    ) %>%
    pivot_longer(
      cols = all_of(classes),
      names_to = "class",
      values_to = "probability"
    ) %>%
    mutate(
      task = task_name,
      feature_space = feature_space,
      model = model_name,
      repeat_id = repeat_id,
      fold = fold_id,
      dataset_role = dataset_role,
      .before = 1
    )
}

# ------------------------------------------------------------------------------
# 11. Main repeated grouped-CV loop
# ------------------------------------------------------------------------------

oof_parts <- list()
unknown_parts <- list()
sepulchral_parts <- list()
fold_metric_parts <- list()
fold_map_parts <- list()
counter <- 1L

for (task_name in names(task_classes)) {
  cat("\n--- TASK:", task_name, "---\n")

  task_dat <- make_task_data(task_name, known_reference)
  classes <- task_classes[[task_name]]

  task_dat <- task_dat %>%
    mutate(task_label = factor(task_label, levels = classes)) %>%
    arrange(record_id)

  fold_map <- make_grouped_stratified_folds(
    task_dat,
    v = CV_FOLDS,
    repeats = CV_REPEATS,
    seed = SEED + match(task_name, names(task_classes)) * 10000L
  ) %>%
    mutate(task = task_name, .before = 1)

  fold_map_parts[[task_name]] <- fold_map

  task_out <- file.path(ROOT_OUT, task_name)
  dir.create(task_out, recursive = TRUE, showWarnings = FALSE)
  write_csv(fold_map, file.path(task_out, "grouped_cv_fold_map.csv"))

  for (space_name in names(feature_sets)) {
    vars <- feature_sets[[space_name]]
    cat("  Feature space:", space_name, "(", length(vars), "requested variables )\n")

    for (r in seq_len(CV_REPEATS)) {
      fm_r <- fold_map %>% filter(repeat_id == r)

      for (f in seq_len(CV_FOLDS)) {
        val_sites <- fm_r %>%
          filter(fold == f) %>%
          pull(site_name)

        train_idx <- which(!task_dat$site_name %in% val_sites)
        val_idx <- which(task_dat$site_name %in% val_sites)

        train <- task_dat[train_idx, , drop = FALSE]
        valid <- task_dat[val_idx, , drop = FALSE]

        pp <- fit_preprocessor(train, vars)

        Xtr_s <- pp$train_scaled
        Xtr_r <- pp$train_raw

        Xva_s <- pp$transform_scaled(valid)
        Xva_r <- pp$transform_raw(valid)

        Xu_s <- pp$transform_scaled(unknown)
        Xu_r <- pp$transform_raw(unknown)

        Xsep_s <- pp$transform_scaled(sepulchral)
        Xsep_r <- pp$transform_raw(sepulchral)

        ytr <- as.character(train$task_label)
        yva <- as.character(valid$task_label)

        bal_idx <- balance_indices(
          ytr,
          seed = SEED +
            match(task_name, names(task_classes)) * 100000L +
            match(space_name, names(feature_sets)) * 10000L +
            r * 100L + f
        )

        Xtr_s_bal <- Xtr_s[bal_idx, , drop = FALSE]
        Xtr_r_bal <- Xtr_r[bal_idx, , drop = FALSE]
        ytr_bal <- ytr[bal_idx]

        for (model_name in model_names) {
          # Fit each fold-trained classifier ONCE, then use that exact fitted
          # classifier for the held-out fold + Unknown + Sepulchral targets.
          n_val <- nrow(valid)
          n_unk <- nrow(unknown)
          n_sep <- nrow(sepulchral)

          Xnew_s <- rbind(Xva_s, Xu_s, Xsep_s)
          Xnew_r <- rbind(Xva_r, Xu_r, Xsep_r)

          Pall <- predict_model(
            model_name,
            Xtr_s_bal, Xtr_r_bal, ytr_bal,
            Xnew_s, Xnew_r,
            classes
          )

          Pva <- Pall[seq_len(n_val), , drop = FALSE]

          met <- calc_metrics(yva, Pva, classes) %>%
            mutate(
              task = task_name,
              feature_space = space_name,
              model = model_name,
              repeat_id = r,
              fold = f,
              n_train_before_balance = length(ytr),
              n_train_after_balance = length(ytr_bal),
              n_features_used = ncol(Xtr_s),
              .before = 1
            )

          fold_metric_parts[[counter]] <- met
          counter <- counter + 1L

          oof_parts[[length(oof_parts) + 1L]] <- prob_to_long(
            Pva,
            valid,
            classes,
            task_name,
            space_name,
            model_name,
            r,
            f,
            "OOF reference"
          )

          if (n_unk > 0) {
            idx_u <- n_val + seq_len(n_unk)
            Pu <- Pall[idx_u, , drop = FALSE]

            unknown_parts[[length(unknown_parts) + 1L]] <- prob_to_long(
              Pu,
              unknown,
              classes,
              task_name,
              space_name,
              model_name,
              r,
              f,
              "Unknown/Unkown"
            )
          }

          if (n_sep > 0) {
            idx_s <- n_val + n_unk + seq_len(n_sep)
            Ps <- Pall[idx_s, , drop = FALSE]

            sepulchral_parts[[length(sepulchral_parts) + 1L]] <- prob_to_long(
              Ps,
              sepulchral,
              classes,
              task_name,
              space_name,
              model_name,
              r,
              f,
              "Sepulchral external"
            )
          }
        }
      }
    }
  }
}

fold_maps <- bind_rows(fold_map_parts)
fold_metrics <- bind_rows(fold_metric_parts)
oof_predictions <- bind_rows(oof_parts)
unknown_predictions <- bind_rows(unknown_parts)
sepulchral_predictions <- bind_rows(sepulchral_parts)

write_csv(fold_maps, file.path(ROOT_OUT, "all_task_grouped_cv_fold_maps.csv"))
write_csv(fold_metrics, file.path(ROOT_OUT, "fold_level_performance.csv"))
write_csv(oof_predictions, file.path(ROOT_OUT, "OOF_reference_probabilities_long.csv"))
write_csv(unknown_predictions, file.path(ROOT_OUT, "Unknown_probabilities_all_fold_models_long.csv"))
write_csv(sepulchral_predictions, file.path(ROOT_OUT, "Sepulchral_probabilities_all_fold_models_long.csv"))

# ------------------------------------------------------------------------------
# 12. OOF ensemble probabilities and pooled validation performance
# ------------------------------------------------------------------------------

# Each reference record has one OOF prediction per CV repeat. Average those five
# OOF probabilities before calculating the principal pooled metrics.
oof_site_ensemble <- oof_predictions %>%
  group_by(
    task, feature_space, model,
    record_id, site_name, task_label, class
  ) %>%
  summarise(
    probability = mean(probability, na.rm = TRUE),
    n_oof_repeats = n(),
    .groups = "drop"
  )

write_csv(
  oof_site_ensemble,
  file.path(ROOT_OUT, "OOF_reference_probabilities_site_ensemble.csv")
)

calc_metrics_from_long <- function(df, task_name) {
  # dplyr::group_modify() stores grouping keys in .y rather than .x, so the
  # task name is passed explicitly from the group key.
  desired <- task_classes[[task_name]]
  classes_present <- unique(as.character(df$class))
  classes <- desired[desired %in% classes_present]

  if (length(classes) < 2) {
    stop("Fewer than two task classes were found while calculating pooled OOF metrics for task: ", task_name)
  }

  truth <- df %>%
    distinct(record_id, task_label) %>%
    arrange(record_id)

  pw <- df %>%
    select(record_id, class, probability) %>%
    pivot_wider(
      names_from = class,
      values_from = probability
    ) %>%
    arrange(record_id)

  missing_prob_cols <- setdiff(classes, names(pw))
  if (length(missing_prob_cols) > 0) {
    stop(
      "Missing probability column(s) in pooled OOF summary for task ",
      task_name, ": ", paste(missing_prob_cols, collapse = ", ")
    )
  }

  P <- as.matrix(pw[, classes, drop = FALSE])
  colnames(P) <- classes

  calc_metrics(
    truth = as.character(truth$task_label),
    P = P,
    classes = classes
  )
}

pooled_performance <- oof_site_ensemble %>%
  group_by(task, feature_space, model) %>%
  group_modify(~ calc_metrics_from_long(.x, task_name = .y$task[[1]])) %>%
  ungroup()

fold_performance_summary <- fold_metrics %>%
  group_by(task, feature_space, model) %>%
  summarise(
    n_assessment_folds = n(),
    BA_mean = mean(balanced_accuracy, na.rm = TRUE),
    BA_sd = sd(balanced_accuracy, na.rm = TRUE),
    AUC_mean = mean(roc_auc, na.rm = TRUE),
    AUC_sd = sd(roc_auc, na.rm = TRUE),
    Precision_mean = mean(precision, na.rm = TRUE),
    Precision_sd = sd(precision, na.rm = TRUE),
    Recall_mean = mean(recall, na.rm = TRUE),
    Recall_sd = sd(recall, na.rm = TRUE),
    Brier_mean = mean(brier, na.rm = TRUE),
    Brier_sd = sd(brier, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  pooled_performance,
  file.path(ROOT_OUT, "pooled_OOF_performance.csv")
)
write_csv(
  fold_performance_summary,
  file.path(ROOT_OUT, "fold_performance_mean_SD.csv")
)

# ------------------------------------------------------------------------------
# 13. Unknown and sepulchral ensemble probabilities
# ------------------------------------------------------------------------------

# Mean across the 15 fold-trained models for each task x feature x algorithm.
unknown_pipeline_ensemble <- unknown_predictions %>%
  group_by(
    task, feature_space, model,
    record_id, site_name, class
  ) %>%
  summarise(
    probability = mean(probability, na.rm = TRUE),
    probability_sd_across_fold_models = sd(probability, na.rm = TRUE),
    n_fold_models = n(),
    .groups = "drop"
  )

sepulchral_pipeline_ensemble <- sepulchral_predictions %>%
  group_by(
    task, feature_space, model,
    record_id, site_name, class
  ) %>%
  summarise(
    probability = mean(probability, na.rm = TRUE),
    probability_sd_across_fold_models = sd(probability, na.rm = TRUE),
    n_fold_models = n(),
    .groups = "drop"
  )

write_csv(
  unknown_pipeline_ensemble,
  file.path(ROOT_OUT, "Unknown_probability_by_task_space_model.csv")
)
write_csv(
  sepulchral_pipeline_ensemble,
  file.path(ROOT_OUT, "Sepulchral_probability_by_task_space_model.csv")
)

# Average the four algorithms within each feature space.
unknown_space_ensemble <- unknown_pipeline_ensemble %>%
  group_by(task, feature_space, record_id, site_name, class) %>%
  summarise(
    probability = mean(probability, na.rm = TRUE),
    model_probability_sd = sd(probability, na.rm = TRUE),
    n_models = n(),
    .groups = "drop"
  )

sepulchral_space_ensemble <- sepulchral_pipeline_ensemble %>%
  group_by(task, feature_space, record_id, site_name, class) %>%
  summarise(
    probability = mean(probability, na.rm = TRUE),
    model_probability_sd = sd(probability, na.rm = TRUE),
    n_models = n(),
    .groups = "drop"
  )

write_csv(
  unknown_space_ensemble,
  file.path(ROOT_OUT, "Unknown_probability_by_task_feature_space.csv")
)
write_csv(
  sepulchral_space_ensemble,
  file.path(ROOT_OUT, "Sepulchral_probability_by_task_feature_space.csv")
)

# ------------------------------------------------------------------------------
# 14. Primary fc_vs_rest unknown affinity table
# ------------------------------------------------------------------------------

unknown_fc_raw <- unknown_predictions %>%
  filter(
    task == "fc_vs_rest",
    class == "Funerary Context"
  )

unknown_fc_space <- unknown_fc_raw %>%
  group_by(record_id, site_name, feature_space) %>%
  summarise(
    P_F = mean(probability, na.rm = TRUE),
    SD_P_F = sd(probability, na.rm = TRUE),
    Decision_F_percent = 100 * mean(probability >= CLASS_THRESHOLD, na.rm = TRUE),
    n_underlying_predictions = n(),
    .groups = "drop"
  )

write_csv(
  unknown_fc_space,
  file.path(ROOT_OUT, "Unknown_fc_vs_rest_feature_space_affinity_long.csv")
)

# Full seven-space probability table.
unknown_fc_7space_wide <- unknown_fc_space %>%
  select(record_id, site_name, feature_space, P_F) %>%
  pivot_wider(
    names_from = feature_space,
    values_from = P_F
  )

# Compact five-space table matching the requested report format.
compact_long <- unknown_fc_space %>%
  filter(feature_space %in% COMPACT_AFFINITY_SPACES)

compact_counts <- compact_long %>%
  group_by(record_id, site_name) %>%
  summarise(
    n_spaces = n(),
    F_spaces_0.60plus = sum(P_F > AFFINITY_HIGH),
    Middle_spaces_0.40to0.60 = sum(P_F >= AFFINITY_LOW & P_F <= AFFINITY_HIGH),
    NF_spaces_below0.40 = sum(P_F < AFFINITY_LOW),
    Hard_F_spaces_P0.50plus = sum(P_F >= CLASS_THRESHOLD),
    Hard_NF_spaces_below0.50 = sum(P_F < CLASS_THRESHOLD),
    Mean_P_F = mean(P_F),
    Min_P_F = min(P_F),
    Max_P_F = max(P_F),
    Range_P_F = max(P_F) - min(P_F),
    .groups = "drop"
  ) %>%
  mutate(
    Interpretation = case_when(
      Hard_F_spaces_P0.50plus == n_spaces &
        Mean_P_F >= 0.60 ~ "robust funerary affinity",

      Hard_NF_spaces_below0.50 == n_spaces &
        Mean_P_F <= 0.40 ~ "robust non-funerary affinity",

      (Max_P_F > 0.60 & Min_P_F < 0.40) |
        Range_P_F >= 0.35 ~ "representation-sensitive",

      Hard_F_spaces_P0.50plus >= ceiling(0.80 * n_spaces) &
        Mean_P_F >= 0.55 ~ "predominantly funerary affinity",

      Hard_NF_spaces_below0.50 >= ceiling(0.80 * n_spaces) &
        Mean_P_F <= 0.45 ~ "predominantly non-funerary affinity",

      TRUE ~ "intermediate / ambiguous"
    )
  )

compact_wide <- compact_long %>%
  select(record_id, site_name, feature_space, P_F) %>%
  pivot_wider(
    names_from = feature_space,
    values_from = P_F
  )

# Overall vote across ALL raw fold-model predictions for this task.
overall_fc_votes <- unknown_fc_raw %>%
  group_by(record_id, site_name) %>%
  summarise(
    Overall_mean_P_F = mean(probability, na.rm = TRUE),
    Overall_SD_P_F = sd(probability, na.rm = TRUE),
    Decision_F_percent_all_models_spaces_folds =
      100 * mean(probability >= CLASS_THRESHOLD, na.rm = TRUE),
    n_total_predictions = n(),
    .groups = "drop"
  )

unknown_fc_compact_report <- compact_wide %>%
  left_join(compact_counts, by = c("record_id", "site_name")) %>%
  left_join(overall_fc_votes, by = c("record_id", "site_name")) %>%
  arrange(desc(Mean_P_F))

write_csv(
  unknown_fc_compact_report,
  file.path(ROOT_OUT, "UNKNOWN_REPORT_TABLE_fc_vs_rest_compact_5spaces.csv")
)

unknown_fc_7space_report <- unknown_fc_7space_wide %>%
  left_join(overall_fc_votes, by = c("record_id", "site_name")) %>%
  arrange(desc(Overall_mean_P_F))

write_csv(
  unknown_fc_7space_report,
  file.path(ROOT_OUT, "UNKNOWN_REPORT_TABLE_fc_vs_rest_all_7spaces.csv")
)

unknown_fc_7space_counts <- unknown_fc_space %>%
  group_by(record_id, site_name) %>%
  summarise(
    n_spaces = n(),
    F_spaces_0.60plus = sum(P_F > AFFINITY_HIGH),
    Middle_spaces_0.40to0.60 = sum(P_F >= AFFINITY_LOW & P_F <= AFFINITY_HIGH),
    NF_spaces_below0.40 = sum(P_F < AFFINITY_LOW),
    Hard_F_spaces_P0.50plus = sum(P_F >= CLASS_THRESHOLD),
    Hard_NF_spaces_below0.50 = sum(P_F < CLASS_THRESHOLD),
    Mean_P_F = mean(P_F),
    Min_P_F = min(P_F),
    Max_P_F = max(P_F),
    Range_P_F = max(P_F) - min(P_F),
    .groups = "drop"
  ) %>%
  arrange(desc(Mean_P_F))

write_csv(
  unknown_fc_7space_counts,
  file.path(ROOT_OUT, "UNKNOWN_fc_vs_rest_7space_vote_and_middle_counts.csv")
)

# ------------------------------------------------------------------------------
# 15. Cannibalism-vs-Funerary unknown table
# ------------------------------------------------------------------------------

unknown_cann_fc <- unknown_predictions %>%
  filter(
    task == "cannibalism_vs_fc",
    class == "Funerary Context"
  ) %>%
  group_by(record_id, site_name, feature_space) %>%
  summarise(
    P_F = mean(probability, na.rm = TRUE),
    Decision_F_percent = 100 * mean(probability >= CLASS_THRESHOLD, na.rm = TRUE),
    .groups = "drop"
  )

unknown_cann_fc_wide <- unknown_cann_fc %>%
  select(record_id, site_name, feature_space, P_F) %>%
  pivot_wider(names_from = feature_space, values_from = P_F) %>%
  arrange(site_name)

write_csv(
  unknown_cann_fc_wide,
  file.path(ROOT_OUT, "UNKNOWN_REPORT_TABLE_cannibalism_vs_fc_P_F.csv")
)

# ------------------------------------------------------------------------------
# 16. Four-class unknown probabilities and consensus
# ------------------------------------------------------------------------------

unknown_4class_space <- unknown_space_ensemble %>%
  filter(task == "4class")

unknown_4class_space_prediction <- unknown_4class_space %>%
  group_by(record_id, site_name, feature_space) %>%
  slice_max(probability, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    record_id, site_name, feature_space,
    predicted_class = class,
    predicted_probability = probability
  )

unknown_4class_overall_probs <- unknown_predictions %>%
  filter(task == "4class") %>%
  group_by(record_id, site_name, class) %>%
  summarise(
    probability = mean(probability, na.rm = TRUE),
    .groups = "drop"
  )

unknown_4class_overall_prediction <- unknown_4class_overall_probs %>%
  group_by(record_id, site_name) %>%
  slice_max(probability, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    record_id, site_name,
    overall_predicted_class = class,
    overall_predicted_probability = probability
  )

unknown_4class_votes <- unknown_4class_space_prediction %>%
  count(record_id, site_name, predicted_class, name = "feature_spaces_voting") %>%
  group_by(record_id, site_name) %>%
  mutate(
    total_feature_spaces = sum(feature_spaces_voting),
    vote_percent = 100 * feature_spaces_voting / total_feature_spaces
  ) %>%
  ungroup()

write_csv(
  unknown_4class_space_prediction,
  file.path(ROOT_OUT, "Unknown_4class_prediction_by_feature_space.csv")
)
write_csv(
  unknown_4class_overall_probs,
  file.path(ROOT_OUT, "Unknown_4class_overall_class_probabilities.csv")
)
write_csv(
  unknown_4class_overall_prediction,
  file.path(ROOT_OUT, "Unknown_4class_overall_prediction.csv")
)
write_csv(
  unknown_4class_votes,
  file.path(ROOT_OUT, "Unknown_4class_feature_space_vote_counts.csv")
)

# ------------------------------------------------------------------------------
# 17. Cross-conformal summaries from ensemble OOF predictions
# ------------------------------------------------------------------------------

# One calibration score per reference record per task x space x model:
# score = 1 - OOF ensemble probability assigned to the true class.
calibration_scores <- oof_site_ensemble %>%
  filter(as.character(task_label) == class) %>%
  transmute(
    task, feature_space, model,
    record_id, site_name,
    true_class = class,
    nonconformity = 1 - probability
  )

write_csv(
  calibration_scores,
  file.path(ROOT_OUT, "cross_conformal_calibration_scores.csv")
)

conformal_for_targets <- function(target_pipeline_probs, calibration, alpha = 0.10) {
  target_pipeline_probs %>%
    mutate(
      nonconformity = 1 - probability,
      calibration_n = pmap_int(
        list(task, feature_space, model, class),
        function(tt, ss, mm, cc) {
          calibration %>%
            filter(
              task == tt,
              feature_space == ss,
              model == mm,
              true_class == cc
            ) %>%
            nrow()
        }
      ),
      conformal_p = pmap_dbl(
        list(task, feature_space, model, class, nonconformity),
        function(tt, ss, mm, cc, snew) {
          z <- calibration %>%
            filter(
              task == tt,
              feature_space == ss,
              model == mm,
              true_class == cc
            ) %>%
            pull(nonconformity)

          if (length(z) == 0) {
            NA_real_
          } else {
            (1 + sum(z >= snew)) / (length(z) + 1)
          }
        }
      ),
      included = is.finite(conformal_p) & conformal_p > alpha
    )
}

unknown_conformal <- conformal_for_targets(
  unknown_pipeline_ensemble,
  calibration_scores,
  alpha = CONFORMAL_ALPHA
)

sepulchral_conformal <- conformal_for_targets(
  sepulchral_pipeline_ensemble,
  calibration_scores,
  alpha = CONFORMAL_ALPHA
)

write_csv(
  unknown_conformal,
  file.path(ROOT_OUT, "Unknown_cross_conformal_by_task_space_model_class.csv")
)
write_csv(
  sepulchral_conformal,
  file.path(ROOT_OUT, "Sepulchral_cross_conformal_by_task_space_model_class.csv")
)

unknown_conformal_summary <- unknown_conformal %>%
  group_by(task, record_id, site_name, class) %>%
  summarise(
    pipelines = n(),
    Included_percent = 100 * mean(included, na.rm = TRUE),
    Median_conformal_p = median(conformal_p, na.rm = TRUE),
    Mean_conformal_p = mean(conformal_p, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  unknown_conformal_summary,
  file.path(ROOT_OUT, "Unknown_cross_conformal_summary.csv")
)

# ------------------------------------------------------------------------------
# 18. Simple probability plots for the primary unknown task
# ------------------------------------------------------------------------------

plot_dat <- unknown_fc_space %>%
  mutate(
    feature_space = factor(feature_space, levels = names(feature_sets)),
    site_name = forcats::fct_reorder(site_name, P_F, .fun = mean)
  )

if (nrow(plot_dat) > 0) {
  gp <- ggplot(
    plot_dat,
    aes(x = P_F, y = site_name, shape = feature_space)
  ) +
    geom_vline(xintercept = 0.5, linetype = 2, colour = "grey45") +
    geom_vline(xintercept = c(AFFINITY_LOW, AFFINITY_HIGH),
               linetype = 3, colour = "grey70") +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      x = "Mean P(F) across 4 algorithms and repeated grouped-CV fold models",
      y = NULL,
      shape = "Feature space",
      title = "Unknown-site funerary-reference affinity",
      subtitle = "fc_vs_rest task; 0.5 is the hard decision threshold"
    ) +
    theme_bw()

  ggsave(
    file.path(ROOT_OUT, "Unknown_fc_vs_rest_feature_space_affinity.png"),
    gp,
    width = 11,
    height = max(7, 0.28 * n_distinct(plot_dat$site_name) + 3),
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 19. Reporting notes
# ------------------------------------------------------------------------------

writeLines(
  c(
    "SUPERVISED ML REPORTING NOTES",
    "",
    paste0("Input: ", INPUT_FILE),
    "Independent CV grouping variable: site_name.",
    paste0("Repeated grouped CV: ", CV_REPEATS, " repeats x ", CV_FOLDS, " folds."),
    "A common fold map is reused across all seven feature spaces within each task.",
    "All preprocessing parameters are estimated from the training fold only.",
    "Minority classes are randomly oversampled within the training fold only.",
    "Validation folds, Unknown/Unkown targets and Sepulchral targets are never oversampled.",
    "",
    "Algorithms:",
    paste0("- Ridge logistic regression: glmnet alpha=0, lambda=", LOGISTIC_LAMBDA),
    "- For the one-variable C-only and Heip-only spaces, a constant-zero dummy column is appended only to satisfy glmnet's >=2-column interface; it contains no information and is not a substantive feature.",
    paste0("- Random forest: ranger, ", RF_TREES, " trees, min.node.size=", RF_MIN_NODE),
    paste0("- RBF-SVM: cost=", SVM_COST, ", gamma=1/p"),
    paste0("- distance-weighted KNN: k=", KNN_K),
    "",
    "Seven feature spaces:",
    paste0("- ", names(feature_sets), " (", lengths(feature_sets), " variables requested)"),
    "",
    "Unknown-site P(F) values are reference-affinity scores under a balanced classification problem.",
    "They should not be described as literal probabilities that an archaeological site was funerary.",
    "The compact report table averages the four algorithms and all fold-trained models within each feature space.",
    paste0("Descriptive affinity bands: < ", AFFINITY_LOW, " = NF-like; ",
           AFFINITY_LOW, "-", AFFINITY_HIGH, " = middle; > ",
           AFFINITY_HIGH, " = F-like."),
    "Hard model decisions use P(F) >= 0.5.",
    "",
    "Cross-conformal prediction:",
    paste0("- alpha = ", CONFORMAL_ALPHA),
    "- calibration uses ensemble OOF probabilities, so no site's own fitted prediction is used as its calibration score.",
    "- conformal inclusion is a set-valued reference-compatibility statement, not a posterior probability."
  ),
  file.path(ROOT_OUT, "REPORTING_NOTES_supervised_pipeline.txt")
)

# ------------------------------------------------------------------------------
# 20. Compact console summary
# ------------------------------------------------------------------------------

cat("\n=== PIPELINE COMPLETE ===\n")
cat("Known reference rows :", nrow(known_reference), "\n")
cat("Unknown rows         :", nrow(unknown), "\n")
cat("Sepulchral rows      :", nrow(sepulchral), "\n")
cat("Tasks                :", paste(names(task_classes), collapse = ", "), "\n")
cat("Feature spaces       :", paste(names(feature_sets), collapse = ", "), "\n")
cat("Algorithms           :", paste(model_names, collapse = ", "), "\n")
cat("CV assessments/pipeline:", CV_REPEATS * CV_FOLDS, "\n")
cat("Output directory     :", ROOT_OUT, "\n\n")

cat("Primary unknown table:\n")
cat(file.path(ROOT_OUT, "UNKNOWN_REPORT_TABLE_fc_vs_rest_compact_5spaces.csv"), "\n")

writeLines(
  paste0("Run finished: ", Sys.time()),
  file.path(ROOT_OUT, "RUN_FINISHED.txt")
)

capture.output(
  sessionInfo(),
  file = file.path(ROOT_OUT, "sessionInfo.txt")
)

# ==============================================================================
# END
# ==============================================================================

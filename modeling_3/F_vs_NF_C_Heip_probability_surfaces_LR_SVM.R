# ==============================================================================
# F vs NF probability surfaces in C x Heip space
# Models: ridge logistic regression + RBF-SVM
#
# DISCRETE-COLOUR VERSION
# -----------------------
# No gradient / continuous colour maps are produced.
#
# Probability zones are fixed at:
#   P(F) < 0.25        = strong NF affinity
#   0.25 <= P(F) < .50 = NF leaning
#   0.50 <= P(F) < .75 = F leaning
#   P(F) >= 0.75       = strong F affinity
#
# Three plot variants are produced for each model and palette:
#   1) reference only: no Unknown sites
#   2) Sima only: 89 -> 90 -> 95 -> 15 -> 24 joined by a spline
#   3) restricted Unknown sites + the Sima spline
#
# Validation:
#   - grouped by site_name
#   - 3 folds x 5 repeats
#   - training-fold median imputation
#   - training-fold standardization
#   - random oversampling in training folds only
#
# The probability surface is the MEAN P(F) across the 15 fold-trained models.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. User settings
# ------------------------------------------------------------------------------

INPUT_FILE <- "tapho_db_sum.csv"
OUT_DIR <- "F_vs_NF_C_Heip_probability_surfaces"

SEED <- 42L
CV_FOLDS <- 3L
CV_REPEATS <- 5L

LOGISTIC_LAMBDA <- 1
SVM_COST <- 1

GRID_N <- 280L
CLASS_THRESHOLD <- 0.50

# Restricted Unknown sites requested for the final panel.
RESTRICTED_UNKNOWN_SITES <- c(
  "Sima de los Huesos 15",
  "Sima de los Huesos 24",
  "Sima de los Huesos 95",
  "Spy cave",
  "Kleine Feldhofer Grotte",
  "Saint-Césaire",
  "Sima de los Huesos 90",
  "Sima de los Huesos 89",
  "Trou de la Naulette",
  "Cova del Parpalló (Parpalló Cave)",
  "Weimar-Ehringsdorf",
  "Peștera Muierii",
  "Bolomor Cave",
  "Bilzingsleben"
)

# Requested Sima sequence.
SIMA_TRAJECTORY <- c(
  "Sima de los Huesos 89",
  "Sima de los Huesos 90",
  "Sima de los Huesos 95",
  "Sima de los Huesos 15",
  "Sima de los Huesos 24"
)

# --------------------------------------------------------------------------
# DISCRETE palettes only.
#
# Background colours are deliberately light so the scatter points remain
# clearly visible. Points use much darker colours than any background zone.
# --------------------------------------------------------------------------

PALETTES <- list(

  # Similar to the original figures, but with stronger separation.
  Archaeology = c(
    "<25% F"   = "#D8CDE8",  # pale purple
    "25-50% F" = "#C8DCE9",  # pale blue
    "50-75% F" = "#CBE6D8",  # pale green
    ">=75% F"  = "#F3DFA0"   # pale gold
  ),

  # Colour-blind-friendly light backgrounds.
  Colourblind = c(
    "<25% F"   = "#D9D9D9",  # light grey
    "25-50% F" = "#BFD7EA",  # light blue
    "50-75% F" = "#DDE8B6",  # yellow-green
    ">=75% F"  = "#F6CF9B"   # light orange
  ),

  # Slightly stronger categorical contrast while still preserving points.
  Contrast = c(
    "<25% F"   = "#D4C2E8",  # purple
    "25-50% F" = "#AFCFE3",  # blue
    "50-75% F" = "#B9DDAF",  # green
    ">=75% F"  = "#F0CF7C"   # gold
  )
)

# Dark point colours for maximum contrast with the light background.
REFERENCE_NF_COLOUR <- "#1565A8"
REFERENCE_NF_FILL   <- "#5AA0CF"

REFERENCE_F_COLOUR  <- "#A95500"
REFERENCE_F_FILL    <- "#E89A49"

UNKNOWN_COLOUR      <- "#176B2C"
SIMA_CURVE_COLOUR   <- "#111111"

# Probability-contour line colours.
CONTOUR_25_COLOUR <- "#5B3A78"
CONTOUR_50_COLOUR <- "#111111"
CONTOUR_75_COLOUR <- "#8A6500"


# ------------------------------------------------------------------------------
# 1. Packages and data
# ------------------------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "glmnet",
  "e1071",
  "ggrepel"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Install these packages and rerun: ",
    paste(missing_packages, collapse = ", ")
  )
}

library(tidyverse)
library(ggrepel)

set.seed(SEED)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!file.exists(INPUT_FILE)) {
  stop(
    "Input file not found: ",
    INPUT_FILE,
    "\nCurrent working directory: ",
    getwd()
  )
}

d <- readr::read_csv(
  INPUT_FILE,
  show_col_types = FALSE
) %>%
  mutate(
    record_id = row_number()
  )

required_columns <- c(
  "site_name",
  "taphonomic_context",
  "Mean_part_abundance",
  "mnau_heip_evenness"
)

missing_columns <- setdiff(
  required_columns,
  names(d)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

C_VAR <- "Mean_part_abundance"
E_VAR <- "mnau_heip_evenness"
VARS <- c(C_VAR, E_VAR)


# ------------------------------------------------------------------------------
# 2. Archaeological labels
# ------------------------------------------------------------------------------

collapse_class4 <- function(x) {
  dplyr::case_when(
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
  mutate(
    class4 = collapse_class4(taphonomic_context)
  )

reference <- d %>%
  filter(!is.na(class4)) %>%
  mutate(
    task_label = if_else(
      class4 == "Funerary Context",
      "Funerary Context",
      "Rest"
    ),

    task_label = factor(
      task_label,
      levels = c(
        "Rest",
        "Funerary Context"
      )
    )
  )

unknown_all <- d %>%
  filter(
    taphonomic_context %in% c(
      "Unknown",
      "Unkown"
    )
  )

missing_requested_unknown <- setdiff(
  RESTRICTED_UNKNOWN_SITES,
  unique(unknown_all$site_name)
)

if (length(missing_requested_unknown) > 0) {
  warning(
    "These requested Unknown sites were not found exactly: ",
    paste(missing_requested_unknown, collapse = ", ")
  )
}

cat("\nReference rows:", nrow(reference), "\n")
cat(
  "Unique reference sites:",
  n_distinct(reference$site_name),
  "\n"
)

cat("\nReference labels:\n")
print(
  reference %>%
    count(task_label)
)

cat(
  "\nUnknown/Unkown rows:",
  nrow(unknown_all),
  "\n"
)


# ------------------------------------------------------------------------------
# 3. Grouped stratified 3-fold x 5-repeat map
# ------------------------------------------------------------------------------

make_grouped_stratified_folds <- function(
  dat,
  v = 3L,
  repeats = 5L,
  seed = 42L
) {

  if (
    any(
      is.na(dat$site_name) |
        dat$site_name == ""
    )
  ) {
    stop(
      "site_name contains missing or blank values."
    )
  }

  group_labels <- dat %>%
    distinct(
      site_name,
      task_label
    )

  inconsistent <- group_labels %>%
    count(
      site_name,
      name = "n_labels"
    ) %>%
    filter(
      n_labels > 1
    )

  if (nrow(inconsistent) > 0) {
    stop(
      "At least one site_name has multiple F-vs-NF labels: ",
      paste(
        inconsistent$site_name,
        collapse = ", "
      )
    )
  }

  maps <- vector(
    "list",
    repeats
  )

  for (r in seq_len(repeats)) {

    set.seed(
      seed + r * 1009L
    )

    maps[[r]] <- group_labels %>%
      group_by(task_label) %>%
      group_modify(
        ~ {
          shuffled <- sample(
            .x$site_name,
            length(.x$site_name),
            replace = FALSE
          )

          tibble(
            site_name = shuffled,
            fold = rep(
              seq_len(v),
              length.out = length(shuffled)
            )
          )
        }
      ) %>%
      ungroup() %>%
      mutate(
        repeat_id = r
      )
  }

  bind_rows(maps)
}

fold_map <- make_grouped_stratified_folds(
  reference,
  v = CV_FOLDS,
  repeats = CV_REPEATS,
  seed = SEED
)

write_csv(
  fold_map,
  file.path(
    OUT_DIR,
    "grouped_cv_fold_map_C_Heip.csv"
  )
)


# ------------------------------------------------------------------------------
# 4. Fold-specific preprocessing
# ------------------------------------------------------------------------------

fit_preprocessor <- function(
  train,
  vars
) {

  medians <- vapply(
    vars,
    function(v) {
      m <- suppressWarnings(
        median(
          train[[v]],
          na.rm = TRUE
        )
      )

      if (is.finite(m)) {
        m
      } else {
        0
      }
    },
    numeric(1)
  )

  impute_raw <- function(dat) {

    x <- as.data.frame(
      dat[, vars, drop = FALSE]
    )

    for (v in vars) {

      bad <- is.na(x[[v]]) |
        !is.finite(x[[v]])

      x[[v]][bad] <- medians[[v]]
    }

    as.matrix(x)
  }

  x_train_raw <- impute_raw(train)

  means <- colMeans(
    x_train_raw
  )

  sds <- apply(
    x_train_raw,
    2,
    sd
  )

  if (
    any(
      !is.finite(sds) |
        sds <= 0
    )
  ) {
    stop(
      "C or Heip has zero variance in a training fold."
    )
  }

  transform_scaled <- function(dat) {

    x <- impute_raw(dat)

    sweep(
      sweep(
        x,
        2,
        means,
        "-"
      ),
      2,
      sds,
      "/"
    )
  }

  list(
    medians = medians,
    means = means,
    sds = sds,
    train_scaled = transform_scaled(train),
    transform_scaled = transform_scaled
  )
}


# ------------------------------------------------------------------------------
# 5. Random oversampling inside training folds only
# ------------------------------------------------------------------------------

balance_indices <- function(
  y,
  seed
) {

  y <- as.character(y)

  classes <- sort(
    unique(y)
  )

  counts <- table(y)

  target_n <- max(counts)

  set.seed(seed)

  idx <- unlist(
    lapply(
      classes,
      function(cl) {

        ii <- which(
          y == cl
        )

        if (
          length(ii) < target_n
        ) {

          c(
            ii,
            sample(
              ii,
              target_n - length(ii),
              replace = TRUE
            )
          )

        } else {

          ii
        }
      }
    ),
    use.names = FALSE
  )

  sample(
    idx,
    length(idx),
    replace = FALSE
  )
}


# ------------------------------------------------------------------------------
# 6. Models
# ------------------------------------------------------------------------------

CLASSES <- c(
  "Rest",
  "Funerary Context"
)

MODELS <- c(
  "Logistic regression",
  "RBF-SVM"
)

align_probability_matrix <- function(
  P,
  classes = CLASSES
) {

  P <- as.matrix(P)

  if (
    is.null(
      colnames(P)
    )
  ) {

    if (
      ncol(P) != length(classes)
    ) {
      stop(
        "Probability matrix has an unexpected number of columns."
      )
    }

    colnames(P) <- classes
  }

  missing <- setdiff(
    classes,
    colnames(P)
  )

  if (
    length(missing) > 0
  ) {
    stop(
      "Probability output is missing class(es): ",
      paste(
        missing,
        collapse = ", "
      )
    )
  }

  P <- P[
    ,
    classes,
    drop = FALSE
  ]

  P[P < 0] <- 0

  rs <- rowSums(P)

  rs[
    !is.finite(rs) |
      rs <= 0
  ] <- 1

  P / rs
}


fit_predict_logistic <- function(
  X_train,
  y_train,
  X_new
) {

  fit <- glmnet::glmnet(
    x = X_train,
    y = factor(
      y_train,
      levels = CLASSES
    ),
    family = "binomial",
    alpha = 0,
    lambda = LOGISTIC_LAMBDA,
    standardize = FALSE
  )

  p_f <- as.numeric(
    predict(
      fit,
      newx = X_new,
      type = "response",
      s = LOGISTIC_LAMBDA
    )
  )

  P <- cbind(
    Rest = 1 - p_f,
    `Funerary Context` = p_f
  )

  align_probability_matrix(P)
}


fit_predict_svm <- function(
  X_train,
  y_train,
  X_new
) {

  fit <- e1071::svm(
    x = X_train,
    y = factor(
      y_train,
      levels = CLASSES
    ),
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

  P <- attr(
    pred,
    "probabilities"
  )

  align_probability_matrix(P)
}


fit_predict_model <- function(
  model_name,
  X_train,
  y_train,
  X_new
) {

  if (
    model_name == "Logistic regression"
  ) {
    return(
      fit_predict_logistic(
        X_train,
        y_train,
        X_new
      )
    )
  }

  if (
    model_name == "RBF-SVM"
  ) {
    return(
      fit_predict_svm(
        X_train,
        y_train,
        X_new
      )
    )
  }

  stop(
    "Unknown model: ",
    model_name
  )
}


# ------------------------------------------------------------------------------
# 7. Metrics
# ------------------------------------------------------------------------------

binary_auc <- function(
  y01,
  score
) {

  ok <- is.finite(score) &
    !is.na(y01)

  y01 <- y01[ok]
  score <- score[ok]

  n1 <- sum(
    y01 == 1
  )

  n0 <- sum(
    y01 == 0
  )

  if (
    n1 == 0 |
      n0 == 0
  ) {
    return(
      NA_real_
    )
  }

  r <- rank(
    score,
    ties.method = "average"
  )

  (
    sum(
      r[y01 == 1]
    ) -
      n1 * (n1 + 1) / 2
  ) /
    (
      n1 * n0
    )
}


calc_binary_metrics <- function(
  truth,
  p_f
) {

  truth <- as.character(truth)

  truth_f <- truth ==
    "Funerary Context"

  pred_f <- p_f >=
    CLASS_THRESHOLD

  TP <- sum(
    truth_f &
      pred_f
  )

  FN <- sum(
    truth_f &
      !pred_f
  )

  TN <- sum(
    !truth_f &
      !pred_f
  )

  FP <- sum(
    !truth_f &
      pred_f
  )

  sensitivity <- if (
    TP + FN > 0
  ) {
    TP / (TP + FN)
  } else {
    NA_real_
  }

  specificity <- if (
    TN + FP > 0
  ) {
    TN / (TN + FP)
  } else {
    NA_real_
  }

  precision <- if (
    TP + FP > 0
  ) {
    TP / (TP + FP)
  } else {
    NA_real_
  }

  recall <- sensitivity

  tibble(
    balanced_accuracy = mean(
      c(
        sensitivity,
        specificity
      ),
      na.rm = TRUE
    ),

    roc_auc = binary_auc(
      as.integer(truth_f),
      p_f
    ),

    precision = precision,

    recall = recall,

    specificity = specificity,

    brier = mean(
      (
        p_f -
          as.integer(truth_f)
      ) ^ 2,
      na.rm = TRUE
    )
  )
}


# ------------------------------------------------------------------------------
# 8. Build the C x Heip prediction grid
# ------------------------------------------------------------------------------

restricted_unknown_raw <- unknown_all %>%
  filter(
    site_name %in%
      RESTRICTED_UNKNOWN_SITES
  )

coords_for_limits <- bind_rows(
  reference %>%
    select(
      all_of(VARS)
    ),

  restricted_unknown_raw %>%
    select(
      all_of(VARS)
    )
)

x_obs <- coords_for_limits[[C_VAR]]
y_obs <- coords_for_limits[[E_VAR]]

x_obs <- x_obs[
  is.finite(x_obs)
]

y_obs <- y_obs[
  is.finite(y_obs)
]

expand_range <- function(
  z,
  fraction = 0.04
) {

  rr <- range(
    z,
    na.rm = TRUE
  )

  span <- diff(rr)

  if (
    !is.finite(span) |
      span == 0
  ) {
    span <- 1
  }

  c(
    rr[1] - fraction * span,
    rr[2] + fraction * span
  )
}

xlim_raw <- expand_range(x_obs)
ylim_raw <- expand_range(y_obs)

if (
  all(
    x_obs >= 0 &
      x_obs <= 1
  )
) {

  xlim_raw <- c(
    max(
      0,
      xlim_raw[1]
    ),
    min(
      1,
      xlim_raw[2]
    )
  )
}

if (
  all(
    y_obs >= 0 &
      y_obs <= 1
  )
) {

  ylim_raw <- c(
    max(
      0,
      ylim_raw[1]
    ),
    min(
      1,
      ylim_raw[2]
    )
  )
}

grid_raw <- tidyr::expand_grid(
  Mean_part_abundance = seq(
    xlim_raw[1],
    xlim_raw[2],
    length.out = GRID_N
  ),

  mnau_heip_evenness = seq(
    ylim_raw[1],
    ylim_raw[2],
    length.out = GRID_N
  )
)


# ------------------------------------------------------------------------------
# 9. Repeated grouped CV
# ------------------------------------------------------------------------------

oof_parts <- list()
unknown_parts <- list()
grid_parts <- list()
fold_metric_parts <- list()

ii_oof <- 1L
ii_unknown <- 1L
ii_grid <- 1L
ii_metric <- 1L

for (
  r in seq_len(
    CV_REPEATS
  )
) {

  cat(
    "\nRepeat",
    r,
    "of",
    CV_REPEATS,
    "\n"
  )

  one_repeat <- fold_map %>%
    filter(
      repeat_id == r
    )

  for (
    f in seq_len(
      CV_FOLDS
    )
  ) {

    cat(
      "  Fold",
      f,
      "of",
      CV_FOLDS,
      "\n"
    )

    valid_sites <- one_repeat %>%
      filter(
        fold == f
      ) %>%
      pull(
        site_name
      )

    train <- reference %>%
      filter(
        !site_name %in%
          valid_sites
      )

    valid <- reference %>%
      filter(
        site_name %in%
          valid_sites
      )

    if (
      any(
        train$site_name %in%
          valid$site_name
      )
    ) {
      stop(
        "Grouped-CV leakage detected."
      )
    }

    pp <- fit_preprocessor(
      train,
      VARS
    )

    X_train <- pp$train_scaled

    X_valid <- pp$transform_scaled(
      valid
    )

    X_unknown <- pp$transform_scaled(
      unknown_all
    )

    X_grid <- pp$transform_scaled(
      grid_raw
    )

    y_train <- as.character(
      train$task_label
    )

    bal_idx <- balance_indices(
      y_train,
      seed =
        SEED +
        r * 1000L +
        f * 100L
    )

    X_train_bal <- X_train[
      bal_idx,
      ,
      drop = FALSE
    ]

    y_train_bal <- y_train[
      bal_idx
    ]

    for (
      model_name in MODELS
    ) {

      P_valid <- fit_predict_model(
        model_name,
        X_train_bal,
        y_train_bal,
        X_valid
      )

      P_unknown <- fit_predict_model(
        model_name,
        X_train_bal,
        y_train_bal,
        X_unknown
      )

      P_grid <- fit_predict_model(
        model_name,
        X_train_bal,
        y_train_bal,
        X_grid
      )

      p_valid_f <- P_valid[
        ,
        "Funerary Context"
      ]

      p_unknown_f <- P_unknown[
        ,
        "Funerary Context"
      ]

      p_grid_f <- P_grid[
        ,
        "Funerary Context"
      ]

      oof_parts[[ii_oof]] <- valid %>%
        transmute(
          model = model_name,
          repeat_id = r,
          fold = f,
          record_id,
          site_name,
          truth = as.character(
            task_label
          ),
          C = .data[[C_VAR]],
          Heip = .data[[E_VAR]],
          P_F = p_valid_f
        )

      ii_oof <- ii_oof + 1L


      unknown_parts[[ii_unknown]] <- unknown_all %>%
        transmute(
          model = model_name,
          repeat_id = r,
          fold = f,
          record_id,
          site_name,
          C = .data[[C_VAR]],
          Heip = .data[[E_VAR]],
          P_F = p_unknown_f
        )

      ii_unknown <- ii_unknown + 1L


      grid_parts[[ii_grid]] <- grid_raw %>%
        transmute(
          model = model_name,
          repeat_id = r,
          fold = f,
          C = .data[[C_VAR]],
          Heip = .data[[E_VAR]],
          P_F = p_grid_f
        )

      ii_grid <- ii_grid + 1L


      fold_metric_parts[[ii_metric]] <-
        calc_binary_metrics(
          truth = valid$task_label,
          p_f = p_valid_f
        ) %>%
        mutate(
          model = model_name,
          repeat_id = r,
          fold = f,
          .before = 1
        )

      ii_metric <- ii_metric + 1L
    }
  }
}

oof_long <- bind_rows(
  oof_parts
)

unknown_long <- bind_rows(
  unknown_parts
)

grid_long <- bind_rows(
  grid_parts
)

fold_metrics <- bind_rows(
  fold_metric_parts
)

write_csv(
  oof_long,
  file.path(
    OUT_DIR,
    "C_Heip_OOF_predictions_long.csv"
  )
)

write_csv(
  unknown_long,
  file.path(
    OUT_DIR,
    "C_Heip_unknown_predictions_15fits_long.csv"
  )
)

write_csv(
  fold_metrics,
  file.path(
    OUT_DIR,
    "C_Heip_fold_metrics.csv"
  )
)


# ------------------------------------------------------------------------------
# 10. Pooled OOF performance
# ------------------------------------------------------------------------------

oof_ensemble <- oof_long %>%
  group_by(
    model,
    record_id,
    site_name,
    truth
  ) %>%
  summarise(
    P_F = mean(
      P_F,
      na.rm = TRUE
    ),
    n_OOF_predictions = n(),
    .groups = "drop"
  )

pooled_performance <- oof_ensemble %>%
  group_by(model) %>%
  group_modify(
    ~ calc_binary_metrics(
      .x$truth,
      .x$P_F
    )
  ) %>%
  ungroup()

fold_performance_summary <- fold_metrics %>%
  group_by(model) %>%
  summarise(
    n_folds = n(),

    BA_mean = mean(
      balanced_accuracy,
      na.rm = TRUE
    ),

    BA_sd = sd(
      balanced_accuracy,
      na.rm = TRUE
    ),

    AUC_mean = mean(
      roc_auc,
      na.rm = TRUE
    ),

    AUC_sd = sd(
      roc_auc,
      na.rm = TRUE
    ),

    Precision_mean = mean(
      precision,
      na.rm = TRUE
    ),

    Precision_sd = sd(
      precision,
      na.rm = TRUE
    ),

    Recall_mean = mean(
      recall,
      na.rm = TRUE
    ),

    Recall_sd = sd(
      recall,
      na.rm = TRUE
    ),

    Brier_mean = mean(
      brier,
      na.rm = TRUE
    ),

    Brier_sd = sd(
      brier,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

write_csv(
  pooled_performance,
  file.path(
    OUT_DIR,
    "C_Heip_pooled_OOF_performance.csv"
  )
)

write_csv(
  fold_performance_summary,
  file.path(
    OUT_DIR,
    "C_Heip_fold_performance_mean_SD.csv"
  )
)


# ------------------------------------------------------------------------------
# 11. Mean probability surface and DISCRETE zones
# ------------------------------------------------------------------------------

grid_mean <- grid_long %>%
  group_by(
    model,
    C,
    Heip
  ) %>%
  summarise(
    Mean_P_F = mean(
      P_F,
      na.rm = TRUE
    ),

    SD_P_F = sd(
      P_F,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    Probability_zone = cut(
      Mean_P_F,

      breaks = c(
        -Inf,
        0.25,
        0.50,
        0.75,
        Inf
      ),

      labels = c(
        "<25% F",
        "25-50% F",
        "50-75% F",
        ">=75% F"
      ),

      right = FALSE
    )
  )

write_csv(
  grid_mean,
  file.path(
    OUT_DIR,
    "C_Heip_probability_surface_grid.csv"
  )
)

# Report which fixed probability bands actually occur for each model.
zone_occupancy <- grid_mean %>%
  count(
    model,
    Probability_zone,
    name = "n_grid_cells"
  ) %>%
  group_by(model) %>%
  mutate(
    percent_of_grid =
      100 *
      n_grid_cells /
      sum(n_grid_cells)
  ) %>%
  ungroup()

write_csv(
  zone_occupancy,
  file.path(
    OUT_DIR,
    "probability_zone_occupancy.csv"
  )
)


# ------------------------------------------------------------------------------
# 12. Unknown-site mean probabilities
# ------------------------------------------------------------------------------

unknown_site_model <- unknown_long %>%
  group_by(
    model,
    repeat_id,
    fold,
    site_name
  ) %>%
  summarise(
    C = mean(
      C,
      na.rm = TRUE
    ),

    Heip = mean(
      Heip,
      na.rm = TRUE
    ),

    P_F = mean(
      P_F,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  group_by(
    model,
    site_name
  ) %>%
  summarise(
    C = mean(
      C,
      na.rm = TRUE
    ),

    Heip = mean(
      Heip,
      na.rm = TRUE
    ),

    Mean_P_F = mean(
      P_F,
      na.rm = TRUE
    ),

    SD_P_F = sd(
      P_F,
      na.rm = TRUE
    ),

    n_fitted_models = n(),

    .groups = "drop"
  ) %>%
  mutate(
    label_pct = paste0(
      site_name,
      " (",
      sprintf(
        "%.0f%%",
        100 * Mean_P_F
      ),
      ")"
    )
  )

write_csv(
  unknown_site_model,
  file.path(
    OUT_DIR,
    "C_Heip_unknown_site_probabilities_LR_SVM.csv"
  )
)

unknown_restricted_model <- unknown_site_model %>%
  filter(
    site_name %in%
      RESTRICTED_UNKNOWN_SITES
  )

write_csv(
  unknown_restricted_model,
  file.path(
    OUT_DIR,
    "C_Heip_unknown_site_probabilities_LR_SVM_restricted.csv"
  )
)


# ------------------------------------------------------------------------------
# 13. Reference plotting data
# ------------------------------------------------------------------------------

reference_plot <- reference %>%
  transmute(
    site_name,

    C = .data[[C_VAR]],

    Heip = .data[[E_VAR]],

    Reference_class = if_else(
      task_label ==
        "Funerary Context",
      "Reference F",
      "Reference NF"
    )
  ) %>%
  filter(
    is.finite(C),
    is.finite(Heip)
  )


# ------------------------------------------------------------------------------
# 14. Sima trajectory and spline
# ------------------------------------------------------------------------------

# Coordinates are independent of the model, so take one unique set of coordinates.
sima_points <- unknown_site_model %>%
  distinct(
    site_name,
    C,
    Heip
  ) %>%
  filter(
    site_name %in%
      SIMA_TRAJECTORY
  ) %>%
  mutate(
    trajectory_order = match(
      site_name,
      SIMA_TRAJECTORY
    )
  ) %>%
  arrange(
    trajectory_order
  )

missing_sima <- setdiff(
  SIMA_TRAJECTORY,
  sima_points$site_name
)

if (
  length(missing_sima) > 0
) {
  warning(
    "Missing Sima trajectory site(s): ",
    paste(
      missing_sima,
      collapse = ", "
    )
  )
}

build_sima_curve <- function(
  sima_df
) {

  if (
    nrow(sima_df) < 3
  ) {
    return(NULL)
  }

  t_in <- sima_df$trajectory_order

  # Parametric spline:
  # one spline for x=C(t) and one for y=Heip(t).
  #
  # method="natural" gives a stable smooth curve without forcing a polynomial
  # through the entire trajectory.
  sx <- stats::spline(
    x = t_in,
    y = sima_df$C,
    n = 250,
    method = "natural"
  )

  sy <- stats::spline(
    x = t_in,
    y = sima_df$Heip,
    n = 250,
    method = "natural"
  )

  tibble(
    trajectory_t = sx$x,
    C = sx$y,
    Heip = sy$y
  )
}

sima_curve_df <- build_sima_curve(
  sima_points
)

write_csv(
  sima_points,
  file.path(
    OUT_DIR,
    "Sima_trajectory_points.csv"
  )
)

if (
  !is.null(
    sima_curve_df
  )
) {
  write_csv(
    sima_curve_df,
    file.path(
      OUT_DIR,
      "Sima_trajectory_spline_curve.csv"
    )
  )
}


# ------------------------------------------------------------------------------
# 15. Plot-variant helpers
# ------------------------------------------------------------------------------

get_unknown_variant <- function(
  model_name,
  variant
) {

  dat <- unknown_site_model %>%
    filter(
      model == model_name
    )

  if (
    variant == "no_unknowns"
  ) {

    return(
      dat[
        0,
        ,
        drop = FALSE
      ]
    )
  }

  if (
    variant == "sima_only"
  ) {

    return(
      dat %>%
        filter(
          site_name %in%
            SIMA_TRAJECTORY
        )
    )
  }

  if (
    variant == "restricted_with_sima"
  ) {

    return(
      dat %>%
        filter(
          site_name %in%
            RESTRICTED_UNKNOWN_SITES
        )
    )
  }

  stop(
    "Unknown plot variant: ",
    variant
  )
}


# ------------------------------------------------------------------------------
# 16. Plot function
# ------------------------------------------------------------------------------

make_probability_plot <- function(
  model_name,
  palette_name,
  variant = c(
    "no_unknowns",
    "sima_only",
    "restricted_with_sima"
  )
) {

  variant <- match.arg(
    variant
  )

  palette_values <- PALETTES[[palette_name]]

  surf <- grid_mean %>%
    filter(
      model == model_name
    )

  unk <- get_unknown_variant(
    model_name,
    variant
  )

  pooled_row <- pooled_performance %>%
    filter(
      model == model_name
    )

  ba_txt <- if (
    nrow(pooled_row) == 1
  ) {
    sprintf(
      "%.3f",
      pooled_row$balanced_accuracy
    )
  } else {
    "NA"
  }

  auc_txt <- if (
    nrow(pooled_row) == 1
  ) {
    sprintf(
      "%.3f",
      pooled_row$roc_auc
    )
  } else {
    "NA"
  }


  # --------------------------------------------------------------------------
  # DISCRETE probability background.
  # --------------------------------------------------------------------------

  p <- ggplot() +

    geom_raster(
      data = surf,
      aes(
        x = C,
        y = Heip,
        fill = Probability_zone
      ),
      interpolate = FALSE
    ) +

    scale_fill_manual(
      values = palette_values,

      breaks = c(
        "<25% F",
        "25-50% F",
        "50-75% F",
        ">=75% F"
      ),

      drop = FALSE,

      name = "Mean P(F)"
    ) +


    # ------------------------------------------------------------------------
    # Exact 25%, 50%, and 75% contours.
    # ------------------------------------------------------------------------

    geom_contour(
      data = surf,
      aes(
        x = C,
        y = Heip,
        z = Mean_P_F
      ),
      breaks = 0.25,
      colour = CONTOUR_25_COLOUR,
      linewidth = 0.65
    ) +

    geom_contour(
      data = surf,
      aes(
        x = C,
        y = Heip,
        z = Mean_P_F
      ),
      breaks = 0.50,
      colour = CONTOUR_50_COLOUR,
      linewidth = 1.0
    ) +

    geom_contour(
      data = surf,
      aes(
        x = C,
        y = Heip,
        z = Mean_P_F
      ),
      breaks = 0.75,
      colour = CONTOUR_75_COLOUR,
      linewidth = 0.65
    ) +


    # ------------------------------------------------------------------------
    # Reference NF: dark-outlined blue circles.
    # ------------------------------------------------------------------------

    geom_point(
      data = reference_plot %>%
        filter(
          Reference_class ==
            "Reference NF"
        ),

      aes(
        x = C,
        y = Heip
      ),

      shape = 21,
      size = 2.8,
      stroke = 0.65,
      colour = REFERENCE_NF_COLOUR,
      fill = REFERENCE_NF_FILL,
      alpha = 0.90
    ) +


    # ------------------------------------------------------------------------
    # Reference F: dark-outlined orange triangles.
    # ------------------------------------------------------------------------

    geom_point(
      data = reference_plot %>%
        filter(
          Reference_class ==
            "Reference F"
        ),

      aes(
        x = C,
        y = Heip
      ),

      shape = 24,
      size = 3.0,
      stroke = 0.65,
      colour = REFERENCE_F_COLOUR,
      fill = REFERENCE_F_FILL,
      alpha = 0.92
    )


  # --------------------------------------------------------------------------
  # Sima spline for Sima-containing variants.
  # --------------------------------------------------------------------------

  if (
    variant %in%
      c(
        "sima_only",
        "restricted_with_sima"
      ) &&
      !is.null(
        sima_curve_df
      )
  ) {

    p <- p +

      geom_path(
        data = sima_curve_df,

        aes(
          x = C,
          y = Heip
        ),

        colour = SIMA_CURVE_COLOUR,
        linewidth = 1.15,
        alpha = 0.95
      )
  }


  # --------------------------------------------------------------------------
  # Unknown sites.
  # --------------------------------------------------------------------------

  if (
    nrow(unk) > 0
  ) {

    p <- p +

      geom_point(
        data = unk,

        aes(
          x = C,
          y = Heip
        ),

        shape = 4,
        size = 4.0,
        stroke = 1.35,
        colour = UNKNOWN_COLOUR
      ) +

      ggrepel::geom_text_repel(
        data = unk,

        aes(
          x = C,
          y = Heip,
          label = label_pct
        ),

        size = 3.0,
        fontface = "plain",
        colour = "#111111",

        min.segment.length = 0,
        segment.colour = "#555555",
        segment.alpha = 0.65,

        box.padding = 0.28,
        point.padding = 0.18,

        max.overlaps = Inf,

        seed = SEED
      )
  }


  # --------------------------------------------------------------------------
  # User-requested titles / labels for ALL images.
  # --------------------------------------------------------------------------

  p +
    coord_cartesian(
      xlim = xlim_raw,
      ylim = ylim_raw,
      expand = FALSE
    ) +

    labs(
      title = paste0(
        "F vs NF probability surface - ",
        model_name
      ),

      subtitle = paste0(
        "Accuracy = ",
        ba_txt,
        ", AUROC = ",
        auc_txt
      ),

      x = "Mean skeletal-part abundance C",

      y = "Heip anatomical evenness",

      caption = paste0(
        "Reference NF = circles; Reference F = triangles."
      )
    ) +

    theme_bw(
      base_size = 12
    ) +

    theme(
      plot.title = element_text(
        face = "bold"
      ),

      plot.subtitle = element_text(
        size = 10.0
      ),

      legend.position = "right",

      panel.grid.minor = element_blank(),

      panel.grid.major = element_line(
        colour = "#E5E5E5",
        linewidth = 0.30
      ),

      legend.key.height = grid::unit(
        0.8,
        "cm"
      ),

      legend.key.width = grid::unit(
        1.0,
        "cm"
      )
    )
}


# ------------------------------------------------------------------------------
# 17. Save all requested discrete-colour versions
# ------------------------------------------------------------------------------

plot_specs <- tidyr::expand_grid(

  model_name = MODELS,

  palette_name = names(
    PALETTES
  ),

  variant = c(
    "no_unknowns",
    "sima_only",
    "restricted_with_sima"
  )
)

output_manifest <- vector(
  "list",
  nrow(plot_specs)
)

for (
  i in seq_len(
    nrow(plot_specs)
  )
) {

  spec <- plot_specs[
    i,
    ,
    drop = FALSE
  ]

  p <- make_probability_plot(
    model_name = spec$model_name,
    palette_name = spec$palette_name,
    variant = spec$variant
  )

  model_stub <- ifelse(
    spec$model_name ==
      "Logistic regression",

    "Logistic_regression",

    "RBF_SVM"
  )

  file_stub <- paste(
    "F_vs_NF",
    "C_Heip",
    model_stub,
    spec$variant,
    spec$palette_name,
    sep = "_"
  )

  png_file <- file.path(
    OUT_DIR,
    paste0(
      file_stub,
      ".png"
    )
  )

  pdf_file <- file.path(
    OUT_DIR,
    paste0(
      file_stub,
      ".pdf"
    )
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = 13,
    height = 8.5,
    dpi = 320
  )

  ggsave(
    filename = pdf_file,
    plot = p,
    width = 13,
    height = 8.5
  )

  output_manifest[[i]] <- tibble(
    model = spec$model_name,
    variant = spec$variant,
    palette = spec$palette_name,
    png_file = png_file,
    pdf_file = pdf_file
  )
}

output_manifest <- bind_rows(
  output_manifest
)

write_csv(
  output_manifest,
  file.path(
    OUT_DIR,
    "plot_output_manifest.csv"
  )
)


# ------------------------------------------------------------------------------
# 18. Restricted Unknown summary
# ------------------------------------------------------------------------------

restricted_summary <- unknown_restricted_model %>%
  select(
    model,
    site_name,
    C,
    Heip,
    Mean_P_F,
    SD_P_F,
    n_fitted_models
  ) %>%
  arrange(
    model,
    desc(
      Mean_P_F
    )
  )

write_csv(
  restricted_summary,
  file.path(
    OUT_DIR,
    "restricted_unknown_site_summary_LR_SVM.csv"
  )
)


# ------------------------------------------------------------------------------
# 19. Console summary
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n"
)

cat(
  "FINISHED\n"
)

cat(
  "============================================================\n\n"
)

cat(
  "Pooled out-of-fold performance:\n"
)

print(
  pooled_performance
)

cat(
  "\nFold performance mean +/- SD:\n"
)

print(
  fold_performance_summary
)

cat(
  "\nProbability-zone occupancy:\n"
)

print(
  zone_occupancy
)

cat(
  "\nRestricted Unknown-site probabilities:\n"
)

print(
  restricted_summary
)

cat(
  "\nOutput directory:\n  ",
  OUT_DIR,
  "\n",
  sep = ""
)

cat(
  "\n18 discrete-colour figures are produced:\n"
)

cat(
  "  2 models x 3 palettes x 3 plot variants\n"
)

cat(
  "\nPlot variants:\n"
)

cat(
  "  no_unknowns\n"
)

cat(
  "  sima_only\n"
)

cat(
  "  restricted_with_sima\n"
)

cat(
  "\nSima spline order:\n"
)

cat(
  "  89 -> 90 -> 95 -> 15 -> 24\n"
)

cat(
  "\nNo gradient-colour figures are produced.\n"
)

cat(
  "\nIf Logistic regression contains fewer than four background colours,\n"
)

cat(
  "check probability_zone_occupancy.csv: this means its predictions do not\n"
)

cat(
  "enter all four fixed 25/50/75 probability intervals.\n"
)

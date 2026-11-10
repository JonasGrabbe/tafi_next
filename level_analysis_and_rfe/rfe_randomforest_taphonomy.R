# =====================================================================
#  Recursive Feature Elimination (RFE) with randomForest
#  Target: taphonomic_context
#
#  Procedure (as specified):
#    * train RF on ALL features, score accuracy
#    * while > 40 features remain  -> drop the 40 least important
#    * once <= 40 (but > 10)       -> drop the 10 least important
#    * once <= 10                  -> train once more, stop & plot
#    * accuracy is the selection / reporting metric
#    * write a step-by-step ranked CSV of eliminated features
# =====================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(readr)
library(randomForest)
library(rlang)

# ---------------------------------------------------------------------
# 0. Metadata columns to remove (everything except the target).
#    taphonomic_context is intentionally NOT here -- it is the label.
# ---------------------------------------------------------------------
meta_cols <- c(
  "site_name", "ind_name", "site_id", "individual_id", "layer_name",
  "count_type", "geological_period", "cultural_period", "culture_period_2",
  "start_date_cal_bp", "end_date_cal_bp", "mis_stage",
  "karstic_system", "karstic_system_2", "open_air_site", "open_air_site_2",
  "region", "country", "mni", "age_category",
  "n_infant", "n_juvenil", "n_subadult", "n_adult", "gender",
  "funerary_context", "discovery_context", "degree_of_excavation",
  "hominin_species", "observation", "reference"
)

# ---------------------------------------------------------------------
# 1. Your task recoder (unchanged) — builds the labelled data frame.
# ---------------------------------------------------------------------
task_recode_data <- function(raw_df,
                             task_name,
                             target_variable = "taphonomic_context") {
  tvsym <- rlang::sym(target_variable)
  if (task_name == "4class") {
    class_levels <- c("Cannibalism", "Carnivore",
                      "Funerary Context", "Geological transport")
    df <- raw_df %>%
      dplyr::filter(.data[[target_variable]] %in% class_levels) %>%
      dplyr::mutate(!!tvsym := factor(.data[[target_variable]], levels = class_levels))
    use_downsample <- FALSE
  } else if (task_name == "cannibalism_vs_fc") {
    class_levels <- c("Cannibalism", "Funerary Context")
    df <- raw_df %>%
      dplyr::filter(.data[[target_variable]] %in% class_levels) %>%
      dplyr::mutate(!!tvsym := factor(.data[[target_variable]], levels = class_levels))
    use_downsample <- TRUE
  } else if (task_name == "fc_vs_rest") {
    rest_classes <- c("Cannibalism", "Carnivore", "Geological transport")
    class_levels <- c("Funerary Context", "rest")
    df <- raw_df %>%
      dplyr::mutate(!!tvsym := dplyr::case_when(
        .data[[target_variable]] %in% rest_classes      ~ "rest",
        .data[[target_variable]] == "Funerary Context"  ~ "Funerary Context",
        TRUE                                            ~ NA_character_)) %>%
      dplyr::filter(!is.na(.data[[target_variable]])) %>%
      dplyr::mutate(!!tvsym := factor(.data[[target_variable]], levels = class_levels))
    use_downsample <- TRUE
  } else {
    stop("Unknown task_name: ", task_name)
  }
  message("Class distribution (", task_name, "):")
  print(table(df[[target_variable]], useNA = "ifany"))
  missing_classes <- setdiff(levels(df[[target_variable]]), unique(df[[target_variable]]))
  if (length(missing_classes))
    warning("Some expected classes missing: ", paste(missing_classes, collapse = ", "))
  list(df = df, class_levels = levels(df[[target_variable]]),
       use_downsample = use_downsample)
}

# ---------------------------------------------------------------------
# 2. Drop meta + non-usable columns; keep target + clean numeric features.
#    - removes metadata and any list-columns / character columns
#    - logicals -> numeric; Inf/NaN -> NA; NA -> column median
#    - drops zero-variance / all-NA predictors
# ---------------------------------------------------------------------
prep_features <- function(df,
                          target_variable = "taphonomic_context",
                          drop_meta       = meta_cols) {

  y <- df[[target_variable]]
  if (!is.factor(y)) y <- factor(y)

  X <- df %>%
    dplyr::select(-dplyr::any_of(setdiff(drop_meta, target_variable))) %>%
    dplyr::select(-dplyr::any_of(target_variable)) %>%
    dplyr::select(where(~ !is.list(.x))) %>%          # drop list-columns
    dplyr::select(where(~ is.numeric(.x) || is.logical(.x))) %>%
    dplyr::mutate(dplyr::across(where(is.logical), as.numeric)) %>%
    dplyr::mutate(dplyr::across(everything(),
                                ~ ifelse(is.nan(.x) | is.infinite(.x), NA_real_, .x)))

  # drop all-NA or zero-variance predictors
  keep <- vapply(X, function(col) {
    v <- col[!is.na(col)]
    length(v) > 0 && stats::var(v) > 0
  }, logical(1))
  X <- X[, keep, drop = FALSE]

  # median-impute remaining NAs
  X <- X %>% dplyr::mutate(dplyr::across(everything(),
              ~ ifelse(is.na(.x), stats::median(.x, na.rm = TRUE), .x)))

  message("Predictors kept after cleaning: ", ncol(X))
  list(X = X, y = y)
}

# ---------------------------------------------------------------------
# 3. Optional k-fold CV accuracy (used only if cv_folds is set).
#    If cv_folds is NULL, OOB accuracy is used instead (fast default).
# ---------------------------------------------------------------------
cv_accuracy <- function(X, y, ntree, folds, seed) {
  set.seed(seed)
  fold_id <- sample(rep_len(seq_len(folds), length(y)))
  accs <- vapply(seq_len(folds), function(f) {
    tr <- fold_id != f; te <- !tr
    if (length(unique(y[tr])) < 2) return(NA_real_)
    m  <- randomForest::randomForest(x = X[tr, , drop = FALSE], y = y[tr], ntree = ntree)
    pr <- predict(m, X[te, , drop = FALSE])
    mean(pr == y[te])
  }, numeric(1))
  mean(accs, na.rm = TRUE)
}

# ---------------------------------------------------------------------
# 4. The RFE loop.
# ---------------------------------------------------------------------
rfe_rf <- function(X, y,
                   big_step        = 40,
                   small_threshold = 40,
                   small_step      = 10,
                   min_features    = 10,
                   ntree           = 500,
                   cv_folds        = NULL,   # NULL -> OOB accuracy
                   guard_overshoot = TRUE,   # avoid skipping the 10-step phase
                   seed            = 42) {

  current    <- names(X)
  history    <- list()   # accuracy at each model size
  eliminated <- list()   # features removed at each step
  step       <- 0L

  repeat {
    step <- step + 1L
    Xc   <- X[, current, drop = FALSE]

    set.seed(seed)
    rf <- randomForest::randomForest(x = Xc, y = y, ntree = ntree)

    # accuracy (OOB by default, else CV)
    acc <- if (is.null(cv_folds)) {
      cm <- rf$confusion[, setdiff(colnames(rf$confusion), "class.error"), drop = FALSE]
      sum(diag(cm)) / sum(cm)
    } else {
      cv_accuracy(Xc, y, ntree = ntree, folds = cv_folds, seed = seed)
    }

    imp <- randomForest::importance(rf, type = 2)          # MeanDecreaseGini
    imp_df <- tibble::tibble(feature    = rownames(imp),
                             importance = as.numeric(imp[, 1])) %>%
              dplyr::arrange(importance)                    # least important first

    history[[step]] <- tibble::tibble(step = step,
                                      n_features = length(current),
                                      accuracy   = acc)
    n <- length(current)
    message(sprintf("step %2d | %3d features | accuracy = %.4f", step, n, acc))

    if (n <= min_features) break                            # final model trained; stop

    k <- if (n > small_threshold) big_step else small_step
    if (guard_overshoot && n > small_threshold)
      k <- min(k, n - small_threshold)                      # land at/above 40, keep 10-phase
    k <- min(k, n - 1L)                                     # never drop every feature

    drop_feats <- imp_df$feature[seq_len(k)]
    eliminated[[step]] <- tibble::tibble(
      step              = step,
      n_features_before = n,
      accuracy_before   = acc,
      feature           = drop_feats,
      importance        = imp_df$importance[seq_len(k)]
    )
    current <- setdiff(current, drop_feats)
  }

  history    <- dplyr::bind_rows(history)
  eliminated <- dplyr::bind_rows(eliminated)
  retained   <- current

  # global elimination ranking: order of removal = ascending importance.
  # rank 1 = eliminated first (least useful); retained features rank last.
  elim_ranked <- eliminated %>%
    dplyr::mutate(status = "eliminated") %>%
    dplyr::arrange(step, importance) %>%
    dplyr::mutate(elimination_rank = dplyr::row_number())

  retained_tbl <- tibble::tibble(
    step              = NA_integer_,
    n_features_before = length(retained),
    accuracy_before   = dplyr::last(history$accuracy),
    feature           = retained,
    importance        = NA_real_,
    status            = "retained",
    elimination_rank  = nrow(elim_ranked) + seq_along(retained)
  )

  ranking <- dplyr::bind_rows(elim_ranked, retained_tbl)

  best <- history %>% dplyr::slice_max(accuracy, n = 1, with_ties = FALSE)

  list(history = history,
       eliminated = eliminated,
       ranking = ranking,
       retained = retained,
       best = best)
}

# ---------------------------------------------------------------------
# 5. Driver: run end-to-end, plot, and write CSV.
#     `tapho_db` is your fully-featured data frame from the loader.
# ---------------------------------------------------------------------
run_rfe <- function(tapho_db,
                    task_name = "4class",
                    out_csv   = "feature_elimination_ranking.csv",
                    out_plot  = "feature_elimination_curve.png",
                    cv_folds  = NULL,
                    seed      = 42) {

  rec  <- task_recode_data(tapho_db, task_name = task_name,
                           target_variable = "taphonomic_context")
  prep <- prep_features(rec$df, target_variable = "taphonomic_context")

  res <- rfe_rf(prep$X, prep$y, cv_folds = cv_folds, seed = seed)

  # ---- plot: accuracy vs number of features (elimination read left->right) ----
  p <- ggplot(res$history, aes(x = n_features, y = accuracy)) +
    geom_line(linewidth = 0.7, colour = "grey40") +
    geom_point(size = 2) +
    geom_point(data = res$best, colour = "firebrick", size = 3.5) +
    geom_text(data = res$best,
              aes(label = sprintf("best: %d feats\nacc %.3f", n_features, accuracy)),
              vjust = -0.8, hjust = 0.5, size = 3, colour = "firebrick") +
    scale_x_reverse() +
    labs(title = paste0("RFE (randomForest) — ", task_name),
         subtitle = "Accuracy as features are recursively eliminated",
         x = "Number of features (eliminated left to right)",
         y = "Accuracy") +
    theme_minimal(base_size = 12)

  ggsave(out_plot, p, width = 8, height = 5, dpi = 150)
  print(p)

  # ---- ranked CSV of eliminated (then retained) features ----
  readr::write_csv(res$ranking, out_csv)
  message("Wrote ranking -> ", out_csv)
  message("Wrote plot    -> ", out_plot)
  message(sprintf("Best subset: %d features at accuracy %.4f",
                  res$best$n_features, res$best$accuracy))

  invisible(res)
}

# ---------------------------------------------------------------------
# 6. Example call
# ---------------------------------------------------------------------
# res <- run_rfe(tapho_db, task_name = "4class")           # OOB accuracy (fast)
# res <- run_rfe(tapho_db, task_name = "fc_vs_rest", cv_folds = 5)   # 5-fold CV
# res$best        # best feature-count / accuracy
# res$retained    # the surviving features

# hyperparam_tuning_core.R
# Hyperparameter tuning pipeline for one (task, variant, level, model)

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(themis)
  library(doParallel)
  library(openxlsx)
  library(rlang)
  library(knitr)
})

# You MUST source("hpa_sdk.R") BEFORE calling tune_single_model_combo()
# because we reuse compute_classification_report_from_confmats() if you later
# want confusion-matrix style summaries. We won't compute that in tuning,
# but we keep the convention.

tidymodels_prefer()


# -----------------------------------------------------------
# 1. Skeleton feature metadata, same idea as before
# -----------------------------------------------------------
skeleton_info_for_level <- function(df) {


    skeleton_df <- data.frame(
      element = df$element[[1]],
      count = df$One_Skeleton_count[[1]]
    )


  list(
    feature_names = df$element[[1]],
    skeleton_df   = skeleton_df,
    n_features    = length(df$element[[1]]), 
    level = df$bone.level[[1]]
  )
}


# -----------------------------------------------------------
# 2. Task-specific class recode (4class / cannibalism_vs_fc / fc_vs_rest)
# -----------------------------------------------------------
task_recode_data <- function(raw_df,
                             task_name,
                             target_variable = "taphonomic_context") {

  tvsym <- rlang::sym(target_variable)

  if (task_name == "4class") {

    class_levels <- c(
      "Cannibalism",
      "Carnivore",
      "Funerary Context",
      "Geological transport"
    )

    df <- raw_df %>%
      dplyr::filter(.data[[target_variable]] %in% class_levels) %>%
      dplyr::mutate(
        !!tvsym := factor(.data[[target_variable]],
                          levels = class_levels)
      )

    message("Class distribution (4class):")
    print(table(df[[target_variable]], useNA = "ifany"))

    # downsampling with themis is typically binary; turn it off
    use_downsample <- FALSE

  } else if (task_name == "cannibalism_vs_fc") {

    class_levels <- c("Cannibalism", "Funerary Context")

    df <- raw_df %>%
      dplyr::filter(.data[[target_variable]] %in% class_levels) %>%
      dplyr::mutate(
        !!tvsym := factor(.data[[target_variable]],
                          levels = class_levels)
      )

    message("Class distribution (cannibalism_vs_fc):")
    print(table(df[[target_variable]], useNA = "ifany"))

    use_downsample <- TRUE

  } else if (task_name == "fc_vs_rest") {

    rest_classes <- c("Cannibalism", "Carnivore", "Geological transport")

    class_levels <- c("Funerary Context", "rest")

    df <- raw_df %>%
      dplyr::mutate(
        !!tvsym := dplyr::case_when(
          .data[[target_variable]] %in% rest_classes ~ "rest",
          .data[[target_variable]] == "Funerary Context" ~ "Funerary Context",
          TRUE ~ NA_character_
        )
      ) %>%
      dplyr::filter(!is.na(.data[[target_variable]])) %>%
      dplyr::mutate(
        !!tvsym := factor(.data[[target_variable]],
                          levels = class_levels)
      )

    message("Class distribution (fc_vs_rest):")
    print(table(df[[target_variable]], useNA = "ifany"))

    use_downsample <- TRUE

  } else {
    stop("Unknown task_name: ", task_name)
  }

  # helpful check
  missing_classes <- setdiff(levels(df[[target_variable]]),
                             unique(df[[target_variable]]))
  if (length(missing_classes)) {
    warning("Some expected classes missing: ",
            paste(missing_classes, collapse = ", "))
  }

  list(
    df              = df,
    class_levels    = levels(df[[target_variable]]),
    use_downsample  = use_downsample
  )
}
# -----------------------------------------------------------
# -. ppasst Build recipe for tuning
# -----------------------------------------------------------
library(recipes)
library(dplyr)
library(rlang)
library(themis)

identity_aug <- function(x) x

default_smote_fun <- function(df, target_col, neighbors, over_ratio) {
  themis::smote(
    df,
    var = target_col,
    k = neighbors,
    over_ratio = over_ratio
  )
}

make_smote_train_sum_id <- function(model_df_train,
  raw_df_id,
  raw_df_sum,
  site_id_col,
  target_col,
  mni_col,
  feature_cols,
  augment_fun = my_augment_fun,
  smote_fun   = default_smote_fun,
  neighbors   = 5,
  over_ratio  = 1) {

# Cast ALL three tables to character first so bind_rows never sees a type clash
model_df_train[[site_id_col]] <- as.character(model_df_train[[site_id_col]])  # <-- add
raw_df_id[[site_id_col]]      <- as.character(raw_df_id[[site_id_col]])
raw_df_sum[[site_id_col]]     <- as.character(raw_df_sum[[site_id_col]])

train_site_ids <- unique(model_df_train[[site_id_col]])   # now always character

target_levels            <- levels(factor(model_df_train[[target_col]]))
raw_df_id[[target_col]]  <- factor(raw_df_id[[target_col]],  levels = target_levels)
raw_df_sum[[target_col]] <- factor(raw_df_sum[[target_col]], levels = target_levels)

dplyr::bind_rows(
dplyr::filter(raw_df_id,  .data[[site_id_col]] %in% train_site_ids),
dplyr::filter(raw_df_sum, .data[[site_id_col]] %in% train_site_ids)
) |>
augment_fun() |>
dplyr::select(dplyr::all_of(c(target_col, mni_col, feature_cols))) |>
smote_fun(target_col = target_col, neighbors = neighbors, over_ratio = over_ratio) |>
dplyr::select(-dplyr::all_of(mni_col))
}

# -----------------------------------------------------------
# 3. Build recipe for tuning
site_id_col <- "site_id"
target_col <- "taphonomic_context"
mni_col <- "mni"
# -----------------------------------------------------------
build_recipe_for_task <- function(df,
  target_variable = "taphonomic_context",
  site_name_col   = "site_name",
  raw_df_id       = NULL,
  raw_df_sum      = NULL,
  site_id_col     = "site_id",
  mni_col         = "mni",
  feature_cols    = NULL,
  use_smote       = FALSE,
  use_downsample  = FALSE) {

rec <- recipes::recipe(
as.formula(paste(target_variable, "~ .")),
data = df
) %>%
recipes::update_role(!!rlang::sym(site_name_col), new_role = "ID") %>%
# Scale BEFORE SMOTE so synthetic points are generated in scaled space.
# If you want raw-space SMOTE, move this step after step_meta_smote.
recipes::step_mutate(
dplyr::across(recipes::all_numeric_predictors(), ~ .x / 100)
)

if (use_smote && !is.null(raw_df_id) && !is.null(raw_df_sum)) {
rec <- rec %>%
step_meta_smote(
raw_df_id    = raw_df_id,
raw_df_sum   = raw_df_sum,
site_id_col  = site_id_col,
target_col   = target_variable,
mni_col      = mni_col,
feature_cols = feature_cols,
neighbors    = 5,
over_ratio   = 3,
skip         = TRUE
)
}

#if (use_downsample) {
#rec <- rec %>%
#themis::step_downsample(
#!!rlang::sym(target_variable),
#under_ratio = 1,
#skip        = TRUE,
#seed        = sample.int(10^5, 1)
#)
#}

rec
}
# -----------------------------------------------------------
# --. paasssgt  Model spec + tunable parameters for each classifier
# -----------------------------------------------------------

step_meta_smote <- function(recipe,
  raw_df_id,
  raw_df_sum,
  site_id_col = "site_id",
  target_col,
  mni_col     = "mni",
  feature_cols,
  augment_fun = identity_aug,
  smote_fun   = default_smote_fun,
  neighbors   = 1,
  over_ratio  = 1,
  skip        = TRUE,
  id          = recipes::rand_id("meta_smote")) {

recipes::add_step(recipe,
structure(
list(raw_df_id    = raw_df_id,
raw_df_sum   = raw_df_sum,
site_id_col  = site_id_col,
target_col   = target_col,
mni_col      = mni_col,
feature_cols = feature_cols,
augment_fun  = augment_fun,
smote_fun    = smote_fun,
neighbors    = neighbors,
over_ratio   = over_ratio,
role         = NA,        # required by recipes internals
trained      = FALSE,
skip         = skip,
id           = id),
class = c("step_meta_smote", "step")
)
)
}

prep.step_meta_smote <- function(x, training, info = NULL, ...) {
x$trained <- TRUE
x
}

bake.step_meta_smote <- function(object, new_data, ...) {
augmented <- make_smote_train_sum_id(
model_df_train = new_data,
raw_df_id      = object$raw_df_id,
raw_df_sum     = object$raw_df_sum,
site_id_col    = object$site_id_col,
target_col     = object$target_col,
mni_col        = object$mni_col,
feature_cols   = object$feature_cols,
augment_fun    = object$augment_fun,
smote_fun      = object$smote_fun,
neighbors      = object$neighbors,
over_ratio     = object$over_ratio
)

# Keep original training rows + synthetic rows
# Drop site_id and mni from new_data to match augmented columns
original <- new_data %>%
dplyr::select(dplyr::all_of(names(augmented)))

dplyr::bind_rows(original, augmented) |>
tibble::as_tibble()
}

print.step_meta_smote <- function(x, ...) {
cat("Meta SMOTE augmentation [",
if (x$trained) "trained" else "untrained", "]\n")
invisible(x)
}

# Register the S3 methods explicitly. Without this, dispatch can silently fail
# on parallel (PSOCK) workers, because the methods are only reached via the
# bake()/prep() generics and tune's global-scan never ships them to workers.
.S3method("prep",  "step_meta_smote", prep.step_meta_smote)
.S3method("bake",  "step_meta_smote", bake.step_meta_smote)
.S3method("print", "step_meta_smote", print.step_meta_smote)

# -----------------------------------------------------------
# 4. Model spec + tunable parameters for each classifier
# -----------------------------------------------------------
get_model_spec <- function(model_name) {

  model_name <- match.arg(
    model_name,
    c("SVM_RBF",
      "SVM_Linear",
      "RandomForest",
      "XGBoost",
      "DecisionTree",
      "KNN",
      "MLP")
  )

  if (model_name == "SVM_RBF") {
    spec <- svm_rbf(
      mode       = "classification",
      cost       = tune(),
      rbf_sigma  = tune()
    ) %>%
      set_engine("kernlab")

  } else if (model_name == "SVM_Linear") {
    spec <- svm_linear(
      mode = "classification",
      cost = tune()
    ) %>%
      set_engine("kernlab")

  } else if (model_name == "RandomForest") {
    spec <- rand_forest(
      mode   = "classification",
      mtry   = tune(),
      min_n  = tune(),
      trees  = tune()
    ) %>%
      set_engine("ranger", importance = "permutation")

  } else if (model_name == "XGBoost") {
    spec <- boost_tree(
      mode            = "classification",
      trees           = tune(),
      tree_depth      = tune(),
      learn_rate      = tune(),
      loss_reduction  = tune(),
      min_n           = tune(),
      sample_size     = tune()
    ) %>%
      set_engine("xgboost")

  } else if (model_name == "DecisionTree") {
    spec <- decision_tree(
      mode            = "classification",
      tree_depth      = tune(),
      min_n           = tune(),
      cost_complexity = tune()
    ) %>%
      set_engine("rpart")

  } else if (model_name == "KNN") {
    spec <- nearest_neighbor(
      mode        = "classification",
      neighbors   = tune(),
      weight_func = tune(),
      dist_power  = tune()
    ) %>%
      set_engine("kknn")

  } else if (model_name == "MLP") {
    spec <- mlp(
      mode         = "classification",
      hidden_units = tune(),
      penalty      = tune(),
      epochs       = tune()
    ) %>%
      set_engine("nnet")

  } else {
    stop("Unsupported model_name: ", model_name)
  }

  spec
}


# -----------------------------------------------------------
# 5. The big function: tune_single_model_combo()
# -----------------------------------------------------------
# Arguments:
# - dataset_path: "pmnau_tb_level_3.csv"
# - level:        3
# - data_variant: "pmnau" or "remnau"
# - task_name:    "4class", "cannibalism_vs_fc", or "fc_vs_rest"
# - model_name:   one of get_model_spec() names above
# - grid_size:    how many hyperparam combos to try
# - n_iterations: how many repeats for repeated v-fold CV (3-fold * repeats)
#
# Output:
# creates "results_tuning_<task>_<variant>_L<level>_<model>/" with:
#   - tuning_results_<model>.csv   (all combos ranked best->worst)
#   - tuning_top10_<model>.md      (markdown table of top10)
#   - tuning_<model>.xlsx          (Excel wb with all_results + top10)
#   - tune_object_<model>.rds      (raw tune_grid result)
#   - plus skeleton + processed data snapshots for traceability

tune_single_model_combo <- function(
  dataset_path,
  level,
  sum, 
  data_variant,
  task_name,
  model_name,
  grid_size     = 25,
  n_iterations  = 5,
  target_variable = "taphonomic_context",
  site_name_col   = "site_id",  # "site_id"
  use_parallel    = FALSE
) {

  level_str <- paste(lvl, collapse = "_")


  cat("\n=====================================================\n")
  cat("TUNING START\n",
      "variant =", data_variant, "\n",
      "level   =", level_str, "\n",
      "task    =", task_name, "\n",
      "model   =", model_name, "\n")
  cat("=====================================================\n")

  # -------------------------------------------------------
  # 0. Prep output dir
  # -------------------------------------------------------
  RESULTS_DIR <- paste0(
    "results_tuning_",
    task_name, "_",
    data_variant, "_L",
    level_str, "_",
    model_name
  )
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

  # -------------------------------------------------------
  # 1. Load & recode data for this task
  # -------------------------------------------------------
  
  source("c:/Users/Usuario/tafi_test/tapho_loader_new.R")
  result        <- tapho.loader.wfl(dataset_path, level = level, data_variant = data_variant )
  raw_df_id  <- result$db_id
  raw_df_sum <- result$db_sum


  raw_df_train <- if (sum) raw_df_sum else raw_df_id
  raw_df_smote <- if (sum) raw_df_id  else raw_df_sum

  task_out <- task_recode_data(
    raw_df         = raw_df_train,
    task_name      = task_name,
    target_variable= target_variable
  )
  

  df             <- task_out$df
  class_levels   <- task_out$class_levels
  # smote always 
  use_downsample <- task_out$use_downsample
  # accuracy always
  is_binary      <- (length(class_levels) == 2)

  task_out_smote <- task_recode_data(
    raw_df         = raw_df_smote,
    task_name      = task_name,
    target_variable= target_variable
  )
  smote_df <- task_out_smote$df

  # snapshot processed data
  saveRDS(
    df,
    file = file.path(
      RESULTS_DIR,
      paste0("processed_", data_variant,
             "_L", level_str, "_", task_name, ".rds")
    )
  )

  # -------------------------------------------------------
  # 2. Save skeleton metadata for this level
  # -------------------------------------------------------
  skel_info <- skeleton_info_for_level(df)    
  readr::write_csv(
    skel_info$skeleton_df,
    file.path(
      RESULTS_DIR,
      paste0("feature_reference_L", level_str, "_",
             data_variant, "_", task_name, ".csv")
    )
  )
  cat("Level", level_str, ": documented",
      skel_info$n_features, "features.\n")

  # -------------------------------------------------------
  # 2/3. Keep only ID, target, and feature columns
  # -------------------------------------------------------

  id_cols <- site_name_col

  feature_cols <- setdiff(
    names(df),
    c(
      id_cols,
      target_variable,

      # metadata / non-feature columns
      names(df)[1:37]
    )
  )

  model_df <- df %>%
    dplyr::select(
      dplyr::all_of(id_cols),
      dplyr::all_of(target_variable),
      dplyr::all_of(feature_cols)
    )
    
    
  # -------------------------------------------------------
  # 3. CV folds
  # -------------------------------------------------------

  set.seed(42)

  cv_folds <- rsample::vfold_cv(
    model_df,
    v       = 3,
    repeats = n_iterations,
    strata  = !!rlang::sym(target_variable)
  )

  saveRDS(
    cv_folds,
    file.path(RESULTS_DIR, "cv_folds_tuning.rds")
  )

  cat("CV folds:", nrow(cv_folds), "resamples.\n")


# -------------------------------------------------------
# 4. Recipe + workflow
# -------------------------------------------------------
# Recode the SMOTE *source* tables into the same classes as the task
# (e.g. fc_vs_rest -> "Funerary Context"/"rest"). Without this, the raw
# tables still hold the original labels, and the factor() coercion inside
# make_smote_train_sum_id() turns every non-"Funerary Context" row into NA,
# leaving 0 observations of "rest" -> themis::smote() errors.
smote_src_id  <- task_recode_data(result$db_id,  task_name, target_variable)$df
smote_src_sum <- task_recode_data(result$db_sum, task_name, target_variable)$df

recipe_obj <- build_recipe_for_task(
  df             = model_df,
  target_variable = target_variable,
  site_name_col  = site_name_col,
  raw_df_id      = smote_src_id,
  raw_df_sum     = smote_src_sum,
  site_id_col    = site_id_col,
  mni_col        = mni_col,
  feature_cols   = feature_cols,
  use_smote      = TRUE,
  use_downsample = use_downsample
)

# -------------------------------------------------------
# 5. Verify baking before wiring into workflow
# -------------------------------------------------------
prepped <- recipes::prep(recipe_obj, training = model_df)

# Should contain original + synthetic rows, no site_name, no mni
augmented_train <- recipes::bake(prepped, new_data = NULL)
cat("Train rows after augmentation:", nrow(augmented_train), "\n")
cat("Columns:", paste(names(augmented_train), collapse = ", "), "\n")
cat("Class balance:\n")
print(table(augmented_train[[target_variable]]))

# Test bake — skip = TRUE means SMOTE is bypassed, original rows only
processed_test <- recipes::bake(prepped, new_data = model_df)
cat("Test rows (should equal nrow(model_df)):", nrow(processed_test), "\n")

# -------------------------------------------------------
# 6. Model spec + workflow
# -------------------------------------------------------
model_spec <- get_model_spec(model_name)

wf <- workflows::workflow() %>%
  workflows::add_recipe(recipe_obj) %>%
  workflows::add_model(model_spec)

  # -------------------------------------------------------
  # 7. Parameter space + grid
  # -------------------------------------------------------
  param_set <- dials::parameters(wf)

  # keep only predictor columns, drop outcome + site_name,
  # and keep only numeric columns (what SVM etc. expect)
  # Keep only real predictor columns for parameter finalization.
  # Drop outcome + ID/site columns, then keep numeric predictors only.
  predictor_df <- model_df %>%
    dplyr::select(
      -dplyr::all_of(target_variable),
      -dplyr::all_of(site_name_col)
    ) %>%
    dplyr::select(where(is.numeric))

  param_set_final <- dials::finalize(
    param_set,
    predictor_df
  )

  grid <- dials::grid_latin_hypercube(
    param_set_final,
    size = grid_size
  )

  cat("Grid size:", nrow(grid), "hyperparameter combos.\n")
  # -------------------------------------------------------
  # 8. Metrics / ranking metric
  # -------------------------------------------------------
  #if (is_binary) {  # always false
  #  metrics <- metric_set(accuracy, f_meas, recall, precision, roc_auc)
  #  ranking_metric <- "roc_auc"
  #} else {
    metrics <- metric_set(accuracy, f_meas, recall, precision)
    ranking_metric <- "accuracy"
  #}

  options(yardstick.event_first = "first")

  # -------------------------------------------------------
  # 9. Parallel backend
  # -------------------------------------------------------
  # NOTE: step_meta_smote is a custom recipe step whose prep/bake S3 methods
  # live in the global environment rather than in a package. On PSOCK workers
  # those methods are not reliably discoverable by UseMethod(), so parallel
  # tuning fails with "no applicable method for 'prep'". We therefore default
  # to SEQUENTIAL execution, which dispatches correctly in the main session.
  # To use parallelism safely you must move step_meta_smote into a small
  # package (see notes), then set use_parallel = TRUE.
  if (isTRUE(use_parallel)) {
    num_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
    cl <- parallel::makeCluster(num_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages({
        library(recipes); library(themis); library(dplyr); library(rlang); library(tibble)
      })
    })
    parallel::clusterExport(
      cl,
      varlist = c("step_meta_smote", "prep.step_meta_smote", "bake.step_meta_smote",
                  "print.step_meta_smote", "make_smote_train_sum_id",
                  "default_smote_fun", "identity_aug"),
      envir = globalenv()
    )
    parallel::clusterCall(cl, function() {
      registerS3method("prep",  "step_meta_smote", prep.step_meta_smote,  envir = globalenv())
      registerS3method("bake",  "step_meta_smote", bake.step_meta_smote,  envir = globalenv())
      registerS3method("print", "step_meta_smote", print.step_meta_smote, envir = globalenv())
      invisible(NULL)
    })
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  # -------------------------------------------------------
  # 10. Run tuning
  # -------------------------------------------------------
  tune_res <- tune_grid(
    object    = wf,
    resamples = cv_folds,
    grid      = grid,
    metrics   = metrics,
    control   = control_grid(
      save_pred     = TRUE,
      parallel_over = "everything",
      verbose       = FALSE
    )
  )

  if (isTRUE(use_parallel)) parallel::stopCluster(cl)

  saveRDS(
    tune_res,
    file = file.path(
      RESULTS_DIR,
      paste0("tune_object_", model_name, ".rds")
    )
  )

  cat("Tuning complete for model:", model_name, "\n")

  # -------------------------------------------------------
  # 11. Rank hyperparameter combos best -> worst
  # -------------------------------------------------------
  best_all <- show_best(
    tune_res,
    metric = ranking_metric,
    n      = Inf
  ) %>%
    arrange(desc(mean)) %>%
    mutate(rank = row_number()) %>%
    select(rank, everything())

  # save full ranking of hyperparams by ranking metric
  csv_full_path <- file.path(
    RESULTS_DIR,
    paste0("tuning_results_", model_name, ".csv")
  )
  readr::write_csv(best_all, csv_full_path)

  # top10 table (markdown)
  top10 <- head(best_all, 10)

  md_table_top10 <- knitr::kable(
    top10,
    format  = "markdown",
    align   = "c",
    caption = paste0(
      "Top 10 hyperparameter combos for ",
      model_name,
      " (task ", task_name,
      ", ", data_variant, " L", level_str, ") ranked by ",
      ranking_metric
    )
  )

  md_top10_path <- file.path(
    RESULTS_DIR,
    paste0("tuning_top10_", model_name, ".md")
  )
  writeLines(md_table_top10, md_top10_path)

  # -------------------------------------------------------
  # 12. Accuracy summary across ALL hyperparameter combos
  #     + accuracy plot
  # -------------------------------------------------------
  acc_all <- collect_metrics(tune_res, summarize = TRUE) %>%
    filter(.metric == "accuracy") %>%
    arrange(desc(mean)) %>%
    mutate(rank_by_accuracy = row_number()) %>%
    select(rank_by_accuracy, .config, mean, std_err, n, dplyr::everything())

  acc_csv_path <- file.path(
    RESULTS_DIR,
    "accuracy_summary_all_combos.csv"
  )
  readr::write_csv(acc_all, acc_csv_path)

  acc_plot_df <- acc_all %>%
    mutate(cfg_label = paste0("cfg", rank_by_accuracy))

  accuracy_plot <- ggplot(acc_plot_df,
                          aes(x = cfg_label, y = mean, fill = cfg_label)) +
    geom_col() +
    geom_errorbar(aes(ymin = mean - std_err,
                      ymax = mean + std_err),
                  width = 0.2) +
    labs(
      title = paste0("Accuracy by hyperparameter combo (",
                     model_name, ", ", task_name,
                     ", ", data_variant, " L", level_str, ")"),
      x = "Hyperparameter combo (ranked by accuracy)",
      y = "Mean CV Accuracy"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  acc_plot_path <- file.path(
    RESULTS_DIR,
    "accuracy_plot_all_combos.png"
  )
  ggsave(
    filename = acc_plot_path,
    plot     = accuracy_plot,
    width    = 8,
    height   = 5,
    dpi      = 300
  )

  # -------------------------------------------------------
  # 13. BEST combo details
  # -------------------------------------------------------
  best_params <- select_best(tune_res, metric = ranking_metric)
  best_config_id <- best_params$.config[1]

  # predictions for best config across CV folds
  preds_all <- collect_predictions(tune_res) %>%
    filter(.config == best_config_id)

  conf_matrices <- preds_all %>%
    group_by(id) %>%
    summarise(
      conf_matrix = list(
        yardstick::conf_mat(
          data     = pick(everything()),
          truth    = !!sym(target_variable),
          estimate = .pred_class
        )
      ),
      .groups = "drop"
    ) %>%
    mutate(id2 = id) %>%
    select(id, id2, conf_matrix)

  # classification report via your helper
  report <- compute_classification_report_from_confmats(conf_matrices)

  pretty_report <- report$raw %>%
    mutate(
      Precision = sprintf("%.3f (%.3f)", precision_mean, precision_sd),
      Recall    = sprintf("%.3f (%.3f)", recall_mean,    recall_sd),
      F1        = sprintf("%.3f (%.3f)", f1_mean,        f1_sd)
    ) %>%
    select(
      Class,
      Precision,
      Recall,
      F1,
      Support = support
    )

  acc_best_cfg <- collect_metrics(tune_res, summarize = FALSE) %>%
    filter(.config == best_config_id, .metric == "accuracy") %>%
    summarise(
      acc_mean = mean(.estimate, na.rm = TRUE),
      acc_sd   = sd(.estimate, na.rm = TRUE),
      n        = n(),
      .groups  = "drop"
    )

  overall_support <- pretty_report %>%
    filter(!(Class %in% c("macro_avg", "weighted_avg"))) %>%
    summarise(Support = sum(Support, na.rm = TRUE)) %>%
    pull(Support)

  accuracy_row <- tibble(
    Class     = "accuracy",
    Precision = "",
    Recall    = "",
    F1        = sprintf("%.3f (%.3f)",
                        acc_best_cfg$acc_mean,
                        acc_best_cfg$acc_sd),
    Support   = overall_support
  )

  pretty_report_final <- bind_rows(
    pretty_report %>% filter(!(Class %in% c("macro_avg", "weighted_avg"))),
    accuracy_row,
    pretty_report %>% filter(Class %in% c("macro_avg", "weighted_avg"))
  )

  # markdown classification report for BEST combo
  md_table_best <- knitr::kable(
    pretty_report_final,
    format  = "markdown",
    align   = "c",
    caption = paste0(
      "Classification report for BEST hyperparameter combo (",
      model_name, ", ", task_name,
      ", ", data_variant, " L", level_str, ")\n",
      "Best .config = ", best_config_id
    )
  )

  best_md_path <- file.path(
    RESULTS_DIR,
    paste0("classification_report_best_", model_name, ".md")
  )
  writeLines(md_table_best, best_md_path)

  # summary metrics table for BEST combo
  summary_best_cfg <- collect_metrics(tune_res, summarize = FALSE) %>%
    filter(.config == best_config_id) %>%
    group_by(.metric) %>%
    summarise(
      mean    = mean(.estimate, na.rm = TRUE),
      std_err = sd(.estimate, na.rm = TRUE) / sqrt(n()),
      std_dev = sd(.estimate, na.rm = TRUE),
      n       = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean))

  readr::write_csv(
    summary_best_cfg,
    file.path(RESULTS_DIR, "summary_metrics_best_config.csv")
  )

  readr::write_csv(
    pretty_report_final,
    file.path(RESULTS_DIR, "classification_report_best_config.csv")
  )

  # -------------------------------------------------------
  # 14. Excel workbook with everything
  # -------------------------------------------------------
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "all_results")
  openxlsx::writeData(wb, "all_results", best_all)

  openxlsx::addWorksheet(wb, "top10")
  openxlsx::writeData(wb, "top10", top10)

  openxlsx::addWorksheet(wb, "accuracy_summary_all")
  openxlsx::writeData(wb, "accuracy_summary_all", acc_all)

  openxlsx::addWorksheet(wb, "best_config_report")
  openxlsx::writeData(wb, "best_config_report", pretty_report_final)

  openxlsx::addWorksheet(wb, "best_config_metrics")
  openxlsx::writeData(wb, "best_config_metrics", summary_best_cfg)

  xlsx_path <- file.path(
    RESULTS_DIR,
    paste0("tuning_", model_name, ".xlsx")
  )
  openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

  # -------------------------------------------------------
  # 15. Prep return object for the driver summary
  # -------------------------------------------------------
  # best score for this model under the ranking metric
  best_metric_row <- summary_best_cfg %>%
    filter(.metric == ranking_metric)

  best_metric_value <- if (nrow(best_metric_row) > 0) {
    best_metric_row$mean[1]
  } else {
    NA_real_
  }

  cat(
    "\nArtifacts saved in:\n  ", RESULTS_DIR, "\n",
    " - tuning_results_", model_name, ".csv (all combos ranked)\n",
    " - accuracy_summary_all_combos.csv\n",
    " - accuracy_plot_all_combos.png\n",
    " - classification_report_best_", model_name, ".md\n",
    " - tuning_", model_name, ".xlsx\n",
    " - tune_object_", model_name, ".rds\n",
    sep = ""
  )

  # return summary for final comparison step
  return(list(
    model_name         = model_name,
    task_name          = task_name,
    data_variant       = data_variant,
    level              = level,
    ranking_metric     = ranking_metric,
    best_metric_value  = best_metric_value,
    best_config_id     = best_config_id,
    results_dir        = RESULTS_DIR,
    summary_best_cfg   = summary_best_cfg
  ))
}

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
skeleton_info_for_level <- function(level) {

  if (level == 1) {
    feature_names <- c(
      "cranium", "mandible", "hyoid","vertebrae_indet", "cervical",
      "thoracic", "lumbar", "sternum", "rib", "clavicle", "scapula",
      "humerus", "radius", "ulna", "carpal", "metacarpal", "hand_phalanx",
      "sacrum", "coccyx", "os_coxae", "femur", "tibia", "fibula",
      "patella", "tarsal", "metatarsal", "foot_phalanx"
    )

    skeleton_df <- data.frame(
      element = feature_names,
      count = c(
        1,1,1,24,7,12,5,1,24,2,2,2,2,2,16,10,28,1,1,2,2,2,2,2,14,10,28
      )
    )

  } else if (level == 2) {

    feature_names <- c(
      "cranium_2", "vertebrae_2", "thorax_2",
      "shoulder_2", "arms_2", "hands_2",
      "pelvic_2", "legs_2", "foot_2"
    )

    skeleton_df <- data.frame(
      element = feature_names,
      count = c(
        2,   # cranium_2
        25,  # vertebrae_2
        25,  # thorax_2
        4,   # shoulder_2
        6,   # arms_2
        54,  # hands_2
        4,   # pelvic_2
        8,   # legs_2
        52   # foot_2
      )
    )

  } else if (level == 3) {

    feature_names <- c(
      "cranium_3","column_3","long_bones_3","hand_3","foot_3"
    )

    skeleton_df <- data.frame(
      element = feature_names,
      count = c(
        2,   # cranium_3
        54,  # column_3
        18,  # long_bones_3
        54,  # hand_3
        52   # foot_3
      )
    )

  } else if (level == 4) {

    feature_names <- c(
      "cranium_4", "mandibula_4", "post_cranium_4"
    )

    skeleton_df <- data.frame(
      element = feature_names,
      count = c(
        1,   # cranium_4
        1,   # mandibula_4
        178  # post_cranium_4
      )
    )

  } else {
    stop("Unsupported level: ", level)
  }

  list(
    feature_names = feature_names,
    skeleton_df   = skeleton_df,
    n_features    = length(feature_names)
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
# 3. Build recipe for tuning
# -----------------------------------------------------------
build_recipe_for_task <- function(df,
                                  target_variable = "taphonomic_context",
                                  site_name_col   = "site_name",
                                  use_downsample  = TRUE) {

  # NOTE: if step_mutate_at() is deprecated in your recipes version,
  # replace with:
  #   step_mutate(across(all_numeric_predictors(), ~ .x / 100))

  rec <- recipe(
    as.formula(paste(target_variable, "~ .")),
    data = df
  ) %>%
    update_role(!!sym(site_name_col), new_role = "ID") %>%
    step_mutate_at(all_numeric_predictors(), fn = ~ . / 100) 

  if (use_downsample) {
    rec <- rec %>%
      themis::step_downsample(
        !!sym(target_variable),
        under_ratio = 1,
        skip        = TRUE,
        seed        = sample.int(10^5, 1)
      )
  }

  rec
}


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
  data_variant,
  task_name,
  model_name,
  grid_size     = 25,
  n_iterations  = 5,
  target_variable = "taphonomic_context",
  site_name_col   = "site_name"
) {

  cat("\n=====================================================\n")
  cat("TUNING START\n",
      "variant =", data_variant, "\n",
      "level   =", level, "\n",
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
    level, "_",
    model_name
  )
  dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

  # -------------------------------------------------------
  # 1. Load & recode data for this task
  # -------------------------------------------------------
  raw_df <- readr::read_csv(dataset_path, show_col_types = FALSE)

  task_out <- task_recode_data(
    raw_df         = raw_df,
    task_name      = task_name,
    target_variable= target_variable
  )

  df             <- task_out$df
  class_levels   <- task_out$class_levels
  use_downsample <- task_out$use_downsample
  is_binary      <- (length(class_levels) == 2)

  # snapshot processed data
  saveRDS(
    df,
    file = file.path(
      RESULTS_DIR,
      paste0("processed_", data_variant,
             "_L", level, "_", task_name, ".rds")
    )
  )

  # -------------------------------------------------------
  # 2. Save skeleton metadata for this level
  # -------------------------------------------------------
  skel_info <- skeleton_info_for_level(level)
  readr::write_csv(
    skel_info$skeleton_df,
    file.path(
      RESULTS_DIR,
      paste0("feature_reference_L", level, "_",
             data_variant, "_", task_name, ".csv")
    )
  )
  cat("Level", level, ": documented",
      skel_info$n_features, "features.\n")

  # -------------------------------------------------------
  # 3. CV folds
  # -------------------------------------------------------
  set.seed(42)
  cv_folds <- vfold_cv(
    df,
    v        = 3,
    repeats  = n_iterations,
    strata   = !!sym(target_variable)
  )
  saveRDS(cv_folds,
          file.path(RESULTS_DIR, "cv_folds_tuning.rds"))
  cat("CV folds:", nrow(cv_folds), "resamples.\n")

  # -------------------------------------------------------
  # 4. Recipe
  # -------------------------------------------------------
  recipe_obj <- build_recipe_for_task(
    df              = df,
    target_variable = target_variable,
    site_name_col   = site_name_col,
    use_downsample  = use_downsample
  )

  # -------------------------------------------------------
  # 5. Model spec with tunable params
  # -------------------------------------------------------
  model_spec <- get_model_spec(model_name)

  # -------------------------------------------------------
  # 6. Workflow
  # -------------------------------------------------------
  wf <- workflow() %>%
    add_recipe(recipe_obj) %>%
    add_model(model_spec)

  # -------------------------------------------------------
  # 7. Parameter space + grid
  # -------------------------------------------------------
  param_set <- parameters(wf)

  # keep only predictor columns, drop outcome + site_name,
  # and keep only numeric columns (what SVM etc. expect)
  predictor_df <- df %>%
    dplyr::select(
      -all_of(target_variable),
      -all_of(site_name_col)
    ) %>%
    dplyr::select(where(is.numeric))
  
  param_set_final <- finalize(param_set, predictor_df)
  
  grid <- grid_latin_hypercube(
    param_set_final,
    size = grid_size
  )
  

  cat("Grid size:", nrow(grid), "hyperparameter combos.\n")

  # -------------------------------------------------------
  # 8. Metrics / ranking metric
  # -------------------------------------------------------
  if (is_binary) {
    metrics <- metric_set(accuracy, f_meas, recall, precision, roc_auc)
    ranking_metric <- "roc_auc"
  } else {
    metrics <- metric_set(accuracy, f_meas, recall, precision)
    ranking_metric <- "accuracy"
  }

  options(yardstick.event_first = "first")

  # -------------------------------------------------------
  # 9. Parallel backend
  # -------------------------------------------------------
  num_cores <- max(1, parallel::detectCores(logical = TRUE) - 1)
  cl <- parallel::makeCluster(num_cores)
  doParallel::registerDoParallel(cl)

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

  parallel::stopCluster(cl)

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
      ", ", data_variant, " L", level, ") ranked by ",
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
                     ", ", data_variant, " L", level, ")"),
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
      ", ", data_variant, " L", level, ")\n",
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

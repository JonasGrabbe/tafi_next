# run_tuning_all.R
# Master tuning sweep across:
#   tasks:       c("4class", "cannibalism_vs_fc", "fc_vs_rest")
#   variants:    c("pmnau", "remnau")
#   levels:      1:4
#   models:      all 7
#
# For each (task, variant, level):
#   - tune each model over its grid
#   - collect each model's best score (AUC if binary, Accuracy if 4-class)
#   - save a final summary + bar plot comparing model winners

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
  library(themis)
  library(doParallel)
  library(openxlsx)
  library(knitr)
  library(rlang)
  library(ggplot2)
})

source("hpa_sdk.R")                # has compute_classification_report_from_confmats()
source("hyperparam_tuning_core.R") # includes skeleton_info_for_level, task_recode_data, build_recipe_for_task,
                                   # get_model_spec, tune_single_model_combo, etc.

tidymodels_prefer()

# -----------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------

tasks        <- c("4class", "cannibalism_vs_fc", "fc_vs_rest")
variants     <- c("pmnau", "remnau")
levels_vec   <- 1:4
models_vec   <- c("SVM_RBF", "SVM_Linear", "RandomForest",
                  "XGBoost", "DecisionTree", "KNN", "MLP")

# repeated CV repeats (3-fold * n_iterations_global)
n_iterations_global <- 30

# per-model hyperparameter grid sizes
# (rough guidance: higher for more complex models / bigger search space)
grid_sizes <- list(
  SVM_RBF      = 1200,  # (cost, rbf_sigma)
  SVM_Linear   = 500,  # (cost)
  RandomForest = 1500,  # (mtry, min_n, trees)
  XGBoost      = 1500,  # (trees, depth, lr, etc.) <- huge space
  DecisionTree = 500,  # (depth, min_n, cost_complexity)
  KNN          = 500,  # (neighbors, weight_func, dist_power)
  MLP          = 1500   # (hidden_units, penalty, epochs) - can behave messy
)

# -----------------------------------------------------------------
# Loop all combos
# -----------------------------------------------------------------

for (task_name in tasks) {
  for (variant in variants) {
    for (lvl in levels_vec) {

      dataset_path <- sprintf("%s_tb_level_%d.csv", variant, lvl)
      cat("\n\n>>> GLOBAL RUN START:",
          "task =", task_name,
          "| variant =", variant,
          "| level =", lvl, "\n")

      # collect per-model best summary for this combo
      combo_best_list <- list()

      for (model_name in models_vec) {

        this_grid_size <- grid_sizes[[model_name]]

        cat("---- tuning model:", model_name,
            "with grid_size =", this_grid_size, "----\n")

        res <- tune_single_model_combo(
          dataset_path    = dataset_path,
          level           = lvl,
          data_variant    = variant,
          task_name       = task_name,
          model_name      = model_name,
          grid_size       = this_grid_size,
          n_iterations    = n_iterations_global,
          target_variable = "taphonomic_context",
          site_name_col   = "site_name"
        )

        combo_best_list[[model_name]] <- res
      }

      # after all 7 models for this (task, variant, level):
      # build summary table of winners
      combo_best_df <- purrr::map_dfr(
        combo_best_list,
        ~ tibble(
            model          = .x$model_name,
            ranking_metric = .x$ranking_metric,
            best_score     = .x$best_metric_value,
            best_config_id = .x$best_config_id,
            results_dir    = .x$results_dir
          )
      )

      # sort by best_score desc
      combo_best_sorted <- combo_best_df %>%
        arrange(desc(best_score))

      # make a summary results dir for this combo
      SUMMARY_DIR <- paste0(
        "results_tuning_summary_",
        task_name, "_",
        variant, "_L",
        lvl
      )
      dir.create(SUMMARY_DIR, showWarnings = FALSE, recursive = TRUE)

      # save CSV summary of the 7 best configs
      summary_csv_path <- file.path(
        SUMMARY_DIR,
        "best_models_summary.csv"
      )
      readr::write_csv(combo_best_sorted, summary_csv_path)

      # bar plot comparing best model from each of the 7
      best_metric_label <- combo_best_sorted$ranking_metric[1]

      best_plot <- ggplot(combo_best_sorted,
                          aes(x = model, y = best_score, fill = model)) +
        geom_col() +
        labs(
          title = paste0("Best ", best_metric_label,
                         " per model (task ", task_name,
                         ", ", variant, " L", lvl, ")"),
          x     = "Model",
          y     = paste0("Best CV ", best_metric_label)
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "none")

      best_plot_path <- file.path(
        SUMMARY_DIR,
        "best_model_scores.png"
      )
      ggsave(
        filename = best_plot_path,
        plot     = best_plot,
        width    = 8,
        height   = 5,
        dpi      = 300
      )

      # final markdown summary of the 7 winners
      final_md_table <- knitr::kable(
        combo_best_sorted,
        format  = "markdown",
        align   = "c",
        caption = paste0(
          "Best tuned performance per model for task ",
          task_name, " (", variant, " L", lvl, "). ",
          "Models ranked by ", best_metric_label, "."
        )
      )

      final_md_path <- file.path(
        SUMMARY_DIR,
        "final_best_models.md"
      )
      writeLines(final_md_table, final_md_path)

      cat("<<< GLOBAL RUN END:",
          "task =", task_name,
          "| variant =", variant,
          "| level =", lvl, "\n",
          "Summary written to: ", SUMMARY_DIR, "\n\n")
    }
  }
}

cat("\nALL DONE.\n")

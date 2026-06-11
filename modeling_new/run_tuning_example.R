# run_tuning_all.R
# Master tuning sweep across:
#   tasks:       c("4class", "cannibalism_vs_fc", "fc_vs_rest")
#   variants:    c("pmnau", "remnau")
#   levels:      random_vectors(n= 500)
#   models:      all 7
#
# For each (task, variant, level):
#   - tune each model over its grid
#   - collect each model's best score (AUC if binary, Accuracy if 4-class)
#   - save a final summary + bar plot comparing model winners
#run <- function(){

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

source("c:/Users/Usuario/tafi_test/modeling/hpa_sdk.R")            # has compute_classification_report_from_confmats()
source("c:/Users/Usuario/tafi_test/modeling/hyperparam_tuning_core.R") # includes skeleton_info_for_level, task_recode_data, build_recipe_for_task,
                                   # get_model_spec, tune_single_model_combo, etc.

tidymodels_prefer()


# -----------------------------------------------------------------
# lvl vec start -> outsource
# -----------------------------------------------------------------
POSITION_MAX <- c(5, 1, 1, 3, 4, 3, 3, 5, 4, 3, 5)

random_vectors <- function(n) {
  # Total number of unique combinations possible
  max_combinations <- prod(POSITION_MAX + 1)
  if (n > max_combinations) stop(paste("n exceeds total unique combinations:", max_combinations))
  
  seen <- list()
  result <- list()
  
  while (length(result) < n) {
    candidate <- mapply(function(mx) sample(1:mx, size = 1), POSITION_MAX)
    key <- paste(candidate, collapse = "-")
    
    if (!key %in% seen) {
      seen <- c(seen, key)
      result <- c(result, list(candidate))
    }
  }
  
  result
}
# -----------------------------------------------------------------
# lvl vec end -> outsource
# -----------------------------------------------------------------

# data path 
excel_path <- "c:/Users/Usuario/tafi_test/hspa.xlsx"

# -----------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------
sum_setting <- TRUE #  or FALSE
tasks        <- c("fc_vs_rest")#, "4class", "cannibalism_vs_fc") #
variants     <-  c("ProportionalRep_MNIbased")  # c("pmnau", "remnau") # other
levels_vec   <- random_vectors(n= 8)  
models_vec   <- c(#"SVM_RBF",
# "SVM_Linear", 
"RandomForest"
#, "XGBoost", "DecisionTree", "KNN", "MLP"
)

# repeated CV repeats (3-fold * n_iterations_global)
n_iterations_global <- 1 #30 # 1

# -----------------------------------------------------------------
# Master accuracy CSV (one row per task/variant/level/model)
# -----------------------------------------------------------------
# This file accumulates the best accuracy of every model at every level.
# It is the input for the two summary plots at the end of the run.
MASTER_CSV <- "level_accuracy_master.csv"

# Names of the 11 skeletal sections, in the SAME ORDER as the level vector
# (POSITION_MAX). EDIT THESE to match your actual anatomy ordering.
SECTION_LABELS <- c(
  "Cranium", "Teeth", "Hyoid", "Vertebrae", "Thorax",
  "Shoulder", "Arm", "Hand", "Pelvis", "Leg", "Foot"
)
stopifnot(length(SECTION_LABELS) == length(POSITION_MAX))

# start each full run from a clean master file
if (file.exists(MASTER_CSV)) file.remove(MASTER_CSV)

# per-model hyperparameter grid sizes
# (rough guidance: higher for more complex models / bigger search space)
grid_sizes <- list(
  SVM_RBF      = 1,#00,  # (cost, rbf_sigma)
  #SVM_Linear   = 200,  # (cost)
  RandomForest = 1#00#,  # (mtry, min_n, trees)
  #XGBoost      = 500,  # (trees, depth, lr, etc.) <- huge space
  #DecisionTree = 500,  # (depth, min_n, cost_complexity)
  #KNN          = 500,  # (neighbors, weight_func, dist_power)
  #MLP          = 500   # (hidden_units, penalty, epochs) - can behave messy
)

# -----------------------------------------------------------------
# Loop all combos
# -----------------------------------------------------------------

for (task_name in tasks) {
  for (variant in variants) {
    for (lvl in levels_vec) {

      dataset_path <-  excel_path   #sprintf("%s_tb_level_%d.csv", variant, lvl)

      lvl_str <- paste(lvl, collapse = "_")

      cat("\n\n>>> GLOBAL RUN START:",
          "task =", task_name,
          "| variant =", variant,
          "| level =", lvl_str, "\n")

      # collect per-model best summary for this combo
      combo_best_list <- list()

      for (model_name in models_vec) {

        this_grid_size <- grid_sizes[[model_name]]

        cat("---- tuning model:", model_name,
            "with grid_size =", this_grid_size, "----\n")

        res <- tune_single_model_combo(
          dataset_path    = dataset_path,
          level           = lvl,
          sum = sum_setting, 
          data_variant    = variant,
          task_name       = task_name,
          model_name      = model_name,
          grid_size       = this_grid_size,
          n_iterations    = n_iterations_global,
          target_variable = "taphonomic_context",
          site_name_col   = "site_id"
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

      # -----------------------------------------------------------
      # Append this (task, variant, level) block to the master CSV.
      # One row per model: level identity + the 11 vector positions +
      # their sum + the model's best accuracy.
      # -----------------------------------------------------------
      master_block <- combo_best_sorted %>%
        transmute(
          task         = task_name,
          variant      = variant,
          level_str    = lvl_str,
          model        = model,
          best_accuracy = best_score
        )

      # add pos1..posN columns (one value of the level vector each) + sum
      for (i in seq_along(lvl)) {
        master_block[[paste0("pos", i)]] <- lvl[i]
      }
      master_block$sum_level <- sum(lvl)

      if (file.exists(MASTER_CSV)) {
        readr::write_csv(master_block, MASTER_CSV, append = TRUE)
      } else {
        readr::write_csv(master_block, MASTER_CSV)
      }
      cat("Appended", nrow(master_block), "rows to", MASTER_CSV, "\n")

      # make a summary results dir for this combo
      SUMMARY_DIR <- paste0(
        "results_tuning_summary_",
        task_name, "_",
        variant, "_L",
        lvl_str
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
                         ", ", variant, " L", lvl_str, ")"),
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
          task_name, " (", variant, " L", lvl_str, "). ",
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
          "| level =", lvl_str, "\n",
          "Summary written to: ", SUMMARY_DIR, "\n\n")
    }
  }
}

cat("\nALL DONE.\n")

# -----------------------------------------------------------------
# Final summary plots from the master CSV
# -----------------------------------------------------------------
source("c:/Users/Usuario/tafi_test/modeling_2/plot_level_accuracy.R")

make_level_plots(
  master_csv     = MASTER_CSV,
  section_labels = SECTION_LABELS,
  out_dir        = "results_level_accuracy_plots",
  per_model      = TRUE,
  combined       = TRUE
)
#}
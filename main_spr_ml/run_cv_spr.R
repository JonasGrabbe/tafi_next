# Repeated cross-validation of skeletal-part representation (SPR).
# Edit the paths and analysis settings below, then source this file in R.
# Results use three stratified folds and prespecified, fixed model settings.

DATA_ROOT <- "c:/Users/jonas.grabbe/tafi_test-2"
MODEL_ROOT <- file.path(DATA_ROOT, "modeling_5")
EXCEL_FILE <- file.path(DATA_ROOT, "hspa.xlsx")
CACHE_ROOT <- file.path(MODEL_ROOT, "ml7_prepared_cache")

SEED <- 42L
CV_FOLDS <- 3L
CV_REPEATS <- 5L

# Eleven anatomical hierarchy levels, one integer per section.
levels_vec <- list(
  c(1,1,1,1,2,3,3,1,1,1,1),
  #c(1,1,1,1,2,2,2,1,1,1,1),
  c(2,1,1,2,2,2,2,2,2,2,2) )

tasks <- c("fc_vs_rest", "cannibalism_vs_fc", "4class")
spaces <- c("abundance", "evenness", "abundance_evenness",
            "remnau", "pmnau", "engineered", "combined")
models_vec <- c("LogisticRegression", "KNN", "SVM_RBF", "RandomForest")

# Aggregate observations remain in the fitted data. Individual and metadata
# observations are interpolation donors only, never additional independent sites.
USE_SMOTE <- FALSE
USE_METADATA_AUG <- FALSE
SMOTE_USE_MNI <- FALSE
SMOTE_NEIGHBORS <- 3L
SYNTHETIC_STRATEGY <- "multiplier" # alternative: "fixed_total"
SYNTHETIC_MULTIPLIER <- 3.0        # common class target: 3 x largest training class
SYNTHETIC_TARGET_PER_CLASS <- 150L # used only when strategy is "fixed_total"
METADATA_MIN_UNIQUE_SITES <- 2L

# One prespecified parameter configuration per classifier, not grid tuning.
FIXED <- list(
  logistic_penalty = 0.01,  # nnet multinomial/binomial L2 weight decay
  knn_neighbors = 5L,
  svm_cost = 1.0,          # e1071 RBF gamma = 1 / number of predictors
  rf_mtry = 5L,           # capped at available predictors
  rf_min_n = 5L,
  rf_trees = 500L
)

RUN_TAG <- paste0("CV_", SYNTHETIC_STRATEGY, "_",
  if (SYNTHETIC_STRATEGY == "multiplier") paste0("x", SYNTHETIC_MULTIPLIER)
  else paste0("n", SYNTHETIC_TARGET_PER_CLASS),
  if (USE_SMOTE) "_interpolated" else "_originalOnly",
  if (USE_METADATA_AUG) "_metadata" else "_noMetadata",
  "_r", CV_REPEATS, "_seed", SEED)
OUTPUT_ROOT <- file.path(MODEL_ROOT, "SPR_publication_results", RUN_TAG)

if (!file.exists(EXCEL_FILE)) stop("Input workbook not found: ", EXCEL_FILE)
if (CV_FOLDS != 3L) stop("This release specifies three CV folds; set CV_FOLDS <- 3L")
setwd(DATA_ROOT) # The existing HSPA loader sources its dependencies relatively.
source(file.path(DATA_ROOT, "tapho_loader_new.R"))
source(file.path(DATA_ROOT, "try_shared_allocator_proposal.R"))
if (USE_METADATA_AUG) {
  source(file.path(MODEL_ROOT, "spr_metadata_augmentation.R"))
}
source(file.path(MODEL_ROOT, "spr_ml_core.R"))
source(file.path(MODEL_ROOT, "spr_common.R"))
foreach::registerDoSEQ() # Custom recipe steps are fitted in the R main process.
set.seed(SEED)
dir.create(CACHE_ROOT, recursive = TRUE, showWarnings = FALSE)

settings <- list(
  experiment = "Repeated stratified site-level CV", seed = SEED,
  cv_folds = CV_FOLDS, cv_repeats = CV_REPEATS, levels = levels_vec,
  tasks = tasks, spaces = spaces, models = models_vec,
  interpolation = USE_SMOTE, metadata_donors = USE_METADATA_AUG,
  smote_use_mni = SMOTE_USE_MNI, neighbors = SMOTE_NEIGHBORS,
  synthetic_strategy = SYNTHETIC_STRATEGY,
  synthetic_multiplier = SYNTHETIC_MULTIPLIER,
  synthetic_target_per_class = SYNTHETIC_TARGET_PER_CLASS,
  metadata_min_unique_sites = METADATA_MIN_UNIQUE_SITES,
  fixed_parameters = FIXED,
  input_workbook = normalizePath(EXCEL_FILE, winslash = "/")
)
source_manifest <- spr_record_run(OUTPUT_ROOT, settings,
                                  DATA_ROOT, MODEL_ROOT, EXCEL_FILE)

fold_parts <- list(); prediction_parts <- list(); metric_parts <- list()
message("SPR repeated CV: ", CV_FOLDS, " folds × ", CV_REPEATS,
        " repeats; seed ", SEED)
for (level in levels_vec) {
  level_id <- paste(level, collapse = "_")
  level_dir <- file.path(OUTPUT_ROOT, paste0("L", level_id))
  dir.create(level_dir, recursive = TRUE, showWarnings = FALSE)
  message("Anatomical level L", level_id)
  tables <- spr_prepare_level(level, DATA_ROOT, MODEL_ROOT, EXCEL_FILE,
                              CACHE_ROOT, source_manifest)
  # Undefined derived indices are explicitly encoded and listed in the output.
  tables <- ml7_resolve_derived_indices(tables,
    audit_file = file.path(level_dir,"derived_index_conventions.csv"))
  spr_write_site_coordinates(tables, level_dir)
  result <- ml7_analyze(
    level_tables = tables, level = level, tasks = tasks, spaces = spaces,
    models = models_vec, repeats = CV_REPEATS, seed = SEED, fixed = FIXED,
    use_smote = USE_SMOTE, use_metadata_aug = USE_METADATA_AUG,
    smote_use_mni = SMOTE_USE_MNI, smote_neighbors = SMOTE_NEIGHBORS,
    synthetic_strategy = SYNTHETIC_STRATEGY,
    synthetic_multiplier = SYNTHETIC_MULTIPLIER,
    synthetic_target_per_class = SYNTHETIC_TARGET_PER_CLASS,
    metadata_min_unique_sites = METADATA_MIN_UNIQUE_SITES,
    output_dir = level_dir)
  fold_parts[[length(fold_parts)+1L]] <- result$fold_metrics
  prediction_parts[[length(prediction_parts)+1L]] <- result$predictions
  metric_parts[[length(metric_parts)+1L]] <- result$summaries
  saveRDS(result, file.path(level_dir,"cv_results.rds"))
}

all_folds <- dplyr::bind_rows(fold_parts)
all_oof <- dplyr::bind_rows(prediction_parts)
all_summary <- dplyr::bind_rows(metric_parts)
readr::write_csv(all_folds, file.path(OUTPUT_ROOT,"metrics_per_cv_fold.csv"))
readr::write_csv(all_oof, file.path(OUTPUT_ROOT,"out_of_fold_predictions.csv"))
readr::write_csv(all_summary, file.path(OUTPUT_ROOT,"metrics_mean_sd_long.csv"))
ml7_publication_table(all_summary, OUTPUT_ROOT, CV_REPEATS, FIXED, SEED,
                      USE_SMOTE, USE_METADATA_AUG)
ml7_plot_heatmaps(all_summary, OUTPUT_ROOT)
message("CV complete. Tables, figures and predictions: ", OUTPUT_ROOT)

## Machine-learning workflow and hyperparameter optimization (pipeline continuation)

ML pipeline implemented in the via R scripts (`hyperparam_tuning_core.R`, `run_tuning_example.R`, with shared helpers in `hpa_sdk.R`). The workflow follows the same overall structure :data variants × anatomical aggregation levels × supervised classifiers under repeated stratified resampling), plus**systematic hyperparameter tuning**, a **larger model sweep**, and **higher-repeat resampling during tuning runs** (configurable).

------

### Input data and feature construction

**Data inputs.** For each anatomical aggregation level (L\in{1,2,3,4}) and each feature variant (V\in{\texttt{pmnau}, \texttt{remnau}}), the pipeline reads a CSV of the form:

[
\texttt{"_tb_level_.csv"}
]

Each dataset includes:

- `taphonomic_context` (target label)
- `site_name` (site identifier; excluded from predictors)
- a set of numeric predictors encoding skeletal-part representation at that aggregation level.

**Anatomical aggregation levels.** The expected feature set for each level is explicitly defined and versioned in `skeleton_info_for_level(level)`. In the current implementation, the feature dimensionality follows:

- **Level 1:** 27 element categories (finest resolution)
- **Level 2:** 9 aggregated regions
- **Level 3:** 5 super-aggregated regions
- **Level 4:** 3 coarse partitions

The function also stores a reference “skeleton count” per feature (number of expected skeletal units). These counts are used upstream for SPR-derived variables and are saved alongside each tuning run to ensure traceability.

**SPR variable derivation (upstream).** The helper functions in `hpa_sdk.R` (`convertMNEtoVariables_*`) formalize how SPR-derived metrics can be computed from an MNE vector and a skeleton reference table (element counts). Across versions, they implement consistent handling of missing values (NAs → 0), safe division for zero-count cases, and computation of standard taphonomic representation measures including:

- ( \text{MNAU}_i = \frac{\text{MNE}_i}{c_i} )
- Percent scaling relative to maximum (e.g., (%\text{MNAU}))
- observed relative representation ( \text{ReMNAU}_i )
- expected relative representation derived from skeleton proportions
- observed/expected ratios (AcReMNAU-style outputs)

In the tuning scripts, these derived metrics are assumed to already be represented in the CSV predictors (the modeling code is agnostic to the exact SPR transform, and treats all numeric columns as candidate predictors).

------

### Classification tasks and label recoding

The pipeline supports three task definitions (`task_recode_data()`):

1. **`4class`**: multiclass classification with four labels
   ({)Cannibalism, Carnivore, Funerary Context, Geological transport(})
2. **`cannibalism_vs_fc`**: binary classification
   ({)Cannibalism vs Funerary Context(})
3. **`fc_vs_rest`**: binary classification
   ({)Funerary Context vs rest(}), where rest combines Cannibalism, Carnivore, Geological transport.

Label recoding is performed deterministically, and class distributions are printed for audit at runtime.

------

### Preprocessing and class balancing

All modeling is implemented using **tidymodels workflows** with a task-specific preprocessing recipe (`build_recipe_for_task()`):

1. **ID handling.** `site_name` is assigned role `"ID"` and excluded from the predictor matrix.
2. **Numeric scaling.** All numeric predictors are rescaled by dividing by 100:
   [
   x' = \frac{x}{100}
   ]
   This standardizes inputs that are naturally on a 0–100 scale (e.g., percent-based SPR metrics), while remaining compatible with kernel methods, distance-based methods, and neural nets.
3. **Downsampling for binary tasks.** For the two binary tasks, the recipe optionally applies `themis::step_downsample()` with `under_ratio = 1`, which undersamples the majority class to match the minority class **within each resampling split** (`skip = TRUE`). This prevents information leakage by ensuring balancing is performed only on training folds.

Downsampling is disabled for the four-class task in the current configuration.

------

### Resampling design (repeated stratified v-fold CV)

Each (task, variant, level) dataset is evaluated using **stratified 3-fold cross-validation** repeated (R) times:

- `v = 3`
- `repeats = n_iterations` (argument passed into `tune_single_model_combo()`)

In the master sweep script (`run_tuning_example.R`), the tuning runs use:

- `n_iterations_global = 30` → (3\times 30 = 90) resamples per tuning call

Stratification is performed on `taphonomic_context` to maintain class proportions across folds. The resampling object is saved (`cv_folds_tuning.rds`) for reproducibility.

------

### Model families and tuned hyperparameters

The updated workflow expands beyond “default settings” by defining **tunable model specifications** (`get_model_spec()`) and optimizing them via grid search.

| Model (tidymodels)            | Engine    | Tuned hyperparameters                                        |
| ----------------------------- | --------- | ------------------------------------------------------------ |
| SVM (RBF) `svm_rbf`           | `kernlab` | `cost`, `rbf_sigma`                                          |
| SVM (linear) `svm_linear`     | `kernlab` | `cost`                                                       |
| Random Forest `rand_forest`   | `ranger`  | `mtry`, `min_n`, `trees`                                     |
| XGBoost `boost_tree`          | `xgboost` | `trees`, `tree_depth`, `learn_rate`, `loss_reduction`, `min_n`, `sample_size` |
| Decision Tree `decision_tree` | `rpart`   | `tree_depth`, `min_n`, `cost_complexity`                     |
| k-NN `nearest_neighbor`       | `kknn`    | `neighbors`, `weight_func`, `dist_power`                     |
| MLP `mlp`                     | `nnet`    | `hidden_units`, `penalty`, `epochs`                          |



------

### Hyperparameter search strategy

For each (task, variant, level, model), the workflow is:

1. **Assemble workflow.** `workflow() + recipe + model_spec`
2. **Define parameter space.** `parameters(wf)`
3. **Finalize data-dependent ranges.** The parameter set is finalized using the predictor matrix (e.g., for `mtry`, the valid range is constrained by the number of predictors):
   - predictors are defined as numeric columns excluding `taphonomic_context` and `site_name`.
4. **Generate grid**

**Grid sizes .** `run_tuning_example.R`:

- SVM_RBF: 1200
- SVM_Linear: 500
- DecisionTree: 500
- KNN: 500
- RandomForest: 1500
- XGBoost: 1500
- MLP: 1500

------

### Metrics, ranking, and model selection

Evaluation metrics are computed within resampling using yardstick:

- **Binary tasks:** `accuracy`, `f_meas`, `recall`, `precision`, `roc_auc`
  **Ranking metric:** `roc_auc`
- **Multiclass task (`4class`):** `accuracy`, `f_meas`, `recall`, `precision`
  **Ranking metric:** `accuracy`

The event level is set consistently (`options(yardstick.event_first = "first")`) for binary metrics.

After tuning (`tune_grid()`), configurations are ranked using `show_best(..., metric = ranking_metric, n = Inf)` and the best configuration is selected with `select_best()`.

To provide an interpretable summary for the best configuration, predictions are collected and confusion matrices are computed per resample; these are converted into a **classification report** using the helper `compute_classification_report_from_confmats()` (`hpa_sdk.R`), returning mean and SD for precision/recall/F1 across folds, plus macro and weighted averages.



------

### Outputs and reproducibility artifacts

Each tuning call creates a dedicated results directory:

[
\texttt{results_tuning___L_/}
]

Saved artifacts include:

- processed, task-filtered dataset snapshot (`processed_*.rds`)
- feature reference table for the level (`feature_reference_*.csv`)
- resampling folds object (`cv_folds_tuning.rds`)
- full ranked hyperparameter results (`tuning_results_*.csv`)
- top-10 configurations (`tuning_top10_*.md`)
- accuracy summary across all configs (`accuracy_summary_all_combos.csv`)
- accuracy bar plot across configs (`accuracy_plot_all_combos.png`)
- best-config classification report (`classification_report_best_*.md` + `.csv`)
- best-config metric summary (`summary_metrics_best_config.csv`)
- full Excel workbook bundling the above (`tuning_*.xlsx`)
- serialized tuning object (`tune_object_*.rds`)

At the end of each (task, variant, level) sweep across all seven models, the master driver script writes a **model-comparison summary** directory:

[
\texttt{results_tuning_summary___L/}
]

containing a CSV and plot comparing the best tuned score per model.

------

### End-to-end workflow summary (algorithmic view)

For each task (T), feature variant (V), and aggregation level (L):

1. Load `"<V>_tb_level_<L>.csv"`
2. Recode labels according to (T)
3. Build preprocessing recipe (scale predictors; downsample if binary)
4. Create stratified repeated CV (3-fold × repeats)
5. For each model family (M):
   - define tunable spec
   - finalize parameter ranges from predictors
   - sample Latin hypercube grid (size depends on (M))
   - run `tune_grid()` on repeated CV
   - select best configuration by AUC (binary) or accuracy (multiclass)
   - export ranked configs + best-config reports
6. Aggregate best tuned score per model into a per-combination summary table and plot

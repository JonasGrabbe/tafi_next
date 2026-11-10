# Publication graphics from completed SPR supervised classification experiments.
# Install this script and spr_ml_figure_core.R in MODEL_ROOT.
# Source this file in R to create all figures and machine-readable tables.

SPR_RUN_FILE <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if (is.null(SPR_RUN_FILE) || !nzchar(SPR_RUN_FILE)) {
  stop('Run this script with source() using its full file path',call.=FALSE)
}
MODEL_ROOT <- dirname(normalizePath(SPR_RUN_FILE,winslash='/',mustWork=TRUE))
ORIGINAL_RESULTS_DIR <- file.path(MODEL_ROOT,'SPR_publication_results')
AUGMENTED_RESULTS_DIR <- file.path(MODEL_ROOT,'SPR_publication_results_data_augment')
FIGURES_OUTPUT_DIR <- file.path(MODEL_ROOT,'SPR_ML_publication_figures')

# Both directories are searched automatically. If augmented outputs were
# written by run_cv_spr.R to SPR_publication_results/CV_*interpolated_metadata*,
# that run is discovered as well. A run's settings file or CV tag determines
# its augmentation status; a directory name alone is fallback evidence.
# If multiple runs of the same type exist, specify the complete CV_* path.

# Available values: augmented, original, augmentation_delta, tradeoff,
#                   confusion, level_delta
SPR_FIGURES_TO_MAKE <- if (exists('SPR_FIGURES_TO_MAKE',inherits=FALSE))
  SPR_FIGURES_TO_MAKE else c('augmented','original','augmentation_delta',
                             'tradeoff','confusion','level_delta')

source(file.path(MODEL_ROOT,'spr_ml_figure_core.R'),local=TRUE)

spr_generate_figures(
  original_root=ORIGINAL_RESULTS_DIR,
  augmented_root=AUGMENTED_RESULTS_DIR,
  output_dir=FIGURES_OUTPUT_DIR,
  figures=SPR_FIGURES_TO_MAKE
)

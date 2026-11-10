# Regenerate only the original publication figure and summary tables.
SPR_FIGURES_TO_MAKE <- 'original'
SPR_PLOT_FILE <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if (is.null(SPR_PLOT_FILE) || !nzchar(SPR_PLOT_FILE))
  stop('Source this figure launcher using its full filename',call.=FALSE)
source(file.path(dirname(normalizePath(SPR_PLOT_FILE,winslash='/',mustWork=TRUE)),
                 'run_all_publication_figures.R'))

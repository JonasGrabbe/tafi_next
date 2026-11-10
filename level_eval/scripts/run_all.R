#!/usr/bin/env Rscript
# Single R entry point: from package root run Rscript scripts/run_all.R.
get_dir <- function() {
  f <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  if(length(f)) return(dirname(normalizePath(sub('^--file=','',f[1]))))
  frames <- sys.frames()
  for (i in rev(seq_along(frames))) if(!is.null(frames[[i]]$ofile))
    return(dirname(normalizePath(frames[[i]]$ofile)))
  if (file.exists(file.path(getwd(),'scripts','01_analyse_Q1.R')))
    return(normalizePath(file.path(getwd(),'scripts')))
  normalizePath(getwd())
}
scripts <- get_dir()
root <- normalizePath(file.path(scripts,'..'),winslash='/',mustWork=TRUE)
options(q1.root=root)
needed <- c('ggplot2','patchwork','sandwich','car','mgcv','officer')
missing <- needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '),
  '\nInstall with: install.packages(c(',paste(sprintf('"%s"',missing),collapse=', '),'))')
set.seed(20260714L)
cat('R pipeline\nRoot: ',root,'\nSeed: 20260714\n',sep='')
for(name in c('01_analyse.R','04_verify_tables.R','02_render_figures.R','03_build_writeup.R')) {
  cat('Running ',name,'...\n',sep='')
  source(file.path(scripts,name),chdir=FALSE,local=FALSE)
}
cat('Completed: tables, 13 PDF+PNG figure pairs, and Markdown+DOCX.\n')

#!/usr/bin/env Rscript
# Verification is written entirely in R; it compares the new fitted tables to
# independently recalculated reference results from the earlier review.
get_root <- function() {
  if (!is.null(getOption('q1.root',NULL))) return(getOption('q1.root'))
  f <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  if(length(f)) return(normalizePath(file.path(dirname(sub('^--file=','',f[1])),'..')))
  normalizePath('..')
}
ROOT <- get_root()
TABS <- file.path(ROOT,'results','tables')
ARCH <- file.path(ROOT,'archive','previously_recalculated_tables')
readfile <- function(f,dir=TABS) utils::read.csv(file.path(dir,f),stringsAsFactors=FALSE,check.names=FALSE)
a <- readfile('01_total_detail_primary.csv')
g <- readfile('02_section_global_HC3_BH.csv')
contr <- readfile('03_section_contrasts_HC3_BH.csv')
m <- readfile('04_section_adjusted_means.csv')
stopifnot(nrow(a)==3L, nrow(g)==27L, nrow(contr)==78L, nrow(m)==105L)
stopifnot(all(table(g$task)==9L), all(table(contr$task)==26L),all(table(m$task)==35L))
checks <- list()
append <- function(family,task,section,level,previous,current) {
  checks[[length(checks)+1L]] <<- data.frame(family=family,task=task,section=section,
                                             level=as.integer(level),reference=previous,
                                             recalculated=current,difference=current-previous)
}
for(i in seq_len(nrow(contr))) {
  z <- contr[i,,drop=FALSE]
  now <- m[m$task==z$task & m$section==z$section & m$level==z$level,'adjusted_mean']
  baseline <- m[m$task==z$task & m$section==z$section & m$level==1L,'adjusted_mean']
  if(length(now)!=1 || length(baseline)!=1) stop('Missing standardized mean for contrast')
  append('internal_contrast_minus_means',z$task,z$section,z$level,z$estimate,now-baseline)
}
if(dir.exists(ARCH)) {
  if(file.exists(file.path(ARCH,'02_total_detail_and_sensitivities.csv'))) {
    old <- readfile('02_total_detail_and_sensitivities.csv',ARCH)
    old <- old[old$outcome=='accuracy_mean' & old$specification=='unadjusted',,drop=FALSE]
    for(i in seq_len(nrow(a))) {
      z <- a[i,,drop=FALSE]
      v <- old$slope_per_level[old$task==z$task]
      if(length(v)==1) append('reference_accuracy_slope',z$task,'',NA,v,z$estimate)
    }
  }
  if(file.exists(file.path(ARCH,'03_global_section_tests.csv'))) {
    old <- readfile('03_global_section_tests.csv',ARCH)
    for(i in seq_len(nrow(g))) {
      z <- g[i,,drop=FALSE]
      v <- old$F_HC3[old$task==z$task & old$section==z$section]
      if(length(v)==1) append('reference_section_global_F',z$task,z$section,NA,v,z$F_hc3)
    }
  }
  if(file.exists(file.path(ARCH,'04_section_contrasts.csv'))) {
    old <- readfile('04_section_contrasts.csv',ARCH)
    for(i in seq_len(nrow(contr))) {
      z <- contr[i,,drop=FALSE]
      v <- old$estimate[old$task==z$task & old$section==z$section & old$level==z$level]
      if(length(v)==1) append('reference_section_contrast',z$task,z$section,z$level,v,z$estimate)
    }
  }
  if(file.exists(file.path(ARCH,'05_adjusted_section_means.csv'))) {
    old <- readfile('05_adjusted_section_means.csv',ARCH)
    for(i in seq_len(nrow(m))) {
      z <- m[i,,drop=FALSE]
      v <- old$adjusted_mean[old$task==z$task & old$section==z$section & old$level==z$level]
      if(length(v)==1) append('reference_adjusted_mean',z$task,z$section,z$level,v,z$adjusted_mean)
    }
  }
}
verified <- do.call(rbind,checks)
utils::write.csv(verified,file.path(TABS,'07_reconciliation_checks.csv'),row.names=FALSE)
maxerror <- max(abs(verified$difference))
message('Reconciliation checks: ',nrow(verified),'; max absolute difference: ',
        format(maxerror,scientific=TRUE,digits=5))
if(maxerror>1e-7) stop('Recomputed analysis does not reproduce its reference to 1e-7. See 07_reconciliation_checks.csv')

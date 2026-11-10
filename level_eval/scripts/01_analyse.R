#!/usr/bin/env Rscript
# Q1: Skeletal representational resolution vs Random Forest performance.
# Reanalysis of already-computed task-by-vector CV summaries; NO model refitting.
# R-only publication pipeline | Version 2.0 (2026-11-08).
# Run: Rscript scripts/run_all.R from the package root (recommended).
# IMPORTANT: CV scores share archaeological observations. Inference is conditional
# on the evaluated vectors/CV design, NOT a test of generalization to new sites.

get_script_dir <- function() {
  if (!is.null(getOption('q1.root', NULL)))
    return(file.path(getOption('q1.root'),'scripts'))
  filearg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(filearg)) return(dirname(normalizePath(sub('^--file=', '', filearg[1]))))
  fr <- sys.frames()
  for (i in rev(seq_along(fr))) if (!is.null(fr[[i]]$ofile))
    return(dirname(normalizePath(fr[[i]]$ofile)))
  normalizePath(getwd())
}
root <- get_script_dir()
args <- commandArgs(trailingOnly = TRUE)
input <- if (length(args) >= 1L && nzchar(args[1]) && is.null(getOption('q1.root'))) args[1] else file.path(getOption('q1.root', file.path(root,'..')), 'data', 'accuracy_f1_RandomForest_by_level.csv')
outdir <- if (length(args) >= 2L && nzchar(args[2]) && is.null(getOption('q1.root'))) args[2] else file.path(getOption('q1.root', file.path(root,'..')), 'results', 'tables')

# Pre-specified design. Change these values only if the experiment differs.
EXPECTED_TASKS <- c('4class', 'cannibalism_vs_fc', 'fc_vs_rest')
EXPECTED_N <- 1000L
REPRESENTATION <- 'ProportionalRep'
DISCRETIZATION <- 'MNIbased'
PARTS <- c('cranium','teeth','hyoid','vertebrae','thorax','shoulder',
           'arms','hands','pelvic','legs','foot')
PART_LABELS <- c('Cranium','Teeth','Hyoid','Vertebrae','Thorax','Shoulder',
                 'Arms','Hands','Pelvis','Legs','Feet')
TASK_LABEL <- c('4class'='Four-class',
                'cannibalism_vs_fc'='Cannibalism vs funerary',
                'fc_vs_rest'='Funerary vs rest')
SEED <- 20260714L
set.seed(SEED)

packages <- c('ggplot2', 'sandwich', 'car', 'mgcv')
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop('Install R packages: ', paste(missing, collapse=', '))
if (!file.exists(input)) stop('Missing input: ', input,
  '\nFirst run collect_scores_by_level.R in the folder containing results_tuning_* directories.')
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

fail <- function(...) stop(paste0(...), call.=FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
fin <- function(x) all(is.finite(x))
write_tab <- function(x, f) utils::write.csv(x, file.path(outdir,f), row.names=FALSE, na='')

raw <- utils::read.csv(input, stringsAsFactors=FALSE, check.names=FALSE)
required <- c('task','level',PARTS,'f1_mean')
if (length(setdiff(required,names(raw)))) fail('Missing columns: ',paste(setdiff(required,names(raw)),collapse=', '))

# The collector's historical f1_mean is ACCURACY, not macro-F1.
# Enforce this rather than silently analysing a different metric.
if ('model' %in% names(raw)) {
  raw <- raw[!is.na(raw$model) & raw$model=='RandomForest',,drop=FALSE]
}
if ('representation' %in% names(raw)) {
  raw <- raw[!is.na(raw$representation) & raw$representation==REPRESENTATION,,drop=FALSE]
}
if ('discretization' %in% names(raw)) {
  raw <- raw[!is.na(raw$discretization) & raw$discretization==DISCRETIZATION,,drop=FALSE]
}
if (!nrow(raw)) fail('No RandomForest rows for ', REPRESENTATION, '/', DISCRETIZATION)
for (col in c(PARTS,'f1_mean',intersect(c('accuracy_mean','macro_f1_mean',
                                       'balanced_accuracy','mtry','trees','min_n'),names(raw)))) {
  raw[[col]] <- num(raw[[col]])
}
if ('accuracy_mean' %in% names(raw)) {
  mismatch <- is.finite(raw$accuracy_mean) & is.finite(raw$f1_mean) &
    abs(raw$accuracy_mean-raw$f1_mean)>1e-6
  if (any(mismatch)) fail('f1_mean differs from accuracy_mean. Confirm outcome semantics.')
}
if (!fin(raw$f1_mean) || any(raw$f1_mean<0 | raw$f1_mean>1))
  fail('Accuracy (f1_mean alias) has missing or out-of-range values.')
for (p in PARTS) if (!fin(raw[[p]]) || any(raw[[p]]<1 | raw[[p]]!=floor(raw[[p]])))
  fail('Invalid non-positive/non-integer skeletal level in ',p)
if (!setequal(unique(raw$task), EXPECTED_TASKS))
  fail('Task labels do not match the three preregistered tasks: ',paste(unique(raw$task),collapse=', '))
if (anyDuplicated(paste(raw$task,raw$level,sep='|')))
  fail('Duplicated task/level combinations after model/representation filtering.')

# Verify filename-encoded resolution matches anatomical level columns.
parse_level <- function(s) {
  vals <- suppressWarnings(as.integer(strsplit(sub('^L','',s),'_',fixed=TRUE)[[1]]))
  if (length(vals)!=length(PARTS) || anyNA(vals)) return(rep(NA_integer_,length(PARTS)))
  vals
}
encoded <- t(vapply(raw$level,parse_level,integer(length(PARTS))))
provided <- as.matrix(raw[,PARTS,drop=FALSE])
if (anyNA(encoded) || any(encoded!=provided)) fail('Level-vector labels disagree with the 11 anatomical columns.')

counts <- as.data.frame(table(factor(raw$task, levels=EXPECTED_TASKS)), stringsAsFactors=FALSE)
names(counts) <- c('task','n_vectors')
counts$n_unique_vectors <- vapply(EXPECTED_TASKS,function(t) length(unique(raw$level[raw$task==t])),integer(1))
counts$expected_n <- EXPECTED_N
counts$complete <- counts$n_vectors==EXPECTED_N & counts$n_unique_vectors==EXPECTED_N
write_tab(counts,'00_sample_integrity.csv')
if (!all(counts$complete)) fail('Expected exactly 1,000 unique level vectors per task; see 00_sample_integrity.csv')
identical_vector_sets <- all(vapply(EXPECTED_TASKS[-1], function(t)
  setequal(raw$level[raw$task==t], raw$level[raw$task==EXPECTED_TASKS[1]]),logical(1)))
if (!identical_vector_sets) warning('Tasks do NOT have identical sets of resolution vectors.')
raw$total_detail <- rowSums(raw[,PARTS,drop=FALSE])
raw$accuracy <- raw$f1_mean

robust <- function(fit) sandwich::vcovHC(fit,type='HC3')
coef_row <- function(fit,term, task, outcome,kind) {
  b <- stats::coef(fit)
  vc <- robust(fit)
  se <- sqrt(diag(vc))
  dfr <- stats::df.residual(fit)
  tval <- unname(b[term]/se[term])
  crit <- stats::qt(.975, dfr)
  data.frame(task=task, outcome=outcome, specification=kind,
    term=term, estimate=unname(b[term]), hc3_se=unname(se[term]),
    ci_low=unname(b[term]-crit*se[term]), ci_high=unname(b[term]+crit*se[term]),
    p=2*stats::pt(-abs(tval),dfr), adjusted_r2=summary(fit)$adj.r.squared)
}

# Nonlinearity check with mgcv (installed as required dependency).
fit_total <- function(d,task,outcome,kind='primary', hps=character()) {
  rhs <- c('total_detail', hps)
  f <- stats::reformulate(rhs,response=outcome)
  lmfit <- stats::lm(f,data=d)
  ans <- coef_row(lmfit,'total_detail',task,outcome,kind)
  if (kind=='primary' && outcome=='accuracy') {
    ans$n <- nrow(d)
    ans$spearman_rho <- suppressWarnings(stats::cor(d$total_detail,d[[outcome]],method='spearman'))
    ans$gam_deviance_explained <- NA_real_
    ans$gam_edf <- NA_real_
    if (requireNamespace('mgcv',quietly=TRUE) && length(unique(d$total_detail))>=6L) {
      s <- mgcv::s
      gam_formula <- stats::reformulate("s(total_detail, k=5, bs='cr')", response=outcome)
      mod <- mgcv::gam(gam_formula,data=d,method='REML')
      ans$gam_deviance_explained <- summary(mod)$dev.expl
      ans$gam_edf <- as.numeric(summary(mod)$s.table[1,'edf'])
    }
  }
  ans
}

# Empirical standardization of categorical additive models.
# Contrasts are vs the lowest observed level of each varying section.
fit_sections <- function(d,task,outcome,kind='primary',hps=character()) {
  varying <- PARTS[vapply(PARTS,function(p) length(unique(d[[p]]))>1,logical(1))]
  md <- d
  for (p in varying) md[[p]] <- factor(md[[p]], levels=sort(unique(md[[p]])))
  fit <- stats::lm(stats::reformulate(c(varying,hps),response=outcome),data=md)
  if (anyNA(stats::coef(fit))) fail('Rank-deficient section model: task=',task,', specification=',kind)
  V <- robust(fit)
  b <- stats::coef(fit)
  dfr <- stats::df.residual(fit)
  crit <- stats::qt(.975,dfr)
  aov_tab <- as.data.frame(car::Anova(fit,type=2,white.adjust='hc3'))
  aov_tab$term <- rownames(aov_tab)
  Fname <- grep('^F',names(aov_tab),value=TRUE)[1]
  pname <- grep('Pr\\(>F\\)',names(aov_tab),value=TRUE)[1]
  dfname <- grep('^Df$',names(aov_tab),value=TRUE)[1]
  global <- data.frame(task=task,outcome=outcome,specification=kind,
      section=aov_tab$term,df1=aov_tab[[dfname]],df2=dfr,
      F_hc3=aov_tab[[Fname]],p=aov_tab[[pname]])
  # Exclude the model residual row: only anatomical predictors form the
  # nine-test family of global section hypotheses.
  global <- global[global$section %in% varying,,drop=FALSE]
  global$q_within_task <- NA_real_
  global$q_within_task[global$section %in% varying] <-
    stats::p.adjust(global$p[global$section %in% varying],method='BH')
  global$partial_eta2_F_approx <- global$F_hc3*global$df1/(global$F_hc3*global$df1+global$df2)
  M <- stats::delete.response(stats::terms(fit))
  mm <- function(part,level) {
    nd <- md
    nd[[part]] <- factor(as.character(level),levels=levels(md[[part]]))
    X <- stats::model.matrix(M, data=nd, contrasts.arg=fit$contrasts)
    colMeans(X[,names(b),drop=FALSE])
  }
  means <- list(); contrasts <- list(); idx <- 0L
  for (sec in varying) {
    levels_sec <- levels(md[[sec]])
    basemat <- mm(sec, levels_sec[1])
    for (lev in levels_sec) {
      a <- mm(sec,lev)
      est <- as.numeric(crossprod(a,b))
      se <- sqrt(max(0,as.numeric(t(a)%*%V%*%a)))
      idx <- idx + 1L
      means[[idx]] <- data.frame(task=task,outcome=outcome,specification=kind,
         section=sec,level=as.integer(lev),adjusted_mean=est,hc3_se=se,
         ci_low=est-crit*se,ci_high=est+crit*se,
         n_vectors_at_level=sum(d[[sec]]==as.integer(lev)))
      if (identical(lev,levels_sec[1])) next
      z <- a-basemat
      c_est <- as.numeric(crossprod(z,b))
      c_se <- sqrt(max(0,as.numeric(t(z)%*%V%*%z)))
      contrasts[[length(contrasts)+1L]] <- data.frame(task=task,outcome=outcome,
         specification=kind, section=sec, level=as.integer(lev),
         reference_level=as.integer(levels_sec[1]), estimate=c_est,
         hc3_se=c_se, ci_low=c_est-crit*c_se,ci_high=c_est+crit*c_se,
         p=if(c_se>0) 2*stats::pt(-abs(c_est/c_se),dfr) else NA_real_)
    }
  }
  con <- do.call(rbind,contrasts)
  con$q_within_task <- stats::p.adjust(con$p,method='BH')
  list(global=global,contrast=con,means=do.call(rbind,means),fit=fit)
}

# RF parameters vary in the collected run summaries. They are included only in
# an explicitly labelled SENSITIVITY analysis, not silently mixed with draft's
# original unadjusted estimates. Conditioning on tuning results is not equivalent
# to nested CV or fixed hyperparameters.
hp_covariates <- function(d) {
  ret <- character(0)
  for (hp in c('mtry','trees','min_n')) {
    if (!(hp %in% names(d))) next
    if (!fin(d[[hp]]) || length(unique(d[[hp]]))<2L) next
    v <- if (hp %in% c('mtry','trees')) log1p(d[[hp]]) else d[[hp]]
    d[[paste0('hp_',hp)]] <- as.numeric(scale(v))
    ret <- c(ret,paste0('hp_',hp))
  }
  list(data=d, terms=ret)
}

primary_total <- list(); primary_global <- list(); primary_contrasts <- list(); primary_means <- list()
sensitivity <- list(); sensitivity_global <- list()
for (task in EXPECTED_TASKS) {
  d <- raw[raw$task==task,,drop=FALSE]
  p <- fit_total(d,task,'accuracy')
  s <- fit_sections(d,task,'accuracy')
  primary_total[[task]] <- p
  primary_global[[task]] <- s$global
  primary_contrasts[[task]] <- s$contrast
  primary_means[[task]] <- s$means

  hp <- hp_covariates(d)
  if (length(hp$terms)) {
    sensitivity[[paste0(task,'_hp')]] <- fit_total(hp$data,task,'accuracy',
                                                   kind='RF_hyperparameter_adjusted',hps=hp$terms)
    sensitivity_global[[task]] <- fit_sections(hp$data,task,'accuracy',
                                               kind='RF_hyperparameter_adjusted',hps=hp$terms)$global
  }
  # Outcome sensitivity with identical unadjusted/mtry-adjusted specifications.
  for (metric in c('macro_f1_mean','balanced_accuracy')) {
    if (!(metric %in% names(d)) || !fin(d[[metric]]) ||
        any(d[[metric]] < 0 | d[[metric]] > 1)) next
    sensitivity[[paste0(task,'_',metric,'_raw')]] <-
      fit_total(d,task,metric,kind='outcome_unadjusted')
    if (length(hp$terms)) {
      sensitivity[[paste0(task,'_',metric,'_hp')]] <-
        fit_total(hp$data,task,metric,kind='RF_hyperparameter_adjusted',hps=hp$terms)
    }
  }
}
total <- do.call(rbind,primary_total)
total$q_across_three_tasks <- stats::p.adjust(total$p,method='BH')
global <- do.call(rbind,primary_global)
contrasts <- do.call(rbind,primary_contrasts)
means <- do.call(rbind,primary_means)
write_tab(total,'01_total_detail_primary.csv')
write_tab(global,'02_section_global_HC3_BH.csv')
write_tab(contrasts,'03_section_contrasts_HC3_BH.csv')
write_tab(means,'04_section_adjusted_means.csv')
if (length(sensitivity)) write_tab(do.call(rbind,sensitivity),'05_sensitivity_total_detail.csv')
if (length(sensitivity_global)) write_tab(do.call(rbind,sensitivity_global),'06_sensitivity_hyperparameter_sections.csv')

# All plotting has been moved to 02_render_figures.R to keep the model-fitting
# stage separate, ensure one inferential model source, and avoid duplicated plots.

manifest <- c(
  paste('Generated UTC:',format(Sys.time(),tz='UTC',usetz=TRUE)),
  paste('Input CSV:', normalizePath(input,winslash='/',mustWork=TRUE)),
  paste('Input MD5:', unname(tools::md5sum(input))),
  paste('Model:', 'RandomForest'),
  paste('Representation:', REPRESENTATION),
  paste('Discretization:', DISCRETIZATION),
  paste('Expected vectors per task:',EXPECTED_N),
  paste('Identical vector sets between tasks:',identical_vector_sets),
  paste('Outcome:', 'f1_mean means CV ACCURACY, not macro-F1'),
  paste('Seed:',SEED),
  paste('Packages:',paste(packages,collapse=', ')),
  'Statistical scope: design-level associations conditional on existing CV summaries; no external validation.',
  'CIs and BH-adjusted tests do NOT account for shared underlying assemblages, CV dependencies, or hyperparameter selection.',
  capture.output(sessionInfo())
)
writeLines(manifest,file.path(outdir,'REPRODUCIBILITY_MANIFEST.txt'))
cat('Q1 analysis complete. Saved to:',normalizePath(outdir),'\n')

# Outputs of this stage are canonical analysis tables. The plot and manuscript
# scripts read these exact tables, avoiding independent second-model refits.

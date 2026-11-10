#!/usr/bin/env Rscript
#Rebuild detailed scientific write-up from the canonical R model tables.
# Both publication-style Markdown and DOCX are generated in R. 
# PDF figures are created by 02_render_figures.R; this script embeds their PNGs.
get_q1_root <- function() {
  opt <- getOption('q1.root',NULL)
  if(!is.null(opt)) return(opt)
  f <- grep('^--file=',commandArgs(FALSE),value=TRUE)
  if(length(f)) return(normalizePath(file.path(dirname(sub('^--file=','',f[1])),'..')))
  normalizePath('..')
}
ROOT <- get_q1_root()
TABS <- file.path(ROOT,'results','tables')
FIGS <- file.path(ROOT,'results','figures')
OUT <- file.path(ROOT,'manuscript')
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
read_table <- function(file) {
  p <- file.path(TABS,file)
  if(!file.exists(p)) stop('Missing canonical table: ',p)
  utils::read.csv(p,stringsAsFactors=FALSE,check.names=FALSE)
}
all <- read_table('01_total_detail_primary.csv')
global <- read_table('02_section_global_HC3_BH.csv')
contrast <- read_table('03_section_contrasts_HC3_BH.csv')
means <- read_table('04_section_adjusted_means.csv')
sens <- read_table('05_sensitivity_total_detail.csv')
TASKS <- c('4class','cannibalism_vs_fc','fc_vs_rest')
TLABS <- c('4class'='Four-class discrimination',
           'cannibalism_vs_fc'='Cannibalism versus funerary context',
           'fc_vs_rest'='Funerary context versus remaining categories')
SECS <- c(cranium='cranium',teeth='teeth',hyoid='hyoid',vertebrae='vertebrae',
          thorax='thorax',shoulder='shoulder',arms='arm',hands='hand',pelvic='pelvis',legs='leg',foot='foot')
fmtpp <- function(z,d=2) sprintf(paste0('%+.',d,'f'),100*z)
fmt <- function(z,d=3) sprintf(paste0('%.',d,'f'),z)
statement <- function(t) {
  x <- all[all$task==t,,drop=FALSE]
  stopifnot(nrow(x)==1L)
  sprintf('%s: accuracy changed by %s percentage points per added summed level (95%% HC3 CI %s to %s pp; adjusted R² = %s; Spearman rho = %s).',
          TLABS[t],fmtpp(x$estimate,3),fmtpp(x$ci_low,3),fmtpp(x$ci_high,3),
          fmt(x$adjusted_r2),fmt(x$spearman_rho))
}
global_sentence <- function(t) {
  x <- global[global$task==t,,drop=FALSE]
  x <- x[x$section %in% names(SECS),,drop=FALSE]
  sig <- x$section[!is.na(x$q_within_task) & x$q_within_task<0.05]
  nos <- x$section[is.na(x$q_within_task) | x$q_within_task>=0.05]
  paste0('Global HC3 Type-II tests: ',length(sig),' of ',nrow(x),
    ' variable anatomical sections pass within-task BH q < 0.05 (',
    paste(SECS[sig],collapse=', '),'). ',
    if(length(nos)) paste0('Not significant: ',paste(SECS[nos],collapse=', '),'.') else '')
}
contrast_sentence <- function(t) {
  x <- contrast[contrast$task==t,,drop=FALSE]
  x <- x[!is.na(x$q_within_task) & x$q_within_task<.05,,drop=FALSE]
  if(!nrow(x)) return('No pairwise comparisons survive BH correction.')
  neg <- x[order(x$estimate),,drop=FALSE]
  pos <- x[order(-x$estimate),,drop=FALSE]
  a <- head(neg[neg$estimate<0,,drop=FALSE],4)
  b <- head(pos[pos$estimate>0,,drop=FALSE],3)
  description <- function(z) paste0(SECS[z$section],' L',z$level,
                                    ' (',fmtpp(z$estimate),' pp; q=',formatC(z$q_within_task,format='g',digits=2),')')
  paste0('Among ',nrow(x),' BH-positive level-to-level-1 comparisons, the largest negative estimates are ',
         paste(vapply(seq_len(nrow(a)),function(i) description(a[i,,drop=FALSE]),character(1)),collapse='; '),
         '. ',if(nrow(b)) paste0('Positive estimates include ',
            paste(vapply(seq_len(nrow(b)),function(i) description(b[i,,drop=FALSE]),character(1)),collapse='; '),'.') else '')
}
sensitivity_sentence <- function(t) {
  x <- sens[sens$task==t & sens$outcome=='accuracy' & sens$specification=='RF_hyperparameter_adjusted',,drop=FALSE]
  if(!nrow(x)) return('mtry-adjusted sensitivity was unavailable.')
  sprintf('After conditioning on varying Random Forest hyperparameters, the summed-detail slope was %s pp per level (95%% HC3 CI %s to %s). Because mtry and representation size are associated, this analysis is diagnostic, not an identified direct effect.',
          fmtpp(x$estimate,3),fmtpp(x$ci_low,3),fmtpp(x$ci_high,3))
}
lines <- c(
  '# Skeletal representational resolution and Random Forest classification performance',
  '',
  '## Detailed methods and results - reproducible R analysis',
  '',
  'Scope: analysis of existing cross-validated model-performance summaries; the Random Forest training and external predictive validation methods are reported separately in the ML chapter.',
  '',
  '### Data and study design',
  '',
  'We analysed 3,000 observations from 1,000 matched unique resolution vectors evaluated in each of three tasks: four-class discrimination, cannibalism versus funerary context, and funerary context versus the remaining classes. Each resolution vector assigns a detail level to cranium, teeth, hyoid, vertebrae, thorax, shoulder, arms, hands, pelvis, legs and feet. Teeth and hyoid were fixed at level 1 and their effects are not identifiable. The outcome is mean three-fold cross-validated accuracy. In the collector CSV the historical `f1_mean` field equals `accuracy_mean` and is not macro-F1.',
  '',
  '### Statistical analyses',
  '',
  'For each task, summed resolution was the sum of the eleven section levels. Primary analyses regressed cross-validated accuracy on summed resolution using OLS and HC3 heteroskedasticity-robust standard errors, with pointwise 95% coefficient intervals. Spearman correlation and a restricted-spline GAM (REML, k=5) provided descriptive monotonicity and nonlinearity checks. Multiple testing across three primary slopes used the Benjamini–Hochberg procedure.',
  '',
  'Mutually adjusted additive linear models treated each varying section as a categorical predictor. HC3 Type-II global tests evaluated each section after adjustment for the others (9 tests per task, within-task BH adjustment). Empirical-standardization means were estimated by assigning a focal section to a given level over all observed configurations of the other sections and averaging fitted predictions. Contrasts compared each level above 1 with level 1 (26 contrasts per task, within-task BH adjustment). Shown confidence intervals are nominal HC3 intervals, not family-wise or BH-adjusted intervals.',
  '',
  'Sensitivity analyses repeated the total-detail regression for true macro-F1 and balanced accuracy, and conditioned on varying RF hyperparameters including log1p(mtry). Neither covariate adjustment nor multiple-comparison correction resolves dependence due to repeated use of the same archaeological observations.',
  '',
  '### Overall summed-detail associations',
  '',
  vapply(TASKS,statement,character(1)),
  '',
  'Across the three tested tasks, the unadjusted total-detail associations are negative, but the amount of explanatory variation is strongly task-dependent. These slopes should not be read as causal effects or as evidence that an arbitrary increase in one skeletal level would achieve the pooled slope.',
  '',
  '### Task-level anatomical inference',
  ''
)
for (t in TASKS) lines <- c(lines,paste0('#### ',TLABS[t]),'',global_sentence(t),'',
                             contrast_sentence(t),'',sensitivity_sentence(t),'')
lines <- c(lines,
  '### Interpretation of the three complementary graph types',
  '',
  '1. Total-detail scatter plus contrast forest: shows the full 1,000-vector empirical distribution and the conditional L>1 versus L1 section estimates. The plot is a representation-design analysis, not a held-out-site learning curve.',
  '2. Effect map (proposed Figure 3): standardized mean accuracy difference from level 1 in percentage points, using a common zero-centered diverging scale across tasks. Bold labels are BH-significant contrasts. Absolute standardized accuracy maps (task-specific sequential scales) are retained in supplementary material; their colors must not be compared across tasks.',
  '3. Raw anatomical response panels: each panel re-uses the SAME 1,000 vector-level accuracy values. Panels are descriptive with repeated outcomes, not eleven independent experiments. The adjusted categorical model provides inferential comparisons.',
  '',
  '### Sensitivity and inferential scope',
  '',
  'Since the level vectors induce differences in model feature representations and optimal mtry, the hyperparameter-adjusted and unadjusted slopes answer different conditional questions. The positive association between summed detail and mtry can create unstable partial slopes, so the adjusted slopes must not be interpreted as direct causal resolution effects. Balanced accuracy and macro-F1 should be reported alongside overall accuracy where they change the interpretation.',
  '',
  'Each of the 1,000 scores per task shares underlying archaeological observations and three-fold cross-validation structure. HC3 coefficients and BH q-values describe uncertainty conditional on the sampled representation vectors and the observed CV summaries; they are not independent-sample biological intervals. The aggregated file cannot independently verify site-grouped CV, leakage controls, per-site predictions or new-assemblage generalizability. These aspects must be cross-referenced to the separate model-methods chapter.',
  '',
  'The analysis is additive and does not establish which full vector is optimal. Section-level level effects are not necessarily stable under untested interactions or unseen combinations. The sample consists of 1,000 configurations rather than an exhaustive skeletal coding universe.',
  '',
  '### Reproduction',
  '',
  'Canonical data: `data/accuracy_f1_RandomForest_by_level.csv`. Execute `Rscript scripts/run_all.R` from the package root. This runs the single canonical R estimation script, creates four publication summary figures and nine original-style reference figures (all in 450-dpi PNG and vector PDF), regenerates this Markdown/DOCX, and writes a session/software manifest. Figure position jitter uses seed 20260714; models are deterministic. No Python scripts are required or included.',
  '',
  'Full supplementary material is supplied as the 26 contrasts per task, nine anatomical-section global tests per task, empirically standardized section-by-level means, sensitivity tables, high-contrast absolute adjusted-accuracy maps, and raw-response panels for all three tasks. See FIGURE_PLAN.md for final captions and allocation.'
)
md <- file.path(OUT,'detailed_methods_results_R.md')
writeLines(lines,md,useBytes=TRUE)
message('Written ',md)

if (!requireNamespace('officer',quietly=TRUE)) {
  stop('Install R package officer to generate DOCX (Markdown was generated): install.packages("officer")')
}

doc <- officer::read_docx()
add <- function(s,style='Normal') {doc <<- officer::body_add_par(doc,value=s,style=style)}
add('Skeletal representational resolution and Random Forest classification performance','heading 1')
add('Detailed methods, fully recalculated results and graphical interpretation','Subtitle')
add('analysis of 3,000 existing task–resolution performance records. Random Forest training and general predictive validation are documented in the separate ML chapter.')
add('Data and design','heading 1')
add('Three tasks were evaluated on the same 1,000 randomly selected skeletal-resolution vectors per task. Each vector specifies levels for 11 anatomical sections. Teeth and hyoid were always level 1 and have no estimable independent effect. All analysed scores are three-fold cross-validated model-performance summaries; no Random Forest models are retrained by this R analysis.')
add('The collector field f1_mean is a historical alias of accuracy_mean. Genuine macro-F1 is macro_f1_mean and is analysed only in labelled sensitivity analyses.')
add('Statistical methods','heading 1')
add('Per-task summed-resolution models use OLS with HC3 robust slope standard errors and nominal 95% intervals, Spearman correlations and a penalized GAM nonlinearity check (REML, k = 5). A separate additive categorical model estimates mutually adjusted anatomical-section effects. Type-II HC3 omnibus tests and 26 higher-level-versus-level-1 contrasts per task use separate within-task Benjamini–Hochberg corrections. Section-level adjusted means use empirical standardization over observed other-section configurations.')
add('Sensitivity analyses evaluate macro-F1, balanced accuracy and the effect of conditioning on varying model hyperparameters. Adjusted RF hyperparameter associations are not interpreted causally.')
add('Total skeletal detail and classification accuracy','heading 1')
for(t in TASKS) add(statement(t))
slopetab <- do.call(rbind,lapply(TASKS,function(t) {
  r <- all[all$task==t,,drop=FALSE]
  data.frame(Task=TLABS[t],`Slope (pp per level)`=fmtpp(r$estimate,3),
             `HC3 95% CI (pp)`=paste0('[',fmtpp(r$ci_low,3),', ',fmtpp(r$ci_high,3),']'),
             `Adj. R2`=fmt(r$adjusted_r2,3),check.names=FALSE)
}))
doc <- officer::body_add_table(doc,slopetab,first_row=TRUE)
add('These regression associations are negative across tasks, but vary in magnitude and should not be generalized to a new assemblage without independent validation.')
add('Anatomical-section results','heading 1')
for(t in TASKS) {
  add(TLABS[t],'heading 2')
  add(global_sentence(t));add(contrast_sentence(t));add(sensitivity_sentence(t))
}
add('Figures: complementary presentation of the same analysis','heading 1')
figspec <- list(
 list('Figure_01_total_detail_three_tasks.png','Figure 1. Three-task observed accuracy and HC3 fitted trends.',6.2,2.18),
 list('Figure_02_all_78_contrasts.png','Figure 2. All 78 conditional section–level contrasts (26 per task).',6.2,4.06),
 list('Figure_03_relative_effect_maps.png','Figure 3. Shared-scale adjusted differences in accuracy versus level 1 (percentage points).',6.2,5.61),
 list('Figure_04_outcome_mtry_sensitivity.png','Supplementary figure: outcome metric and RF mtry sensitivity.',6.2,2.25),
 list('Reference_1_total_and_adjusted_fc_vs_rest.png','Supplementary example A. Total detail and section contrasts for funerary versus rest.',5.2,6.2),
 list('Reference_2_absolute_adjusted_map_4class.png','Supplementary example B. Absolute adjusted accuracy by section and level for the four-class task.',6.2,3.26),
 list('Reference_3_raw_section_panels_4class.png','Supplementary example C. Raw section-panel responses for the four-class task.',6.2,4.30)
)
for(f in figspec) {
  file <- file.path(FIGS,f[[1]])
  if (!file.exists(file)) stop('Figure missing before building DOCX: ',file)
  # Landscape original-style panels are intentionally presented as a representative
  # set; every figure for all three tasks is separately archived in figures/.
  if (grepl('Reference_1',f[[1]],fixed=TRUE)) doc <- officer::body_add_break(doc)
  doc <- officer::body_add_img(doc,src=file,width=f[[3]],height=f[[4]])
  add(f[[2]])
}
add('Scientific limitations and scope','heading 1')
add('The 1,000 vector-level summaries per task are not 1,000 independent archaeological samples. HC3 and BH correction do not account for shared site observations, correlated CV folds, choice of tuned hyperparameters, or performance on unseen assemblages. The present input cannot reconstruct per-site predictions or audit group-based resampling. Those aspects are specified in and must be cross-referenced to the ML chapter.')
add('Categorical models are additive; standardized marginal cells cannot be combined to claim an optimum untested vector.')
add('Reproducibility','heading 1')
add('Execute Rscript scripts/run_all.R. All generated numerical tables and 13 figure pairs derive from the canonical CSV, with seed 20260714 controlling plot jitter. The files results/tables/REPRODUCIBILITY_MANIFEST.txt and results/figures/FIGURE_README.txt document the source and environment. All scripts in the package are R files. See README.md for dependencies and caveats.')
output <- file.path(OUT,'detailed_methods_results_R.docx')
print(doc,target=output)
message('Written ',output)

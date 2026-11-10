#!/usr/bin/env Rscript
# Reads the canonical CSV plus fitted tables from 01_analyse.R.
# No statistical estimates are taken from images, and no RF models are fitted.
# Seed used only for deterministic horizontal point jitter: 20260714.

find_root <- function() {
  opt <- getOption('q1.root', NULL)
  if (!is.null(opt)) return(opt)
  args <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  if (length(args)) return(normalizePath(file.path(dirname(sub('^--file=', '', args[1])), '..'), mustWork=TRUE))
  normalizePath('..', mustWork=TRUE)
}
ROOT <- find_root()
INDIR <- file.path(ROOT, 'data')
TABDIR <- file.path(ROOT, 'results', 'tables')
FIGDIR <- file.path(ROOT, 'results', 'figures')
dir.create(FIGDIR, recursive=TRUE, showWarnings=FALSE)
for (pkg in c('ggplot2','patchwork','sandwich')) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop('Install R package: ',pkg)
}
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
SEED <- 20260714L
set.seed(SEED)
PARTS <- c('cranium','teeth','hyoid','vertebrae','thorax','shoulder','arms','hands','pelvic','legs','foot')
VARIABLE <- setdiff(PARTS, c('teeth','hyoid'))
NAME <- c(cranium='Cranium',teeth='Teeth',hyoid='Hyoid',vertebrae='Vertebrae',
          thorax='Thorax',shoulder='Shoulder',arms='Arms',hands='Hands',pelvic='Pelvis',
          legs='Legs',foot='Feet')
TASKS <- c('4class','cannibalism_vs_fc','fc_vs_rest')
TASK_TITLE <- c('4class'='Four-class','cannibalism_vs_fc'='Cannibalism vs funerary',
                'fc_vs_rest'='Funerary vs rest')
INK <- '#1E2D38'
BLUE <- '#23618A'
NEG <- '#AF5531'
NEUTRAL <- '#F7F7F5'
GRID <- '#E8EDF0'

rd <- function(f) {
  path <- file.path(TABDIR, f)
  if (!file.exists(path)) stop('Missing fitted table ',path,'; run 01_analyse.R first')
  utils::read.csv(path, check.names=FALSE, stringsAsFactors=FALSE)
}
raw_path <- file.path(INDIR, 'accuracy_f1_RandomForest_by_level.csv')
if (!file.exists(raw_path)) stop('Missing canonical CSV: ',raw_path)
raw <- utils::read.csv(raw_path, check.names=FALSE, stringsAsFactors=FALSE)
if (all(c('model','representation','discretization') %in% names(raw))) {
  raw <- raw[raw$model=='RandomForest' & raw$representation=='ProportionalRep' &
             raw$discretization=='MNIbased', , drop=FALSE]
}
stopifnot(nrow(raw)==3000L, all(table(raw$task)==1000L))
raw$accuracy <- as.numeric(raw$accuracy_mean)
stopifnot(all(is.finite(raw$accuracy)), all(abs(raw$accuracy-as.numeric(raw$f1_mean)) < 1e-6))
raw$total_detail <- rowSums(raw[,PARTS,drop=FALSE])
raw$task_label <- factor(TASK_TITLE[raw$task], levels=unname(TASK_TITLE[TASKS]))

total <- rd('01_total_detail_primary.csv')
ctr <- rd('03_section_contrasts_HC3_BH.csv')
means <- rd('04_section_adjusted_means.csv')
sens <- rd('05_sensitivity_total_detail.csv')
stopifnot(nrow(total)==3L, nrow(ctr)==78L, nrow(means)==105L)

save_pair <- function(p, stem, width, height) {
  pdf_dev <- if (capabilities('cairo')) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(FIGDIR,paste0(stem,'.pdf')),plot=p,width=width,height=height,
                  units='in',bg='white',device=pdf_dev,limitsize=FALSE)
  ggplot2::ggsave(file.path(FIGDIR,paste0(stem,'.png')),plot=p,width=width,height=height,
                  units='in',dpi=450,bg='white',limitsize=FALSE)
  message('Saved ',stem,' as vector PDF and 450 dpi PNG')
}
publication_theme <- function(base=10.5) {
  theme_minimal(base_size=base) +
    theme(text=element_text(color=INK),
          plot.title=element_text(face='bold',size=rel(1.24),margin=margin(b=5)),
          plot.subtitle=element_text(color='#485660',size=rel(.94),margin=margin(b=10)),
          plot.caption=element_text(color='#56616A',size=rel(.76),hjust=0,margin=margin(t=10)),
          panel.grid.major.x=element_line(color=GRID,linewidth=.26),
          panel.grid.major.y=element_blank(),
          panel.grid.minor=element_blank(),
          strip.background=element_rect(fill='#EFF2F4',color=NA),
          strip.text=element_text(face='bold',size=rel(.99),margin=margin(t=6,b=6)),
          axis.title=element_text(face='plain'),
          plot.margin=margin(10,14,10,10),
          legend.title=element_text(face='bold'),legend.position='right')
}
# Pointwise robust HC3 confidence band for the task-specific fitted mean.
robust_fit_line <- function(d) {
  fit <- stats::lm(accuracy ~ total_detail, data=d)
  V <- sandwich::vcovHC(fit, type='HC3')
  xx <- seq(min(d$total_detail), max(d$total_detail), length.out=181)
  X <- cbind('(Intercept)'=1, total_detail=xx)
  mu <- as.numeric(X %*% stats::coef(fit))
  se <- sqrt(pmax(0,rowSums((X %*% V)*X)))
  c95 <- stats::qt(.975, stats::df.residual(fit))
  data.frame(total_detail=xx, fit=mu, low=mu-c95*se, high=mu+c95*se)
}
line_data <- do.call(rbind,lapply(TASKS,function(t) {
  z <- robust_fit_line(raw[raw$task==t,,drop=FALSE])
  z$task <- t
  z$task_label <- TASK_TITLE[t]
  z
}))
line_data$task_label <- factor(line_data$task_label,levels=unname(TASK_TITLE[TASKS]))

# FIGURE 1 (main paper). Keep task-specific y ranges explicit in caption.
strip_names <- setNames(vapply(TASKS, function(t) {
  x <- total[total$task==t,,drop=FALSE]
  sprintf('%s\n%+.2f pp / level · adjusted R² %.3f',TASK_TITLE[t],100*x$estimate,x$adjusted_r2)
},character(1)),unname(TASK_TITLE[TASKS]))
main1 <- ggplot(raw,aes(total_detail,accuracy)) +
  geom_point(color='#253541',alpha=.25,size=.95,
             position=position_jitter(width=.13,height=0,seed=SEED)) +
  geom_ribbon(data=line_data,aes(x=total_detail,ymin=low,ymax=high),inherit.aes=FALSE,
              fill=BLUE,alpha=.14) +
  geom_line(data=line_data,aes(x=total_detail,y=fit),inherit.aes=FALSE,
            color=BLUE,linewidth=.95) +
  facet_wrap(~task_label,labeller=labeller(task_label=as_labeller(strip_names)),
             nrow=1,scales='free_y') +
  scale_y_continuous(labels=function(x) sprintf('%.0f%%',x*100)) +
  scale_x_continuous(breaks=seq(10,35,5)) +
  labs(title='Skeletal representational detail and cross-validated accuracy',
       subtitle='One point per tested resolution vector (1,000 per task); line = OLS fit; ribbon = pointwise HC3 95% CI',
       x='Total resolution (sum of 11 anatomical detail levels)',y='Cross-validated accuracy',
       caption='Y-axis ranges differ by task for legibility. Scores reuse the same archaeological data; uncertainty is conditional on the evaluated vectors.') +
  publication_theme(12) + theme(panel.grid.major.y=element_line(color=GRID,linewidth=.25),
                                  panel.spacing.x=grid::unit(1.2,'lines'))
save_pair(main1,'Figure_01_total_detail_three_tasks',15.4,5.6)

# Shared forest-plot rows and section dividers.
ct <- ctr
ct$section <- factor(ct$section,levels=VARIABLE)
ct <- ct[order(ct$section,ct$level),,drop=FALSE]
term_order <- unique(paste0(NAME[as.character(ct$section)],': L',ct$level,' vs L1'))
ct$term <- factor(paste0(NAME[as.character(ct$section)],': L',ct$level,' vs L1'),
                  levels=rev(term_order))
ct$task_label <- factor(TASK_TITLE[ct$task],levels=unname(TASK_TITLE[TASKS]))
ct$bh_sig <- is.finite(ct$q_within_task) & ct$q_within_task < .05
sizes <- vapply(VARIABLE,function(s) sum(ct$section==s & ct$task=='4class'),numeric(1))
dividers <- length(term_order)-cumsum(sizes)
dividers <- dividers[dividers>0]+.5

# FIGURE 2 (main paper). Significance is separate from CI crossing 0.
main2 <- ggplot(ct,aes(x=estimate*100,y=term)) +
  geom_hline(yintercept=dividers,color='#E9EDF0',linewidth=.35) +
  geom_vline(xintercept=0,linetype='22',linewidth=.65,color='#6A7581') +
  geom_segment(aes(x=ci_low*100,xend=ci_high*100,yend=term),
               color='#495861',linewidth=.55) +
  geom_point(aes(fill=bh_sig),shape=21,color=BLUE,size=2.1,stroke=.6) +
  scale_fill_manual(values=c('TRUE'=BLUE,'FALSE'='white'),breaks=c(TRUE,FALSE),
                    labels=c('BH q < 0.05','Not BH-significant'),name=NULL) +
  scale_x_continuous(breaks=seq(-3.0,1.0,1.0)) +
  facet_wrap(~task_label,nrow=1) +
  labs(title='Mutually adjusted anatomical-section contrasts',
       subtitle='All 78 higher-level-versus-level-1 contrasts; change in accuracy in percentage points',
       x='Adjusted accuracy difference (percentage points)',y=NULL,
       caption='Horizontal bars = nominal HC3 95% CIs, not multiplicity-adjusted. Filled circles pass BH correction over 26 contrasts within each task.') +
  publication_theme(11) + theme(legend.position='bottom',panel.grid.major.y=element_blank(),
                                axis.text.y=element_text(size=8.4),panel.spacing.x=grid::unit(1.2,'lines'))
save_pair(main2,'Figure_02_all_78_contrasts',16.5,10.8)

# FIGURE 3 (main paper): differences from a common level-1 baseline.
# Use a SHARED zero-centered divergence palette; absolute accuracy would obscure
# effects because between-task baselines differ more than within-task effects.
base <- means[means$level==1,c('task','section','adjusted_mean')]
names(base)[3] <- 'baseline'
effect <- merge(means,base,by=c('task','section'),sort=FALSE)
effect$change_pp <- 100*(effect$adjusted_mean-effect$baseline)
q <- ct[,c('task','section','level','q_within_task')]
effect <- merge(effect,q,by=c('task','section','level'),all.x=TRUE,sort=FALSE)
effect$significant <- !is.na(effect$q_within_task) & effect$q_within_task<.05
effect$section_label <- factor(NAME[as.character(effect$section)],levels=unname(NAME[VARIABLE]))
effect$task_label <- factor(TASK_TITLE[effect$task],levels=unname(TASK_TITLE[TASKS]))
effect$level_factor <- factor(effect$level,levels=1:5)
effect$cell <- ifelse(effect$level==1,'ref',sprintf('%+.2f',effect$change_pp))
effect$label_weight <- ifelse(effect$significant,'bold','plain')
lim <- ceiling(max(abs(effect$change_pp))/0.5)*0.5
# Text switches to white for high-saturation cells; baseline shown neutral.
effect$text_color <- ifelse(abs(effect$change_pp) >= .62*lim,'white',INK)
main3 <- ggplot(effect,aes(section_label,level_factor,fill=change_pp)) +
  geom_tile(width=.97,height=.95,color='white',linewidth=.55) +
  geom_text(aes(label=cell,color=text_color,fontface=label_weight),
            size=3.65,show.legend=FALSE) +
  scale_color_identity() +
  scale_fill_gradient2(low=NEG,mid=NEUTRAL,high=BLUE,midpoint=0,
                       limits=c(-lim,lim),breaks=pretty(c(-lim,lim),n=7),name='Δ accuracy\n(pp)') +
  facet_grid(task_label ~ .,switch='y') +
  scale_x_discrete(drop=FALSE) +
  labs(title='Resolution effects are specific to skeletal section and classification task',
       subtitle='Adjusted mean accuracy difference from level 1; identical zero-centered color limits across all tasks',
       x='Anatomical section',y='Detail level',
       caption='“ref” = level 1 (0 pp); bold numerical values have within-task BH q < 0.05. Blank = unavailable level. Estimates arise from additive models, not observed joint optima.') +
  publication_theme(11.5) +
  theme(strip.placement='outside',strip.text.y.left=element_text(angle=0,face='bold',size=11),
        strip.background=element_rect(fill='#F0F2F4',color=NA),
        axis.text.x=element_text(angle=35,hjust=1,size=10.0),
        panel.grid=element_blank(),panel.spacing.y=grid::unit(.8,'lines'),
        legend.key.height=grid::unit(2.3,'cm'))
save_pair(main3,'Figure_03_relative_effect_maps',11.6,10.5)

# FIGURE 4 supplementary sensitivity: the mtry-adjusted accuracy estimate is
# one diagnostic; macro F1 / balanced accuracy are unadjusted alternative outcomes.
sensdat <- do.call(rbind,lapply(TASKS,function(t) {
  p <- total[total$task==t,,drop=FALSE]
  z <- sens[sens$task==t,,drop=FALSE]
  make <- function(metric, spec,est,low,high) data.frame(task=t,metric=metric,
    specification=spec,estimate=est,ci_low=low,ci_high=high)
  ans <- make('Accuracy','Primary',p$estimate,p$ci_low,p$ci_high)
  for (i in seq_len(nrow(z))) {
    r <- z[i,,drop=FALSE]
    if(r$outcome=='accuracy' && r$specification=='RF_hyperparameter_adjusted') {
      ans <- rbind(ans,make('Accuracy','mtry adjusted',r$estimate,r$ci_low,r$ci_high))
    } else if(r$outcome=='macro_f1_mean') {
      ans <- rbind(ans,make('Macro-F1','Primary',r$estimate,r$ci_low,r$ci_high))
    } else if(r$outcome=='balanced_accuracy') {
      ans <- rbind(ans,make('Balanced accuracy','Primary',r$estimate,r$ci_low,r$ci_high))
    }
  }
  ans
}))
sensdat$task_label <- factor(TASK_TITLE[sensdat$task],levels=unname(TASK_TITLE[TASKS]))
sensdat$metric <- factor(sensdat$metric,levels=c('Balanced accuracy','Macro-F1','Accuracy'))
sensdat$y <- as.numeric(sensdat$metric)+ifelse(sensdat$specification=='mtry adjusted',-.13,.13)
sensdat$specification <- factor(sensdat$specification,levels=c('Primary','mtry adjusted'))
sensitivity_plot <- ggplot(sensdat,aes(x=estimate*100,y=y,color=specification)) +
  geom_vline(xintercept=0,linetype='22',color='#6C7881',linewidth=.6) +
  geom_segment(aes(x=ci_low*100,xend=ci_high*100,yend=y),linewidth=.68) +
  geom_point(aes(shape=specification),size=2.65,stroke=.72) +
  scale_y_continuous(breaks=1:3,labels=c('Balanced accuracy','Macro-F1','Accuracy'),
                     limits=c(.55,3.45)) +
  scale_color_manual(values=c('Primary'=BLUE,'mtry adjusted'=NEG),name=NULL) +
  scale_shape_manual(values=c('Primary'=16,'mtry adjusted'=17),name=NULL) +
  facet_wrap(~task_label,nrow=1) +
  labs(title='Sensitivity to outcome metric and Random Forest tuning',
       subtitle='Slope per added summed resolution level; conditional HC3 intervals',
       x='Change in model-performance score (percentage points per level)',y=NULL,
       caption='Only accuracy receives an mtry-adjusted sensitivity estimate; conditioning on mtry does not isolate a causal effect of resolution.') +
  publication_theme(11) + theme(legend.position='bottom',panel.grid.major.y=element_blank())
save_pair(sensitivity_plot,'Figure_04_outcome_mtry_sensitivity',14.6,5.3)

# SUPPLEMENT: absolute adjusted means, separate panel/scale PER TASK.
# Task-specific scales enhance within-task contrast but cannot be interpreted as
# common colors across tasks. The main paper uses the shared-scale delta map.
for(t in TASKS) {
  m <- means[means$task==t,,drop=FALSE]
  m$section_label <- factor(NAME[m$section],levels=unname(NAME[VARIABLE]))
  m$level_factor <- factor(m$level,levels=1:5)
  m$cell <- sprintf('%.3f\nn=%d',m$adjusted_mean,m$n_vectors_at_level)
  lo <- min(m$adjusted_mean)
  hi <- max(m$adjusted_mean)
  # High-contrast sequential palette that uses the available within-task range.
  m$text_col <- ifelse(m$adjusted_mean >= lo+.67*(hi-lo),'white',INK)
  p <- ggplot(m,aes(section_label,level_factor,fill=adjusted_mean)) +
    geom_tile(width=.97,height=.95,color='white',linewidth=.55) +
    geom_text(aes(label=cell,color=text_col),size=3.0,lineheight=1.0,show.legend=FALSE) +
    scale_color_identity() +
    scale_fill_gradientn(colours=c('#EDF4F9','#B5D2E6','#73A7C9','#3977A5','#163E65'),
                         limits=c(lo,hi),name='Adjusted\naccuracy',
                         labels=function(x) sprintf('%.3f',x)) +
    labs(title=paste0('Absolute adjusted accuracy by anatomical resolution — ',TASK_TITLE[t]),
         subtitle=sprintf('Within-task scale: %.3f–%.3f; figures have independent color ranges',lo,hi),
         x='Anatomical section',y='Detail level',
         caption='Cells are model-standardized means; n counts sampled vectors at the level, not independent sites. Teeth/hyoid excluded (fixed at level 1).') +
    publication_theme(11) +
    theme(axis.text.x=element_text(angle=35,hjust=1),panel.grid=element_blank(),
          legend.key.height=grid::unit(2.0,'cm'))
  save_pair(p,paste0('Reference_2_absolute_adjusted_map_',t),11.7,6.2)
}

# SUPPLEMENT: raw-response panels (all 11 sections), descriptive only.
for(t in TASKS) {
  d <- raw[raw$task==t,,drop=FALSE]
  rawlong <- do.call(rbind,lapply(PARTS,function(sec) data.frame(
    section=sec,level=as.integer(d[[sec]]),accuracy=d$accuracy)))
  rawlong$section_label <- factor(NAME[rawlong$section],levels=unname(NAME[PARTS]))
  rawlong$level_factor <- factor(rawlong$level,levels=1:5)
  slices <- split(rawlong,interaction(rawlong$section,rawlong$level,drop=TRUE))
  sumdat <- do.call(rbind,lapply(slices,function(z) {
    n <- nrow(z); se <- stats::sd(z$accuracy)/sqrt(n)
    ci <- if(n>1) stats::qt(.975,df=n-1)*se else 0
    data.frame(section=z$section[1],level=z$level[1],n=n,
               mu=mean(z$accuracy),low=mean(z$accuracy)-ci,high=mean(z$accuracy)+ci)
  }))
  sumdat$section_label <- factor(NAME[sumdat$section],levels=unname(NAME[PARTS]))
  sumdat$level_factor <- factor(sumdat$level,levels=1:5)
  p <- ggplot(rawlong,aes(level_factor,accuracy)) +
    geom_point(position=position_jitter(width=.14,height=0,seed=SEED),
               alpha=.14,size=.57,color='#2C333A') +
    geom_line(data=sumdat,aes(x=level_factor,y=mu,group=1),inherit.aes=FALSE,
              color=BLUE,linewidth=.65) +
    geom_pointrange(data=sumdat,aes(x=level_factor,y=mu,ymin=low,ymax=high),
                    inherit.aes=FALSE,color=BLUE,linewidth=.52,size=.29) +
    geom_text(data=sumdat,aes(x=level_factor,y=Inf,label=paste0('n=',n)),
              inherit.aes=FALSE,color='#68757D',size=2.2,vjust=1.22) +
    facet_wrap(~section_label,ncol=4,scales='free_x') +
    scale_y_continuous(labels=function(x) sprintf('%.0f%%',100*x)) +
    scale_x_discrete(drop=TRUE) +
    labs(title=paste0('Observed section-level response distributions — ',TASK_TITLE[t]),
         subtitle='Jittered observations are the same 1,000 vector accuracies reclassified by each skeletal section',
         x='Section detail level',y='Cross-validated accuracy',
         caption='Blue means and 95% bars are descriptive, not independent sample inference. Each section panel reuses the same task-level observations.') +
    publication_theme(10.5) +
    theme(strip.text=element_text(face='bold',size=9.5),
          panel.grid.major.y=element_line(color=GRID,linewidth=.25))
  save_pair(p,paste0('Reference_3_raw_section_panels_',t),14.0,9.8)
}

# Original combined figure type, retained for one-to-one comparability.
for(t in TASKS) {
  d <- raw[raw$task==t,,drop=FALSE]
  band <- line_data[line_data$task==t,,drop=FALSE]
  ct_t <- ct[ct$task==t,,drop=FALSE]
  a <- ggplot(d,aes(total_detail,accuracy)) +
    geom_point(alpha=.23,color=INK,size=.72,position=position_jitter(width=.13,height=0,seed=SEED)) +
    geom_ribbon(data=band,aes(x=total_detail,ymin=low,ymax=high),inherit.aes=FALSE,
                fill=BLUE,alpha=.14) +
    geom_line(data=band,aes(x=total_detail,y=fit),inherit.aes=FALSE,color=BLUE,linewidth=.8) +
    labs(title='A. Total detail and accuracy',x='Summed anatomical resolution',y='Cross-validated accuracy') +
    publication_theme(10.5) + theme(panel.grid.major.y=element_line(color=GRID))
  b <- ggplot(ct_t,aes(x=estimate*100,y=term)) +
    geom_hline(yintercept=dividers,color=GRID,linewidth=.32) +
    geom_vline(xintercept=0,linetype='22',color='#77818A') +
    geom_segment(aes(x=ci_low*100,xend=ci_high*100,yend=term),linewidth=.5,color='#4D5B64') +
    geom_point(aes(fill=bh_sig),shape=21,size=2,color=BLUE) +
    scale_fill_manual(values=c('TRUE'=BLUE,'FALSE'='white'),guide='none') +
    labs(title='B. Conditional anatomical contrasts',
         x='Change in accuracy (percentage points)',y=NULL,
         caption='Dots = adjusted differences vs level 1; bars = nominal HC3 95% CIs; filled = BH q < 0.05.') +
    publication_theme(9.8) + theme(axis.text.y=element_text(size=8),panel.grid.major.y=element_blank())
  both <- a/b + plot_layout(heights=c(1,1.6)) +
    plot_annotation(title=paste0('Resolution and classifier performance — ',TASK_TITLE[t]),
                    theme=theme(plot.title=element_text(face='bold',size=16)))
  save_pair(both,paste0('Reference_1_total_and_adjusted_',t),12.5,14.5)
}
writeLines(c('R publication figure suite (version 3)',
             paste('Input:',normalizePath(raw_path,winslash='/')),
             paste('Seed for plotting:',SEED),
             paste('Shared effect palette midpoint: 0 pp; symmetric limits:',sprintf('±%.2f pp',lim)),
             'Main figures: Figure_01, Figure_02, Figure_03.',
             'Supplementary: Figure_04 sensitivity, Reference_2 absolute maps (task-specific scales),',
             'Reference_3 raw panels; Reference_1 combined originals for provenance.',
             'Each PNG is exported at 450 dpi and each PDF is vector.',
             'No fitted statistics are inferred from rendered graphics.'),
           file.path(FIGDIR,'FIGURE_README.txt'))
message('R publication figures complete: 13 figure pairs.')

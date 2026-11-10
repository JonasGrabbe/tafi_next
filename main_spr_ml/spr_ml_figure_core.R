# Supervised SPR classification: reproducible publication figures and tables.
# Requires exported cross-validation results; model fitting is not repeated.
# Source this file through run_all_publication_figures.R or an individual driver.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})

spr_levels <- c(
  '2_1_1_2_2_2_2_2_2_2_2' = 'Standard',
  '1_1_1_1_2_3_3_1_1_1_1' = 'Performance-guided'
)
spr_tasks <- c(
  '4class' = 'Four-class (chance BA 25%)',
  'cannibalism_vs_fc' = 'Cannibalism vs Funerary (chance BA 50%)',
  'fc_vs_rest' = 'Funerary vs Rest (chance BA 50%)'
)
spr_spaces <- c(
  'abundance' = 'Mean part\nabundance',
  'evenness' = 'Heip\nevenness',
  'abundance_evenness' = 'Abundance\n+ evenness',
  'remnau' = 'ReMNAU',
  'pmnau' = 'MNI-normalized\nrepresentation',
  'engineered' = 'Engineered',
  'combined' = 'Combined'
)
spr_models <- c(
  'LogisticRegression' = 'Logistic regression',
  'KNN' = 'KNN',
  'SVM_RBF' = 'RBF-SVM',
  'RandomForest' = 'Random forest'
)
spr_order <- function(dat) {
  dat |>
    mutate(
      level_lab = factor(unname(spr_levels[as.character(level)]),
                         levels = unname(spr_levels)),
      task_lab = factor(unname(spr_tasks[as.character(task)]),
                        levels = unname(spr_tasks)),
      space_lab = factor(unname(spr_spaces[as.character(space)]),
                         levels = unname(spr_spaces)),
      model_lab = factor(unname(spr_models[as.character(model)]),
                         levels = rev(unname(spr_models)))
    )
}

spr_require <- function(ok, description) {
  if (!isTRUE(ok)) stop(description, call. = FALSE)
}

# The CV runner writes each experiment into a CV_<settings> subdirectory.
# Discover completed runs from their actual exports or archived per-level RDS.
spr_candidate_runs <- function(model_root, extra_roots=character()) {
  roots <- unique(c(extra_roots,
    file.path(model_root, 'SPR_publication_results'),
    file.path(model_root, 'SPR_publication_results_data_augment')))
  roots <- roots[dir.exists(roots)]
  all_files <- character()
  for (d in roots) {
    found <- list.files(d, recursive=TRUE, full.names=TRUE,
      pattern='^(metrics_per_cv_fold\\.csv|cv_results\\.rds)$',
      ignore.case=FALSE)
    all_files <- c(all_files, found)
  }
  all_files <- unique(normalizePath(all_files,winslash='/',mustWork=TRUE))
  all_files <- all_files[!grepl('/SPR_ML_publication_figures/', all_files,
                                     fixed=TRUE)]
  if (!length(all_files)) return(data.frame())
  locations <- vapply(all_files, function(path) {
    run <- dirname(path)
    if (basename(path) == 'cv_results.rds' && grepl('^L[0-9_]+$',basename(run))) {
      run <- dirname(run)
    }
    run
  }, character(1))
  roots_found <- unique(locations)
  tab <- lapply(roots_found, function(run) {
    settings <- file.path(run, 'analysis_settings.rds')
    mode <- NA_character_
    if (file.exists(settings)) {
      config <- tryCatch(readRDS(settings),error=function(e) NULL)
      if (is.list(config) && !is.null(config$interpolation) &&
          !is.null(config$metadata_donors)) {
        if (identical(config$interpolation,FALSE) &&
            identical(config$metadata_donors,FALSE)) mode <- 'Original-only'
        if (identical(config$interpolation,TRUE) &&
            identical(config$metadata_donors,TRUE)) mode <- 'Augmented'
      }
    }
    if (is.na(mode)) {
      key <- basename(run)
      if (grepl('originalOnly.*noMetadata',key,ignore.case=TRUE))
        mode <- 'Original-only'
      else if (grepl('interpolated.*metadata',key,ignore.case=TRUE))
        mode <- 'Augmented'
      else if (grepl('SPR_publication_results_data_augment',run,fixed=TRUE))
        mode <- 'Augmented (directory convention)'
    }
    data.frame(run=run,mode=mode,
      aggregated_csv=file.exists(file.path(run,'metrics_per_cv_fold.csv')),
      per_level_rds=length(list.files(run,pattern='^cv_results\\.rds$',
                           recursive=TRUE))>0L,
      stringsAsFactors=FALSE)
  })
  do.call(rbind,tab)
}

spr_resolve_run <- function(model_root,requested_dir,condition) {
  candidates <- spr_candidate_runs(model_root, requested_dir)
  if (!nrow(candidates)) {
    stop('No exported metrics_per_cv_fold.csv or per-level cv_results.rds ',
         'found under the two SPR results roots. Run the CV experiment to ',
         'completion, or supply the folder containing its exports.\n',
         'Checked: ',file.path(model_root,'SPR_publication_results'), '\n',
         '         ',file.path(model_root,'SPR_publication_results_data_augment'),
         call.=FALSE)
  }
  message('Discovered completed/recoverable CV runs:\n',
          paste(sprintf('  %s: %s [pooled CSV=%s; per-level RDS=%s]',
              ifelse(is.na(candidates$mode),'UNCLASSIFIED',candidates$mode),
              candidates$run,candidates$aggregated_csv,candidates$per_level_rds),
                collapse='\n'))
  exact <- candidates$run == normalizePath(requested_dir,winslash='/',
                                           mustWork=FALSE)
  valid <- !is.na(candidates$mode) &
      (candidates$mode == condition |
       (condition=='Augmented' & candidates$mode=='Augmented (directory convention)'))
  # A supplied exact run directory takes precedence if its settings agree.
  if (sum(exact & valid)==1L) valid <- valid & exact
  # If several candidates meet the condition, prefer one under its requested root.
  req_prefix <- paste0(normalizePath(requested_dir,winslash='/',mustWork=FALSE),'/')
  local <- valid & startsWith(candidates$run,req_prefix)
  if (sum(local)==1L) valid <- local
  matched <- candidates[valid,,drop=FALSE]
  if (nrow(matched)==0L) {
    stop('Could not identify a ',condition,' experiment. Check the list above. ',
         'The output RUN_TAG identifies original-only as originalOnly_noMetadata ',
         'and augmented as interpolated_metadata. ',
         'You can set the desired experiment folder explicitly.',call.=FALSE)
  }
  if (nrow(matched)>1L) {
    stop('Multiple ',condition,' runs are present. Set the corresponding ',
         'RESULTS_DIR to one specific CV_* folder:\n',
         paste(matched$run,collapse='\n'), call.=FALSE)
  }
  message('Selected ',condition,' run: ',matched$run)
  matched$run[[1]]
}

spr_ensure_pooled_results <- function(run_dir, condition, staging_root) {
  merged <- file.path(run_dir,'metrics_per_cv_fold.csv')
  if (file.exists(merged) && file.exists(file.path(run_dir,'out_of_fold_predictions.csv')))
    return(run_dir)
  # A run interrupted after completing both anatomical levels may have
  # cv_results.rds for each level but not the root-level combined exports.
  rds_files <- list.files(run_dir,recursive=TRUE,full.names=TRUE,
                          pattern='^cv_results\\.rds$')
  rds_files <- rds_files[grepl('/L[0-9_]+/cv_results\\.rds$',
                  gsub('\\\\','/',rds_files))]
  if (!length(rds_files))
    stop('Run lacks pooled metrics_per_cv_fold.csv and per-level cv_results.rds: ',
         run_dir,call.=FALSE)
  parts <- lapply(rds_files,readRDS)
  needed <- c('fold_metrics','predictions')
  if (!all(vapply(parts,function(x) is.list(x) &&
              all(needed %in% names(x)),logical(1)))) {
    stop('Per-level cv_results.rds do not contain fold_metrics and predictions: ',
         run_dir, call.=FALSE)
  }
  new_folds <- dplyr::bind_rows(lapply(parts, `[[`, 'fold_metrics'))
  new_oof <- dplyr::bind_rows(lapply(parts, `[[`, 'predictions'))
  spr_require(nrow(new_folds)==2520L && nrow(new_oof)==78400L,
    paste('Archived level files are incomplete:',nrow(new_folds),'fold rows,',
          nrow(new_oof),'predictions; expected 2520 and 78400.'))
  staging <- file.path(staging_root,if (condition=='Original-only')
                      '_reconstructed_original' else '_reconstructed_augmented')
  dir.create(staging,recursive=TRUE,showWarnings=FALSE)
  exported_folds <- file.path(staging,'metrics_per_cv_fold.csv')
  exported_oof <- file.path(staging,'out_of_fold_predictions.csv')
  readr::write_csv(new_folds, exported_folds)
  readr::write_csv(new_oof, exported_oof)
  spr_require(file.exists(exported_folds) && file.exists(exported_oof),
              paste('Could not write reconstructed CSV exports to', staging))
  writeLines(c(paste('Reconstructed from',run_dir),
    paste(normalizePath(rds_files,winslash='/',mustWork=TRUE),
          unname(tools::md5sum(rds_files)),sep=' | ')),
    file.path(staging,'source_rds_provenance.txt'))
  message('Recovered full result exports from ',length(rds_files),
          ' level RDS files into ',staging)
  staging
}

spr_find_file <- function(root, filename, required=TRUE) {
  spr_require(dir.exists(root),
              paste('Experiment directory does not exist:', root))
  # Look first for the exact exported CSV. Reconstructed exports are deliberately
  # stored under SPR_ML_publication_figures/_reconstructed_<condition>/.
  # Their location is valid because spr_ensure_pooled_results() returns that
  # directory explicitly. Never exclude files based on the output-root name.
  direct <- file.path(root, filename)
  if (file.exists(direct)) {
    return(normalizePath(direct, winslash='/', mustWork=TRUE))
  }
  # Recursive fallback supports a user-specified nested experiment directory.
  files <- list.files(root, recursive=TRUE, full.names=TRUE,
                      include.dirs=FALSE)
  files <- files[basename(files)==filename]
  if (!length(files)) {
    if (required) {
      available <- list.files(root, recursive=FALSE, full.names=FALSE)
      stop('Missing ', filename, ' below ', root,
           '. Available entries: ', paste(available, collapse=', '),
           '. Check the selected CV run and the reconstructed CSV folder.',
           call.=FALSE)
    }
    return(NULL)
  }
  if (length(files)>1L) {
    stop('Multiple copies of ',filename,' under ',root,
         '; select one run folder:\n',paste(files,collapse='\n'),
         call.=FALSE)
  }
  normalizePath(files[[1]], winslash='/', mustWork=TRUE)
}

spr_load_experiment <- function(root, condition, load_predictions=TRUE) {
  fold_path <- spr_find_file(root,'metrics_per_cv_fold.csv')
  folds <- read_csv(fold_path,show_col_types=FALSE,progress=FALSE)
  required <- c('level','task','space','model','fold_index','n_features',
                'n_assess','balanced_accuracy','auroc','brier','precision',
                'recall','f1','accuracy')
  spr_require(all(required %in% names(folds)),
              paste('Missing columns in',fold_path,':',
                    paste(setdiff(required,names(folds)),collapse=', ')))
  folds <- folds |> mutate(
    level=as.character(level), task=as.character(task),
    space=as.character(space), model=as.character(model),
    fold_index=as.integer(fold_index), condition=condition
  )
  spr_require(all(is.finite(folds$balanced_accuracy)) &&
              all(is.finite(folds$auroc)) && all(is.finite(folds$brier)),
              paste('Invalid finite metric in',fold_path))
  spr_require(nrow(folds)==2520L,
              paste('Expected 2,520 fold observations in',fold_path,'; got',nrow(folds)))
  spr_require(!anyDuplicated(folds[c('level','task','space','model','fold_index')]),
              paste('Duplicate fold-metric keys in',fold_path))
  spr_require(identical(sort(unique(folds$level)),sort(names(spr_levels))) &&
              identical(sort(unique(folds$task)),sort(names(spr_tasks))) &&
              identical(sort(unique(folds$space)),sort(names(spr_spaces))) &&
              identical(sort(unique(folds$model)),sort(names(spr_models))),
              paste('Experimental levels, tasks, spaces or models differ in',fold_path))
  summary <- folds |>
    group_by(level,task,space,model) |>
    summarise(
      n_features=first(n_features), n_folds=n(),
      BA=mean(balanced_accuracy)*100,
      BA_sd=sd(balanced_accuracy)*100,
      AUROC=mean(auroc)*100,
      AUROC_sd=sd(auroc)*100,
      Brier=mean(brier), Brier_sd=sd(brier),
      F1=mean(f1)*100, F1_sd=sd(f1)*100,
      Accuracy=mean(accuracy)*100, Accuracy_sd=sd(accuracy)*100,
      .groups='drop'
    )
  spr_require(nrow(summary)==168L && all(summary$n_folds==15L),
              paste('Expected 168 configurations × 15 folds in',fold_path))
  predictions <- NULL
  if (load_predictions) {
    oof_path <- spr_find_file(root,'out_of_fold_predictions.csv')
    predictions <- read_csv(oof_path,
      col_types=cols(site_id=col_character(),.default=col_guess()),
      show_col_types=FALSE,progress=FALSE)
    needed <- c('site_id','truth','fold_index','level','task','space','model')
    spr_require(all(needed %in% names(predictions)),
      paste('Missing OOF columns in',oof_path))
    predictions <- predictions |>
      mutate(across(all_of(c('site_id','truth','level','task','space','model')),
                    as.character), fold_index=as.integer(fold_index),
             condition=condition)
    spr_require(nrow(predictions)==78400L,
      paste('Expected 78,400 OOF predictions in',oof_path,'; got',nrow(predictions)))
    spr_require(!anyDuplicated(predictions[c('site_id','level','task','space',
                                            'model','fold_index')]),
      paste('Duplicate site-fold predictions in',oof_path))
  }
  list(folds=folds, summary=summary, oof=predictions,
       file=fold_path, directory=root)
}

spr_pair_results <- function(original,augmented,output_dir) {
  key <- c('level','task','space','model')
  fkey <- c(key,'fold_index')
  a <- original$folds |> arrange(across(all_of(fkey)))
  b <- augmented$folds |> arrange(across(all_of(fkey)))
  spr_require(identical(a[fkey], b[fkey]),
              'Original-only and augmented folds are not matched.')
  spr_require(identical(a$n_assess,b$n_assess),
              'Assessment sample counts differ between corresponding folds.')
  spr_require(identical(a$n_features,b$n_features),
              'Predictor dimensions differ between augmentation conditions.')
  if (!is.null(original$oof) && !is.null(augmented$oof)) {
    ok <- c('site_id','truth','level','task','space','model','fold_index')
    p <- original$oof |> arrange(across(all_of(ok)))
    q <- augmented$oof |> arrange(across(all_of(ok)))
    spr_require(identical(p[ok],q[ok]),
                'Original-only and augmented OOF site/fold membership differs.')
  }
  # Full precision is retained; plots show rounded display labels only.
  paired <- original$summary |>
    inner_join(augmented$summary,by=key,suffix=c('_original','_aug')) |>
    mutate(
      delta_BA_pp=BA_aug-BA_original,
      delta_AUROC_pp=AUROC_aug-AUROC_original,
      delta_Brier=Brier_aug-Brier_original,
      delta_F1_pp=F1_aug-F1_original
    )
  spr_require(nrow(paired)==168L &&
              all(paired$n_features_original==paired$n_features_aug),
              'Incomplete configuration pairing or changed predictor dimensions.')
  write_csv(paired,file.path(output_dir,'All_168_Paired_Configurations.csv'))
  write_csv(original$summary,file.path(output_dir,'OriginalOnly_168_Configurations.csv'))
  write_csv(augmented$summary,file.path(output_dir,'Augmented_168_Configurations.csv'))
  summary <- paired |> group_by(task,level) |>
    summarise(n_configurations=n(),
              mean_delta_BA_pp=mean(delta_BA_pp),
              median_delta_BA_pp=median(delta_BA_pp),
              n_augmented_better=sum(delta_BA_pp>1e-10),
              n_augmented_worse=sum(delta_BA_pp<(-1e-10)),
              mean_delta_AUROC_pp=mean(delta_AUROC_pp),
              mean_delta_Brier=mean(delta_Brier), .groups='drop')
  get_best <- function(x,suffix) {
    x |> group_by(task,level) |>
      slice_max(order_by=BA,n=1,with_ties=FALSE) |>
      ungroup() |>
      select(task,level,model,space,BA,BA_sd) |>
      rename_with(~paste0(.x,suffix), .cols=c('model','space','BA','BA_sd'))
  }
  summary <- summary |>
    left_join(get_best(original$summary,'_original'),by=c('task','level')) |>
    left_join(get_best(augmented$summary,'_aug'),by=c('task','level')) |>
    mutate(task=factor(task,levels=names(spr_tasks)),
           level=factor(level,levels=names(spr_levels))) |>
    arrange(task,level)
  write_csv(summary,file.path(output_dir,'Main_Table_Augmentation_Summary_FullPrecision.csv'))
  manuscript_table <- summary |>
    transmute(Task=unname(c('Four-class','Cannibalism vs Funerary','Funerary vs Rest')[as.integer(task)]),
      Level=ifelse(as.character(level)==names(spr_levels)[1], 'Standard', 'Guided'),
      `Best original-only BA`=sprintf('%.1f ± %.1f%%',BA_original,BA_sd_original),
      `Best augmented BA`=sprintf('%.1f ± %.1f%%',BA_aug,BA_sd_aug),
      `Mean paired delta BA (pp)`=sprintf('%+.2f',mean_delta_BA_pp))
  write_csv(manuscript_table,file.path(output_dir,'Main_Table_ReadyForWord.csv'))
  list(paired=paired,summary=summary,manuscript_table=manuscript_table)
}

spr_theme <- function(base_size=10.5) {
  theme_minimal(base_size=base_size,base_family='sans') +
    theme(
      plot.title=element_text(size=17,face='bold',hjust=.5,margin=margin(b=6)),
      plot.subtitle=element_text(size=10.3,hjust=.5,colour='#48515b',margin=margin(b=12)),
      plot.caption=element_text(size=9,colour='#4d5660',hjust=.5,margin=margin(t=12)),
      strip.text=element_text(face='bold',size=10.5),
      strip.background=element_rect(fill='#f1f3f5',colour='#d5dbe1',linewidth=.3),
      panel.grid=element_blank(),
      axis.text.x=element_text(size=9),
      axis.text.y=element_text(size=9),
      legend.title=element_text(size=9.4),
      legend.text=element_text(size=8.8),
      plot.margin=margin(t=15,r=15,b=12,l=15)
    )
}

spr_save_plot <- function(plot,stem,output_dir,width,height) {
  png_file <- file.path(output_dir,paste0(stem,'.png'))
  pdf_file <- file.path(output_dir,paste0(stem,'.pdf'))
  svg_file <- file.path(output_dir,paste0(stem,'.svg'))
  ggsave(png_file,plot=plot,width=width,height=height,units='in',dpi=300,
         bg='white',limitsize=FALSE)
  if (capabilities('cairo')) {
    ggsave(pdf_file,plot=plot,width=width,height=height,units='in',
           device=grDevices::cairo_pdf,bg='white',limitsize=FALSE)
  } else {
    ggsave(pdf_file,plot=plot,width=width,height=height,units='in',
           device=grDevices::pdf,bg='white',limitsize=FALSE)
  }
  if (requireNamespace('svglite',quietly=TRUE)) {
    ggsave(svg_file,plot=plot,width=width,height=height,units='in',
           device=svglite::svglite,bg='white',limitsize=FALSE)
  } else {
    ggsave(svg_file,plot=plot,width=width,height=height,units='in',
           device=grDevices::svg,bg='white',limitsize=FALSE)
  }
  message('Saved: ',normalizePath(pdf_file,winslash='/',mustWork=TRUE))
  invisible(plot)
}

spr_plot_heatmap <- function(tab,stem,output_dir,mode=c('BA','delta','level_delta')) {
  mode <- match.arg(mode)
  lab_col <- if (mode=='BA') 'BA' else if (mode=='delta') 'delta_BA_pp' else 'delta_level_BA_pp'
  spr_require(lab_col %in% names(tab),paste('Missing metric column:',lab_col))
  pdat <- spr_order(tab)
  pdat$z <- as.numeric(pdat[[lab_col]])
  spr_require(all(is.finite(pdat$z)),'Invalid heatmap values')
  pdat$label <- if (mode=='BA') sprintf('%.1f',pdat$z) else sprintf('%+.1f',pdat$z)
  if (mode=='BA') {
    pdat$label_colour <- ifelse(pdat$z>=56 & pdat$z<=82,'white','black')
    fill <- scale_fill_gradientn(
      colours=c('#440154','#3b528b','#21918c','#5ec962','#fde725'),
      limits=c(25,95),oob=scales::squish,
      name='Mean balanced\naccuracy (%)')
  } else {
    pdat$label_colour <- ifelse(abs(pdat$z)>=8.5,'white','black')
    fill <- scale_fill_gradientn(
      colours=c('#b2182b','#ef8a62','#f7f7f7','#67a9cf','#2166ac'),
      values=scales::rescale(c(-17,-8.5,0,8.5,17)),
      limits=c(-17,17),oob=scales::squish,
      name='Change in BA\n(percentage points)')
  }
  title <- switch(mode,
    BA='Cross-validated balanced accuracy across skeletal-representation feature spaces',
    delta='Effect of training-fold synthetic augmentation on balanced accuracy',
    level_delta='Effect of anatomical-resolution choice on balanced accuracy')
  subtitle <- switch(mode,
    BA='Mean of 15 held-out folds; three classification tasks and two anatomical resolutions',
    delta='Augmented minus original-only training; matched assessment sites and folds',
    level_delta='Performance-guided minus standard anatomical resolution; augmented experiment')
  # The level-difference figure is one column (3 tasks x 1 level comparison).
  if (mode=='level_delta') {
    panel <- facet_grid(rows=vars(task_lab),switch='y')
  } else {
    panel <- facet_grid(rows=vars(task_lab),cols=vars(level_lab),switch='y')
  }
  plot <- ggplot(pdat,aes(x=space_lab,y=model_lab,fill=z)) +
    geom_tile(colour='white',linewidth=1.25) +
    geom_text(aes(label=label,colour=label_colour),size=3.0,
              fontface='plain',show.legend=FALSE) +
    scale_colour_identity() + fill + panel +
    scale_x_discrete(drop=FALSE,position='bottom') +
    scale_y_discrete(drop=FALSE) +
    labs(title=title,subtitle=subtitle,
      x=NULL,y=NULL,
      caption=if (mode=='BA')
        'Cells are means of 15 held-out fold scores; fold SD and other metrics are reported in the tables.'
      else 'Paired differences are descriptive; configurations reuse archaeological reference sites.') +
    spr_theme() +
    theme(axis.text.x=element_text(size=8.4,margin=margin(t=5)),
          strip.placement='outside',strip.text.y.left=element_text(angle=90),
          panel.spacing=unit(11,'pt'),
          legend.position='right',
          axis.ticks=element_blank(),
          panel.border=element_rect(colour='#dce2e8',fill=NA,linewidth=.3))
  spr_save_plot(plot,stem,output_dir,
                width=if (mode=='level_delta') 11.7 else 17.6,
                height=if (mode=='level_delta') 10.8 else 11.5)
}

spr_plot_tradeoff <- function(paired,output_dir) {
  d <- paired |>
    mutate(task_lab=factor(unname(spr_tasks[task]),levels=unname(spr_tasks)),
           level_lab=factor(unname(spr_levels[level]),levels=unname(spr_levels)),
           model_lab=factor(unname(spr_models[model]),levels=unname(spr_models)))
  pal <- c('Logistic regression'='#c06d3d','KNN'='#587399',
           'RBF-SVM'='#467d71','Random forest'='#8565a8')
  plot <- ggplot(d,aes(x=delta_BA_pp,y=delta_Brier,colour=model_lab,shape=level_lab)) +
    geom_hline(yintercept=0,linetype='dashed',colour='#63707b',linewidth=.45) +
    geom_vline(xintercept=0,linetype='dashed',colour='#63707b',linewidth=.45) +
    geom_point(size=2.25,alpha=.78) +
    facet_wrap(vars(task_lab),nrow=1,scales='free_y') +
    scale_colour_manual(values=pal,name='Classifier') +
    scale_shape_manual(values=c('Standard'=16,'Performance-guided'=17),name='Resolution') +
    labs(title='Effects of synthetic augmentation on discrimination and probability scoring',
       subtitle='Each point represents one feature-space × classifier × anatomical-resolution comparison',
       x='Change in balanced accuracy (percentage points)',
       y='Change in Brier score (augmented − original-only)',
       caption='Rightward change improves balanced accuracy; upward change worsens Brier score. Binary and multiclass Brier scales differ.') +
    spr_theme() +
    theme(legend.position='bottom',panel.grid.major=element_line(colour='#e5e9ee',linewidth=.25))
  spr_save_plot(plot,'Supplement_Figure_Paired_BA_Brier_Change',output_dir,16.5,6)
}

spr_confusion_data <- function(oof,condition,task,level,space,model) {
  classes <- switch(task,
    '4class'=c('Funerary Context','Cannibalism','Carnivore','Geological transport'),
    'cannibalism_vs_fc'=c('Funerary Context','Cannibalism'),
    'fc_vs_rest'=c('Funerary Context','rest'))
  cols <- paste0('P_',gsub(' ','.',classes,fixed=TRUE))
  spr_require(all(cols %in% names(oof)),
              paste('Missing probabilities:',paste(setdiff(cols,names(oof)),collapse=',')))
  d <- oof |>
    filter(.data$task==task,.data$level==level,
           .data$space==space,.data$model==model)
  spr_require(nrow(d)>0L,'No OOF predictions for the requested confusion panel')
  d <- d |>
    group_by(site_id) |>
    summarise(truth=first(truth),n_assess=n(),
       across(all_of(cols),mean), .groups='drop')
  spr_require(all(d$n_assess==5L),
              'Every site must have five held-out predictions for the confusion display')
  probs <- as.matrix(d[,cols,drop=FALSE])
  spr_require(all(is.finite(probs)),'Missing class probabilities in OOF file')
  d$predicted <- classes[max.col(probs,ties.method='first')]
  count <- as.data.frame(table(
    truth=factor(d$truth,levels=classes),
    predicted=factor(d$predicted,levels=classes)),stringsAsFactors=FALSE)
  names(count)[3] <- 'n_sites'
  count <- count |>
    group_by(truth) |>
    mutate(row_percentage=100*n_sites/sum(n_sites)) |>
    ungroup() |>
    mutate(task=task,level=level,space=space,model=model,
           condition=condition)
  count
}

spr_plot_confusions <- function(original_oof, augmented_oof,output_dir) {
  level <- names(spr_levels)[1]
  selected <- data.frame(
    task=names(spr_tasks),
    space=c('abundance_evenness','pmnau','abundance'),
    model=c('RandomForest','SVM_RBF','SVM_RBF'),
    stringsAsFactors=FALSE)
  panels <- list()
  for (i in seq_len(nrow(selected))) {
    for (j in 1:2) {
      condition <- c('Original-only','Augmented')[[j]]
      d <- if(j==1) original_oof else augmented_oof
      panels[[length(panels)+1L]] <- spr_confusion_data(
        d,condition,selected$task[[i]],level,
        selected$space[[i]],selected$model[[i]])
    }
  }
  dat <- bind_rows(panels)
  write_csv(dat,file.path(output_dir,'Supplement_Table_Confusion_SiteAveraged.csv'))
  short <- c('Funerary Context'='Funerary','Cannibalism'='Cannibalism',
             'Carnivore'='Carnivore','Geological transport'='Geological',
             'rest'='Rest')
  dat <- dat |> mutate(
    task_name=unname(spr_tasks[task]),
    panel=factor(paste(condition,task_name,sep=' • '),
        levels=c(paste('Original-only',unname(spr_tasks),sep=' • '),
                 paste('Augmented',unname(spr_tasks),sep=' • '))),
    truth_name=unname(short[as.character(truth)]),
    predicted_name=unname(short[as.character(predicted)]),
    cell_label=sprintf('%d\n%.1f%%',n_sites,row_percentage)
  )
  # Facet-specific class labels; drop absent classes from the binary task panels.
  dat$truth_name <- factor(dat$truth_name,levels=rev(unname(short)))
  dat$predicted_name <- factor(dat$predicted_name,levels=unname(short))
  plot <- ggplot(dat,aes(x=predicted_name,y=truth_name,fill=row_percentage)) +
    geom_tile(colour='white',linewidth=1.1) +
    geom_text(aes(label=cell_label,colour=ifelse(row_percentage>55,'white','black')),
              show.legend=FALSE,size=3.2,lineheight=.95) +
    scale_colour_identity() +
    scale_fill_gradientn(colours=c('#f7fbff','#6baed6','#08519c'),
                         limits=c(0,100),name='Within-class\npercentage') +
    facet_wrap(vars(panel),nrow=2,scales='free') +
    scale_x_discrete(drop=TRUE) + scale_y_discrete(drop=TRUE) +
    labs(title='Site-level confusion matrices from repeated out-of-fold predictions',
      subtitle='Standard anatomical resolution; identical classifier/feature configurations in both conditions',
      x='Predicted reference category',y='Observed reference category',
      caption='Each archaeological site contributes once per panel; predictions use its mean probability across five held-out assessments.') +
    spr_theme(base_size=11) +
    theme(axis.text.x=element_text(angle=32,hjust=1),
          panel.spacing=unit(19,'pt'),legend.position='right',
          panel.border=element_rect(fill=NA,colour='#dce2e8',linewidth=.4))
  spr_save_plot(plot,'Supplement_Figure_Confusion_Original_vs_Augmented',
                output_dir,16.5,12.2)
}

spr_level_difference <- function(augmented_summary,output_dir) {
  standard <- names(spr_levels)[1]; guided <- names(spr_levels)[2]
  a <- augmented_summary |>
    filter(level==standard) |>
    select(task,space,model,BA_standard=BA)
  b <- augmented_summary |>
    filter(level==guided) |>
    select(task,space,model,BA_guided=BA)
  d <- inner_join(a,b,by=c('task','space','model')) |>
    mutate(level=standard,delta_level_BA_pp=BA_guided-BA_standard)
  spr_require(nrow(d)==84L,'Expected 84 paired model/feature/task comparisons between levels')
  write_csv(d,file.path(output_dir,'Augmented_GuidedMinusStandard_BA.csv'))
  spr_plot_heatmap(d,'Supplement_Figure_GuidedMinusStandard_Augmented_BA',
                   output_dir,mode='level_delta')
}

spr_generate_figures <- function(original_root,augmented_root,output_dir,
                                 figures=c('augmented','original','augmentation_delta',
                                           'tradeoff','confusion','level_delta')) {
  allowed <- c('augmented','original','augmentation_delta',
               'tradeoff','confusion','level_delta')
  spr_require(all(figures %in% allowed),paste('Unknown requested figure:',
            paste(setdiff(figures,allowed),collapse=', ')))
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  model_root <- dirname(normalizePath(output_dir,winslash='/',mustWork=FALSE))
  orig_dir <- spr_resolve_run(model_root, original_root, 'Original-only')
  aug_dir <- spr_resolve_run(model_root, augmented_root, 'Augmented')
  spr_require(orig_dir != aug_dir,
              'Original-only and augmented runs resolve to the same directory')
  orig_dir <- spr_ensure_pooled_results(orig_dir,'Original-only',output_dir)
  aug_dir <- spr_ensure_pooled_results(aug_dir,'Augmented',output_dir)
  orig <- spr_load_experiment(orig_dir,'Original-only',
                              load_predictions='confusion' %in% figures)
  aug <- spr_load_experiment(aug_dir,'Augmented',
                             load_predictions='confusion' %in% figures)
  result <- spr_pair_results(orig,aug,output_dir)
  if ('augmented' %in% figures)
    spr_plot_heatmap(aug$summary,
      'Figure_ML_BalancedAccuracy_Augmented_3Tasks_2Levels',output_dir,'BA')
  if ('original' %in% figures)
    spr_plot_heatmap(orig$summary,
      'Supplement_Figure_BalancedAccuracy_OriginalOnly',output_dir,'BA')
  if ('augmentation_delta' %in% figures)
    spr_plot_heatmap(result$paired,
      'Main_Figure_Augmentation_Delta_BalancedAccuracy',output_dir,'delta')
  if ('tradeoff' %in% figures)
    spr_plot_tradeoff(result$paired,output_dir)
  if ('confusion' %in% figures)
    spr_plot_confusions(orig$oof,aug$oof,output_dir)
  if ('level_delta' %in% figures)
    spr_level_difference(aug$summary,output_dir)
  # This manifest allows figure creation to be reproduced from exact input files.
  writeLines(c(
    paste0('Generated: ',format(Sys.time(),'%Y-%m-%d %H:%M:%S %z')),
    paste0('Original-only fold source: ',orig$file),
    paste0('Original-only MD5: ',unname(tools::md5sum(orig$file))),
    paste0('Augmented fold source: ',aug$file),
    paste0('Augmented MD5: ',unname(tools::md5sum(aug$file))),
    paste0('Figures: ',paste(figures,collapse=', ')),
    '', capture.output(sessionInfo())
  ),con=file.path(output_dir,'publication_figures_reproducibility.txt'))
  message('Complete: ',normalizePath(output_dir,winslash='/',mustWork=TRUE))
  invisible(result)
}

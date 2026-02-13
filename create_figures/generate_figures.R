####################################################################
#
# Create plots
#
# data_path: path to RData file containing v (output from effect_size master() function)
# output_dir: directory to save plots
# plot_output_style: 'shiny' or 'manuscript'
# - shiny: save plots for the shiny app
# - manuscript: save plots concatenated across all studies for manuscript figures
# all_effect_size_types: vector of effect size types to plot (default: c('d', 'r_sq'))
# all_motion: vector of motion correction types to plot (default: c('none', 'regression', 'threshold', 'threshold2'))
# all_pooling: vector of pooling types to plot (default: c('net', 'none'))
# - net: network-based pooling
# - none: no pooling (voxel-wise or edge-wise)
# plot_multi: whether to plot multivariate results (default: FALSE)
# all_plot_combination_styles: vector of plot combination styles to plot (default: c('single', 'meta'))
# - single: individual study plots
# - meta: meta-analysis plots 
# all_grouping_var: vector of grouping variables for meta-analysis (default: c('category'))
# all_manuscript_plot_types: vector of plot types to use for manuscript (default: c('simci')
# - other options: 'spatial', 'spatial_pow_thr', 'spatial_pow_n_thr', 'density', 'density_binned', 'power', 'power_n','power_binned','power_n_binned'
# make_plots: whether to make plots (default: TRUE)
# save_plots: whether to save (export) plots (default: FALSE)


####################################################################

generate_figures <- function(data_path, output_dir, plot_output_style, all_effect_size_types = c('d', 'r_sq'), all_motion = c('none', 'regression', 'threshold', 'threshold2'), 
  all_pooling = c('net', 'none'), plot_multi = FALSE, all_plot_combination_styles = c('single', 'meta'), all_grouping_var = c('category'), all_manuscript_plot_types = c('simci'),
  make_plots = TRUE, save_plots = TRUE) {

  ## Load libraries
  library(ggpubr) # requires svglite
  library(devtools) # for installing/loading utils
  library(BrainEffeX.utils)

  # load data (v)
  # if(!exists("v")) {
  #   load(data_path)
  # }
  load(data_path)
  
  # reorder_by_grouping_var <- TRUE # for single plots - force to always be true
  grouping_var_order <- NULL
  grouping_var_order$orig_stat_type <- c('r','t2','t')
  grouping_var_order$category <- c('cognitive','cognitive (task)','psychiatric','biometric','sex (demographic)','age (demographic)')

  # special settings for shiny vs. manuscript

  if (plot_output_style == 'shiny') {
    all_plot_types <- c('simci-spatial')
    if ('overlapping' %in% all_plot_combination_styles) {
      warning("Overlapping plots not supported in shiny mode. Removing.")
      all_plot_combination_styles <- all_plot_combination_styles[all_plot_combination_styles != 'overlapping']
    }
    # if (save_logs == TRUE) {
    # warning("Logs not used in shiny. Setting save_logs to FALSE.")
    save_logs <- FALSE
    # }
    add_plt_description <- TRUE # text at bottom of screen
    use_minimal_title <- FALSE
    
  } else if (plot_output_style == 'manuscript') {
    all_plot_types <- all_manuscript_plot_types
    all_plot_combination_styles <- all_plot_combination_styles
    save_logs = TRUE
    add_plt_description <- FALSE # text at bottom of screen
    use_minimal_title <- TRUE # cleaner title for manuscript
  }

  ## Loop over plot types and styles

  for (pooling in all_pooling) {
    for (motion in all_motion) {
      for (grouping_var in all_grouping_var) {
        for (effect_size_type in all_effect_size_types) {
          for (plot_combination_style in all_plot_combination_styles) {
            for (plot_type in all_plot_types) {
              
              print(paste0('Doing plot_combination_style: ', plot_combination_style, ' | plot_type: ', plot_type, ' | pooling: ', pooling, ' | motion: ', motion, ' | grouping var: ', grouping_var, ' | effect size: ', effect_size_type))
              
              
              ## Set up strings
              
              combo_name <- paste0('pooling.', pooling, '.motion.', motion, '.mv.none')
              mv_combo_basename <- paste0('pooling.', pooling, '.motion.', motion, '.mv.multi')
              
              ## Run meta-analysis, if specified
              
              # TODO: move to combine_gl with other more intensive processing / stat estimates & run all relevantmeta beforehand
              run_meta <- FALSE
              try_meta_file <- FALSE
              
              if (plot_combination_style == 'meta') {
                
                # setup string and filename
                meta_str <- paste0('meta_',grouping_var)
                meta_fn_dir <- system.file("meta/", package = "BrainEffeX.utils") # TODO: set this somewhere else
                meta_fn <- file.path(meta_fn_dir, "v.RData")
                
                # check if this var contains the meta for this grouping var / combo
                if (!(meta_str %in% names(v))) { # check if this grouping var already exists in data
                  try_meta_file <- TRUE
                } else if (!(combo_name %in% names(v[[meta_str]]$data[[1]]))) { # check if this combo_name has been run
                  try_meta_file <- TRUE
                }
              }
              
            
              if (make_plots) {
                
                ## Set up unique identifiers for each plot
                
                plot_info__idx <- list() # each row = list of study(s) in data to include in each plot
                # for single plots: each row = 1 entry per study to index into v$data
                # for meta-analysis plots: each row = 1 entry per category to index v[[meta_str]]$data
                # for overlapping plots: each row = list of indices per group (x map type) to index into v$data
                
                plot_info__grouping_var <- list() # each row = grouping variable (same value repeated for each plot)
                plot_info__group_level <- list() # each row = level within grouping variable
                plot_info__ref <- list() # each row = ref(s) used for a study or grouping variable
                
                if (plot_combination_style == 'single') {  # name by study
                  
                  all_study_names <- names(v$data)
                  
                  for (i in 1:length(v$data)) {
                    plot_info__idx[[all_study_names[[i]]]] <- i
                    plot_info__grouping_var[[all_study_names[[i]]]] <- "none"  # overwrite any other grouping var if doing single plots
                    plot_info__group_level[[all_study_names[[i]]]] <- NA
                    plot_info__ref[[all_study_names[[i]]]] <- v$study$ref[i]
                  }
                  cat("Warning: grouping_var set to 'none' for single plots\n")
                  
                  
                  # sort all rows of plot_info__idx by specified group order (studies in the same group next to each other)
                  # if (reorder_by_grouping_var) {
                  idx_reorder <- NULL
                  for (j in 1:length(grouping_var_order[[grouping_var]])) {
                    idx_reorder <- c(idx_reorder, which(v$study[[grouping_var]] == grouping_var_order[[grouping_var]][j]))
                  }
                  # idx_reorder <- c(which(v$study$orig_stat_type == 'r'), which(v$study$orig_stat_type == 't2'), which(v$study$orig_stat_type == 't'))
                  plot_info__idx <- plot_info__idx[idx_reorder]
                  plot_info__grouping_var <- plot_info__grouping_var[idx_reorder]
                  plot_info__group_level <- plot_info__group_level[idx_reorder]
                  plot_info__ref <- plot_info__ref[idx_reorder]
                  # }
                  
                  
                } else if (plot_combination_style == 'meta') { # name by average of grouping var
                  
                  for (i in 1:length(v[[meta_str]]$data)) {
                    plot_info__idx[[names(v[[meta_str]]$data)[[i]]]] <- i
                    plot_info__grouping_var[[names(v[[meta_str]]$data)[[i]]]] <- grouping_var
                    plot_info__group_level[[names(v[[meta_str]]$data)[[i]]]] <- v[[meta_str]]$study$group_level[i]
                    plot_info__ref[[names(v[[meta_str]]$data)[[i]]]] <- v[[meta_str]]$study$ref[i]
                  }
                  
                  # sort all rows of plot_info__idx by specified group order (studies in the same group next to each other)
                  idx_reorder <- order(match(v[[meta_str]]$study$group_level,grouping_var_order[[grouping_var]]))
                  plot_info__idx <- plot_info__idx[idx_reorder]
                  plot_info__grouping_var <- plot_info__grouping_var[idx_reorder]
                  plot_info__group_level <- plot_info__group_level[idx_reorder]
                  plot_info__ref <- plot_info__ref[idx_reorder]
                  
                  
                } else if (plot_combination_style == 'overlapping') { # overlapping individual plots
                  
                  if (grouping_var == 'category') {
                    study_group_name <- v$study$category
                  } else if (grouping_var == 'orig_stat_type') {
                    study_group_name <- v$study$orig_stat_type
                  }
                  
                  # sort existing groups by specified group order
                  all_group_names <- grouping_var_order[[grouping_var]][grouping_var_order[[grouping_var]] %in% unique(study_group_name)]
                  all_map_types <- unique(v$study$map_type)
                  
                  for (this_group_name in all_group_names) {
                    for (this_map_type in all_map_types) {
                      idx <- which(study_group_name == this_group_name & v$study$map_type == this_map_type)
                      plot_info__idx[[paste0(this_group_name, '.', this_map_type)]] <- idx
                      plot_info__grouping_var[[paste0(this_group_name, '.', this_map_type)]] <- grouping_var
                      plot_info__group_level[[paste0(this_group_name, '.', this_map_type)]] <- this_group_name
                      plot_info__ref[[paste0(this_group_name, '.', this_map_type)]] <- unique(v$study$ref[idx])
                    }
                  }
                }
                
                plot_info <- data.frame(
                  idx = I(plot_info__idx),
                  grouping_var = unlist(plot_info__grouping_var),
                  group_level = unlist(plot_info__group_level),
                  ref = I(plot_info__ref),
                  row.names = names(plot_info__idx),
                  stringsAsFactors = FALSE
                )
                rm(plot_info__idx, plot_info__grouping_var, plot_info__group_level, plot_info__ref)
                ## Print plot info
                
                
                ## Make Plots
                
                panel_list <- list() # list of panels
                panel_list_2 <- list() # list of panels
                
                panel_list <-list()
                log_list <- list() # list of logs
                
                for (i in 1:length(plot_info$idx)) { # loop over panels - this_study_or_group is the name of the group or study
                  
                  this_study_or_group <- rownames(plot_info)[i]
                  this_plot_info <- plot_info[this_study_or_group,]
                  
                  pd_list <- list() # list of plot info for single panel
                  pd_list_2 <- list() # add'l list of plot info for plot type to add to pd_list
                  ld_list <- list() # list of log info for single panel
                  
                  n_studies_in_pd_list <- 1
                  
                  # 1. Prep
                  
                  for (j in plot_info$idx[[i]]) {
                    
                    # change metadata based on whether using meta-analysis
                    
                    if (plot_combination_style == 'meta') {
                      
                      data <- v[[meta_str]]$data[[j]]
                      study_details <- list()
                      brain_masks <- v[[meta_str]]$brain_masks[[j]]$pooling.none.motion.none.mv.none # TODO: this is because we explicitly set this for meta but not for single studies - assuming motion type shouldn't affect the mask and always using an external mask for pooling
                      
                    } else {
                      
                      data <- v$data[[j]]
                      study_details <- v$study[j, ]
                      brain_masks <- v$brain_masks[[j]]
                    }
                    
                    if (plot_multi) {
                      mv_combo_name <- names(v$data[[j]])[grepl(mv_combo_basename,names(v$data[[j]]))]
                      combo_name <- mv_combo_name
                    }
                    
                    if (combo_name %in% names(data)) { # if combo_name exists in data (e.g., not all studies have net)
                      if (any(!is.na(data[[combo_name]][[effect_size_type]])) > 0) {  # data is not just NA
                        
                        # prep
                        
                        
                        if (plot_type == 'simci-spatial') { # shiny
                          
                          pd_list[[n_studies_in_pd_list]] <- prep_data_for_plot(data = data, study_details = study_details, combo_name = combo_name, mv_combo_name = mv_combo_basename, estimate = effect_size_type, plot_info = this_plot_info)
                          # pd_list_2[[n_studies_in_pd_list]] <- prep_data_for_spatial_plot(data = data, brain_masks = brain_masks, study_details = study_details, combo_name = combo_name, mv_combo_name = mv_combo_basename, estimate = effect_size_type, plot_info = this_plot_info)
                          pd_list_2[[n_studies_in_pd_list]] <- prep_data_for_plot(data = data, study_details = study_details, combo_name = combo_name, mv_combo_name = mv_combo_basename, estimate = effect_size_type, plot_info = this_plot_info, prep_spatial = TRUE, brain_masks = brain_masks)
                          
                        } else if (grepl('spatial', plot_type)) {
                          
                          # TODO: finally combined prep for spatial with general prep. test more thoroughly for orig spatial and remove
                          # pd_list[[n_studies_in_pd_list]] <- prep_data_for_spatial_plot(data = data, brain_masks = brain_masks, study_details = study_details, combo_name = combo_name, mv_combo_name = mv_combo_basename, estimate = effect_size_type, plot_info = this_plot_info)
                          pd_list[[n_studies_in_pd_list]] <- prep_data_for_plot(data = data, study_details = study_details, combo_name = combo_name, mv_combo_name = mv_combo_basename, estimate = effect_size_type, plot_info = this_plot_info, prep_spatial = TRUE, brain_masks = brain_masks)
                          
                        } else {
                          
                          pd_list[[n_studies_in_pd_list]] <- prep_data_for_plot(data = data, study_details = study_details, combo_name = combo_name, mv_combo_name = mv_combo_basename, estimate = effect_size_type, plot_info = this_plot_info)
                          
                        }
                        
                        ld_list[[n_studies_in_pd_list]] <- get_summary_info(pd_list[[n_studies_in_pd_list]]$study_details, pd_list[[n_studies_in_pd_list]]$extra_study_details)
                        n_studies_in_pd_list <- n_studies_in_pd_list + 1
                        
                      }
                    }
                  }
                  
                  
                  # 2. Plot & Log
                  
                  if (length(pd_list) > 0) { # plot only if pd_list isn't empty
                    
                    # set up plot
                    if (length(ld_list) > 1) {
                      log_list[[i]] <- combine_summary_info(ld_list)
                    } else {
                      log_list[[i]] <- ld_list[[1]]
                    }
                    
                    if (plot_type == 'simci-spatial') {
                      
                      panel_list[[i]] <- create_plots(pd_list, plot_type = 'simci', effect_type = effect_size_type, do_multivariate = plot_multi, add_description = add_plt_description, do_minimal_title = use_minimal_title, log_list[[i]])
                      panel_list_2[[i]] <- create_plots(pd_list_2, plot_type = 'spatial', effect_type = effect_size_type, do_multivariate = plot_multi,add_description = add_plt_description, do_minimal_title = use_minimal_title, log_list[[i]])
                      
                    } else {
                      
                      panel_list[[i]] <- create_plots(pd_list, plot_type = plot_type, effect_type = effect_size_type, do_multivariate = plot_multi,add_description = add_plt_description, do_minimal_title = use_minimal_title, log_list[[i]])
                      
                    }
                    
                  }
                  
                }
                
                
                
                # General plot parameters
                # TODO: figure out what we want to set up here vs. to pass or set up in
                # create_plots, which gets passed to plot_sim_ci, etc.
                # Should at least set all panel / canvas dimensions here
                
                pp <- list()
                pp$width_per_panel <- 7 
                pp$height_per_panel <- 6 
                pp$res <- 100
                pp$units <- "in"
                pp$title_size <- 20
                if (plot_type == 'simci-spatial') { # shiny
                  pp$ncol <- 2
                  pp$nrow <- 1
                } else {
                  # remove nulls
                  panel_list <- panel_list[!sapply(panel_list,is.null)]
                  if (plot_type == 'power' | plot_type == 'power_n') {
                    pp$ncol <- length(panel_list)
                    pp$nrow <- 1
                  } else {
                    pp$ncol <- 1
                    pp$nrow <- length(panel_list)
                  }
                }
                if (plot_output_style == 'manuscript') {
                  do_minimal_title <- TRUE
                } else {
                  do_minimal_title <- FALSE
                }
                
                # prep: set up dir and file names
                
                if (save_plots) {
                  if (plot_combination_style == 'meta') {
                    grouping_var_str <- paste0('_', grouping_var)
                  } else {
                    grouping_var_str <- ''
                  }
                  
                  if (plot_multi) {
                    multi_str <- 'mv_'
                  } else {
                    multi_str <- ''
                  }
                  
                  if (plot_output_style == 'manuscript') {
                    grouping_var_str <- paste0('_', plot_type)
                  } else { 
                    grouping_var_str <- paste0('_', plot_type, '_', effect_size_type)
                  }
                  
                  out_basename <- paste0(output_dir, plot_output_style, '/', effect_size_type, '/motion_', motion, '/pooling_', pooling, '/', multi_str, plot_combination_style, grouping_var_str)
                  if (plot_type == 'simci-spatial') { # use out_basename as dir, otherwise use as basename for concat plots
                    out_basename <- paste0(out_basename,'/')
                  }
                  actual_dir <- sub("/[^/]*$", "", out_basename)
                  if (!dir.exists(actual_dir)) {
                    dir.create(actual_dir, recursive = TRUE)
                  }
                }
                
                
                # plot multi panels
                
                if (plot_type == 'simci-spatial') {
                  
                  for (i in 1:length(panel_list)) {
                    
                    t <- list(panel_list[[i]], panel_list_2[[i]])
                    t_master_title <- t[[1]]$labels$title
                    t[[1]]$labels$title <- ""
                    t[[2]]$labels$title <- ""
                    
                    multi_plot <- ggarrange(plotlist = t, ncol = pp$ncol, nrow = pp$nrow)
                    multi_plot <- annotate_figure(multi_plot,
                                                  top = text_grob(t_master_title, face = "bold", size = pp$title_size))
                    
                    if (save_plots) {
                      
                      study_name <- tolower(names(plot_info$idx)[[i]])
                      
                      plot_fn <- paste0(out_basename, study_name, '.png')
                      ggsave(plot_fn, plot = multi_plot, width = pp$width_per_panel * pp$ncol, height = pp$height_per_panel * pp$nrow, units = pp$units, dpi = pp$res, bg = "white", device = "png")
                      
                      if (save_logs) { # TODO: some logs are identical - see about saving only single
                        log_fn <- paste0(out_basename, study_name,'.txt')
                        writeLines(unlist(lapply(log_list, function(x) c(x$title_text, x$bottom_text, ""))), log_fn)
                      }
                      
                    }
                    
                  }
                  
                } else {
                  
                  multi_plot <- ggarrange(plotlist = panel_list, ncol=pp$ncol, nrow=pp$nrow)
                  
                  if (save_plots) {
                    
                    plot_fn <- tolower(paste0(out_basename, '.png'))
                    
                    cat("Saving plots to...\n", plot_fn, "\n", sep = "")
                    
                    ggsave(plot_fn, plot = multi_plot, width = pp$width_per_panel * pp$ncol, height = pp$height_per_panel * pp$nrow, units = pp$units, dpi = pp$res, bg = "white", device = "png", limitsize = FALSE)
                    
                    if (save_logs) { # TODO: some logs are identical - see about saving only single
                      log_fn <- paste0(out_basename, '.txt')
                      writeLines(unlist(lapply(log_list, function(x) c(x$title_text, x$bottom_text, ""))), log_fn)
                    }
                    
                  }
                  
                }
                
              }
              
              
              ## close loop over plot types and styles
            } # effect_size_type
          } # pooling
        } # plot_combination_style
      } # plot_type
    } # motion
  } # grouping_var
  }
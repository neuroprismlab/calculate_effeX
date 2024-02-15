# plot histograms, CI ribbons, and make logs

library(ggridges) # ridgeplot
library(ggplot2)
library(scales) # to get those ggplot colors
library(ggcorrplot)
library(dplyr)
#library(viridis)
#library(hrbrthemes)
library(knitr) # saving tables
library(kableExtra) # saving tables

# Sort by magnitude of d - both for plot and for imputation
d_long__sorted <- d_long
d_long__sorted <- arrange(d_long__sorted, study, d_long__by_study)
d_long__sorted$idx_within_study <- ave(d_long__sorted$d_long__by_study,
                    d_long__sorted$study,
                    FUN = seq_along)

# Impute missing CI data (Note: requires CI to be sorted first for locf)
# Rationale: sometimes d.ci has difficulty estimating CI when alpha and d are both very small, especially small alpha (e.g., d.ci(0.01,100,100,alpha=1e-6) works but d.ci(0.01,100,100,alpha=1e-7) yields NA) - sometimes it works again close to d=0
# TODO: move this all this above when CI's are first being defined. Important: will need to sort, which is currently not done above and certainly not in the list of list. Will prob help to have spatial IDs (e.g., this is edge #1 according to the Shen atlas)

# 1. Replace any NA in ci_ub with d+(d-ci_lb) and in ci_lb with d-(ci_ub-d) (both 2*d-ci_xx) - TODO: move earlier, right after calculated
d_long__sorted$ci_ub[is.na(d_long__sorted$ci_ub)] <- 2*d_long__sorted$d_long__by_study[is.na(d_long__sorted$ci_ub)] - d_long__sorted$ci_lb[is.na(d_long__sorted$ci_ub)]
d_long__sorted$ci_lb[is.na(d_long__sorted$ci_lb)] <- 2*d_long__sorted$d_long__by_study[is.na(d_long__sorted$ci_lb)] - d_long__sorted$ci_ub[is.na(d_long__sorted$ci_lb)]
d_long__sorted$ci_height <- na.locf(d_long__sorted$ci_ub-d_long__sorted$ci_lb)
d_long__sorted$ci_ub[is.na(d_long__sorted$ci_ub)] <- d_long__sorted$d_long__by_study[is.na(d_long__sorted$ci_ub)] + d_long__sorted$ci_height[is.na(d_long__sorted$ci_ub)]/2
d_long__sorted$ci_lb[is.na(d_long__sorted$ci_lb)] <- d_long__sorted$d_long__by_study[is.na(d_long__sorted$ci_lb)] - d_long__sorted$ci_height[is.na(d_long__sorted$ci_lb)]/2

# 2. Replace any NA for both ci_lb and ci_ub with last observation carried forward (na.locf) - this is a good conservative estimate since CI height ~ d, and these d's to be filled in are all closer to 0, so actual CI height would be narrower
d_long__sorted$ci_contains_zero <- ifelse(d_long__sorted$ci_lb<0 & d_long__sorted$ci_ub>0,TRUE,FALSE)
d_long__sorted$ci_contains_zero__group <- rleid(d_long__sorted$ci_contains_zero)
d_long__sorted$smallest_ci_bound <- ifelse(abs(d_long__sorted$ci_lb)<abs(d_long__sorted$ci_ub),d_long__sorted$ci_lb,d_long__sorted$ci_ub)
d_long__sorted$most_conservative_effect <- (!d_long__sorted$ci_contains_zero) * d_long__sorted$smallest_ci_bound



######## Descriptive Logs ########

# TODO: add option to save
if (do_log) {

  ## Basic study info
  kable(study_info[2:8], row.names=F) %>%
    kable_styling(bootstrap_options = "striped") %>%
    kable_styling(full_width = F) %>%
    save_kable(file = study_dataset_info_filename) # TODO: save should be an option

  ## Summary of most conservative effects
  ## Max most_conservative_effect, grouped by statistic

  # Percent of edges for which ci_contains_zero==1, grouped by statistic
  perc_ci_contains_zero__by_statistic <- setNames(aggregate(data = d_long__sorted, ci_contains_zero ~ statistic, FUN = function(x) 100*mean(x)), c('statistic','percent_ci_contains_zero'))
  max_conservative_effect__by_statistic <- setNames(aggregate(data = d_long__sorted, most_conservative_effect ~ statistic, FUN = max), c('statistic','max_conservative_effect'))
  conservative_effect_summary__by_statistic <- merge(perc_ci_contains_zero__by_statistic, max_conservative_effect__by_statistic, by='statistic')

  ## power for each row of most_conservative_effect__by_statistic
  library("pwr")
  conservative_effect_summary__by_statistic$power_n100__max_effect_by_statistic<-list()
  conservative_effect_summary__by_statistic$n_needed__max_effect_by_statistic<-list()
  for (i in 1:dim(max_conservative_effect__by_statistic)[1]) {
    if (max_conservative_effect__by_statistic$statistic[i]=='r' | max_conservative_effect__by_statistic$statistic[i]=='t2') {
      this_type<-"two.sample"
    } else {
      this_type<-"paired"
    }
    conservative_effect_summary__by_statistic$power_n100__max_effect_by_statistic[[i]] <- pwr.t.test(d=max_conservative_effect__by_statistic$max_conservative_effect[i],sig.level=0.05,n=100,type=this_type)$power
    conservative_effect_summary__by_statistic$n_needed__max_effect_by_statistic[[i]] <- pwr.t.test(d=max_conservative_effect__by_statistic$max_conservative_effect[i],sig.level=0.05,power=0.8,type=this_type)$n
  }
  # power__by_statistic <- matrix(NA, nrow = dim(max_conservative_effect__by_statistic)[1], ncol = 1)
  # for (i in 1:dim(max_conservative_effect__by_statistic)[1]) {
  #   power__by_statistic[i,1] <- pwr.t.test(d=max_conservative_effect__by_statistic$most_conservative_effect[i],sig.level=0.05,power=0.8,type="two.sample")$power
  # }

  # reorder rows c("r","t2","t","t(act)")
  conservative_effect_summary__by_statistic$statistic <- factor(conservative_effect_summary__by_statistic$statistic, levels = c("r","t2","t","t (act)"))
  conservative_effect_summary__by_statistic__sorted <- conservative_effect_summary__by_statistic[order(conservative_effect_summary__by_statistic$statistic),]
  sink(file = conservative_effect_summary_by_statistic_filename)
  print(conservative_effect_summary__by_statistic__sorted)
  sink(NULL)

  


  ## Percent above d-thresholds

  # Understand sample sizes for conventional small, medium, and large effect sizes
    # note: for 2-sample, n is number in *each* group
  conventional_d_thresh__n <- matrix(NA, nrow = length(conventional_d_thresh), ncol = length(ttest_types)) 
  for (i in 1:dim(conventional_d_thresh__n)[1]) { # large, medium, and small sample sizes (corresponds with expected small, medium, and large effect sizes)
    for (j in 1:dim(conventional_d_thresh__n)[2]) { # one-sample and two-sample t-tests
      conventional_d_thresh__n[i,j]<-pwr.t.test(d=conventional_d_thresh[i],sig.level=0.05,power=0.8,type=ttest_types[j])$n
    }
  }
  kable(conventional_d_thresh__n, row.names=F) %>%
    kable_styling(bootstrap_options = "striped") %>%
    kable_styling(full_width = F) %>%
    save_kable(file = n_for_conventional_dthresh_filename) # TODO: save should be an option

  # Record percent above d-thresholds
  d_perc_gt_thresh1__by_study <- aggregate(d_long__by_study ~ statistic, data = d_long, FUN = function(x) 100*mean(abs(x) > conventional_d_thresh[1]))
  d_perc_gt_thresh2__by_study <- aggregate(d_long__by_study ~ statistic, data = d_long, FUN = function(x) 100*mean(abs(x) > conventional_d_thresh[2]))
  names(d_perc_gt_thresh1__by_study)[2] <- "percent_gt_dthresh"
  names(d_perc_gt_thresh2__by_study)[2] <- "percent_gt_dthresh"

  kable(d_perc_gt_thresh1__by_study, row.names=F) %>%
    kable_styling(bootstrap_options = "striped") %>%
    kable_styling(full_width = F) %>%
    save_kable(file = d_perc_gt_thresh1_filename) # TODO: save should be an option

  kable(d_perc_gt_thresh2__by_study, row.names=F) %>%
    kable_styling(bootstrap_options = "striped") %>%
    kable_styling(full_width = F) %>%
    save_kable(file = d_perc_gt_thresh2_filename) # TODO: save should be an option

  # stop("pausing here")
  ## Model fits and parameters
  # mdl__full_filename <- paste0(results_dir, 'mdl__full', skip_str__vis, ".txt")
  # # write.csv(summary(mdl__full_filename), file=mdl__full_filename)
  # sink(file = mdl__full_filename)
  # print(summary(mdl__full))
  # sink(NULL)

  sink(file = mdl__by_stat_filename)
  print(summary(mdl__by_stat))
  sink(NULL)

  sink(file = mdl__by_stat__within_pairwise_filename)
  print(summary(mdl__by_stat__within_pairwise))
  sink(NULL)
  
  sink(file = mdl__by_stat__between_pairwise_filename)
  print(summary(mdl__by_stat__between_pairwise))
  sink(NULL)

  
  
  
  
  # d_gt_thresh1 <- abs(d_long$d_long__by_study) > conventional_d_thresh[1]
  

  # create d thresholds for small-, medium-, and large-N corresponding effect sizes
  # results: One-sample: for small d: n=198; medium: n=33; large: n=14
  # results: Two-sample: for small d: n=393; medium: n=64; large: n=26 (n here is PER GROUP)
  # my_large_n <- 2000
  # d_thresh <- matrix(NA, nrow = 3, ncol = 2)
  # for (i in 1:3) { # large, medium, and small sample sizes (corresponds with expected small, medium, and large effect sizes)
  #   for (j in 1:2) { # one-sample and two-sample t-tests
  #     d_thresh[i,j]<-pwr.t.test(n=my_large_n/(10^(i-1)),sig.level=0.05,power=0.8,type=c("one.sample","two.sample")[j])$d
  #   }
  # }

  # effect size magnitudes
  # d_perc_gt_thresh__by_study <- aggregate(d_long__by_study ~ study, data = d_long, FUN = function(x) 100*mean(x > 0.2))
  # d_perc_gt_thresh__by_stat <- aggregate(d_long__by_study ~ statistic, data = d_long, FUN = function(x) 100*mean(x > 0.2))

  # for different thresholds per one- or two-sample t-test
  # d_perc_gt_thresh1_a <- aggregate(d_long__by_study ~ statistic, data = d_long, FUN = function(x) 100*mean(abs(x) > conventional_d_thresh[1,1]))
  # d_perc_gt_thresh1_b <- aggregate(d_long__by_study ~ statistic, data = d_long, FUN = function(x) 100*mean(abs(x) > conventional_d_thresh[1,2]))
  # d_perc_gt_thresh1 <- d_perc_gt_thresh1_a
  # d_perc_gt_thresh1$d_long__by_study[1] <- d_perc_gt_thresh1_b$d_long__by_study[1]
  # d_perc_gt_thresh1$d_long__by_study[4] <- d_perc_gt_thresh1_b$d_long__by_study[4]
  # names(d_perc_gt_thresh1)[2] <- "percent_gt_thresh"

  # d_perc_gt_custom_sm__by_stat <- aggregate(d_long__by_study ~ statistic, data = d_long, FUN = function(x) {
  #   if (unique(d_long$statistic) == "t") {
  #     100*mean(x > d_thresh[1,1])
  #   } else {
  #     100*mean(x > d_thresh[1,2])
  #   }
  # })
}



######## Visualization ########

if (make_figs) {

  ## Set up for plotting

  # Get unique colors for each dataset - note: will need to assign colors again if ever plot with a different order of datasets
  unique_cols_by_dataset <- hue_pal()(length(unique(study_info$dataset))) # for hist and forest
  unique_cols_by_stat<- hue_pal(l=100,h = c(50, 150))(4) # for hist #hue_pal(l = 90)(9) is lightness 90 and 9 squares
  unique_cols_by_stat<- unique_cols_by_stat[c(2,3,3,2)] # weirdly only this histogram switches colors (1,4,2,3), so let's reorder - goes 1st, then 4th, then 2nd, then 3rd
  unique_datasets <- unique(mean_d__by_study$dataset) # for forest, both color and legends # TODO - this and the next can be cleaned up, in comparison with the above - get this from the study_info instead?
  color_by_dataset <- unique_cols_by_dataset[match(mean_d__by_study$dataset, unique_datasets)] # for forest 
  # unique_studies <- unique(d_long$study) # a little different, since starting with the full (d_long) than summary dataset (mean_d__by_study)
  # color_by_study <- unique_cols_by_dataset[match(d_long$study, unique_studies)]

  # show_col(unique_cols_by_stat) # for testing

  # Make factor for the study and statistic, so we can group the studies by study in hist, CI, etc. - TODO: prob move this above
  d_long$study_factor <- factor(d_long$study, levels=unique(d_long$study))
  d_long$statistic_factor <- factor(d_long$statistic, levels=unique(d_long$statistic))
  # this is the reverse version - maybe this was used for histogram
#   d_long$study_factor <- factor(d_long$study, levels=rev(unique(d_long$study)))
#   d_long$statistic_factor <- factor(d_long$statistic, levels=rev(unique(d_long$statistic)))


  ## CI ribbons

  # Downsample for plotting

  d_long__sorted__for_plot <- d_long__sorted
  # plot every: 10th row for studies with <10000 n_edges, 100th row for studies with >=10000 n_edges, and 1000th row for studies with >100000 n_edges
  d_long__sorted__for_plot$plot_idx_within_study <- d_long__sorted__for_plot$idx_within_study
  d_long__sorted__for_plot$plot_idx_within_study <- d_long__sorted__for_plot$plot_idx_within_study %% 100 == 0 # this restarts the counting for each study
  idx_small_n_edges_studies <- which(d_long__sorted__for_plot$study %in% study_info$study[study_info$n_edges < 10000])
  d_long__sorted__for_plot$plot_idx_within_study[idx_small_n_edges_studies] <- d_long__sorted__for_plot$idx_within_study[idx_small_n_edges_studies] %% 10 == 0
  idx_large_n_edges_studies <- which(d_long__sorted__for_plot$study %in% study_info$study[study_info$n_edges > 100000])
  d_long__sorted__for_plot$plot_idx_within_study[idx_large_n_edges_studies] <- d_long__sorted__for_plot$idx_within_study[idx_large_n_edges_studies] %% 1000 == 0
  # also plot the first and last entry for each study
  d_long__sorted__for_plot$plot_idx_within_study[d_long__sorted__for_plot$idx_within_study == 1] <- TRUE
  d_long__sorted__for_plot$plot_idx_within_study[d_long__sorted__for_plot$idx_within_study[-1] - d_long__sorted__for_plot$idx_within_study[-length(d_long__sorted__for_plot$idx_within_study)] < 1] <- TRUE
  ci_contains_zero_transition_idx<-which(d_long__sorted__for_plot$ci_contains_zero[-1] != d_long__sorted__for_plot$ci_contains_zero[-length(d_long__sorted__for_plot$ci_contains_zero)])
  # also keep the row before and after the transition into CI containing 0 - prevents gap between the red and blue CI ribbons in the plot
  d_long__sorted__for_plot$plot_idx_within_study[ci_contains_zero_transition_idx] <- TRUE
  d_long__sorted__for_plot$plot_idx_within_study[ci_contains_zero_transition_idx+1] <- TRUE
  d_long__sorted__for_plot <- d_long__sorted__for_plot[d_long__sorted__for_plot$plot_idx_within_study,]
#   d_long__sorted__for_plot2<-d_long__sorted__for_plot[d_long__sorted__for_plot$plot_idx_within_study,] # write into separate plotting variable to be safe - TODO: when confident, apply in place instead


  # TEMPORARY - remove IMAGEN
  d_long__sorted__for_plot <- d_long__sorted__for_plot[d_long__sorted__for_plot$dataset != "IMAGEN",]

  # Plot ascending d with CI ribbons
  # https://thenode.biologists.com/visualizing-data-with-r-ggplot2/education/
  # make separate plots via facets: https://www3.nd.edu/~steve/computing_with_data/13_Facets/facets.html
  # a little faster with shape = "." - https://stackoverflow.com/questions/56972234/improve-performance-for-facet-grid-plot-on-big-data
  # changing colors for ribbon: https://stackoverflow.com/questions/50217794/change-colors-of-a-single-geom-ribbon-depending-on-variable
  # getting the aesthetics right: https://stackoverflow.com/questions/54648389/ggplotly-with-geom-ribbon-grouping
  # background: https://r-graphics.org/recipe-appearance-hide-gridlines https://www.geeksforgeeks.org/remove-grid-and-background-from-plot-using-ggplot2-in-r/
  # spacing: https://stackoverflow.com/questions/32426951/in-ggplot2-and-facet-wrap-how-to-remove-all-margins-and-padding-yet-keep-strip
  
  my_minimal_theme <- function(){ 
    theme_light() %+replace%
    theme(
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      strip.background = element_blank(), # facet titles
      strip.text = element_text(color = "black",margin = margin(0,0,0.04,0, "cm")),
      panel.spacing=unit(0,'npc'),
      panel.margin.y = unit(-0.5, "lines"),
      legend.position="bottom"
    )
  }

  
  ci_singleplot_width <- 2.8
  ci_singleplot_height <- 1.265
  # ci_singleplot_small_width <- 2.8
  ci_singleplot_small_height <- 0.633
  ci_ncol <- 1
  ci_nplots <- length(unique(d_long__sorted__for_plot$study))
  ci_multiplot_height <- ci_singleplot_height * (ceiling(ci_nplots/ci_ncol)+1)
  ci_multiplot_width <- ci_singleplot_width * ci_ncol
  ci_multiplot_small_height <- ci_singleplot_small_height * (ceiling(ci_nplots/ci_ncol)+1)
  # ci_multiplot_small_width <- ci_singleplot_small_width * ci_ncol



#   d_long__sorted__for_plot__r <- d_long__sorted__for_plot[d_long__sorted__for_plot$statistic=="r",]
  ggplot(d_long__sorted__for_plot,  aes(x= idx_within_study, y = d_long__by_study, group=ci_contains_zero__group, colour=ci_contains_zero, fill=ci_contains_zero)) +
    geom_line(linewidth=0.2) +
    geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub), alpha=0.2, linetype=0) +
    coord_cartesian(ylim = ylim_ci_plot) +
    facet_wrap(~study_factor, ncol=ci_ncol, scales = "free_x") + 
    my_minimal_theme() 
  ggsave(d_ci_ribbon_filename, width = ci_multiplot_width, height = ci_multiplot_height, units = "in", dpi = 300, limitsize = FALSE)

  # replot smaller (special params: ylim, height, ybreaks)
  ggplot(d_long__sorted__for_plot,  aes(x= idx_within_study, y = d_long__by_study, group=ci_contains_zero__group, colour=ci_contains_zero, fill=ci_contains_zero)) +
    geom_line(linewidth=0.2) +
    geom_ribbon(aes(ymin = ci_lb, ymax = ci_ub), alpha=0.2, linetype=0) +
    coord_cartesian(ylim = ylim_ci_plot_smaller) +
    facet_wrap(~study_factor, ncol=ci_ncol, scales = "free_x") + 
    scale_y_continuous(breaks = breaks_ci_plot_smaller) +
    my_minimal_theme()
  ggsave(d_ci_ribbon_small_filename, width = ci_multiplot_width, height = ci_multiplot_small_height, units = "in", dpi = 300, limitsize = FALSE)


  # Plot most conservative effect size estimate
  ggplot(d_long__sorted__for_plot,  aes(x= idx_within_study, y = most_conservative_effect, colour=ci_contains_zero,)) +
    geom_line(linewidth=0.5) +
    coord_cartesian(ylim = ylim_ci_plot) +
    facet_wrap(~study_factor, ncol=1, scales = "free_x") +
    my_minimal_theme()
  ggsave(d_conservative_filename, width = ci_multiplot_width, height = ci_multiplot_height, units = "in", dpi = 300, limitsize = FALSE)

  # replot smaller
  ggplot(d_long__sorted__for_plot,  aes(x= idx_within_study, y = most_conservative_effect, colour=ci_contains_zero,)) +
    geom_line(linewidth=0.5) +
    coord_cartesian(ylim = ylim_ci_plot_smaller) +
    facet_wrap(~study_factor, ncol=1, scales = "free_x") +
    scale_y_continuous(breaks = breaks_ci_plot_smaller) +
    my_minimal_theme()
  ggsave(d_conservative_small_filename, width = ci_multiplot_width, height = ci_multiplot_small_height, units = "in", dpi = 300, limitsize = FALSE)



  ## Histograms

  # Histograms by study, colored by dataset, grouped by statistic (we reordered the studies by statistics already)
  ggplot(d_long,  aes(x = d_long__by_study, y = study_factor, fill = dataset)) +
    geom_density_ridges(scale = 5) +
    theme_ridges() + 
    theme(legend.position = "none") +
    xlim(-.5,.5) + 
    scale_fill_manual(values=unique_cols_by_dataset) # I know this is the same palette as default, but keeping in case want to change again
  ggsave(hist_by_study_filename, width = 8.5, height = 10.5, units = "in", dpi = 300)

  # Histograms by statistic
  #  ggplot(d_long,  aes(x = d_long__by_study, y = statistic, fill = statistic)) +
  ggplot(d_long,  aes(x = d_long__by_study, y = statistic_factor, fill = statistic)) +
    geom_density_ridges(scale = 2.1) +
    theme_ridges() + 
    theme(legend.position = "none") +
    xlim(-1,1) +
    scale_fill_manual(values=unique_cols_by_stat) # I know this is the same palette as default, but keeping in case want to change again
  ggsave(hist_by_stat_filename, width = 8.5, height = 10.5, units = "in", dpi = 300)


  ## Forest plots

  # Forest plot
  # forest__full_filename <- paste0(results_dir, 'forest__full', skip_str__vis,forest_img_suffix)
  # pdf(forest__full_filename,width=forest_img_width,height=forest_img_height,pointsize=forest_img_pointsize)
  # forest(mdl__full, slab = mean_d__by_study$study, colout=color_by_dataset) # paste(mean_d__by_study$study, as.character(mean_d__by_study$year), sep = ", "))
  # legend("top", legend = unique_datasets, fill=unique_cols_by_dataset, horiz = TRUE,box.lwd=box_border_thickness, cex = legend_closeness)
  # dev.off()

  # Forest plot by statistic
  pdf(forest__by_stat_filename,width=forest_img_width,height=forest_img_height,pointsize=forest_img_pointsize)
  forest(mdl__by_stat, colout=color_by_dataset, slab = mean_d__by_study$study) #, order=mean_d__by_study$statistic)
  legend("top", legend = unique_datasets, fill=unique_cols_by_dataset, horiz = TRUE,box.lwd=box_border_thickness, cex = legend_closeness)
  dev.off()

  # Funnel plots
  pdf(funnel__by_stat_filename,width=forest_img_width,height=forest_img_height,pointsize=forest_img_pointsize)
  funnel(mdl__by_stat) # probably meaningless since the esz was essentially calculated with variance using abs
  dev.off()








  ## other saving options:
  # kable(mdl__by_stat, row.names=F) %>%
  #   kable_styling(bootstrap_options = "striped") %>%
  #   kable_styling(full_width = F) %>%
  #   save_kable(file = mdl__by_stat_filename) # TODO: save should be an option

  # cat("For stat params:", "\n")
  # cat(mdl__by_stat, "\n")

  # forest__within_by_stat_filename <- paste0(results_dir, 'forest__within_only', skip_str__vis, forest_img_suffix)
  # pdf(forest__within_by_stat_filename,width=forest_img_width,height=forest_img_height,pointsize=forest_img_pointsize)
  # forest(mdl__by_stat__within, colout=color_by_dataset[this_idx__within], slab = mean_d__by_study$study[this_idx__within]) #,order=mean_d__by_study$statistic[this_idx])
  # legend("top", legend = unique(mean_d__by_study$dataset[this_idx__within]), fill=unique(color_by_dataset[this_idx__within]), horiz = TRUE,box.lwd=box_border_thickness, cex = legend_closeness)
  # dev.off()

  # mdl__by_stat__within_filename <- paste0(results_dir, 'mdl_by_stat__within_only', skip_str__vis, ".txt")
  # sink(file = mdl__by_stat__within_filename)
  # mdl__by_stat__within
  # sink(NULL)


  # forest_stats_filename=paste0(results_dir,"metafor_stats.txt") # TODO: put in setparams
  # sink(forest_stats_filename)

  # mdl__by_mapping__within_filename <- paste0(results_dir, 'forest__within_only', skip_str__vis, ".png")
  # kable(mdl__by_mapping__within, row.names=F) %>%
  #   kable_styling(bootstrap_options = "striped") %>%
  #   kable_styling(full_width = F) %>%
  #   save_kable(file = mdl__by_mapping__within_filename) # TODO: save should be an option

  # cat("By statistic, for one-sample only:", "\n")
  # cat(mdl__by_mapping__within, "\n")

  # sink()


  ## Scatterplot for sample size vs. effect size - IN PROGRESS

  # # effect of sample size
  # variance_d <- aggregate(d_long__by_study ~ study, data = d_long, FUN = var)

  # plot(study_info$n_subs,variance_d$d_long__by_study,col = as.factor(study_info$statistic)) #),xlim=c(0,2000)) # ,xlim=c(0,100000)
  # filter_stat="t2" # r, t, t2, t (act)
  # ids_stat=study_info$statistic == filter_stat
  # points(study_info$n_subs[ids_stat],variance_d$d_long__by_study[ids_stat], col = "red", pch = 16)

  # legend("topright", legend = unique(study_info$statistic), col = unique(as.factor(study_info$statistic)), pch = 1, title = "Group")


  ## Boxplot
  # boxplot(d ~ statistic, data = mean_d__by_study[this_idx,])

  # old saving:
  # meta_img_height <- 1500
  # meta_img_width <- 700
  # sc=3
  # png(forest__full_filename,width=450, height=600, res=600)
  # dev.copy(png, forest__full_filename, width = meta_img_width, height = meta_img_height)
  # dev.off()




  # # Summarize
  # thresh <- c(0.2,0.5,0.8)
  # for (i in thresh) {
  #   perc_gt_thresh <- sum(abs(d_long$d) > i)/length(d_long$d) * 100
  #   sprintf('Percent effects above d=%0.1f: %0.2f%%', i, perc_gt_thresh)
  # }

}


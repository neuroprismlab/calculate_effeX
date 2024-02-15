# format data for R plotting and make a few plots
# study 1, for testing: "ABCD.fc.r.rest.age" 

# library(R.matlab) # reads matlab file formats into R - https://bookdown.org/content/b298e479-b1ab-49fa-b83d-a57c2b034d49/distributions.html#ridgeline-chart

library(pwr)
library(metafor)
library(esc) # for effect size estimation
library(psych) # for CI
library(multcomp) # for multiple comparisons

library(magrittr)
library(data.table) # for rleid for helping CIs imputing
library(zoo) # for na.locf for imputing missing CIs

library(dplyr) # for "rename"

######## Setup ########

## Set Params
source("setparams.R") # as of 050223, only to set effect_maps_filename and some plot params

## Create Data Frames
# Note that this first part relies on keeping study_info and d_master in the same order

# Read variables from matfile
load_data <- readRDS(effect_maps_filename)
d_master <- load_data$d
# ci_lb <- load_data$ci_lb # TODO: append this and ci_ub to d_master so don't have to manipulate each of these independently
# ci_ub <- load_data$ci_ub

# study <- load_data$study$name
n_subs <- t(data.frame(load_data$n))
n_subs__g1 <- t(data.frame(load_data$n1))
n_subs__g2 <- t(data.frame(load_data$n2))
n_groups <- t(data.frame(load_data$n_groups))


# create study_info to store study descriptive info
study_info <- data.frame(n_subs, n_groups)
n_subs__by_group <- data.frame(n_subs__g1, n_subs__g2)
study_info <- merge(study_info, n_subs__by_group, by = "row.names", all.x = TRUE) # careful - default behavior is to reorder rows, and even re-specifying sort=False will still put some things out of order
row.names(study_info) <- study_info$Row.names
colnames(study_info)[1] <- "study"
rm(n_subs, n_subs__g1, n_subs__g2, n_groups)

# alphabetize d_master and ci's to match study_info (needed for renaming rows later)
d_master <- d_master[order(names(d_master))]
# ci_lb <- ci_lb[order(names(ci_lb))]
# ci_ub <- ci_ub[order(names(ci_ub))]

# skip the too small sample hbn*evt and hbn*ppvt studies
if (filter__vis) {
    studies_to_remove <- grepl("hbn.*evt", study_info$study)
    studies_to_remove <- studies_to_remove | grepl("hbn.*ppvt", study_info$study)
    d_master <- d_master[!studies_to_remove]
    # ci_lb <- ci_lb[!studies_to_remove]
    # ci_ub <- ci_ub[!studies_to_remove]
    study_info <- study_info[!studies_to_remove,]
}

# get study fullnames and extract meta-data from titles
study_name_parts <- data.frame(matrix(unlist(strsplit(study_info$study, "_")), nrow=length(study_info$study), byrow=T))
colnames(study_name_parts) <- c("dataset","type","statistic","condition1","condition2")

# reassign ambiguous "d" in statistic and new study name
study_newname <- study_info$study
ids_d=study_name_parts$statistic=="d"
ids_act=study_name_parts$type=="act"
study_name_parts$statistic[ids_d] <- ifelse(study_name_parts$type[ids_d] == "act", "t (act)", "t")
study_newname[ids_d]=gsub(".d.","_t_",study_info$study[ids_d])
study_info$study<-study_newname

# update row names for both study_info and d_master (which is in the same order as study_info) and bind new name parts to study_info
row.names(study_info)<-study_newname
names(d_master)<-study_newname
# names(ci_lb)<-study_newname
# names(ci_ub)<-study_newname
study_info <- cbind(study_info[1],study_name_parts,study_info[2:ncol(study_info)])
rm(study_name_parts, study_newname)
# TODO: rename all this stuff earlier, ideally in combine_and_standardize_effects.R - does this still apply?

# Reorder both d_master and study_info based on statistic: r, t2, t, t (act)
idx_reorder <- order(match(study_info$statistic, desired_stat_order))
d_master<-d_master[idx_reorder]
# ci_lb<-ci_lb[idx_reorder]
# ci_ub<-ci_ub[idx_reorder]
study_info<-study_info[idx_reorder, ]

# Get data dimensions and make uniform (mainly matrices)
study_n_variables <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(study_n_variables) <- c("n_edges", "n_nodes")

for (this_study in study_info$study) {
  if (! study_info[this_study,]$statistic=="t (act)") {
    
    n_edges <- length(d_master[[this_study]])

    # convert square to lower triangle if necessary
    sqrt_n_edges <- sqrt(n_edges)
    if (sqrt_n_edges %% 1==0) {
      # TODO: double check getting same upper triangular indices as the other matlab matrices

      triu_msk<-upper.tri(matrix(1, nrow=sqrt_n_edges, ncol=sqrt_n_edges))

      tmp <- d_master[[this_study]]
      tmp <- matrix(tmp, nrow=sqrt_n_edges)
      d_master[[this_study]] <- matrix(tmp[triu_msk])

      # tmp2 <- ci_lb[[this_study]]
      # tmp2 <- matrix(tmp2, nrow=sqrt_n_edges)
      # ci_lb[[this_study]] <- matrix(tmp2[triu_msk])

      # tmp3 <- ci_ub[[this_study]]
      # tmp3 <- matrix(tmp3, nrow=sqrt_n_edges)
      # ci_ub[[this_study]] <- matrix(tmp3[triu_msk])

      n_edges <- length(d_master[[this_study]])
    }

    n_nodes <- ((-1+sqrt(1+8*n_edges))/2)+1

    # convert HCP matrix correlations from lower to upper triangle, for consistency with all else
    if (study_info[this_study,]$dataset=="hcp" && study_info[this_study,]$statistic=="r"){ 

      triu_msk<-upper.tri(matrix(1, nrow=n_nodes, ncol=n_nodes))
      tril_msk<-lower.tri(matrix(1, nrow=n_nodes, ncol=n_nodes))

      tmp<-tril_msk
      tmp[tril_msk]<-d_master[[this_study]]
      tmp<-t(tmp)
      d_master[[this_study]] <- matrix(tmp[triu_msk])

      # tmp2<-tril_msk
      # tmp2[tril_msk]<-ci_lb[[this_study]]
      # tmp2<-t(tmp2)
      # ci_lb[[this_study]] <- tmp2[triu_msk]

      # tmp3<-tril_msk
      # tmp3[tril_msk]<-ci_ub[[this_study]]
      # tmp3<-t(tmp3)
      # ci_ub[[this_study]] <- tmp3[triu_msk]
    }
  } else {
    n_edges <- length(d_master[[this_study]]) # this is really n_voxels
    n_nodes <- NaN
  }
  
  # if matrix and first dim is 1, reshape to so second dim is 1
  # TODO: this should really be done in the initial loading (combine script) - tmp$stats[[1]] loads stuff as 1 x n_edges - maybe really everything should be made into a vector
  if (is.matrix(d_master[[this_study]]) && dim(d_master[[this_study]])[1]==1) {
    d_master[[this_study]] <- matrix(d_master[[this_study]], ncol=1)
  }

  study_n_variables[this_study,]<-c(n_edges, n_nodes)
}

# Bind all study info
study_info <- cbind(study_info, study_n_variables)
rm(study_n_variables, n_edges, n_nodes, sqrt_n_edges, tmp)





######## Calculate simultaneous CI (e.g., multiple testing corrected) ########
# 3.5 min to run
# TODO: IN PROGRESS - only one output for hcp_fc_t_* and pnc_fc_t2_malerest_femalerest - check the class used for apply

ci_lb2 <- list()
ci_ub2 <- list()
# ci_lb2 <- matrix(NA, nrow=nrow(ci_lb), ncol=1)
# ci_lb2<-ci_lb # TODO: replace, this isn't a good initialization and will prob remove prev calc ci's anyways
# ci_ub2<-ci_ub # TODO: replace, same as above

for (this_study in study_info$study) { # testing: this_study<-names(ci_lb)[1]
  this_alpha_corr<-this_alpha/study_info[this_study,]$n_edges
  if (study_info[this_study,]$statistic=="t" | study_info[this_study,]$statistic=="t (act)") {
    if (is.matrix(d_master[[this_study]])) {
      ci_tmp <- apply(d_master[[this_study]],1,function(x) d.ci(x,n1=study_info[this_study,]$n_subs,alpha=this_alpha_corr))
    } else { # vector
      ci_tmp <- sapply(d_master[[this_study]], function(x) d.ci(x, n1=study_info[this_study,]$n_subs, alpha=this_alpha_corr))
    }
    ci_lb2[[this_study]] <- ci_tmp[1,]
    ci_ub2[[this_study]] <- ci_tmp[3,]

    # ci_lb[[this_study]] <- tmp$ci[1,] / sqrt(n[[this_study]])
    # ci_ub[[this_study]] <- tmp$ci[2,] / sqrt(n[[this_study]])
  } else if (study_info[this_study,]$statistic=="t2") {
    # TODO: test
    # ci_tmp<-d.ci(d[[this_study]],n1=n1[[this_study]],n2=n2[[this_study]][1],alpha=.05)
    ci_tmp<-apply(d_master[[this_study]],1,function(x) d.ci(x,n1=study_info[this_study,]$n_subs__g1,
                                            study_info[this_study,]$n_subs__g2,alpha=this_alpha_corr))
    ci_lb2[[this_study]] <- ci_tmp[1,]
    ci_ub2[[this_study]] <- ci_tmp[3,]

    # ci_lb[[this_study]] <- tmp$ci[1,] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])
    # ci_ub[[this_study]] <- tmp$ci[2,] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])
  } else if (study_info[this_study,]$statistic=="r") {
      r_tmp <- d_master[[this_study]] / (sqrt(d_master[[this_study]]^2 + num_sdx_r2d^2)) # confirmed num_sdx_r2d^2 instead of 4 by trial and error

      # TODO: compare against d CI as if converted to t2 (need to make some group size assumption)
      z_95 <- qnorm(1 - this_alpha_corr) # e.g., 0.05 = 1.96
      r_ci_lb <- tanh(atanh(r_tmp) - z_95 / sqrt(study_info[this_study,]$n_subs - 3)) # TESTING: n[[this_study]][1,1]
      r_ci_ub <- tanh(atanh(r_tmp) + z_95 / sqrt(study_info[this_study,]$n_subs - 3))
      # https://stats.stackexchange.com/questions/109861/exact-central-confidence-interval-for-a-correlation
      # http://faculty.washington.edu/gloftus/P317-318/Useful_Information/r_to_z/PearsonrCIs.pdf 
      # TODO: look into the "3" adjustment
      ci_lb2[[this_study]] <- num_sdx_r2d * r_ci_lb / (1 - r_ci_lb ^ 2) ^ (1/2)
      ci_ub2[[this_study]] <- num_sdx_r2d * r_ci_ub / (1 - r_ci_ub ^ 2) ^ (1/2)
  }
}

## Unwrap data (i.e., make one big vector)
# Purpose: need to assign group per variable for group-based plots, e.g., histograms
# also append relevant study data to each variable
d_long__by_study<-setNames(unlist(d_master, use.names=F), rep(names(d_master), lengths(d_master))) # using original study name
names(d_master)<-study_info$statistic # statistic var
d_long__by_stat<-setNames(unlist(d_master, use.names=F), rep(names(d_master), lengths(d_master)))
names(d_master)<-study_info$dataset # dataset var
d_long__by_dataset<-setNames(unlist(d_master, use.names=F), rep(names(d_master), lengths(d_master)))
names(ci_lb2)<-study_info$study # dataset var
ci_lb_long__by_study <- setNames(unlist(ci_lb2, use.names=F),rep(names(ci_lb2), lengths(ci_lb2)))
names(ci_ub2)<-study_info$study # dataset var
ci_ub_long__by_study <- setNames(unlist(ci_ub2, use.names=F),rep(names(ci_ub2), lengths(ci_ub2)))


names(d_master)<-study_info$study # revert to original study name

# put it all together in one long data frame
# TODO: why doesn't this rename first col to "d" instead of "d_long__by_study"?
d_long<-data.frame(study=names(d_long__by_study),
  d=data.frame(d_long__by_study),statistic=names(d_long__by_stat),
  dataset=names(d_long__by_dataset),ci_lb=ci_lb_long__by_study,
  ci_ub=ci_ub_long__by_study) 
# save(d_long, file = "d_long.rds")

# want: study name - d, dataset, data_type, study_type, conditions, sample_size, atlas_size
# conditions: condition1, condition2
# atlas size: n_nodes, n_edges
# sample size: n_subs, n_groups, n_subs__g1, n_subs__g2


######## Per study per variable CI's ########

# TODO: use pvalues
# # Calculate CI's for each study
# for (this_stat in desired_stat_order) {
#   this_idx=d_long$statistic==this_stat
#   if (this_stat=="t" | this_stat=="t (act)") {
#     d_long$CI_lb[this_idx] <- escalc(measure="SMCC", di=d, ni=n_subs, data=mean_d__by_study[this_idx,])$ci.lb
#   } else if (this_stat=="t2") {
#     d_long$CI_lb[this_idx] <- escalc(measure="SMD", di=d, n1i=n_subs__g1, n2i=n_subs__g2, data=mean_d__by_study[this_idx,])$ci.lb
#   } else if (this_stat=="r") {

#     d_old<-mean_d__by_study$d[this_idx]
#     rpb=d_old/sqrt(d_old^2+4) # convert back to r
#     r_to_d <- esc_rpb(r = rpb,      # point-biserial correlation
#             totaln = mean_d__by_study$n_subs[this_idx]*2,
#             es.type = "d")    
#     d_long$CI_lb[this_idx] <- r_to_d$ci.lb
#   }
# }




######## Meta-Analysis ########

## 1.1 Typical Effect Size: Get "typical study effect size" (mean magnitude) for each study

mean_d__by_study <- aggregate(d_long__by_study ~ study, data = d_long, FUN = function(x) mean(abs(x)))
mean_d__by_study <- rename(mean_d__by_study, d = d_long__by_study)
# append study_info
mean_d__by_study <- merge(mean_d__by_study, study_info, by = "study", all.x = TRUE,sort=FALSE)

# reorder both d_master and study info based on statistic: r, t2, t, t (act)
idx_reorder <- order(match(mean_d__by_study$statistic, desired_stat_order)) # desired order defined above
mean_d__by_study<-mean_d__by_study[idx_reorder, ]

## 1.2 Calculate variance for each study

# calculate sampling variance from hedge's g using escalc - TODO: but SMD accounts for the bias, which we don't need to do here - double check variance estimates
# for two-sample, "SMD" uses vtype="LS" (Hedges 1982) to calculate variance by default
# for one-sample, "SMCC" uses di or ti as input, and ni as sample size
for (this_stat in desired_stat_order) {
  this_idx=mean_d__by_study$statistic==this_stat
  if (this_stat=="t" | this_stat=="t (act)") {
    mean_d__by_study$var[this_idx] <- escalc(measure="SMCC", di=d, ni=n_subs, data=mean_d__by_study[this_idx,])$vi
  } else if (this_stat=="t2") {
    mean_d__by_study$var[this_idx] <- escalc(measure="SMD", di=d, n1i=n_subs__g1, n2i=n_subs__g2, data=mean_d__by_study[this_idx,])$vi
  } else if (this_stat=="r") {

    d_old<-mean_d__by_study$d[this_idx]
    rpb=d_old/sqrt(d_old^2+4) # convert back to r
    r_to_d <- esc_rpb(r = rpb,      # point-biserial correlation
            totaln = mean_d__by_study$n_subs[this_idx]*2,
            es.type = "d")    
    mean_d__by_study$var[this_idx] <- r_to_d$var

    # d_new <- r_to_d$es # for comparison - yay, they're identical
    # diff <- d_old - d_new # only differences in the 18th place
  }
}

## 1.3 Fit meta-analysis model

# add non-brain codes
nonbrain_code_filename <- paste0(local_effect_map_dir, 'non_brain_measure_code.csv') # TODO: put in setparams - also need to update for latest data
nonbrain_code_df <- read.csv(nonbrain_code_filename)
colnames(nonbrain_code_df)[1]<-"study" # weird import of col name
mean_d__by_study$nonbrain_code <- nonbrain_code_df$code[match(mean_d__by_study$study, nonbrain_code_df$study)]


# Define within and between stat type indices
this_idx__between=(mean_d__by_study$statistic=="r" | mean_d__by_study$statistic=="t2")
this_idx__within=(mean_d__by_study$statistic=="t" | mean_d__by_study$statistic=="t (act)")

 # TODO: make sure these next 2 lines don't change plotting or anything
mean_d__by_study$statistic<-as.factor(mean_d__by_study$statistic)
mean_d__by_study$nonbrain_code<-relevel(as.factor(mean_d__by_study$nonbrain_code),ref="demographic (age)")
# relevel(as.factor(statistic),ref="t2")

# Run meta-analyses (all, by stat, by stat - within, by stat - between)
mdl__full <- rma(d, var, mods=~as.factor(statistic), data=mean_d__by_study)
mdl__by_stat <- rma(d, var, mods=~as.factor(statistic)-1, data=mean_d__by_study)
mdl__by_stat__within_pairwise <- rma(d, var, mods=~as.factor(statistic), data=mean_d__by_study[this_idx__within,])
mdl__by_stat__between_pairwise <- rma(d, var, mods=~as.factor(statistic), data=mean_d__by_study[this_idx__between,])

# mdl__by_nonbrain_code__between_pairwise <- rma(d, var, mods=~as.factor(nonbrain_code)*as.factor(statistic), data=mean_d__by_study[this_idx__between,])
# mdl__by_nonbrain_code__between_pairwise <- rma(d, var, mods=~relevel(as.factor(nonbrain_code),ref="demographic (age)")*relevel(as.factor(statistic),ref="t2"), data=mean_d__by_study[this_idx__between,])

# pairwise comparison between categories: 1. estimates for each category (no intercept -> no contrasts), 2. all pairwise contrasts via glht
mdl__by_stat_and_nonbrain_code__no_int__between__for_pairwise <- rma(d, var, mods = ~nonbrain_code + statistic-1, data=mean_d__by_study[this_idx__between,])
pairwise_summary__by_stat_and_nonbrain_code__no_int__between <- summary(glht(mdl__by_stat_and_nonbrain_code__no_int__between__for_pairwise, linfct=cbind(contrMat(rep(1,7), type="Tukey"))), test=adjusted("none"))
# mdl__by_stat_and_nonbrain_code__between_pairwise <- rma(d, var, mods = ~nonbrain_code + statistic, data=mean_d__by_study[this_idx__between,])


# nesting... in progress
# https://stats.stackexchange.com/questions/116659/mixed-effects-meta-regression-with-nested-random-effects-in-metafor-vs-mixed-mod
# mdl__full <- rma.mv(d, var, mods=~as.factor(statistic)-1, data=mean_d__by_study, method="REML", random=~1|dataset/statistic)
# mean_d__by_study$dataset_by_stat <- paste(mean_d__by_study$dataset,mean_d__by_study$statistic, sep = "-")

# mdl__by_stat <- rma.mv(d, var, mods=~as.factor(statistic)-1, data=mean_d__by_study, method="REML", random=~1|dataset/statistic)
# mdl__by_stat__within_pairwise <- rma.mv(d, var, mods=~as.factor(statistic), data=mean_d__by_study[this_idx__within,], method="REML", random=~1|dataset/statistic)
# mdl__by_stat__between_pairwise <- rma.mv(d, var, mods=~as.factor(statistic), data=mean_d__by_study[this_idx__between,], method="REML", random=~1|dataset/statistic)

# mdl__by_stat__between_pairwise <- rma(d, var, mods=~relevel(as.factor(statistic),ref="t"), data=mean_d__by_study[this_idx__between,])

## Meta-regression
# # How is d influenced by: study type? sample size? atlas size?
# # one.way <- aov(abs(d) ~ study_type, data = d_long)
# # summary(one.way)


## 2.1 Edge-level meta (Shen atlas): use d_long but append nonbrain code for meta
if (do_meta_all_edges) {

  # reorder stuff that needs reordering
  # TODO: integrate better

  # d_master_orig <- d_master
  d_master__reordered <- d_master

  # TODO: bandaids, but this fix was already added to combine_and_standardize_data.R - remove when ready to re-run that script
  # TODO: IMPORTANT - still need to change study name for plots (-> pnc_fc_t2_malerest_femalerest)
  d_master__reordered$pnc_fc_t2_malerest_femalerest <- -d_master__reordered$pnc_fc_t2_malerest_femalerest
  # rename pnc_fc_t2_malerest_femalerest to pnc_fc_t2_femalerest_malerest
  # d_master <- d_master %>% rename(pnc_fc_t2_femalerest_malerest = pnc_fc_t2_malerest_femalerest)
  
  mat_tmp<-list()
  this_n_nodes<-268 # hard-coded for Shen
  this_n_edges<-35778 # hard-coded for Shen
  triumask <- upper.tri(matrix(1, this_n_nodes, this_n_nodes))

  # only for Shen atlas
  for (this_study in study_info$study) {
      if (study_info[this_study,"n_edges"]==this_n_edges) {

          # reorder if needed
          if (grepl("hcp", this_study)) { # HCP is special - usually reordered
            mat_tmp<-d_master__reordered[[this_study]] 
            if (grepl("t2", this_study)) { # TODO: temporary fix
              mat_tmp <- structure_data(mat_tmp,triangleside="lower")[[1]]
              mat_tmp <- as.numeric(mat_tmp[triumask])
            }
            
          } else {
            this_d<-d_master__reordered[[this_study]]
            mat_tmp<-triumask
            mat_tmp[triumask] <- this_d
            # reorder_matrix_by_atlas
            mat_tmp <- mat_tmp + t(mat_tmp)
            mat_tmp <- reorder_matrix_by_atlas(mat_tmp)
            mat_tmp <- mat_tmp[triumask]
          }

          d_master__reordered[[this_study]]<-mat_tmp

      }
  }

  # subset shen fc data from relevant variables 
  studies_shen<-study_info$study[study_info$n_edges==35778]

  d_master__reordered__shenfc<-d_master__reordered[names(d_master__reordered) %in% studies_shen]
  study_info__shenfc<-study_info[study_info$study %in% studies_shen,]
  
  d_long__by_study__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc))) # using original study name
  names(d_master__reordered)<-study_info$statistic # statistic var
  d_long__by_stat__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  names(d_master__reordered)<-study_info$dataset # dataset var
  d_long__by_dataset__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  # NO CI because I have not yet reordered the CI by atlas - TODO

  names(d_master__reordered)<-study_info$study # revert to original study name

  # put it all together in one long data frame
  # TODO: why doesn't this rename first col to "d" instead of "d_long__by_study"?
  d_long__shenfc<-data.frame(study=names(d_long__by_study__shenfc),
    d=data.frame(d_long__by_study__shenfc),statistic=names(d_long__by_stat__shenfc),
    dataset=names(d_long__by_dataset__shenfc)) 


  # Define n_subs, n_subs__g1, n_subs__g2, nonbrain_code
  # superficial stuff to just get those variables in there
  mean_d__by_study__shenfc<-mean_d__by_study[mean_d__by_study$study %in% studies_shen,]
  this_idx__between__shenfc<-this_idx__between[mean_d__by_study$study %in% studies_shen]
  this_idx__within__shenfc<-this_idx__within[mean_d__by_study$study %in% studies_shen]
  names(d_master__reordered__shenfc)<-mean_d__by_study__shenfc$n_subs # dataset var
  d_long__n_subs__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  names(d_master__reordered__shenfc)<-mean_d__by_study__shenfc$n_subs__g1 # dataset var
  d_long__n_subs__g1__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  names(d_master__reordered__shenfc)<-mean_d__by_study__shenfc$n_subs__g2 # dataset var
  d_long__n_subs__g2__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  names(d_master__reordered__shenfc)<-mean_d__by_study__shenfc$nonbrain_code # dataset var
  d_long__by_nonbrain_code__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  names(d_master__reordered__shenfc)<-mean_d__by_study__shenfc$statistic # dataset var
  d_long__statistic__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  # add numeric index to each sub-list in d_master__reordered__shenfc
  d_long__edge_ids__shenfc<-unlist(lapply(d_master__reordered__shenfc, function(x) seq_along(x)), use.names=F)

  # Define within and between stat type indices
  # names(d_master__reordered__shenfc)<-this_idx__between__shenfc # dataset var
  # this_idx__between__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))
  # names(d_master__reordered__shenfc)<-this_idx__within__shenfc # dataset var
  # this_idx__within__by_edge__shenfc<-setNames(unlist(d_master__reordered__shenfc, use.names=F), rep(names(d_master__reordered__shenfc), lengths(d_master__reordered__shenfc)))

  names(d_master__reordered__shenfc)<-study_info__shenfc$study # revert to original study name

  # add to long data frame
  d_long__shenfc$n_subs<-as.numeric(names(d_long__n_subs__shenfc))
  d_long__shenfc$n_subs__g1<-as.numeric(names(d_long__n_subs__g1__shenfc))
  d_long__shenfc$n_subs__g2<-as.numeric(names(d_long__n_subs__g2__shenfc))
  d_long__shenfc$nonbrain_code<-relevel(as.factor(names(d_long__by_nonbrain_code__shenfc)),ref="demographic (age)")
  d_long__shenfc$statistic<-as.factor(names(d_long__statistic__shenfc))
  d_long__shenfc$edge_ids<-d_long__edge_ids__shenfc

  # TODO: double check releveling
  # mean_d__by_study$nonbrain_code<-relevel(as.factor(mean_d__by_study$nonbrain_code),ref="demographic (age)") # TODO: make sure this doesn't change plotting or anything




  ## 2.2 Edge-level meta: Calculate variance for each edge

  # calculate sampling variance from hedge's g using escalc - TODO: but SMD accounts for the bias, which we don't need to do here - double check variance estimates
  # for two-sample, "SMD" uses vtype="LS" (Hedges 1982) to calculate variance by default
  # for one-sample, "SMCC" uses di or ti as input, and ni as sample size

  # one_sample_idx<-(mean_d__by_study$statistic=="t" | mean_d__by_study$statistic=="t (act)")
  # d_long__shenfc$var[one_sample_idx]<-escalc(measure="SMCC", di=d_long__by_study, ni=n_subs, data=d_long__shenfc[one_sample_idx,])$vi
  # d_long__shenfc$var[!one_sample_idx]<-escalc(measure="SMCC", di=d_long__by_study, ni=n_subs, data=d_long__shenfc[! one_sample_idx,])$vi

  for (this_stat in desired_stat_order) {
    this_idx=d_long__shenfc$statistic==this_stat
    if (any(this_idx)) {
      if (this_stat=="t" | this_stat=="t (act)") {
        d_long__shenfc$var[this_idx] <- escalc(measure="SMCC", di=d_long__by_study__shenfc, ni=n_subs, data=d_long__shenfc[this_idx,])$vi
      } else if (this_stat=="t2") {
        d_long__shenfc$var[this_idx] <- escalc(measure="SMD", di=d_long__by_study__shenfc, n1i=n_subs__g1, n2i=n_subs__g2, data=d_long__shenfc[this_idx,])$vi
      } else if (this_stat=="r") {

        d_old<-d_long__shenfc$d_long__by_study[this_idx]
        rpb=d_old/sqrt(d_old^2+4) # convert back to r
        r_to_d <- esc_rpb(r = rpb,      # point-biserial correlation
                totaln = d_long__shenfc$n_subs[this_idx]*2,
                es.type = "d")    
        d_long__shenfc$var[this_idx] <- r_to_d$var
      }
    }
  }

  ## 2.3 Edge-level meta: Fit meta-analysis model
  # Timing: 6-10 min to complete (10 min for rma.mv, 6 min for rma.uni)

  # Run meta-analyses (all, by stat, by stat - within, by stat - between)
  # TODO: why as.factor for the latter ones but not the first? And same for above?
  # for non-convergence: https://stackoverflow.com/questions/68817204/why-did-the-fisher-scoring-algorithm-not-converge-after-adjusting

  # mdl__full__each_shenedge<-list()
  # mdl__by_shenedge__by_stat<-list()
  # mdl__by_shenedge__by_stat__within_pairwise<-list()
  # mdl__by_shenedge__by_stat__between_pairwise<-list()
  mdl__stat_cat_dataset__each_shenedge<-list()
    # mdl__stat_cat_dataset__each_shenedge_old<-mdl__stat_cat_dataset__each_shenedge


  start_time <- Sys.time()
  for (this_edge_id in unique(d_long__shenfc$edge_ids)) {
    this_edge_idx <- d_long__shenfc$edge_ids==this_edge_id
    d_long__this_edge <- d_long__shenfc[this_edge_idx,]
    
    # fit models
    # mdl__full__each_shenedge[[this_edge_id]]<-rma(data=d_long__this_edge, yi=d_long__by_study, var, mods=~as.factor(statistic))
    
    # old?
    # # mdl__by_shenedge__by_stat[[this_edge_id]]<-rma(yi=d_long__by_study, var, mods=~as.factor(statistic)-1, data=d_long__this_edge)
    # # mdl__by_shenedge__by_stat__within_pairwise[[this_edge_id]]<-rma(yi=d_long__by_study, var, mods=~as.factor(statistic), data=d_long__this_edge[this_idx__within__shenfc,])
    # mdl__by_shenedge__by_stat__between_pairwise[[this_edge_id]]<-rma(yi=d_long__by_study, var, mods=~as.factor(statistic), data=d_long__this_edge[this_idx__between__shenfc,],control=list(stepadj=0.5, maxiter=10000))

    # note: "between" model needed more iterations/different steps to converge
    # pairwise comparison between categories: 1. estimates for each category (no intercept -> no contrasts), 2. all pairwise contrasts via glht
    mdl__stat_cat_dataset__each_shenedge[[this_edge_id]]<-rma.mv(data=d_long__this_edge, yi=d_long__by_study__shenfc, var, mods = ~nonbrain_code-1, random=~1|dataset)
    # mdl__stat_cat_dataset__each_shenedge[[this_edge_id]]<-rma.mv(data=d_long__this_edge, yi=d_long__by_study, var, mods = ~nonbrain_code + n_subs, random=~1|dataset) # n_subs association is expected to be an artifact of smaller n -> larger variance
    # mdl__stat_cat_dataset__each_shenedge[[this_edge_id]]<-rma.mv(data=d_long__this_edge, yi=d_long__by_study, var, mods = ~nonbrain_code + statistic + n_subs, random=~1|dataset)

  }
  end_time <- Sys.time()
  time_by_edge__mdl1 <- end_time - start_time

  # mdl__by_stat_and_nonbrain_code__no_int__between_pairwise__by_edge<-list()
  # start_time <- Sys.time()
  # for (this_edge_id in unique(d_long__shenfc$edge_ids)) {
  #   this_edge_idx <- d_long__shenfc$edge_ids==this_edge_id
  #   d_long__this_edge <- d_long__shenfc[this_edge_idx,]
  #   mdl__by_stat_and_nonbrain_code__no_int__between_pairwise__by_edge[[this_edge_id]] <- rma(data=d_long__this_edge, d_long__by_study__shenfc, var, mods = ~nonbrain_code + statistic-1)
  #   # # summary(glht(mdl__by_stat_and_nonbrain_code__no_int__between_pairwise__by_edge, linfct=cbind(contrMat(rep(1,7), type="Tukey"))), test=adjusted("none"))
  #   # # # mdl__by_stat_and_nonbrain_code__between_pairwise__by_edge[[this_edge_id]] <- rma(d, var, mods = ~nonbrain_code + statistic, data=mean_d__by_study[this_idx__between,])
  # }
  # end_time <- Sys.time()
  # time_by_edge__mdl2 <- end_time - start_time

  # move stuff into a matrix
  cat_beta_matrix<-matrix(NA, nrow=35778, ncol=7)
  for (this_edge_id in unique(d_long__shenfc$edge_ids)) {
    # sig__stat_cat_dataset__each_shenedge[[this_edge_id]]<-mdl__stat_cat_dataset__each_shenedge[[this_edge_id]]$pval*35778<0.05
    cat_beta_matrix[this_edge_id,]<-mdl__stat_cat_dataset__each_shenedge[[this_edge_id]]$beta
  }
  colnames(cat_beta_matrix)<-rownames(mdl__stat_cat_dataset__each_shenedge[[this_edge_id]]$beta)

  # get: point estimate + CI (from SE) for each stat x study -> edge
  # pval_uncorr__t_vs_tact<-sapply(mdl__full__each_shenedge, function(x) x$pval)
  # pval_corr__t_vs_tact<-p.adjust(pval_uncorr[2,], method="fdr")

  # NO within (t vs. t(act) - t(act)) - does not exist for FC
  # pval_uncorr__t_vs_tact<-sapply(mdl__by_shenedge__by_stat__within_pairwise, function(x) x$pval)
  # pval_corr__t_vs_tact<-p.adjust(pval_uncorr__t_vs_tact, method="fdr")
  # proportion_sig__t_vs_tact<-sum(pval_corr__t_vs_tact<0.05)/length(pval_corr__t_vs_tact<0.05)

  # pvals
  pval_uncorr__r_v_t2 <- sapply(mdl__by_shenedge__by_stat__between_pairwise, function(x) x$pval)
  pval_corr__r_v_t2 <- p.adjust(pval_uncorr__r_v_t2, method="fdr")
  proportion_sig__r_v_t2 <- sum(pval_corr__r_v_t2<0.05)/length(pval_corr__r_v_t2<0.05)

}




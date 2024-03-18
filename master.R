##############################################
#
# Summarize effect maps
#
# Prerequisites:
#   1. Calculate group effect maps
#       - Instructions for estimating effect maps: https://docs.google.com/document/d/1Sj2nC_4VocOEzOOEx6GokruxYy4Zp1vwsG72gBUSQf0/edit#
#   2. Transfer group effect maps to Discovery
#       - run: transfer_orig_maps.sh
#       - moves data from MRRC -> Discovery
#         (i.e., MRRC /data_dustin/store3/training/effect_size/ -> Discovery /work/neuroprism/effect_size/)
#
# Input: Data naming conventions - separate studies:
#   - File names:
#       (dataset)_(act|fc)_(t|t2|r)_(var/group1)_(var/group2)
#       Ex 1 (FC ttest):            hcp_fc_t_malerest_femalerest.mat
#       Ex 2 (activation ttest):    hcp_act_t_emotion_rest.mat
#   - Variables:
#       correlation: r: mx1,  p: mx1, std_x: 1xm
#       t-test:      p: 1xm,  ci: 2xm,  stats.tstat: 1xm, n: 1    if ttest2, instead of n: n1, n2
#   - Study script names: (matfile name).m
#
##############################################

# Libraries
source("setparams.R") # USER-DEFINED: review for user-defined parameters
source("clean_data.R")
#source("convert_to_d.R")
#source("estimate_simci.R")
# source("effeX_vis.R") # rely on Hallee's visualizations

# Clean effect maps and combine into single data frame
# TODO: looks like variables are inherited from setparams.R - may not be the best idea - figure out some separation
# access e.g., with cleaned_data$effect_map
cleaned_data <- clean_data()

# # Convert effect maps to d
d <- calc_d(cleaned_data$study, cleaned_data$effect_map)

# # Estimate simci
ci_sim <- estimate_simci(d, cleaned_data$study)

# ci_sim contains all data, with d, and sim_ci
# cleaned_data$study contains all study attributes

# # Visualize effect maps
# effeX_vis(d, ci_sim, cleaned_data$study)


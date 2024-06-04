##############################################
#
# Estimate simultaneous confidence intervals
#
# In: 
# - cleaned_data$study: table with study attributes
# - d: list with effect maps, including d
#
# Output: 
# - effect_map: list of effect maps including simultaneous confidence intervals, stored as attribute "ci"
#
##############################################

estimate_simci <- function(d, study, alpha = 0.05, num_sdx_r2d = 2, output_file = '/work/neuroprism/effect_size/output/estimate_simci.RData') {

for (i in 1:length(d)) {
    print(i)
    if (study$orig_stat_type[i] == "t") {
        # load data
        this_d <- d[[i]]$d
        this_alpha_corrected <- alpha / length(this_d)
        this_n <- d[[i]]$n[1]
        this_n_groups <- 1
        this_n_vars <- length(this_d)

        # calculate sim CI
        ci_tmp <- sapply(this_d, function(x) d.ci(x, n1 = this_n, alpha = this_alpha_corrected))
        ci_lb <- ci_tmp[1,]
        ci_ub <- ci_tmp[3,]

        # add sim CI to d
        d[[i]]$sim_ci_lb <- ci_lb
        d[[i]]$sim_ci_ub <- ci_ub
        }

        # check if effect map is a correlation by checking if study$orig_stat_type is equal to "r"
        else if (study$orig_stat_type[i] == "r") { 
        # load data
        this_r <- d[[i]]$orig_stat
        this_alpha_corrected <- alpha / length(this_r)
        this_n <- d[[i]]$n[1]
        
        z_95 <- qnorm(1 - this_alpha_corrected) # e.g., 0.05 = 1.96

        r_ci_lb <- tanh(atanh(this_r) - z_95 / sqrt(this_n-3))
        r_ci_ub <- tanh(atanh(this_r) + z_95 / sqrt(this_n-3))

        d_ci_lb <- num_sdx_r2d * r_ci_lb / (1 - r_ci_lb ^ 2) ^ (1/2)
        d_ci_ub <- num_sdx_r2d * r_ci_ub / (1 - r_ci_ub ^ 2) ^ (1/2)

        # add sim CI to d
        d[[i]]$sim_ci_lb <- d_ci_lb
        d[[i]]$sim_ci_ub <- d_ci_ub
        }

        # check if effect map is a d value by checking if study$orig_stat_type is equal to "d"
        else if (study$orig_stat_type[i] == "d") {
        # load data
        this_d <- d[[i]]$d
        this_alpha_corrected <- alpha / length(this_d)
        this_n <- d[[i]]$n[1]
        this_n_groups <- 1
        this_n_vars <- length(this_d)

        # calculate sim CI
        ci_tmp <- sapply(this_d, function(x) d.ci(x, n1 = this_n, alpha = this_alpha_corrected))
        ci_lb <- ci_tmp[1,]
        ci_ub <- ci_tmp[3,]

        # add sim CI to d
        d[[i]]$sim_ci_lb <- ci_lb
        d[[i]]$sim_ci_ub <- ci_ub
        }

        # check if effect map is a t2 value by checking if study$orig_stat_type is equal to "t2"
        else if (study$orig_stat_type[i] == "t2") {
        # load data
        this_d <- d[[i]]$d
        this_alpha_corrected <- alpha / length(this_d)
        this_n1 <- d[[i]]$n1[1]
        this_n2 <- d[[i]]$n2[1]
        this_n <- this_n1 + this_n2
        this_n_groups <- 2
        this_n_vars <- length(this_d)

        # calculate sim CI
        ci_tmp <- sapply(this_d, function(x) d.ci(x, n1 = this_n1, n2 = this_n2, alpha = this_alpha_corrected))
        ci_lb <- ci_tmp[1,]
        ci_ub <- ci_tmp[3,]

        # add sim CI to d
        d[[i]]$sim_ci_lb <- ci_lb
        d[[i]]$sim_ci_ub <- ci_ub
        }

        else {
        print("Error: could not calculate simultaneous CI. Check that orig_stat_type is one of: r, t, d, t2.")
        }
    }


     # for HCP activation and FC studies, switch orig_stat_type to t from d TODO: this is a temporary fix
    for (i in 1:length(d)) {
      study_name <- study$name[i]
      if ((study$dataset[i] == "HCP") & (study$orig_stat_type[i] == "d")) {
        study$orig_stat_type[i] <- "t"
      }
    }

    # save results
    save(effect_map = d, study = study, file = output_file)

    # return d
    return(list(d, study))

    
}

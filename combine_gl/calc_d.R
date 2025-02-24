##############################################
#
# Convert effect maps to d
#
# Input: 
# - cleaned_d_maps$study: table with study attributes
# - cleaned_d_maps$data: list of lists with effect maps
# - output_basename: output basename prefix
# - alpha: significance threshold
# - num_sdx_r2d: for r stat: number of standard deviations in X to use for Maya's r-to-d conversion
#
# Output: 
# - d_maps: list with effect maps converted to d, stored as attribute "d"
#
##############################################

# Convert effect maps to d

calc_d <- function(study, d_maps, output_dir, output_basename = 'd_maps', alpha = 0.05, num_sdx_r2d = 2) {
  
  output_path = file.path(output_dir, paste0(output_basename, '_', Sys.Date(), '.RData'))
  
  for (i in 1:length(d_maps)) { # unique studies

    effect_name <- names(d_maps)[i]
    study_idx <- which(study$name == effect_name)
    
    cat('study:', effect_name, '\n')

    for (t in names(d_maps[[i]])) { # unique tests within studies (e.g., w/w-o motion correction)

      cat('   test:', t, '\n')
      
      stat <- d_maps[[i]][[t]]$stat
      alpha_corrected <- alpha / length(stat)
      
      # calculate effect sizes & confidence intervals
      
      result <- calculate_effect_size(stat, study$orig_stat_type[[i]], d_maps, i, t, num_sdx_r2d, alpha_corrected)
      
      d_maps[[i]][[t]]$d <- result$d
      d_maps[[i]][[t]]$sim_ci_lb <- result$ci[1,]
      d_maps[[i]][[t]]$sim_ci_ub <- result$ci[2,]
      d_maps[[i]][[t]]$r_sq <- result$r_sq
      d_maps[[i]][[t]]$r_sq_sim_ci_lb <- result$r_sq_ci[1,]
      d_maps[[i]][[t]]$r_sq_sim_ci_ub <- result$r_sq_ci[2,]
 
      
      # if motion regression, also calculate effects for full residualization case
      
      if (d_maps[[i]][[t]]$motion.method == "regression") {

          stat.fullres <- d_maps[[i]][[t]]$stat.fullres
          result.fullres <- calculate_effect_size(stat.fullres, study$orig_stat_type[[i]], d_maps, i, t, num_sdx_r2d, alpha_corrected)

          d_maps[[i]][[t]]$d.fullres <- result.fullres$d
          d_maps[[i]][[t]]$sim_ci_lb.fullres <- result.fullres$ci[1,]
          d_maps[[i]][[t]]$sim_ci_ub.fullres <- result.fullres$ci[2,]
          d_maps[[i]][[t]]$r_sq.fullres <- result.fullres$r_sq
          d_maps[[i]][[t]]$r_sq_sim_ci_lb.fullres <- result.fullres$r_sq_ci[1,]
          d_maps[[i]][[t]]$r_sq_sim_ci_ub.fullres <- result.fullres$r_sq_ci[2,]
 
      }

    }
  }
  
  save(d_maps, file = output_path)
  return(d_maps)

}




####### Stats Functions ########

d_ci <- function(d, n1, n2 = NULL, alpha = 0.05) {
    if (is.null(n2)) { # one-sample
        n <- n1
        se <- sqrt(1 / n + (d^2 / (2 * n)))
        df <- n - 1
    } else { # two-sample
        se <- sqrt((n1 + n2) / (n1 * n2) + (d^2 / (2 * (n1 + n2))))
        df <- n1 + n2 - 2
    }

    t_crit <- qt(1 - alpha / 2, df = df)
    lower_bound <- d - t_crit * se
    upper_bound <- d + t_crit * se

    return(list(lower_bound, upper_bound))
}


r_sq_ci <- function(r_sq, n, alpha = 0.05) {
  
  r_sq_se <- sqrt((4 * r_sq * (1 - r_sq)^2 * (n - 2)^2) / ((n^2 - 1) * (n + 3)))

  z_crit <- qnorm(1 - alpha) # e.g., 0.05 = 1.96
  
  lower_bound <- r_sq - z_crit * r_sq_se
  upper_bound <- r_sq + z_crit * r_sq_se
  
  return(list(lower_bound, upper_bound))
}

r_ci <- function(r, n, alpha = 0.05) {
  
  z_95 <- qnorm(1 - alpha) # e.g., 0.05 = 1.96
  
  lower_bound <- tanh(atanh(r) - z_95 / sqrt(n-3))
  upper_bound <- tanh(atanh(r) + z_95 / sqrt(n-3))
  
  return(list(lower_bound, upper_bound))
}

#d_ci__from_r_ci <- function(r_ci_lower, r_ci_upper, num_sdx_r2d) {
#
#    lower_bound <- num_sdx_r2d * r_ci_lower / (1 - r_ci_lower ^ 2) ^ (1/2)
#    upper_bound <- num_sdx_r2d * r_ci_upper / (1 - r_ci_upper ^ 2) ^ (1/2)
#
#    return(list(lower_bound, upper_bound))
#}


d_ci__from_r <- function(r, n, num_sdx_r2d, alpha = 0.05) {

    z_95 <- qnorm(1 - alpha) # e.g., 0.05 = 1.96
      
    r_ci_lower <- tanh(atanh(r) - z_95 / sqrt(n-3))
    r_ci_upper <- tanh(atanh(r) + z_95 / sqrt(n-3))

    lower_bound <- num_sdx_r2d * r_ci_lower / (1 - r_ci_lower ^ 2) ^ (1/2)
    upper_bound <- num_sdx_r2d * r_ci_upper / (1 - r_ci_upper ^ 2) ^ (1/2)

    return(list(lower_bound, upper_bound))
}


calculate_effect_size <- function(stat, stat_type, d_maps, i, t, num_sdx_r2d, alpha_corrected) {
    
  switch(stat_type,
         "r" = {
           d <- num_sdx_r2d * stat / ((1 - stat^2) ^ (1/2))
           ci <- sapply(stat, function(x) d_ci__from_r(x, n = d_maps[[i]][[t]]$n[1], num_sdx_r2d = num_sdx_r2d, alpha = alpha_corrected))
           
           r_sq <- stat^2
           r_sq_ci <- sapply(r_sq, function(x) r_sq_ci(x, n = d_maps[[i]][[t]]$n[1], alpha = alpha_corrected))
           #r_sq_ci <- apply(r_ci, c(1, 2), function(x) unlist(x)^2) # TODO: doing list-wise because sapply returns matrix of lists - simplify


         },
         
         "t2" = {
           d <- stat * sqrt(1/d_maps[[i]][[t]]$n1[1] + 1/d_maps[[i]][[t]]$n2[1])
           if (!is.null(d_maps[[i]][[t]]$n1[1])) { # TODO: catch this issue elsewhere then remove
             ci <- sapply(d, function(x) d_ci(x, n1 = d_maps[[i]][[t]]$n1[1], n2 = d_maps[[i]][[t]]$n2[1], alpha = alpha_corrected))
             
             r <- stat / sqrt(stat^2 + (d_maps[[i]][[t]]$n1[1] + d_maps[[i]][[t]]$n2[1] - 2))
             r_sq <- r^2
             r_sq_ci <- sapply(r_sq, function(x) r_sq_ci(x, n = d_maps[[i]][[t]]$n[1], alpha = alpha_corrected))
             #r_sq_ci <- apply(r_ci, c(1, 2), function(x) unlist(x)^2) # TODO: doing list-wise because sapply returns matrix of lists - simplify
           }
         },
         
         "t" = {
           d <- stat / sqrt(d_maps[[i]][[t]]$n[1])
           ci <- sapply(d, function(x) d_ci(x, n1 = d_maps[[i]][[t]]$n[1], alpha = alpha_corrected))
           
           r_sq <- matrix(NaN)
           r_sq_ci <- matrix(c(NaN, NaN), ncol = 1)
         }
  )
  
  
  # quick cleanup - convert empty results to NaN (e.g., we do not calculate multiple regression analog of multi_t so will be empty)
  if (length(d)==0) {
    d <- matrix(NaN)
    ci <- matrix(c(NaN, NaN), ncol = 1)
    r_sq <- matrix(NaN)
    r_sq_ci <- matrix(c(NaN, NaN), ncol = 1)
  }
  
  return(list(d = d, ci = ci, r_sq = r_sq, r_sq_ci = r_sq_ci))
}


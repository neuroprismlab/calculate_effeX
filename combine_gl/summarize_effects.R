##############################################
#
# Summarize effect size data#
# Calculate confidence intervals and simultaneous confidence intervals
# TODO: convert to R^2
#
##############################################

summarize_effects <- function(study, effect_map) {
        
    source("setparams.R") # TODO: remove here because set in master, and just pass necessary variables??

    if (file.exists(effect_maps_filename)) {
        overwrite <- readline(prompt = paste0("Cohen's d coefficient data file ", effect_maps_filename, " already exists. Re-convert data to d? (y/n)"))
        if (overwrite == 'y') {
            convert_to_d <- 1
        } else {
            convert_to_d <- 0
        }
    } else {
        convert_to_d <- 1
    }

    if (convert_to_d) {

        r<-NULL
        t<-NULL
        t2<-NULL
        d<-NULL
        ci_lb<-NULL
        ci_ub<-NULL
        n<-NULL
        n1<-NULL
        n2<-NULL
        n_groups<-NULL

        for (s in 1:length(study$basefile)) {
        
            this_study <- study$name[s]
            this_orig_stat_type <- study$orig_stat_type[s]
            this_filename <- paste(study$folder[s], '/', study$basefile[s], sep = "")

            # Convert to Cohen's d and estimate CI - todo: move CI part to "summarize" script
            tryCatch({
                switch(this_orig_stat_type,
                    # for CI, note that alpha=0.05 is the default for matlab t test
                    "d" = {
                        # ci_tmp <- d.ci(d[[this_study]],n1=n[[this_study]][1],alpha=.05)
                        ci_tmp <- sapply(d[[this_study]], function(x) d.ci(x, n1=n[[this_study]][1], alpha=this_alpha))
                        ci_lb[[this_study]] <- ci_tmp[1,]
                        ci_ub[[this_study]] <- ci_tmp[3,]
                    },
                    "t" = {
                        # TODO: test
                        # ci_tmp<-d.ci(d[[this_study]],n1=n[[this_study]][1],alpha=.05)
                        ci_tmp<-apply(d[[this_study]],2,function(x) d.ci(x,n1=n[[this_study]][1],alpha=this_alpha))
                        ci_lb[[this_study]] <- ci_tmp[1,]
                        ci_ub[[this_study]] <- ci_tmp[3,]

                        # ci_lb[[this_study]] <- tmp$ci[1,] / sqrt(n[[this_study]])
                        # ci_ub[[this_study]] <- tmp$ci[2,] / sqrt(n[[this_study]])
                    },
                    "t2" = {
                        # also uses n1/n2

                        # TODO: test
                        # ci_tmp<-d.ci(d[[this_study]],n1=n1[[this_study]],n2=n2[[this_study]][1],alpha=.05)
                        ci_tmp<-apply(d[[this_study]],2,function(x) d.ci(x,n1=n1[[this_study]],n2=n2[[this_study]][1],alpha=this_alpha))
                        ci_lb[[this_study]] <- ci_tmp[1,]
                        ci_ub[[this_study]] <- ci_tmp[3,]

                        # ci_lb[[this_study]] <- tmp$ci[1,] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])
                        # ci_ub[[this_study]] <- tmp$ci[2,] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])
                    },
                    "r" = {
                        # TODO: compare against d CI as if converted to t2
                        z_95 <- qnorm(1 - this_alpha) # e.g., 0.05 = 1.96
                        r_ci_lb <- tanh(atanh(r[[this_study]]) - z_95 / sqrt(n[[this_study]][1] - 3))
                        r_ci_ub <- tanh(atanh(r[[this_study]]) + z_95 / sqrt(n[[this_study]][1] - 3))
                        # https://stats.stackexchange.com/questions/109861/exact-central-confidence-interval-for-a-correlation
                        # http://faculty.washington.edu/gloftus/P317-318/Useful_Information/r_to_z/PearsonrCIs.pdf 
                        # TODO: look into the "3" adjustment
                        ci_lb[[this_study]] <- num_sdx_r2d * r_ci_lb / (1 - r_ci_lb ^ 2) ^ (1/2)
                        ci_ub[[this_study]] <- num_sdx_r2d * r_ci_ub / (1 - r_ci_ub ^ 2) ^ (1/2)
                    })
            }, error = function(e) {
                if (e$message == "no item called ‘d’") {
                errfnames <- paste(names(tmp), collapse = " ")
                msg <- sprintf("Could not find required variable for %s statistic.\nExisting variables include: %s.\nEither file name specifies the wrong statistic type (common for t->t2) or required variables were not saved in the file.", this_orig_stat_type, errfnames)
                stop(msg)
                } else {
                stop(e)
                }
            })
        }

        saveRDS(list( ci_lb=ci_lb, ci_ub=ci_ub), effect_maps_filename)

        return(list(ci_lb=ci_lb, ci_ub=ci_ub), effect_maps_filename)
    }
}

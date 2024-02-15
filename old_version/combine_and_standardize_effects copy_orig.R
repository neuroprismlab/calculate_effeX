##############################################
#
# Load effect maps and convert all to d
#
# Prerequisites:
#   1. Effect maps previously calculated
#       - Instructions for estimating effect maps: https://docs.google.com/document/d/1Sj2nC_4VocOEzOOEx6GokruxYy4Zp1vwsG72gBUSQf0/edit#
#   2. Data mounted locally: sshfs smn33@10.168.128.189:/data_dustin/store3/training/effect_size/ mnt/
#       - Original data dir: /data_dustin/store3/training/effect_size/
# 
# Next step:
#   - summarize_effects
#
#
# Data conventions:
#
# File names:
#   (dataset)_(act|fc)_(t|t2|r)_(var/group1)_(var/group2)
#   Ex 1 (FC ttest):            hcp_fc_t_malerest_femalerest.mat
#   Ex 2 (activation ttest):    hcp_act_t_emotion_rest.mat
#
# Variables:
#   correlation: r: mx1,  p: mx1, std_x: 1xm
#   t-test:      p: 1xm,  ci: 2xm,  stats.tstat: 1xm, n: 1    if ttest2, instead of n: n1, n2
# 
# Study script names: (matfile name).m
#
##############################################

library(R.matlab) # for reading effect maps from .mat files
library(oro.nifti)
library(psych) # for CI

source("setparams.R")

if (file.exists(effect_maps_precursor_filename)) { # check if data exists
  # prompt user whether to overwrite existing data
  overwrite <- readline(prompt = paste0("Local pre-loaded data file ", effect_maps_precursor_filename, " already exists. Re-load data from scratch? (y/n)"))
  if (overwrite == 'y') {
    combine_maps <- 1
  } else {
    combine_maps <- 0
  }
} else {
  combine_maps <- 1
}

if (file.exists(effect_maps_filename)) { # check if data exists
  # prompt user whether to overwrite existing data
    overwrite <- readline(prompt = paste0("Cohen's d coeffient data file ", effect_maps_filename, " already exists. Convert data to d? (y/n)"))
    if (overwrite == 'y') {
        convert_to_d <- 1
    } else {
        convert_to_d <- 0
    }
} else {
    convert_to_d <- 1
}

# Combine effects across studies and convert all to d
if (combine_maps) {

    # check connected
    if (!dir.exists(data_dir)) {
        stop(paste0('Data directory ', data_dir, ' not found. Check the mount point and try again.'))
    }

    # get filenames
    study <- data.frame()
    for (i in 1:length(exts)) {
        fn <- list.files(path = data_dir, pattern = exts[i], full.names = TRUE)
        fn <- grep(fn,pattern='__helper', invert=TRUE, value=TRUE)
        bn <- data.frame(basefile = basename(fn))
        study <- rbind(study, bn)
    }
    study$folder <- data_dir
    
    # parse filenames: (dataset)_(act|fc)_(t|t2|r)_(var/group1)_(var/group2)
    colnames(study)[1] <- "basefile" # replace name with more descriptive "basefile"
    study$basefile <- gsub("\\.nii\\.gz", ".niigz", study$basefile) # special replace for nii to make same number of elements - will be undone below
    idx_pnc <- grep("pnc", study$basefile)
    study$basefile[idx_pnc] <- gsub("_rt", "xrt", study$basefile[idx_pnc]) # special filename preprocessing for pnc 
    study$basefile[idx_pnc] <- gsub("_rc", "xrc", study$basefile[idx_pnc]) # special filename preprocessing for pnc 
    study$basefile[idx_pnc] <- gsub("_acc", "xacc", study$basefile[idx_pnc]) # special filename preprocessing for pnc 
    
    tmp <- strsplit(study$basefile, "\\.") # split filenames
    study$name <- sapply(tmp, "[", 1)
    study$ext <- sapply(tmp, "[", 2)
    
    tmp <- strsplit(study$name, "_") # split filename components
    study$dataset <- sapply(tmp, "[", 1)
    study$map_type <- sapply(tmp, "[", 2)
    study$stat_type <- sapply(tmp, "[", 3)
    study$var1 <- sapply(tmp, "[", 4)
    study$var2 <- sapply(tmp, "[", 5)
    
    # undo special replace for nii to make same number of elements
    study$basefile <- gsub("\\.niigz", ".nii.gz", study$basefile)
    # study$name <- gsub("\\.niigz", ".nii.gz", study$name)
    study$ext <- gsub("\\.niigz", ".nii.gz", study$ext)
    study$basefile[idx_pnc] <- gsub("xrt", "_rt", study$basefile[idx_pnc]) # special filename preprocessing for pnc 
    study$basefile[idx_pnc] <- gsub("xrc", "_rc", study$basefile[idx_pnc]) # special filename preprocessing for pnc 
    study$basefile[idx_pnc] <- gsub("xacc", "_acc", study$basefile[idx_pnc]) # special filename preprocessing for pnc 
    
    if (skip_nii) {
        files_to_remove <- grepl("nii", study$ext)
        study <- study[!files_to_remove,]
    }
  


    ## Load & convert to Cohen's D ##
    # Note that d is "uncorrected" bc uses sample averages, not population estimates
    # Nice article about Cohen's d: https://www.frontiersin.org/articles/10.3389/fpsyg.2013.00863/full
    #
    # Correlation: using Maya's r-to-d conversion
    # E.g., for 2 SD_X: d reflects the fit-line difference between estimates 2 std apart in X
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7906439
    #       https://rdrr.io/cran/MetaUtility/man/r_to_d.html
    #
    # more on effect size: https://psychology.okstate.edu/faculty/jgrice/psyc3214/Ferguson_EffectSizes_2009.pdf
    # compare with: large d=0.8, large r=0.5

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
        this_stat_type <- study$stat_type[s]
        this_filename <- paste(study$folder[s], '/', study$basefile[s], sep = "")
        
        if (study$ext[s] == 'mat') {
            tmp <- readMat(this_filename)
            if (grepl("hbn_fc_r_rest_nihtoolbox", this_filename)) {
                # special for nih toolbox since multiple outputs saved in here (card, flanker, list, process)--TODO: break apart original file (prob good to do for all "special cases" to avoid special parsing in script)
                    tmp$r <- tmp$card.r
                    # tmp$p <- tmp$card.p # not used
                    # tmp$n is okay - already saved as tmp$n, and same for all, as expected
            }
        } else if (grepl('nii', study$ext[s])) {
            # activation maps stored in nifti format
            # note that if skip_nii=1, there will be no nii files
            tmp_3D <- oro.nifti::readNIfTI(this_filename)
            msk <- tmp_3D!= 0
            tmp[[this_stat_type]] <- tmp_3D[msk]
            rm(tmp_3D)

            # get n saved in separate helper file
            this_helper_filename <- gsub("\\.nii(.*)", "__helper\\.mat", this_filename)
            tmp2 <- readMat(this_helper_filename)
            tmp$n <- tmp2$n

        }

        # temporary fixes - flip signs for (rest-task -> task-rest) and (male-female -> female-male)
        
        # rest-task -> task-rest
        if (grepl('*[td]_rest_*', this_filename)) {
            this_study <- sub("rest_(.*)","\\1_rest",this_study)
            
            study$name[s] <- this_study
            tmp2 <- strsplit(study$name[s], "_") # split filename components
            study$var1[s] <- sapply(tmp2, "[", 4)
            study$var2[s] <- sapply(tmp2, "[", 5)

            tmp$stats[[1]] <- -tmp$stats[[1]] 
        }

        # male-female -> female-male
        if (grepl('*2_malerest_femalerest*', this_filename)) {
            this_study <- sub("malerest_femalerest","femalerest_malerest",this_study)
            
            study$name[s] <- this_study
            tmp2 <- strsplit(study$name[s], "_") # split filename components
            study$var1[s] <- sapply(tmp2, "[", 4)
            study$var2[s] <- sapply(tmp2, "[", 5)

            tmp$stats[[1]] <- -tmp$stats[[1]] 
        }

        # Convert to Cohen's d and estimate CI - todo: move CI part to "summarize" script
        tryCatch({
            switch(this_stat_type,
                # for CI, note that alpha=0.05 is the default for matlab t test
                "d" = {
                    # TODO: still need to characterize original dcoeff according to stat_type (usually t2)
                    d[[this_study]] <- tmp$d
                    tryCatch({
                    n[[this_study]] <- tmp$n
                    }, error = function(e) {
                    stop(sprintf("Probably need to define n for this study: %s", this_study))
                    })
                    n_groups[[this_study]] <- 1

                    # ci_tmp <- d.ci(d[[this_study]],n1=n[[this_study]][1],alpha=.05)
                    ci_tmp <- sapply(d[[this_study]], function(x) d.ci(x, n1=n[[this_study]][1], alpha=this_alpha))
                    ci_lb[[this_study]] <- ci_tmp[1,]
                    ci_ub[[this_study]] <- ci_tmp[3,]
                },
                "t" = {
                    # TODO: why using double brackets here??
                    # also uses n
                    # t[[this_study]] <- tmp$stats$tstat # NO IDEA why this doesn't work anymore, but now there are no fieldnames from the read variable
                    t[[this_study]] <- tmp$stats[[1]]
                    n[[this_study]] <- tmp$n[1]
                    n_groups[[this_study]] <- 1
                    d[[this_study]] <- t[[this_study]] / sqrt(n[[this_study]])
                    
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
                    # t2[[this_study]] <- tmp$stats$tstat # same as above - this should (and previously did) work
                    t2[[this_study]] <- tmp$stats[[1]]
                    n[[this_study]] <- tmp$n1[1] + tmp$n2[1]
                    n1[[this_study]] <- tmp$n1[1] # get also individual
                    n2[[this_study]] <- tmp$n2[1]
                    n_groups[[this_study]] <- 2
                    d[[this_study]] <- t2[[this_study]] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])

                    # TODO: test
                    # ci_tmp<-d.ci(d[[this_study]],n1=n1[[this_study]],n2=n2[[this_study]][1],alpha=.05)
                    ci_tmp<-apply(d[[this_study]],2,function(x) d.ci(x,n1=n1[[this_study]],n2=n2[[this_study]][1],alpha=this_alpha))
                    ci_lb[[this_study]] <- ci_tmp[1,]
                    ci_ub[[this_study]] <- ci_tmp[3,]

                    # ci_lb[[this_study]] <- tmp$ci[1,] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])
                    # ci_ub[[this_study]] <- tmp$ci[2,] * sqrt(1/n1[[this_study]] + 1/n2[[this_study]])
                },
                "r" = {
                    r[[this_study]] <- tmp$r
                    d[[this_study]] <- num_sdx_r2d * r[[this_study]] / (1 - r[[this_study]] ^ 2) ^ (1/2)
                    n[[this_study]] <- tmp$n
                    n_groups[[this_study]] <- 1

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
            msg <- sprintf("Could not find required variable for %s statistic.\nExisting variables include: %s.\nEither file name specifies the wrong statistic type (common for t->t2) or required variables were not saved in the file.", this_stat_type, errfnames)
            stop(msg)
            } else {
            stop(e)
            }
        })
    }

    saveRDS(list(d=d, study=study, n=n, n1=n1, n2=n2, n_groups=n_groups, ci_lb=ci_lb, ci_ub=ci_ub), effect_maps_filename)
}
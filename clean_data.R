##############################################
#
# Clean effect maps and aggregate data and study info
#
# In: Data naming conventions - separate studies:
#   - File names:
#       (dataset)_(act|fc)_(t|t2|r)_(var/group1)_(var/group2)
#       Ex 1 (FC ttest):            hcp_fc_t_malerest_femalerest.mat
#       Ex 2 (activation ttest):    hcp_act_t_emotion_rest.mat # HALLEE: can we standardize which of var1 and var2 is the "task" and which is the "rest" when it's a task vs. rest comparison?
#   - Variables:
#       correlation: r: mx1,  p: mx1, std_x: 1xm # HALLEE: what does mx1 and 1xm mean?
#       t-test:      p: 1xm,  ci: 2xm,  stats.tstat: 1xm, n: 1    if ttest2, instead of n: n1, n2 # TODO: ask Steph, do we want ci here or do we calculate them later now?
#   - Study script names: (matfile name).m
#
# Output: Data attributes - combined studies:
#   - "Study" variable attributes:
#       - basefile, folder, name, ext, dataset, map_type, orig_stat_type, var1, var2
#        # TODO: add phenotype code
#        # HALLEE: should sample size be in this study variable instead of effect map variable?
#   - "Effect map" variable attributes:
#       - matfile variables:
#           - for correlation: r, p, std_x # HALLEE: add n too (unless I'm missing something)
#           - for t-test:      p, ci, stats.tstat, n (if ttest2: n1, n2)
#
# TODO: move all renaming, sign flipping, and triangle->full correction operations to idiosync script where we create a new file for each effect map. Ex: rename r_rest_cognitive_rt -> r_rest_rt, transfer r_rest_cognitive_rt to old/
#
##############################################

# for testing, here are some parameters:
data_dir <- '/work/neuroprism/effect_size/data/individual_studies/'
exts <- c('mat$', 'nii$', 'nii.gz$')
req_fields <- list(
    d = c("d", "n"), # only pre-calculated for t activation map niftis
    t = c("stats", "n"),
    t2 = c("stats", "n1", "n2"),
    r = c("r", "n")
)
output_file <- '/work/neuroprism/effect_size/data/combined_studies/output/effect_maps_clean.RData'

# TODO: address CBCL studies

clean_data <- function(data_dir = '/work/neuroprism/effect_size/data/individual_studies/', exts = c('mat$', 'nii$', 'nii.gz$'), skip_nii = FALSE, skip = NULL, testing = FALSE, req_fields = list(d = c("d", "n"), t = c("stats", "n"), t2 = c("stats", "n1", "n2"), r = c("r", "n")), output_file = '/work/neuroprism/effect_size/effect_maps_clean.RData') {
    
    
    # Libraries
    # check if libraries are installed, and install if not
    # R.matlab for reading .mat effect map files
    # oro.nifti for reading .nii effect map files
    # tidyverse for pipe
    for (pkg in c("R.matlab", "Rcpp", "oro.nifti", "tidyverse", "neurobase")) {
        if (!require(pkg, character.only = TRUE)) {
            install.packages(pkg, dependencies = TRUE)
            library(pkg)
            if (!require(pkg, character.only = TRUE)) {
                stop(paste("Package", pkg, "not found and could not be installed."))
            }
        }
    }

    library(R.matlab)
    library(oro.nifti)
    library(tidyverse)
    library(Rcpp)
    library(neurobase)
    
    # check connected
    if (!dir.exists(data_dir)) {
        stop(paste0('Data directory ', data_dir, ' not found. Check the mount point and try again.'))
    }
    
    # get filenames
    file_paths <- list.files(path = data_dir, pattern = paste(exts, collapse = "|"), full.names = TRUE)
    file_paths <- grep(file_paths, pattern = '__helper', invert = TRUE, value = TRUE) # exclude helper files
    
    study <- data.frame(basefile = basename(file_paths), folder = rep(data_dir, length(file_paths)))

    # TODO: idiosyncratic
    # temporary replace for nii.gz and pnc tasks and hcp_ep so will be treated as single field - will be undone below
    study$basefile <- gsub("\\.nii\\.gz", ".niigz", study$basefile) # nii.gz
    
    idx_pnc <- grep("pnc", study$basefile) # pnc: _rt, _rc, _acc
    pnc_patterns_to_temp_mask <- c("_rt", "_rc", "_acc")
    pnc_patterns_undo <- c()
    
    for (pattern in pnc_patterns_to_temp_mask) {
        replacement <- sub("_", "x", pattern)
        study$basefile[idx_pnc] <- gsub(pattern, replacement, study$basefile[idx_pnc])
        pnc_patterns_undo <- c(pnc_patterns_undo, replacement)
    }
    
    idx_hcp_ep <- grep("hcp_ep", study$basefile) # TODO: change this once we figure out what the ep means 
    hcp_ep_patterns_to_temp_mask <- c("_ep")
    hcp_ep_patterns_undo <- c()
    
    for (pattern in hcp_ep_patterns_to_temp_mask) {
        replacement <- sub("_", "x", pattern)
        study$basefile[idx_hcp_ep] <- gsub(pattern, replacement, study$basefile[idx_hcp_ep])
        hcp_ep_patterns_undo <- c(hcp_ep_patterns_undo, replacement)
    } 
    
    # remove nii files to skip
    if (skip_nii) {
        files_to_remove <- grepl("nii", study$basefile)
        study <- study[!files_to_remove,]
    }
    
    # get names of nih toolbox tasks
    idx_nih <- grep("nih", study$basefile)
    nih_tmp <- readMat(paste0(data_dir, study$basefile[idx_nih]))
    idx = grep("\\.r", attributes(nih_tmp)$names)
    nih_names <- attributes(nih_tmp)$names[idx]
    nih_names <- gsub("\\.r", "", nih_names)
    
    # get fields from filenames
    
    study <- study %>%
        mutate(name = str_split(basefile, "\\.") %>% map_chr(1),
               ext = str_split(basefile, "\\.") %>% map_chr(2),
               dataset = str_split(name, "_") %>% map_chr(1),
               map_type = str_split(name, "_") %>% map_chr(2),
               orig_stat_type = str_split(name, "_") %>% map_chr(3),
               var1 = str_split(name, "_") %>% map_chr(4),
               var2 = str_split(name, "_") %>% map_chr(5))
    # TODO: account for extra details like in this study: ABCD_fc_r_rest_cbcl_scr_syn_somatic_t_FU1
    # right now only the first 5 splits are used, but we should account for the possibility of more splits
    
    # TODO: idiosyncratic
    # add new rows at the bottom of study for each nih task
    for (i in 1:length(nih_names)) {
        study <- rbind(study, study[idx_nih,]) # copy nih toolbox row to the bottom of study
        study$name[dim(study)[1]] <- paste0(study$name[idx_nih], nih_names[i]) # change study name to include nih task
        study$var2[dim(study)[1]] <- nih_names[i] # change var2 to the name of the nih task
    }
    
    # undo temporary replace for nii.gz and pnc tasks
    fields_to_update <- c("basefile", "ext")
    for (field in fields_to_update) {
        study[[field]] <- gsub("\\.niigz", ".nii.gz", study[[field]])
    }
    
    for (pattern in pnc_patterns_undo) {
        replacement <- sub("^x", "_", pattern)
        # re-find the index for the pnc studies since they may have changed from removing other studies
        idx_pnc <- grep("pnc", study$basefile)
        study$basefile[idx_pnc] <- gsub(pattern, replacement, study$basefile[idx_pnc])
    }
    
    for (pattern in hcp_ep_patterns_undo) {
        replacement <- sub("^x", "_", pattern)
        # re-find the index for the hcp_ep studies since they may have changed from removing other studies
        idx_hcp_ep <- grep("hcpxep", study$basefile)
        study$basefile[idx_hcp_ep] <- gsub(pattern, replacement, study$basefile[idx_hcp_ep])
    } # TODO: make sure this is working properly and that in the end study and effect_maps both have proper form of study name
    
    
    # remove other files to skip
    # if skip = NULL, then no files are skipped
    # if skip = c("hbn") for example, then all files containing "hbn" are skipped
    if (!is.null(skip)) {
        files_to_remove <- grepl(paste(skip, collapse = "|"), study$basefile)
        study <- study[!files_to_remove,]
    }
    
    
    if (testing) {   # only use a couple "hard" files
        # choose files with "_rc" and "_malerest_femalerest"
        files_to_keep <- grepl("_rc", study$basefile) | grepl("_malerest_femalerest", study$basefile)
        # skip files containing "v1" for now
        # TODO: update for those
        files_to_keep <- files_to_keep & !grepl("v2", study$basefile)
        study <- study[files_to_keep,]
    }
    
    # TODO: SLIM studies have r instead of t2 in their .mat files
    # steph said it's okay to keep them as r values named t2 values (see meeting notes)
    
    effect_map <- list()
    for (s in 1:length(study$basefile)) {
        
        this_study <- study$name[s]
        this_orig_stat_type <- study$orig_stat_type[s]
        this_filename <- paste(study$folder[s], '/', study$basefile[s], sep = "")
        
        # print this_study as a progress update for debugging
        print(this_study)
        
        # load data into tmp (format: .mat or .nii.gz)
        
        if (study$ext[s] == 'mat') {
            tmp <- readMat(this_filename)
        } else if (grepl('nii', study$ext[s])) {
            # activation maps stored in nifti format
            tmp_3D <- oro.nifti::readNIfTI(this_filename)
            msk <- tmp_3D!= 0
            tmp[[this_orig_stat_type]] <- tmp_3D[msk]
            rm(tmp_3D)
            
            # get n saved in separate helper file
            this_helper_filename <- gsub("\\.nii(.*)", "__helper\\.mat", this_filename)
            tmp2 <- readMat(this_helper_filename)
            tmp$n <- tmp2$n
        }
        
        # check that all fields are loaded as needed for each orig_stat_type
        
        tryCatch({
            if (grepl("nihtoolbox", this_study)) {
                # if nih toolbox, continue on for now and we will check for required fields later
                # TODO: check for required nihtoolbox fields later in the code once parsed
            }
            else{
                stopifnot(all(req_fields[[this_orig_stat_type]] %in% names(tmp)))
            }
        }, error = function(e) {
            existing_fields <- paste(names(tmp), collapse = " ")
            msg <- sprintf("Could not find required variable for %s statistic.\nNeed to provide: %s.\nBut existing variables only include: %s.\nEither file name specifies the wrong statistic type (common for t->t2) or required variables were not saved in the file.", this_orig_stat_type, req_fields[[this_orig_stat_type]], existing_fields)
            stop(msg)
        })
        
        # save stat (r,d,t,t2) to "orig_stat" - we already store the orig_stat_type in the study dataframe for reference
        #   for d and r, orig_stat data has been previously defined in tmp$d and tmp$r
        #   for t and t2, it has been defined in tmp$stats[[1]]
        
        if ("stats" %in% names(tmp)) {
            tmp$orig_stat <- tmp$stats[[1]] # what are stats[[2]] and stats[[3]]?
        } else {
            if (grepl("nihtoolbox", this_study)) {
                task <- study$var2[s]
                tmp$orig_stat <- tmp[[paste0(task, ".r")]]
                for (other_task in nih_names[nih_names != task]) {
                    tmp[[paste0(other_task, ".r")]] <- NULL
                    tmp[[paste0("std.y.", other_task)]] <- NULL
                    tmp[[paste0(other_task, ".p")]] <- NULL
                }
                tmp[[paste0(task, ".r")]] <- NULL
                tmp$p <- tmp[[paste0(task, ".p")]]
                tmp[[paste0(task, ".p")]] <- NULL
                tmp$std.Y <- tmp[[paste0("std.y.", task)]]
                tmp[[paste0("std.y.", task)]] <- NULL
            }
            else {
                tmp$orig_stat <- tmp[[this_orig_stat_type]]
                tmp[[this_orig_stat_type]] <- NULL
            }
        }
        
        # if n, n1, or n2 is null, remove that field
        if ("n" %in% names(tmp) && is.nan(tmp$n)) {
            tmp$n <- NULL
        }
        if ("n1" %in% names(tmp) && is.nan(tmp$n1)) {
            tmp$n1 <- NULL
        }
        if ("n2" %in% names(tmp) && is.nan(tmp$n2)) {
            tmp$n2 <- NULL
        }

        # for each FC study, use the triangle function to make sure the data is in the correct format:
        if (grepl("FC", toupper(this_study)) & (study$var2[s] != "nihtoolbox")) {
            map_path <- "/work/neuroprism/effect_size/data/helper_data/map268_subnetwork.csv"
            if (grepl("UKB", toupper(this_study))) {
                map_path <- NA
            } # TODO: add other maps as we have them!
            tmp$orig_stat <- square_to_triangle(drop(tmp$orig_stat), map_path, show_plot = FALSE)
        }
        
        
        # TODO: maybe change this part so that it allows for nih toolbox?
        
        effect_map[[this_study]] <- tmp
        
    }
    
    # Temporary fixes: flip signs for (rest-task -> task-rest) and (male-female -> female-male)
    # (moved up to here, so that we can use the updated study dataframe when we make effect_maps)
    
    # TODO: testing
    # TODO: do this in a separate script
    
    studies_to_flip <- grepl("[td]_rest_", study$name) | grepl("_malerest_femalerest", study$name)
    # patterns <- c('*[td]_rest_*', '*_malerest_femalerest*') 
    # replacements <- c("\\1_rest", "femalerest_malerest")
    
    # patterns <- c('*[td]_rest_*', '*2_malerest_femalerest*')
    # replacements <- c("\\1_rest", "femalerest_malerest")
    
    for (i in seq_along(studies_to_flip)) {
        if (studies_to_flip[i]) {
            # Update names
            old_name <- study$name[i]
            new_name <- paste0(study$dataset[i], "_", study$map_type[i], "_", study$orig_stat_type[i], "_", study$var2[i], "_", study$var1[i])
            study$name[i] <- new_name
            names(effect_map)[names(effect_map) %in% old_name] <- new_name
            
            # Assign to var1 and var2
            tmp2 <- strsplit(study$name[i], "_")
            study$var1[i] <- sapply(tmp2, "[", 4)
            study$var2[i] <- sapply(tmp2, "[", 5)
            
            # if switching malerest with femalerest, need to switch n1 and n2 as well
            if (grepl("_malerest_femalerest", old_name)) {
                tmp_n1 <- effect_map[[new_name]]$n1
                effect_map[[new_name]]$n1 <- effect_map[[new_name]]$n2
                effect_map[[new_name]]$n2 <- tmp_n1
            } # still need to test this! #TODO
            
            # flip orig_stat sign
            effect_map[[study$name[i]]]$orig_stat <- -effect_map[[study$name[i]]]$orig_stat
        }
    }
    
    # can remove the nihtoolbox studies that have been parsed now
    study <- study[-idx_nih,]

    # remove original nih toolbox study from effect_map
    effect_map[idx_nih] <- NULL
    
    # check to make sure nihtoolbox studies have all required fields in efect_map
    tryCatch({
        nih_idx <- grep("nihtoolbox", study$name)
        for (i in nih_idx) {
            this_study <- study$name[i]
            task <- study$var2[i]
            stopifnot(all(c("n", "std.Y", "std.X", "orig_stat", "p") %in% names(effect_map[[this_study]])))
        }
    }, error = function(e) {
        existing_fields <- paste(names(effect_map[[this_study]]), collapse = " ")
        msg <- sprintf("Could not find required variable for nihtoolbox statistic.\nNeed to provide: %s.\nBut existing variables only include: %s.\nEither file name specifies the wrong statistic type (common for t->t2) or required variables were not saved in the file.", c("n", "std.Y", "std.X", "orig_stat", "p"), existing_fields)
        stop(msg) # TODO: update this error message to be correct for nih toolbox check here
    })
    
    # TODO: right now the nihtoolbox fix only works for if there's one nihtoolbox study
    #       just need to make some for loops to address this!
    
    
    # add HCP activation studies (saved as .nii.gz files with abnormal naming conventions currently)
    # TODO: update this when we have the correct naming conventions
    for (task in c("emotion", "gambling", "relational", "social", "wm")) {
        task_l = task
        task_u = toupper(task_l) # upper-case
        
        n = case_when(
            task_l == "emotion" ~ readMat(paste0(data_dir, 'hcp_act_d_emotion_rest__helper.mat'))$n[1], # TODO: check these numbers since they differ from the fc sample sizes??
            task_l == "gambling" ~ readMat(paste0(data_dir, 'hcp_act_d_gambling_rest__helper.mat'))$n[1],
            task_l == "relational" ~ readMat(paste0(data_dir, 'hcp_act_d_relational_rest__helper.mat'))$n[1],
            task_l == "social" ~ readMat(paste0(data_dir, 'hcp_act_d_social_rest__helper.mat'))$n[1],
            task_l == "wm" ~ readMat(paste0(data_dir, 'hcp_act_d_wm_rest__helper.mat'))$n[1]
        )
        
        for (i in 1:length(task_l)) {
            
            # load nifti
            data_files <- list.files(path = substr(data_dir, 1, nchar(data_dir)-1), full.names = TRUE)
            pattern <- paste0(task_l[i], ".*\\.nii\\.gz")
            matching_file <- grep(pattern, data_files, value = TRUE)
            nifti <- readnii(matching_file)
            
            nifti_data <- img_data(nifti)
            
            list <- as.list(nifti_data)
            
            # remove zero values (likely those outside of the brain, but could theoretically include some in the brain)
            # NEW: don't remove zero values here, remove them when plotting instead because if we remove here then it removes some in the brain and results in 
            # uneven lengths of the effect maps, then can't average across maps when grouping by stat or category
            # nonzero_list <- list[list!=0]
            
            # convert to numeric
            nonzero_numeric <- unlist(list)
            
            # create new list entry to add to cleaned_data$effect_map
            new_data <- list(d = nonzero_numeric, 
                             orig_stat = nonzero_numeric,
                             n = n)
            
            new_study_data <- data.frame(
                basefile = basename(matching_file),
                folder = dirname(matching_file),
                name = paste0("hcp_", "act", "_", "d", "_", task_l, "_rest"),
                ext = ".nii.gz",
                dataset = "hcp",
                map_type = "act",
                orig_stat_type = "d",
                var1 = task_l,
                var2 = "rest")
            
            # combine new data with cleaned_data
            effect_map[[length(effect_map) + 1]] <- new_data
            effect_map <- `names<-`(effect_map, `[<-`(names(effect_map), (length(effect_map)), new_study_data$name))
            study <- rbind(study, new_study_data)
            
        }
    }
    
    # create a group activation map for HCP that removes any voxels that are zero across all HCP tasks:
    # create array of zeros with dimensions of d from one of the HCP tasks by the number of tasks
    all_hcp_act_d <- array(0, dim = c(length(effect_map[[grep("hcp_act", names(effect_map))[[1]]]]$d), length(grep("hcp_act", names(effect_map)))))

    count <- 1
    for (task in c("emotion", "gambling", "relational", "social", "wm")) {
        # get the d values from this task
        d_values <- effect_map[[grep(paste0("hcp_act_d_", task, "_rest"), names(effect_map))]]$d
        # add d_values to the mask
        all_hcp_act_d[,count] <- d_values
        count <- count + 1
    }

    # create a binary mask with 1s for voxels that are nonzero in any of the tasks
    mask <- apply(all_hcp_act_d, 1, function(x) any(x != 0))

    # for each task, remove voxels that are zero in all tasks
    for (task in c("emotion", "gambling", "relational", "social", "wm")) {
        # get the d values from this task
        d_values <- effect_map[[grep(paste0("hcp_act_d_", task, "_rest"), names(effect_map))]]$d
        # get the orig_stat values from this task
        orig_stat_values <- effect_map[[grep(paste0("hcp_act_d_", task, "_rest"), names(effect_map))]]$orig_stat
        # remove zero values from d
        effect_map[[grep(paste0("hcp_act_d_", task, "_rest"), names(effect_map))]]$d <- d_values[mask]
        # remove zero values from orig_stat
        effect_map[[grep(paste0("hcp_act_d_", task, "_rest"), names(effect_map))]]$orig_stat <- orig_stat_values[mask]
    }


    
    # clean up for consistency
    # make all datasets, map types, and names capitalized
    study$dataset <- toupper(study$dataset)
    study$map_type <- toupper(study$map_type)
    study$name <- toupper(study$name)


    ## Save study and effect_maps
    
    save(study, effect_map, file = output_file)
    
    return(list(study = study, effect_map = effect_map))
    
    # TODO: make sure we use "orig_stat" and "orig_stat_type" instead of "stat" and "stat_type" in the other scripts
    
    # TODO: check that we have triangular matrices not full matrices for all studies
    
}

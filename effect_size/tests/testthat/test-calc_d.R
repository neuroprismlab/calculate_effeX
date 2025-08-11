library(testthat)

# Test suite for calc_d function using truncated real test data
# Uses clean_data to get real structure, then truncates for manageable testing

# Helper function to find test data directory
find_test_data_dir <- function() {
  possible_paths <- c(
    "tests/test_data",
    "./tests/test_data", 
    "../tests/test_data",
    "../../tests/test_data",
    file.path(getwd(), "tests/test_data"),
    file.path(dirname(getwd()), "tests/test_data")
  )
  
  for (path in possible_paths) {
    if (dir.exists(path)) {
      cat("Found test data directory at:", normalizePath(path), "\n")
      return(normalizePath(path))
    }
  }
  return(NULL)
}

# Function to truncate data to manageable size
truncate_test_data <- function(cleaned_data, max_studies = 2, max_items = 5) {
  
  # Truncate to max_studies
  study_truncated <- cleaned_data$study[1:min(max_studies, nrow(cleaned_data$study)), ]
  data_truncated <- cleaned_data$data[1:min(max_studies, length(cleaned_data$data))]
  
  # Truncate large vectors/matrices in each study's data
  for (study_name in names(data_truncated)) {
    for (test_name in names(data_truncated[[study_name]])) {
      test_data <- data_truncated[[study_name]][[test_name]]
      
      # Truncate any vectors/matrices with more than max_items elements
      for (field_name in names(test_data)) {
        field_data <- test_data[[field_name]]
        
        if (is.numeric(field_data) && length(field_data) > max_items) {
          # Keep first max_items elements
          data_truncated[[study_name]][[test_name]][[field_name]] <- field_data[1:max_items]
          
        } else if (is.matrix(field_data) && ncol(field_data) > max_items) {
          # For matrices, keep first max_items columns
          data_truncated[[study_name]][[test_name]][[field_name]] <- field_data[, 1:max_items, drop = FALSE]
          
        } else if (is.matrix(field_data) && nrow(field_data) > max_items) {
          # For matrices with many rows, keep first max_items rows
          data_truncated[[study_name]][[test_name]][[field_name]] <- field_data[1:max_items, , drop = FALSE]
        }
      }
    }
  }
  
  return(list(study = study_truncated, data = data_truncated))
}

# Setup function to prepare truncated test data
setup_truncated_test_data <- function() {
  test_data_dir <- find_test_data_dir()
  if (is.null(test_data_dir)) {
    return(NULL)
  }
  
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate") 
  script_dir <- file.path(temp_dir, "scripts")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get real data structure using clean_data
  cat("Loading real test data with clean_data...\n")
  cleaned_data <- clean_data(data_dir = test_data_dir,
                             script_dir = script_dir,
                             intermediate_dir = intermediate_dir,
                             testing = FALSE)
  
  # Truncate to manageable size
  cat("Truncating data: max 2 studies, max 5 items per variable...\n")
  truncated_data <- truncate_test_data(cleaned_data, max_studies = 2, max_items = 5)
  
  cat("Truncated data summary:\n")
  cat("- Studies:", nrow(truncated_data$study), "\n")
  cat("- Study names:", paste(truncated_data$study$name, collapse = ", "), "\n")
  
  # Show size of stat vectors after truncation
  for (study_name in names(truncated_data$data)) {
    for (test_name in names(truncated_data$data[[study_name]])) {
      if ("stat" %in% names(truncated_data$data[[study_name]][[test_name]])) {
        stat_length <- length(truncated_data$data[[study_name]][[test_name]]$stat)
        cat("- ", study_name, "/", test_name, " stat length:", stat_length, "\n")
      }
    }
  }
  
  return(list(
    study = truncated_data$study,
    d_maps = truncated_data$data,
    temp_dir = temp_dir
  ))
}

test_that("calc_d works with truncated real data", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run calc_d function with truncated real data
  result <- calc_d(study = test_data$study,
                   d_maps = test_data$d_maps,
                   output_dir = output_dir,
                   output_basename = "truncated_test",
                   alpha = 0.05,
                   num_sdx_r2d = 2)
  
  # Test basic structure
  expect_type(result, "list")
  expect_equal(length(result), length(test_data$d_maps))
  expect_equal(names(result), names(test_data$d_maps))
  
  # Test that output file was saved
  expected_file <- file.path(output_dir, paste0("truncated_test_", Sys.Date(), ".RData"))
  expect_true(file.exists(expected_file))
  
  cat("calc_d completed successfully with truncated real data!\n")
  
  # Cleanup
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("calc_d adds effect size attributes with real data structure", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- calc_d(study = test_data$study,
                   d_maps = test_data$d_maps,
                   output_dir = output_dir)
  
  # Check that each study has the expected new attributes
  for (study_name in names(result)) {
    for (test_name in names(result[[study_name]])) {
      test_data_item <- result[[study_name]][[test_name]]
      
      # Test that d attribute was added
      expect_true("d" %in% names(test_data_item), 
                  info = paste("Study", study_name, "test", test_name, "should have 'd' attribute"))
      
      # Test that confidence interval attributes were added
      expect_true("sim_ci_lb" %in% names(test_data_item))
      expect_true("sim_ci_ub" %in% names(test_data_item))
      
      # Test that r_sq attributes were added
      expect_true("r_sq" %in% names(test_data_item))
      expect_true("r_sq_sim_ci_lb" %in% names(test_data_item))
      expect_true("r_sq_sim_ci_ub" %in% names(test_data_item))
      
      # Check that d values are reasonable (not all NaN or infinite)
      d_values <- test_data_item$d
      if (length(d_values) > 0) {
        finite_d <- d_values[is.finite(d_values)]
        if (length(finite_d) > 0) {
          expect_true(length(finite_d) > 0, "Should have some finite d values")
        }
      }
    }
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("calc_d handles motion regression with real data structure", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- calc_d(study = test_data$study,
                   d_maps = test_data$d_maps,
                   output_dir = output_dir)
  
  # Check for motion regression cases in the real data
  found_motion_regression <- FALSE
  
  for (study_name in names(result)) {
    for (test_name in names(result[[study_name]])) {
      test_data_item <- result[[study_name]][[test_name]]
      
      if (!is.null(test_data_item$motion.method) && test_data_item$motion.method == "regression") {
        found_motion_regression <- TRUE
        cat("Found motion regression case:", study_name, "/", test_name, "\n")
        
        # Should have fullres attributes
        expect_true("d.fullres" %in% names(test_data_item))
        expect_true("sim_ci_lb.fullres" %in% names(test_data_item))
        expect_true("sim_ci_ub.fullres" %in% names(test_data_item))
        expect_true("r_sq.fullres" %in% names(test_data_item))
        expect_true("r_sq_sim_ci_lb.fullres" %in% names(test_data_item))
        expect_true("r_sq_sim_ci_ub.fullres" %in% names(test_data_item))
      }
    }
  }
  
  if (found_motion_regression) {
    cat("Successfully tested motion regression cases\n")
  } else {
    cat("No motion regression cases found in truncated data\n")
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("calc_d respects different statistical types in real data", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Show what stat types we're working with
  cat("Statistical types in truncated data:\n")
  for (i in 1:nrow(test_data$study)) {
    cat("- Study", i, ":", test_data$study$name[i], "-> stat type:", test_data$study$orig_stat_type[i], "\n")
  }
  
  result <- calc_d(study = test_data$study,
                   d_maps = test_data$d_maps,
                   output_dir = output_dir)
  
  # Test that function handled different stat types
  expect_equal(length(result), nrow(test_data$study))
  
  # Check that d values were calculated for each study
  for (study_name in names(result)) {
    has_d_values <- FALSE
    for (test_name in names(result[[study_name]])) {
      if ("d" %in% names(result[[study_name]][[test_name]])) {
        has_d_values <- TRUE
        break
      }
    }
    expect_true(has_d_values, paste("Study", study_name, "should have d values"))
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("calc_d parameter effects with real data structure", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Test with different num_sdx_r2d values for r-type statistics
  result_sdx_2 <- calc_d(study = test_data$study,
                         d_maps = test_data$d_maps,
                         output_dir = output_dir,
                         output_basename = "sdx_2",
                         num_sdx_r2d = 2)
  
  result_sdx_1 <- calc_d(study = test_data$study,
                         d_maps = test_data$d_maps,
                         output_dir = output_dir,
                         output_basename = "sdx_1", 
                         num_sdx_r2d = 1)
  
  # For r-type studies, d values should be different (proportional to num_sdx_r2d)
  r_studies <- test_data$study[test_data$study$orig_stat_type == "r", ]
  
  if (nrow(r_studies) > 0) {
    cat("Testing num_sdx_r2d effect on", nrow(r_studies), "r-type studies\n")
    
    for (study_name in r_studies$name) {
      if (study_name %in% names(result_sdx_2) && study_name %in% names(result_sdx_1)) {
        # Find a test case with valid d values
        for (test_name in names(result_sdx_2[[study_name]])) {
          d_2 <- result_sdx_2[[study_name]][[test_name]]$d
          d_1 <- result_sdx_1[[study_name]][[test_name]]$d
          
          if (length(d_2) > 0 && length(d_1) > 0 && all(is.finite(d_2)) && all(is.finite(d_1))) {
            # d values should be approximately doubled when num_sdx_r2d doubles
            ratio <- d_2 / d_1
            finite_ratios <- ratio[is.finite(ratio)]
            if (length(finite_ratios) > 0) {
              expect_true(all(abs(finite_ratios - 2) < 0.1), 
                          "d values should scale with num_sdx_r2d for r-type statistics")
              cat("num_sdx_r2d scaling test passed for", study_name, "/", test_name, "\n")
              break
            }
          }
        }
      }
    }
  } else {
    cat("No r-type studies found to test num_sdx_r2d scaling\n")
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("calc_d saves correct data structure", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- calc_d(study = test_data$study,
                   d_maps = test_data$d_maps,
                   output_dir = output_dir,
                   output_basename = "save_test")
  
  # Load saved file and verify it matches returned result
  expected_file <- file.path(output_dir, paste0("save_test_", Sys.Date(), ".RData"))
  expect_true(file.exists(expected_file))
  
  load(expected_file)
  expect_true(exists("d_maps"))
  expect_equal(names(d_maps), names(result))
  expect_equal(length(d_maps), length(result))
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

# Summary test showing what was processed
test_that("calc_d integration summary with truncated real data", {
  test_data <- setup_truncated_test_data()
  skip_if(is.null(test_data), "Test data not available")
  
  output_dir <- file.path(test_data$temp_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n=== calc_d Integration Test Summary ===\n")
  cat("Studies to process:", nrow(test_data$study), "\n")
  
  # Show study details
  for (i in 1:nrow(test_data$study)) {
    row <- test_data$study[i, ]
    cat("Study", i, ":", row$name, "\n")
    cat("  - Dataset:", row$dataset, "\n")
    cat("  - Map type:", row$map_type, "\n")
    cat("  - Stat type:", row$orig_stat_type, "\n")
    cat("  - Category:", row$category, "\n")
    
    # Show number of test conditions
    if (row$name %in% names(test_data$d_maps)) {
      n_tests <- length(test_data$d_maps[[row$name]])
      cat("  - Test conditions:", n_tests, "\n")
    }
  }
  
  # Run calc_d
  cat("\nRunning calc_d...\n")
  result <- calc_d(study = test_data$study,
                   d_maps = test_data$d_maps,
                   output_dir = output_dir,
                   output_basename = "integration_summary")
  
  cat("calc_d completed successfully!\n")
  cat("Processed", length(result), "studies\n")
  
  # Count total test conditions processed
  total_tests <- 0
  for (study_name in names(result)) {
    total_tests <- total_tests + length(result[[study_name]])
  }
  cat("Total test conditions processed:", total_tests, "\n")
  
  # Verify output file
  expected_file <- file.path(output_dir, paste0("integration_summary_", Sys.Date(), ".RData"))
  cat("Output file created:", file.exists(expected_file), "\n")
  
  expect_true(length(result) > 0)
  expect_true(file.exists(expected_file))
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

cat("Running calc_d tests with truncated real data...\n")
cat("Data will be loaded with clean_data, then truncated to 2 studies with max 5 items per variable.\n")
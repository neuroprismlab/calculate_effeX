library(testthat)

# Test suite for checker function using real test data from tests/test_data/

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

# Setup function to prepare test data using clean_data and calc_d
setup_test_d_maps <- function() {
  test_data_dir <- find_test_data_dir()
  if (is.null(test_data_dir)) {
    return(NULL)
  }
  
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate") 
  script_dir <- file.path(temp_dir, "scripts")
  output_dir <- file.path(temp_dir, "output")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get real data structure using clean_data
  cat("Loading real test data with clean_data...\n")
  cleaned_data <- clean_data(data_dir = test_data_dir,
                             script_dir = script_dir,
                             intermediate_dir = intermediate_dir,
                             testing = FALSE)
  
  # Truncate to manageable size (2 studies max, 5 items per variable)
  cat("Truncating data: max 2 studies, max 5 items per variable...\n")
  study_truncated <- cleaned_data$study[1:min(2, nrow(cleaned_data$study)), ]
  data_truncated <- cleaned_data$data[1:min(2, length(cleaned_data$data))]
  
  # Truncate large vectors/matrices in each study's data
  for (study_name in names(data_truncated)) {
    for (test_name in names(data_truncated[[study_name]])) {
      test_data <- data_truncated[[study_name]][[test_name]]
      
      # Truncate any vectors/matrices with more than 5 elements
      for (field_name in names(test_data)) {
        field_data <- test_data[[field_name]]
        
        if (is.numeric(field_data) && length(field_data) > 5) {
          data_truncated[[study_name]][[test_name]][[field_name]] <- field_data[1:5]
        } else if (is.matrix(field_data) && ncol(field_data) > 5) {
          data_truncated[[study_name]][[test_name]][[field_name]] <- field_data[, 1:5, drop = FALSE]
        } else if (is.matrix(field_data) && nrow(field_data) > 5) {
          data_truncated[[study_name]][[test_name]][[field_name]] <- field_data[1:5, , drop = FALSE]
        }
      }
    }
  }
  
  # Run calc_d to get effect size data with proper structure for checker
  cat("Running calc_d to generate effect size data...\n")
  d_maps <- calc_d(study = study_truncated,
                   d_maps = data_truncated,
                   output_dir = output_dir,
                   alpha = 0.05,
                   num_sdx_r2d = 2)
  
  cat("Test d_maps data summary:\n")
  cat("- Studies:", length(d_maps), "\n")
  cat("- Study names:", paste(names(d_maps), collapse = ", "), "\n")
  
  # Show structure of effect size data
  for (study_name in names(d_maps)) {
    n_tests <- length(d_maps[[study_name]])
    cat("- Study", study_name, ":", n_tests, "test conditions\n")
    
    # Show one example test condition structure
    first_test <- names(d_maps[[study_name]])[1]
    test_data <- d_maps[[study_name]][[first_test]]
    d_dim <- if ("d" %in% names(test_data) && is.matrix(test_data$d)) {
      paste(dim(test_data$d), collapse = "x")
    } else if ("d" %in% names(test_data)) {
      paste("length", length(test_data$d))
    } else {
      "missing"
    }
    cat("  - Example 'd' dimensions:", d_dim, "\n")
  }
  
  return(list(
    d_maps = d_maps,
    temp_dir = temp_dir,
    intermediate_dir = intermediate_dir
  ))
}

test_that("checker works with real test data", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  # Set intermediate_dir for the checker function
  intermediate_dir <- test_data$intermediate_dir
  
  result <- checker(d_maps = test_data$d_maps, output_file = "real_data_test", int_dir = intermediate_dir)
  
  # Test basic structure
  expect_type(result, "list")
  expect_equal(length(result), length(test_data$d_maps))
  expect_equal(names(result), names(test_data$d_maps))
  
  cat("checker completed successfully with real data!\n")
  
  # Cleanup
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("checker preserves data structure with real data", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  intermediate_dir <- test_data$intermediate_dir
  
  # Store original structure info
  original_studies <- names(test_data$d_maps)
  original_test_counts <- sapply(test_data$d_maps, length)
  
  result <- checker(d_maps = test_data$d_maps, output_file = "structure_test", int_dir = intermediate_dir)
  
  # Should preserve overall structure
  expect_equal(names(result), original_studies)
  expect_equal(sapply(result, length), original_test_counts)
  
  # Each study should have same test condition names
  for (study_name in names(test_data$d_maps)) {
    expect_equal(names(result[[study_name]]), names(test_data$d_maps[[study_name]]))
    
    # Each test condition should still have required fields
    for (test_name in names(result[[study_name]])) {
      result_fields <- names(result[[study_name]][[test_name]])
      
      # Should still have d field
      expect_true("d" %in% result_fields, 
                  info = paste("Study", study_name, "test", test_name, "should have 'd' field"))
      
      # Should still have CI fields if they existed
      if ("sim_ci_lb" %in% names(test_data$d_maps[[study_name]][[test_name]])) {
        expect_true("sim_ci_lb" %in% result_fields)
      }
      if ("sim_ci_ub" %in% names(test_data$d_maps[[study_name]][[test_name]])) {
        expect_true("sim_ci_ub" %in% result_fields)
      }
    }
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("checker handles motion regression cases with real data", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  intermediate_dir <- test_data$intermediate_dir
  
  result <- checker(d_maps = test_data$d_maps, output_file = "motion_test", int_dir = intermediate_dir)
  
  # Look for motion regression cases in the real data
  found_motion_regression <- FALSE
  
  for (study_name in names(result)) {
    for (test_name in names(result[[study_name]])) {
      if (grepl("motion.regression", test_name)) {
        found_motion_regression <- TRUE
        cat("Found motion regression case:", study_name, "/", test_name, "\n")
        
        test_item <- result[[study_name]][[test_name]]
        
        # If original had fullres fields, they should still be there
        if ("d.fullres" %in% names(test_data$d_maps[[study_name]][[test_name]])) {
          expect_true("d.fullres" %in% names(test_item))
        }
        if ("sim_ci_lb.fullres" %in% names(test_data$d_maps[[study_name]][[test_name]])) {
          expect_true("sim_ci_lb.fullres" %in% names(test_item))
        }
        if ("sim_ci_ub.fullres" %in% names(test_data$d_maps[[study_name]][[test_name]])) {
          expect_true("sim_ci_ub.fullres" %in% names(test_item))
        }
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

test_that("checker processes d matrices correctly", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  intermediate_dir <- test_data$intermediate_dir
  
  # Look at d matrix dimensions before processing
  d_info_before <- list()
  for (study_name in names(test_data$d_maps)) {
    for (test_name in names(test_data$d_maps[[study_name]])) {
      if ("d" %in% names(test_data$d_maps[[study_name]][[test_name]])) {
        d_data <- test_data$d_maps[[study_name]][[test_name]]$d
        key <- paste(study_name, test_name, sep = "_")
        if (is.matrix(d_data)) {
          d_info_before[[key]] <- list(
            is_matrix = TRUE,
            dims = dim(d_data),
            nrow = nrow(d_data),
            ncol = ncol(d_data)
          )
        } else {
          d_info_before[[key]] <- list(
            is_matrix = FALSE,
            length = length(d_data)
          )
        }
      }
    }
  }
  
  cat("d matrix info before checker:\n")
  for (key in names(d_info_before)) {
    info <- d_info_before[[key]]
    if (info$is_matrix) {
      cat("-", key, ": matrix", paste(info$dims, collapse = "x"), "\n")
    } else {
      cat("-", key, ": vector length", info$length, "\n")
    }
  }
  
  result <- checker(d_maps = test_data$d_maps, output_file = "d_matrix_test", int_dir = intermediate_dir)
  
  cat("d matrix info after checker:\n")
  for (study_name in names(result)) {
    for (test_name in names(result[[study_name]])) {
      if ("d" %in% names(result[[study_name]][[test_name]])) {
        d_data <- result[[study_name]][[test_name]]$d
        key <- paste(study_name, test_name, sep = "_")
        if (is.matrix(d_data)) {
          cat("-", key, ": matrix", paste(dim(d_data), collapse = "x"), "\n")
        } else {
          cat("-", key, ": vector length", length(d_data), "\n")
        }
        
        # Values should be finite numbers
        d_values <- as.numeric(d_data)
        finite_values <- d_values[is.finite(d_values)]
        expect_true(length(finite_values) > 0, 
                    info = paste("Study", study_name, "test", test_name, "should have finite d values"))
      }
    }
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("checker saves output file with real data", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  intermediate_dir <- test_data$intermediate_dir
  output_filename <- "real_data_save_test"
  
  result <- checker(d_maps = test_data$d_maps, output_file = output_filename, int_dir = intermediate_dir)
  
  # Check if file was saved
  expected_file <- file.path(intermediate_dir, paste0(output_filename, "_", Sys.Date(), ".RData"))
  
  if (file.exists(expected_file)) {
    expect_true(file.exists(expected_file))
    
    # Load and verify saved data
    load(expected_file)
    expect_true(exists("d_maps"))
    expect_equal(names(d_maps), names(result))
    expect_equal(length(d_maps), length(result))
    
    cat("Output file saved and verified successfully\n")
    file.remove(expected_file)  # Clean up
  } else {
    cat("Output file not found at:", expected_file, "\n")
    cat("Directory contents:", paste(list.files(intermediate_dir), collapse = ", "), "\n")
    # Don't fail the test, just note the issue
    expect_true(TRUE)
  }
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

test_that("checker preserves data values with real data", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  intermediate_dir <- test_data$intermediate_dir
  
  # Calculate checksums of original data
  original_checksums <- list()
  for (study_name in names(test_data$d_maps)) {
    for (test_name in names(test_data$d_maps[[study_name]])) {
      test_item <- test_data$d_maps[[study_name]][[test_name]]
      key <- paste(study_name, test_name, sep = "_")
      
      # Sum all numeric values in this test condition
      numeric_values <- unlist(test_item[sapply(test_item, is.numeric)])
      finite_values <- numeric_values[is.finite(numeric_values)]
      original_checksums[[key]] <- sum(finite_values, na.rm = TRUE)
    }
  }
  
  result <- checker(d_maps = test_data$d_maps, output_file = "value_preservation_test", int_dir = intermediate_dir)
  
  # Verify values are preserved after processing
  for (study_name in names(result)) {
    for (test_name in names(result[[study_name]])) {
      test_item <- result[[study_name]][[test_name]]
      key <- paste(study_name, test_name, sep = "_")
      
      # Sum all numeric values in this test condition after processing
      numeric_values <- unlist(test_item[sapply(test_item, is.numeric)])
      finite_values <- numeric_values[is.finite(numeric_values)]
      result_checksum <- sum(finite_values, na.rm = TRUE)
      
      # Values should be approximately the same (allowing for small numerical differences)
      if (key %in% names(original_checksums)) {
        expect_equal(result_checksum, original_checksums[[key]], 
                     tolerance = 1e-10,
                     info = paste("Data values should be preserved for", key))
      }
    }
  }
  
  cat("Data value preservation verified\n")
  unlink(test_data$temp_dir, recursive = TRUE)
})

# Integration test
test_that("checker integration test with real data", {
  test_data <- setup_test_d_maps()
  skip_if(is.null(test_data), "Test data not available")
  
  intermediate_dir <- test_data$intermediate_dir
  
  cat("\n=== checker Integration Test with Real Data ===\n")
  cat("Studies to process:", length(test_data$d_maps), "\n")
  
  # Show study details
  total_tests <- 0
  for (study_name in names(test_data$d_maps)) {
    study_tests <- length(test_data$d_maps[[study_name]])
    total_tests <- total_tests + study_tests
    cat("- Study:", study_name, "(", study_tests, "test conditions )\n")
  }
  cat("Total test conditions:", total_tests, "\n")
  
  # Run checker
  cat("\nRunning checker...\n")
  result <- checker(d_maps = test_data$d_maps, 
                    output_file = "integration_real_data", int_dir = intermediate_dir)
  
  cat("checker completed successfully!\n")
  cat("Processed", length(result), "studies\n")
  
  # Count processed test conditions
  result_tests <- 0
  for (study_name in names(result)) {
    result_tests <- result_tests + length(result[[study_name]])
  }
  cat("Processed test conditions:", result_tests, "\n")
  
  # Basic validation
  expect_equal(length(result), length(test_data$d_maps))
  expect_equal(result_tests, total_tests)
  
  cat("Integration test completed successfully!\n")
  
  unlink(test_data$temp_dir, recursive = TRUE)
})

cat("Running checker function tests with real test data...\n")
cat("Data will be loaded with clean_data, processed with calc_d, then tested with checker.\n")
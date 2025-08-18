library(testthat)

# Test suite for clean_data function using real test data
# Looks for test data files in tests/test_data/ directory

# Helper function to find test data directory
find_test_data_dir <- function() {
  # Possible locations for test data
  possible_paths <- c(
    "tests/test_data",           # From project root
    "./tests/test_data",         # Current directory
    "../tests/test_data",        # One level up
    "../../tests/test_data",     # Two levels up
    file.path(getwd(), "tests/test_data"),  # Absolute from current working directory
    file.path(dirname(getwd()), "tests/test_data")  # From parent directory
  )
  
  for (path in possible_paths) {
    if (dir.exists(path)) {
      cat("Found test data directory at:", normalizePath(path), "\n")
      return(normalizePath(path))
    }
  }
  
  # If not found, show current working directory and what we tried
  cat("Current working directory:", getwd(), "\n")
  cat("Searched for test data in:\n")
  cat(paste("  -", possible_paths, collapse = "\n"), "\n")
  cat("Please ensure your test data directory exists in one of these locations.\n")
  
  return(NULL)
}

test_that("clean_data returns correct structure with real data", {
  # Find test data directory
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  # Skip test if test data directory doesn't exist
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  # Create required directories
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run the function
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Test main structure
  expect_type(result, "list")
  expect_named(result, c("study", "data", "brain_masks"))
  expect_s3_class(result$study, "data.frame")
  expect_type(result$data, "list")
  expect_type(result$brain_masks, "list")
  
  # Test that data was actually processed
  expect_gt(nrow(result$study), 0, "Should process at least one file")
  expect_gt(length(result$data), 0, "Should have stat_maps data")
  expect_gt(length(result$brain_masks), 0, "Should have brain_masks data")
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data study dataframe has correct structure", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Test study dataframe columns
  expected_cols <- c("basefile", "folder", "name", "ext", "dataset", "map_type",
                     "orig_stat_type", "test_component_1", "test_component_2", 
                     "category", "ref")
  expect_named(result$study, expected_cols)
  
  # Test data types and basic validation
  expect_true(all(result$study$ext == ".mat"))
  expect_true(all(result$study$folder == test_data_dir))
  expect_true(all(!is.na(result$study$basefile)))
  expect_true(all(!is.na(result$study$name)))
  expect_true(all(!is.na(result$study$dataset)))
  expect_true(all(!is.na(result$study$map_type)))
  
  # Test that ref assignment logic works
  expect_true(all(result$study$ref %in% c("voxel", "ukb_55", "shen_268")))
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data processes stat_maps correctly", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Test that stat_maps has entries for each processed file
  expect_equal(length(result$data), nrow(result$study))
  
  # Test that stat_maps names match study names
  expect_setequal(names(result$data), result$study$name)
  
  # Test that each stat_map entry contains data
  for (name in names(result$data)) {
    expect_type(result$data[[name]], "list")
  }
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data processes brain_masks correctly", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Test that brain_masks has entries for each processed file
  expect_equal(length(result$brain_masks), nrow(result$study))
  
  # Test that brain_masks names match study names
  expect_setequal(names(result$brain_masks), result$study$name)
  
  # Test that each brain_mask entry has a mask
  for (name in names(result$brain_masks)) {
    expect_true("mask" %in% names(result$brain_masks[[name]]), 
                info = paste("brain_masks entry for", name, "should have a mask"))
  }
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data ref assignment logic works correctly", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Test ref assignment logic based on map_type and dataset
  for (i in 1:nrow(result$study)) {
    row <- result$study[i, ]
    
    if (row$map_type == "act") {
      expect_equal(row$ref, "voxel", 
                   info = paste("Activation studies should have ref='voxel'"))
    } else if (row$dataset == "ukb") {
      expect_equal(row$ref, "ukb_55", 
                   info = paste("UKB datasets should have ref='ukb_55'"))
    } else {
      expect_equal(row$ref, "shen_268", 
                   info = paste("Other cases should have ref='shen_268'"))
    }
  }
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data testing mode saves file correctly", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Run with testing = TRUE
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = TRUE)
  
  # Check if file was saved
  expected_file <- file.path(intermediate_dir, paste0("clean_data_", Sys.Date(), ".RData"))
  expect_true(file.exists(expected_file), "RData file should be saved when testing=TRUE")
  
  # Load and verify saved data
  load(expected_file)
  expect_true(exists("study"), "Saved file should contain 'study' object")
  expect_true(exists("stat_maps"), "Saved file should contain 'stat_maps' object")
  expect_true(exists("brain_masks"), "Saved file should contain 'brain_masks' object")
  
  # Verify saved data matches returned data
  expect_equal(study, result$study)
  expect_equal(stat_maps, result$data)
  expect_equal(brain_masks, result$brain_masks)
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data custom output_file parameter works", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  custom_output <- "custom_test_data"
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       output_file = custom_output,
                       intermediate_dir = intermediate_dir,
                       testing = TRUE)
  
  # Check if custom named file was created
  expected_file <- file.path(intermediate_dir, paste0(custom_output, "_", Sys.Date(), ".RData"))
  expect_true(file.exists(expected_file), "Custom named RData file should be created")
  
  unlink(temp_dir, recursive = TRUE)
})

test_that("clean_data handles test_component_2 assignment", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Test that test_component_2 is either NA or a valid string
  for (i in 1:nrow(result$study)) {
    component_2 <- result$study$test_component_2[i]
    expect_true(is.na(component_2) || is.character(component_2),
                info = "test_component_2 should be NA or character")
  }
  
  # Test that test_component_1 is always present and not NA
  expect_true(all(!is.na(result$study$test_component_1)),
              "test_component_1 should never be NA")
  expect_true(all(is.character(result$study$test_component_1)),
              "test_component_1 should be character")
  
  unlink(temp_dir, recursive = TRUE)
})

# Summary test that provides overview of what was processed
test_that("clean_data processes test data files summary", {
  test_data_dir <- find_test_data_dir()
  temp_dir <- tempdir()
  intermediate_dir <- file.path(temp_dir, "intermediate")
  script_dir <- file.path(temp_dir, "scripts")
  
  skip_if(is.null(test_data_dir), "Test data directory not found")
  
  # List actual .mat files in test directory
  mat_files <- list.files(test_data_dir, pattern = "\\.mat$", full.names = FALSE)
  cat("\nTest data files found:\n")
  cat(paste(mat_files, collapse = "\n"), "\n")
  
  dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
  
  result <- clean_data(data_dir = test_data_dir, 
                       script_dir = script_dir,
                       intermediate_dir = intermediate_dir,
                       testing = FALSE)
  
  # Verify we processed the expected number of files
  expect_equal(nrow(result$study), length(mat_files))
  
  # Print summary for manual verification
  cat("\nProcessed files summary:\n")
  print(result$study[, c("name", "dataset", "map_type", "category", "ref")])
  
  unlink(temp_dir, recursive = TRUE)
})

cat("Running clean_data function tests with real test data...\n")
cat("Make sure your test data files are in: tests/test_data/\n")
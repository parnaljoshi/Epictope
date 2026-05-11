#!/usr/bin/env Rscript
# Resource monitoring script for EpicTope Shiny app
# Usage: Rscript monitor_resources.R [output_file]

library(jsonlite)

# Function to get system memory info
get_memory_usage <- function() {
  if (Sys.info()["sysname"] == "Windows") {
    # Windows memory check
    cmd_output <- system('wmic OS get TotalVisibleMemorySize,FreePhysicalMemory /format:csv', 
                        intern = TRUE, ignore.stderr = TRUE)
    
    if (length(cmd_output) > 2) {
      data_line <- cmd_output[3]  # Skip headers
      parts <- strsplit(data_line, ",")[[1]]
      if (length(parts) >= 3) {
        free_kb <- as.numeric(parts[2])
        total_kb <- as.numeric(parts[3])
        used_kb <- total_kb - free_kb
        
        return(list(
          total_gb = round(total_kb / 1024 / 1024, 2),
          used_gb = round(used_kb / 1024 / 1024, 2),
          free_gb = round(free_kb / 1024 / 1024, 2),
          percent_used = round(used_kb / total_kb * 100, 1)
        ))
      }
    }
  } else {
    # Linux/Unix memory check
    mem_info <- readLines("/proc/meminfo")
    total_line <- grep("MemTotal:", mem_info, value = TRUE)
    free_line <- grep("MemAvailable:", mem_info, value = TRUE)
    
    total_kb <- as.numeric(gsub(".*?([0-9]+).*", "\\1", total_line))
    free_kb <- as.numeric(gsub(".*?([0-9]+).*", "\\1", free_line))
    used_kb <- total_kb - free_kb
    
    return(list(
      total_gb = round(total_kb / 1024 / 1024, 2),
      used_gb = round(used_kb / 1024 / 1024, 2),
      free_gb = round(free_kb / 1024 / 1024, 2),
      percent_used = round(used_kb / total_kb * 100, 1)
    ))
  }
  
  return(list(total_gb = NA, used_gb = NA, free_gb = NA, percent_used = NA))
}

# Function to get R process memory usage
get_r_memory_usage <- function() {
  # Try multiple methods to get R memory usage
  r_memory <- list()
  
  # Method 1: pryr package (if available)
  if (requireNamespace("pryr", quietly = TRUE)) {
    r_memory$pryr_mb <- round(as.numeric(pryr::mem_used()) / 1024^2, 2)
  }
  
  # Method 2: gc() output
  gc_info <- gc()
  r_memory$gc_used_mb <- round(sum(gc_info[, "used"] * c(8, 8)) / 1024, 2)  # 8 bytes per unit
  r_memory$gc_max_mb <- round(sum(gc_info[, "max used"] * c(8, 8)) / 1024, 2)
  
  # Method 3: object.size of global environment (rough estimate)
  if (requireNamespace("utils", quietly = TRUE)) {
    env_size <- object.size(globalenv())
    r_memory$env_mb <- round(as.numeric(env_size) / 1024^2, 2)
  }
  
  return(r_memory)
}

# Function to get disk usage
get_disk_usage <- function(path = ".") {
  if (Sys.info()["sysname"] == "Windows") {
    # Windows disk check
    drive <- substr(normalizePath(path), 1, 2)
    cmd_output <- system(paste0('dir ', drive, ' /-c'), intern = TRUE, ignore.stderr = TRUE)
    
    # Parse the last line for free space info
    if (length(cmd_output) > 0) {
      last_lines <- tail(cmd_output, 3)
      free_line <- grep("bytes free", last_lines, value = TRUE)
      if (length(free_line) > 0) {
        # Extract numbers from the free space line
        numbers <- regmatches(free_line, gregexpr("[0-9,]+", free_line))[[1]]
        if (length(numbers) >= 2) {
          free_bytes <- as.numeric(gsub(",", "", numbers[length(numbers)]))
          total_bytes <- as.numeric(gsub(",", "", numbers[length(numbers)-1]))
          used_bytes <- total_bytes - free_bytes
          
          return(list(
            total_gb = round(total_bytes / 1024^3, 2),
            used_gb = round(used_bytes / 1024^3, 2),
            free_gb = round(free_bytes / 1024^3, 2),
            percent_used = round(used_bytes / total_bytes * 100, 1)
          ))
        }
      }
    }
  } else {
    # Linux/Unix disk check
    df_output <- system(paste("df", path), intern = TRUE)
    if (length(df_output) > 1) {
      data_line <- df_output[2]
      parts <- strsplit(trimws(data_line), "\\s+")[[1]]
      if (length(parts) >= 4) {
        total_kb <- as.numeric(parts[2])
        used_kb <- as.numeric(parts[3])
        free_kb <- as.numeric(parts[4])
        
        return(list(
          total_gb = round(total_kb / 1024^2, 2),
          used_gb = round(used_kb / 1024^2, 2),
          free_gb = round(free_kb / 1024^2, 2),
          percent_used = round(used_kb / total_kb * 100, 1)
        ))
      }
    }
  }
  
  return(list(total_gb = NA, used_gb = NA, free_gb = NA, percent_used = NA))
}

# Function to get detailed EpicTope data sizes
get_epictope_data_usage <- function() {
  data_usage <- list()
  
  # Check CDS database sizes
  cds_path <- "data/CDS"
  if (dir.exists(cds_path)) {
    cds_files <- list.files(cds_path, full.names = TRUE, recursive = TRUE)
    cds_sizes <- sapply(cds_files, file.size)
    data_usage$cds_mb <- round(sum(cds_sizes, na.rm = TRUE) / 1024^2, 2)
    data_usage$cds_files <- length(cds_files)
  }
  
  # Check models directory
  models_path <- "data/models"
  if (dir.exists(models_path)) {
    model_files <- list.files(models_path, full.names = TRUE, recursive = TRUE)
    model_sizes <- sapply(model_files, file.size)
    data_usage$models_mb <- round(sum(model_sizes, na.rm = TRUE) / 1024^2, 2)
    data_usage$models_files <- length(model_files)
  }
  
  # Check outputs directory
  outputs_path <- "outputs"
  if (dir.exists(outputs_path)) {
    output_files <- list.files(outputs_path, full.names = TRUE, recursive = TRUE)
    output_sizes <- sapply(output_files, file.size)
    data_usage$outputs_mb <- round(sum(output_sizes, na.rm = TRUE) / 1024^2, 2)
    data_usage$outputs_files <- length(output_files)
  }
  
  return(data_usage)
}

# Function to monitor during analysis
monitor_analysis <- function(output_file = "resource_monitor.json", interval_seconds = 5) {
  cat("Starting resource monitoring (Ctrl+C to stop)...\n")
  cat("Logging to:", output_file, "\n")
  
  monitor_data <- list()
  start_time <- Sys.time()
  
  # Signal handler for clean exit
  tryCatch({
    while (TRUE) {
      timestamp <- Sys.time()
      
      # Collect all metrics
      metrics <- list(
        timestamp = as.character(timestamp),
        elapsed_minutes = round(as.numeric(difftime(timestamp, start_time, units = "mins")), 2),
        system_memory = get_memory_usage(),
        r_memory = get_r_memory_usage(),
        disk_usage = get_disk_usage(),
        epictope_data = get_epictope_data_usage()
      )
      
      # Add to monitor data
      monitor_data[[length(monitor_data) + 1]] <- metrics
      
      # Print current status
      cat(sprintf("[%s] RAM: %.1f%% (%.1f/%.1f GB), R: %.0f MB, Disk: %.1f%% (%.1f/%.1f GB)\n",
                  format(timestamp, "%H:%M:%S"),
                  metrics$system_memory$percent_used,
                  metrics$system_memory$used_gb,
                  metrics$system_memory$total_gb,
                  metrics$r_memory$gc_used_mb,
                  metrics$disk_usage$percent_used,
                  metrics$disk_usage$used_gb,
                  metrics$disk_usage$total_gb))
      
      # Save to file periodically
      if (length(monitor_data) %% 10 == 0) {
        writeLines(toJSON(monitor_data, pretty = TRUE), output_file)
      }
      
      Sys.sleep(interval_seconds)
    }
  }, interrupt = function(e) {
    cat("\nMonitoring stopped. Saving final data...\n")
    writeLines(toJSON(monitor_data, pretty = TRUE), output_file)
    cat("Data saved to:", output_file, "\n")
  })
}

# Function to analyze monitoring results
analyze_monitoring_results <- function(results_file = "resource_monitor.json") {
  if (!file.exists(results_file)) {
    cat("Results file not found:", results_file, "\n")
    return(NULL)
  }
  
  data <- fromJSON(results_file)
  
  cat("=== EpicTope Resource Usage Analysis ===\n\n")
  
  # Memory analysis
  ram_usage <- sapply(data, function(x) x$system_memory$used_gb)
  r_memory <- sapply(data, function(x) x$r_memory$gc_used_mb)
  
  cat("Memory Usage Summary:\n")
  cat("- Peak system RAM:", max(ram_usage, na.rm = TRUE), "GB\n")
  cat("- Average system RAM:", round(mean(ram_usage, na.rm = TRUE), 2), "GB\n")
  cat("- Peak R process:", max(r_memory, na.rm = TRUE), "MB\n")
  cat("- Average R process:", round(mean(r_memory, na.rm = TRUE), 2), "MB\n\n")
  
  # Disk analysis
  if (length(data) > 0 && !is.null(data[[1]]$epictope_data)) {
    last_data <- data[[length(data)]]$epictope_data
    cat("EpicTope Data Storage:\n")
    cat("- CDS databases:", last_data$cds_mb, "MB (", last_data$cds_files, "files)\n")
    cat("- Model cache:", last_data$models_mb, "MB (", last_data$models_files, "files)\n")
    cat("- Output files:", last_data$outputs_mb, "MB (", last_data$outputs_files, "files)\n")
    cat("- Total EpicTope data:", round(last_data$cds_mb + last_data$models_mb + last_data$outputs_mb, 2), "MB\n\n")
  }
  
  # Recommendations
  ram_recommendation <- max(ram_usage, na.rm = TRUE) * 1.5  # 50% buffer
  disk_recommendation <- if (length(data) > 0 && !is.null(data[[1]]$epictope_data)) {
    last_data <- data[[length(data)]]$epictope_data
    (last_data$cds_mb + last_data$models_mb) / 1024 * 2  # 100% buffer for growth
  } else {
    10  # Default 10GB if no data
  }
  
  cat("Server Provisioning Recommendations:\n")
  cat("- RAM:", ceiling(ram_recommendation), "GB (including 50% buffer)\n")
  cat("- Disk space:", ceiling(disk_recommendation), "GB (including growth buffer)\n")
  cat("- CPU: 4+ cores recommended for concurrent users\n")
  
  return(list(
    peak_ram_gb = max(ram_usage, na.rm = TRUE),
    avg_ram_gb = mean(ram_usage, na.rm = TRUE),
    peak_r_mb = max(r_memory, na.rm = TRUE),
    recommended_ram_gb = ceiling(ram_recommendation),
    recommended_disk_gb = ceiling(disk_recommendation)
  ))
}

# Main execution
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) > 0 && args[1] == "analyze") {
    # Analyze existing results
    results_file <- if (length(args) > 1) args[2] else "resource_monitor.json"
    analyze_monitoring_results(results_file)
  } else {
    # Start monitoring
    output_file <- if (length(args) > 0) args[1] else "resource_monitor.json"
    monitor_analysis(output_file)
  }
}

# Export functions for use in other scripts
if (interactive()) {
  cat("Resource monitoring functions loaded. Use:\n")
  cat("- get_memory_usage() for current system memory\n")
  cat("- get_r_memory_usage() for R process memory\n")
  cat("- get_disk_usage() for disk space\n")
  cat("- monitor_analysis() to start continuous monitoring\n")
  cat("- analyze_monitoring_results() to analyze log files\n")
}

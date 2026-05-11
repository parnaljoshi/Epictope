library(shiny)
library(epictope)
library(ggplot2)
library(shinycssloaders)
library(reshape2)
library(dplyr)
library(msaR)
library(plotly)
library(bio3d)
library(shinydashboard)
library(shinyBS)
library(DT)
library(shinyWidgets)
library(shinyjs)
library(jsonlite)
library(httr)
library(scales)
library(httr)

# Helper function for null coalescing
`%||%` <- function(a, b) if(is.null(a) || is.na(a)) b else a

# Resource monitoring functions for server provisioning
get_current_memory_usage <- function() {
  tryCatch({
    # Try pryr package first
    if (requireNamespace("pryr", quietly = TRUE)) {
      mem_mb <- round(as.numeric(pryr::mem_used()) / 1024^2, 2)
      return(paste0(mem_mb, " MB"))
    }
    
    # Fallback to gc()
    gc_info <- gc()
    mem_mb <- round(sum(gc_info[, "used"] * c(8, 8)) / 1024, 2)
    return(paste0(mem_mb, " MB (gc)"))
  }, error = function(e) {
    return("Unable to measure")
  })
}

get_system_memory_info <- function() {
  tryCatch({
    if (Sys.info()["sysname"] == "Windows") {
      # Windows memory check using PowerShell
      cmd <- 'powershell "Get-WmiObject -Class Win32_OperatingSystem | Select-Object TotalVisibleMemorySize,FreePhysicalMemory"'
      result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)
      return("Check Task Manager for system memory")
    } else {
      # Linux memory check
      if (file.exists("/proc/meminfo")) {
        mem_info <- readLines("/proc/meminfo")
        total_line <- grep("MemTotal:", mem_info, value = TRUE)
        avail_line <- grep("MemAvailable:", mem_info, value = TRUE)
        
        total_kb <- as.numeric(gsub(".*?([0-9]+).*", "\\1", total_line))
        avail_kb <- as.numeric(gsub(".*?([0-9]+).*", "\\1", avail_line))
        used_kb <- total_kb - avail_kb
        
        total_gb <- round(total_kb / 1024^2, 2)
        used_gb <- round(used_kb / 1024^2, 2)
        percent <- round(used_kb / total_kb * 100, 1)
        
        return(paste0(used_gb, "/", total_gb, " GB (", percent, "%)"))
      }
    }
    return("System info unavailable")
  }, error = function(e) {
    return("Error reading system memory")
  })
}

log_resource_usage <- function(stage, uniprot_id = NULL) {
  timestamp <- Sys.time()
  r_memory <- get_current_memory_usage()
  sys_memory <- get_system_memory_info()
  
  log_entry <- paste0(
    "[", format(timestamp, "%Y-%m-%d %H:%M:%S"), "] ",
    "Stage: ", stage,
    if (!is.null(uniprot_id)) paste0(" | UniProt: ", uniprot_id) else "",
    " | R Memory: ", r_memory,
    " | System Memory: ", sys_memory
  )
  
  # Log to console
  cat(log_entry, "\n")
  
  # Log to file
  log_file <- "resource_usage.log"
  tryCatch({
    write(log_entry, file = log_file, append = TRUE)
  }, error = function(e) {
    # Ignore file write errors
  })
  
  return(log_entry)
}

# Fetch Pfam domains from InterPro
fetch_pfam_domains <- function(uniprot_id) {
  tryCatch({
    # Use the correct InterPro API URL structure
    url <- paste0("https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/UniProt/", uniprot_id, "/")
    
    # Send request using the working framework
    response <- httr::GET(url, httr::add_headers(Accept = "application/json"))
    
    # Check status and parse
    if (httr::status_code(response) != 200) {
      message("InterPro API returned status: ", httr::status_code(response))
      return(NULL)
    }
    
    json_text <- httr::content(response, as = "text", encoding = "UTF-8")
    json_data <- jsonlite::fromJSON(json_text, simplifyVector = FALSE)
    
    if (is.null(json_data$results) || length(json_data$results) == 0) {
      message("No Pfam domains found for ", uniprot_id)
      return(NULL)
    }
    
    # Extract domain information using the working parsing approach
    domains_list <- list()
    
    for (result in json_data$results) {
      pfam_id <- result$metadata$accession %||% "Unknown"
      pfam_name <- result$metadata$name %||% "Unknown Domain"
      
      for (protein in result$proteins) {
        for (location in protein$entry_protein_locations) {
          for (fragment in location$fragments) {
            start <- as.numeric(fragment$start)
            end <- as.numeric(fragment$end)
            
            domains_list[[length(domains_list) + 1]] <- data.frame(
              domain_id = pfam_id,
              domain_name = pfam_name,
              start = start,
              end = end,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    
    if (length(domains_list) == 0) {
      message("No valid domain locations found for ", uniprot_id)
      return(NULL)
    }
    
    result <- do.call(rbind, domains_list)
    message("Found ", nrow(result), " Pfam domains for ", uniprot_id)
    return(result)
    
  }, error = function(e) {
    message("Error fetching Pfam domains: ", e$message)
    return(NULL)
  })
}

# Fetch CATH domains from TED database
fetch_cath_domains_ted <- function(uniprot_id) {
  tryCatch({
    # TED API v1 endpoint
    url <- paste0("https://ted.cathdb.info/api/v1/uniprot/summary/", uniprot_id)
    
    response <- httr::GET(url, 
                         httr::accept("application/json"),
                         httr::add_headers("User-Agent" = "R/EpicTope Domain Fetcher"),
                         httr::timeout(15))
    
    if (httr::status_code(response) != 200) {
      message("TED API returned status: ", httr::status_code(response))
      return(NULL)
    }
    
    # Use simplifyVector = FALSE to preserve list structure
    json_data <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"), 
                                   simplifyVector = FALSE)
    
    # Parse TED domain response
    domains_list <- list()
    
    if (!is.null(json_data$data) && length(json_data$data) > 0) {
      for (i in seq_along(json_data$data)) {
        domain <- json_data$data[[i]]
        
        ted_id <- domain$ted_id %||% paste0("TED_", i)
        cath_label <- domain$cath_label %||% "Unknown"
        chopping <- domain$chopping %||% ""
        consensus_level <- domain$consensus_level %||% ""
        
        # Parse the chopping field to extract start-end positions
        # Chopping can be like "14-131" or "140-145_169-176_262-454" (multiple segments)
        if (chopping != "") {
          segments <- strsplit(chopping, "_")[[1]]
          
          for (segment in segments) {
            if (grepl("^\\d+-\\d+$", segment)) {
              positions <- strsplit(segment, "-")[[1]]
              start_pos <- as.numeric(positions[1])
              end_pos <- as.numeric(positions[2])
              
              # Create a clean domain name from CATH classification only
              if (length(segments) > 1) {
                # For multi-segment domains, show CATH classification with segment number
                domain_name <- paste0(cath_label, " (seg", which(segments == segment), ")")
              } else {
                # For single-segment domains, show just the CATH classification
                domain_name <- cath_label
              }
              
              domains_list[[length(domains_list) + 1]] <- data.frame(
                domain_id = ted_id,
                domain_name = domain_name,
                start = start_pos,
                end = end_pos,
                database = "TED",
                stringsAsFactors = FALSE
              )
            }
          }
        }
      }
    }
    
    if (length(domains_list) > 0) {
      result <- do.call(rbind, domains_list)
      # Count unique domains, not segments
      unique_domains <- length(unique(result$domain_id))
      total_segments <- nrow(result)
      message("Found ", unique_domains, " CATH domains (", total_segments, " segments) for ", uniprot_id)
      return(result)
    }
    
    return(NULL)
    
  }, error = function(e) {
    message("Error fetching CATH domains: ", e$message)
    return(NULL)
  })
}

ui <- dashboardPage(
  dashboardHeader(
    title = tagList(
      span(class = "logo-lg", 
           tags$i(class = "fa fa-dna", style = "margin-right: 10px;"),
           "EpicTope"),
      span(class = "logo-mini", tags$i(class = "fa fa-dna"))
    ),
    titleWidth = 300
  ),
  dashboardSidebar(
    width = 300,
    div(class = "sidebar-content",
        div(class = "input-section",
            h4("Input Parameters", class = "section-header"),
            div(class = "form-group-custom",
                textInput("uniprot_id", 
                         label = div(icon("id-card"), "UniProt ID:"), 
                         placeholder = "Enter UniProt Accession (e.g., Q9W7E7)"),
                bsTooltip("uniprot_id", "Enter a valid UniProt accession", 
                         placement = "right", trigger = "hover")
            ),
            
            # Advanced Parameters - Collapsible Section
            div(class = "advanced-params-section",
                div(class = "advanced-params-header",
                    actionButton("toggle_advanced", 
                               label = tagList(icon("cog"), "Advanced Parameters"),
                               class = "btn btn-link advanced-toggle",
                               style = "color: #3498db; font-weight: 500; padding: 8px 0; border: none; background: none; text-decoration: none;"),
                    tags$small("(Species Selection)", style = "color: #7f8c8d; margin-left: 10px;")
                ),
                
                conditionalPanel(
                    condition = "input.toggle_advanced % 2 == 1",
                    div(class = "advanced-params-content", style = "margin-top: 15px;",
                        div(class = "species-selection-header",
                            h5(icon("dna"), "Species Selection for Multiple Sequence Alignment", 
                               style = "color: #2c3e50; margin-bottom: 10px; font-weight: 600;"),
                            p(class = "species-help-text", 
                              "Select at least 2 species to include in the phylogenetic analysis. All species are selected by default for comprehensive analysis.",
                              style = "font-size: 13px; color: #7f8c8d; margin-bottom: 15px; line-height: 1.4;")
                        ),
                        div(class = "species-selection-container",
                            div(class = "species-controls", style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px; padding: 0 5px;",
                                div(class = "species-counter", 
                                    span("📊 ", style = "margin-right: 6px;"),
                                    textOutput("species_count_display", inline = TRUE),
                                    style = "font-weight: 500; color: #3498db; font-size: 14px; display: flex; align-items: center;"
                                ),
                                div(class = "species-buttons",
                                    actionButton("select_all_species", 
                                               label = tagList(icon("check-square"), "Select All"), 
                                               class = "btn btn-outline-success btn-sm", 
                                               style = "font-size: 12px;")
                                )
                            ),
                            div(class = "species-grid",
                                checkboxGroupInput("selected_species",
                                                 label = NULL,
                                                 choices = list(
                                                     "🐄 Cattle (Bos taurus)" = "bos_taurus",
                                                     "🐕 Dog (Canis lupus familiaris)" = "canis_lupus_familiaris", 
                                                     "🐔 Chicken (Gallus gallus)" = "gallus_gallus",
                                                     "👤 Human (Homo sapiens)" = "homo_sapiens",
                                                     "🐭 Mouse (Mus musculus)" = "mus_musculus",
                                                     "🐟 Pufferfish (Takifugu rubripes)" = "takifugu_rubripes",
                                                     "🐸 Frog (Xenopus tropicalis)" = "xenopus_tropicalis"
                                                 ),
                                                 selected = c("bos_taurus", "canis_lupus_familiaris", "gallus_gallus", 
                                                            "homo_sapiens", "mus_musculus", "takifugu_rubripes", "xenopus_tropicalis"),
                                                 inline = FALSE)
                            )
                        )
                    )
                )
            )
        ),
        
        # Run Analysis Button - Prominent placement
        div(class = "run-analysis-section", style = "margin: 25px 0; text-align: center;",
            actionBttn("run_analysis", 
                      "Run EpicTope", 
                      style = "jelly",
                      color = "primary",
                      size = "lg",
                      class = "run-btn-main")
        ),
        
        div(class = "action-section",
            h4("Additional Actions", class = "section-header"),
            div(class = "action-buttons",
                actionBttn("clear_all", 
                          "Clear All", 
                          icon = icon("trash"),
                          style = "jelly", 
                          color = "danger",
                          size = "md"),
                br(),
                uiOutput("download_ui")
            )
        ),
        div(class = "info-section",
            h4("Information", class = "section-header"),
            div(class = "info-box",
                p("EpicTope identifies which amino acid positions are suitable for epitope tagging by combining the following protein features:"),
                #p("Features include:"),
                tags$ul(
                  tags$li("Sequence conservation"),
                  tags$li("Secondary structure"),
                  tags$li("Relative solvent accessibility"),
                  tags$li("Disordered binding regions"),
                  #tags$li("Disorder prediction")
                )
            )
        )
    )
  ),
  dashboardBody(
    useShinyjs(),
    tags$head(
      tags$script(src = "https://unpkg.com/molstar@3.42.0/build/viewer/molstar.js"),
      tags$link(rel = "stylesheet", href = "https://unpkg.com/molstar@3.42.0/build/viewer/molstar.css"),
      tags$script(src = "https://unpkg.com/3dmol@latest/build/3Dmol-min.js"),  # Fallback 3D viewer
      tags$style(HTML("
        /* Enhanced Professional Styling */
        @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
        
        body {
          font-family: 'Inter', sans-serif;
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        }
        
        .main-header .navbar {
          background: linear-gradient(90deg, #2c3e50 0%, #34495e 100%);
          border-bottom: 3px solid #3498db;
        }
        
        .main-header .logo {
          background: linear-gradient(90deg, #3498db 0%, #2980b9 100%);
          color: white;
          font-weight: 600;
        }
        
        .sidebar {
          background: linear-gradient(180deg, #2c3e50 0%, #34495e 100%);
        }
        
        .sidebar-content { 
          padding: 20px; 
          color: #ecf0f1;
        }
        
        .section-header {
          color: #3498db;
          font-weight: 600;
          border-bottom: 2px solid #3498db;
          padding-bottom: 8px;
          margin-bottom: 15px;
          font-size: 16px;
        }
        
        .form-group-custom {
          margin-bottom: 20px;
        }
        
        .form-group-custom label {
          color: #ecf0f1;
          font-weight: 500;
          margin-bottom: 8px;
        }
        
        .upload-section {
          background: rgba(52, 73, 94, 0.3);
          padding: 15px;
          border-radius: 8px;
          border: 1px solid #3498db;
          margin: 10px 0;
        }
        
        .species-selection-container {
          background: linear-gradient(145deg, #f8f9fa, #e9ecef);
          border: 2px solid #dee2e6;
          border-radius: 12px;
          padding: 20px;
          box-shadow: 0 4px 6px rgba(0, 0, 0, 0.07);
          transition: all 0.3s ease;
        }
        
        .species-selection-container:hover {
          box-shadow: 0 6px 12px rgba(0, 0, 0, 0.1);
          border-color: #3498db;
        }
        
        .species-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
          gap: 8px;
          margin-top: 10px;
        }
        
        .species-grid .checkbox {
          margin-bottom: 0;
          padding: 8px 12px;
          background: white;
          border-radius: 8px;
          border: 1px solid #e0e6ed;
          transition: all 0.2s ease;
        }
        
        .species-grid .checkbox:hover {
          background: #f0f7ff;
          border-color: #3498db;
          transform: translateY(-1px);
        }
        
        .species-grid .checkbox input[type='checkbox']:checked + span {
          color: #2c3e50;
          font-weight: 600;
        }
        
        .species-grid label {
          font-size: 14px;
          font-weight: 500;
          color: #495057;
          cursor: pointer;
          margin: 0;
          display: flex;
          align-items: center;
          width: 100%;
          padding: 0;
        }
        
        .species-grid input[type='checkbox'] {
          margin-right: 10px;
          transform: scale(1.2);
          accent-color: #3498db;
        }
        
        .species-counter {
          display: flex;
          align-items: center;
          font-weight: 500;
          color: #3498db;
        }
        
        .advanced-params-section {
          margin: 20px 0;
          border: 1px solid #e0e6ed;
          border-radius: 8px;
          background: #f8f9fa;
        }
        
        .advanced-params-header {
          padding: 12px 16px;
          border-bottom: 1px solid #e0e6ed;
          background: linear-gradient(145deg, #f1f3f4, #e9ecef);
          border-radius: 8px 8px 0 0;
          display: flex;
          align-items: center;
        }
        
        .advanced-toggle {
          font-size: 14px;
          transition: all 0.3s ease;
        }
        
        .advanced-toggle:hover {
          color: #2980b9 !important;
          text-decoration: none;
        }
        
        .advanced-params-content {
          padding: 20px;
          background: white;
          border-radius: 0 0 8px 8px;
          border-top: 1px solid #e0e6ed;
        }
        
        .run-analysis-section {
          background: linear-gradient(145deg, #3498db, #2980b9);
          border-radius: 12px;
          padding: 20px;
          box-shadow: 0 4px 8px rgba(52, 152, 219, 0.3);
          border: 2px solid #2980b9;
        }
        
        .run-btn-main {
          font-size: 16px !important;
          font-weight: 600 !important;
          padding: 12px 30px !important;
          box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
          transition: all 0.3s ease;
        }
        
        .run-btn-main:hover {
          transform: translateY(-2px);
          box-shadow: 0 6px 12px rgba(0, 0, 0, 0.15);
        }
        
        .action-buttons { 
          display: flex; 
          flex-direction: column; 
          gap: 15px; 
        }
        
        .info-section {
          margin-top: 30px;
          padding-top: 20px;
          border-top: 1px solid #34495e;
        }
        
        .info-box {
          background: rgba(52, 152, 219, 0.1);
          padding: 15px;
          border-radius: 8px;
          border-left: 4px solid #3498db;
          font-size: 13px;
          line-height: 1.5;
        }
        
        .info-box p { margin-bottom: 10px; }
        .info-box ul { margin-left: 15px; }
        .info-box li { margin-bottom: 5px; }
        
        .content-wrapper { 
          background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
          min-height: 100vh;
        }
        
        .box { 
          border-radius: 10px; 
          box-shadow: 0 4px 15px rgba(0,0,0,0.1);
          border: none;
          transition: transform 0.2s ease, box-shadow 0.2s ease;
        }
        
        .box:hover {
          transform: translateY(-2px);
          box-shadow: 0 8px 25px rgba(0,0,0,0.15);
        }
        
        .box-header {
          background: linear-gradient(90deg, #3498db 0%, #2980b9 100%);
          color: white;
          font-weight: 600;
          border-radius: 10px 10px 0 0;
        }
        
        .box-header .box-title {
          font-size: 18px;
          font-weight: 600;
        }
        
        #status { 
          margin-bottom: 15px; 
          font-weight: 600;
          padding: 15px;
          border-radius: 8px;
          background: #f8f9fa;
          border-left: 4px solid #28a745;
        }
        
        .sequence-box { 
          font-family: 'Courier New', monospace; 
          word-break: break-all; 
          overflow-y: auto; 
          max-height: 200px;
          background: #f8f9fa;
          padding: 15px;
          border-radius: 8px;
          border: 1px solid #dee2e6;
          line-height: 1.4;
          font-size: 12px;
        }
        
        .nav-tabs-custom > .nav-tabs {
          border-bottom: 2px solid #3498db;
          background: white;
          border-radius: 8px 8px 0 0;
        }
        
        .nav-tabs-custom > .nav-tabs > li.active > a {
          background: #3498db;
          color: white;
          font-weight: 600;
        }
        
        .progress-container {
          position: fixed;
          top: 0;
          left: 0;
          width: 100%;
          z-index: 9999;
          background: rgba(0,0,0,0.8);
          height: 100vh;
          display: none;
        }
        
        .summary-card {
          background: white;
          border-radius: 10px;
          padding: 20px;
          margin-bottom: 20px;
          box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        
        .metric-box {
          text-align: center;
          padding: 15px;
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          color: white;
          border-radius: 8px;
          margin-bottom: 10px;
        }
        
        .metric-value {
          font-size: 24px;
          font-weight: 700;
          display: block;
        }
        
        .metric-label {
          font-size: 12px;
          opacity: 0.9;
          text-transform: uppercase;
          letter-spacing: 0.5px;
        }
        
        .dataTables_wrapper {
          margin-top: 20px;
        }
        
        .loading-spinner {
          text-align: center;
          padding: 50px;
          color: #3498db;
          font-size: 18px;
        }
        
        .alert-info {
          background: linear-gradient(90deg, #d1ecf1 0%, #bee5eb 100%);
          border-color: #3498db;
          color: #0c5460;
          border-radius: 8px;
          border-left: 4px solid #3498db;
        }
        
        /* MSA Viewer Styles */
        .msa-row {
          font-family: 'Courier New', monospace;
          font-size: 12px;
          line-height: 1.4;
          margin-bottom: 2px;
          white-space: nowrap;
          overflow: visible;
        }
        
        .msa-block {
          margin-bottom: 15px;
          border-bottom: 1px solid #eee;
          padding-bottom: 10px;
          overflow-x: visible;
        }
        
        .msa-ruler {
          margin-bottom: 5px;
          font-family: 'Courier New', monospace;
          font-size: 10px;
          overflow: visible;
        }
        
        .msa-ruler-label {
          font-weight: bold;
          color: #2c3e50;
        }
        
        .msa-ruler-tick {
          color: #bbb;
        }
        
        .msa-sequence {
          display: inline-block;
          margin-left: 5px;
        }
        
        .msa-label {
          display: inline-block;
          font-weight: bold;
          color: #2c3e50;
          font-size: 11px;
          vertical-align: top;
        }
        
        .msa-position {
          display: inline-block;
          width: 12px;
          text-align: center;
          padding: 1px 2px;
          margin: 0 1px;
          border-radius: 2px;
          transition: all 0.2s ease;
          font-weight: bold;
          font-size: 11px;
          vertical-align: top;
        }
        
        .msa-position.highlighted {
          background-color: #e74c3c !important;
          color: white !important;
          font-weight: bold;
          transform: scale(1.2);
          box-shadow: 0 4px 8px rgba(231, 76, 60, 0.6), 0 0 0 3px rgba(231, 76, 60, 0.3);
          z-index: 100;
          position: relative;
          border: 3px solid #c0392b !important;
          border-radius: 4px;
          transition: all 0.2s ease;
          animation: highlight-pulse 1s ease-in-out;
        }
        
        @keyframes highlight-pulse {
          0% { transform: scale(1); }
          50% { transform: scale(1.3); }
          100% { transform: scale(1.2); }
        }
        
        /* Responsive MSA adjustments */
        @media (max-width: 768px) {
          .msa-row {
            font-size: 10px;
          }
          
          .msa-position {
            width: 10px;
            font-size: 9px;
            padding: 1px;
            margin: 0;
          }
          
          .msa-label {
            font-size: 9px;
          }
        }
        
        @media (max-width: 480px) {
          .msa-row {
            font-size: 9px;
          }
          
          .msa-position {
            width: 8px;
            font-size: 8px;
            padding: 0;
          }
          
          .msa-label {
            font-size: 8px;
          }
          
          .msa-ruler-label, .msa-ruler-tick {
            font-size: 7px;
          }
        }
        
        /* Hover information panel */
        .hover-info {
          position: absolute;
          top: 10px;
          right: 10px;
          background: rgba(0,0,0,0.8);
          color: white;
          padding: 10px;
          border-radius: 5px;
          font-size: 12px;
          display: none;
        }
        
        /* Mol* Integration Styles */
        #structure-viewer .msp-plugin {
          border-radius: 8px;
          overflow: hidden;
        }
        
        #structure-viewer .msp-layout-main {
          background: #ffffff;
        }
        
        #structure-controls {
          pointer-events: auto;
        }
        
        #structure-controls .btn {
          margin-left: 5px;
          opacity: 0.9;
          transition: opacity 0.3s ease;
        }
        
        #structure-controls .btn:hover {
          opacity: 1;
        }
        
        /* Ensure Mol* canvas fills container */
        #structure-viewer canvas {
          width: 100% !important;
          height: 100% !important;
        }
      "))
    ),
    
    # JavaScript for interactive features
    tags$script(HTML("
      // Global variables for interactive features
      var currentHighlightedPosition = null;
      var msaData = null;
      var structureData = null;
      
      // Function to highlight position across all visualizations
      function highlightPosition(position) {
        if (currentHighlightedPosition) {
          // Remove previous highlights
          $('.msa-position').removeClass('highlighted');
        }
        
        if (position) {
          // Add new highlights
          $('.msa-position[data-position=\"' + position + '\"]').addClass('highlighted');
          
          // Update structure highlight (placeholder for now)
          console.log('Highlighting position ' + position + ' in structure');
          
          currentHighlightedPosition = position;
        }
      }
      
      // Function to clear all highlights
      function clearHighlights() {
        $('.msa-position').removeClass('highlighted');
        currentHighlightedPosition = null;
        console.log('All highlights cleared');
      }
      
      // Species selection validation
      function validateSpeciesSelection() {
        const checkboxes = document.querySelectorAll('input[name=\"selected_species\"]:checked');
        const count = checkboxes.length;
        
        // Update the visual counter
        const counterElement = document.querySelector('.species-counter');
        if (counterElement) {
          const textElement = counterElement.querySelector('#species_count_display');
          if (textElement) {
            textElement.textContent = count + ' of 7 species selected';
          }
        }
        
        return count >= 2;
      }
      
      // Prevent unchecking when only 2 species remain
      $(document).on('change', 'input[name=\"selected_species\"]', function() {
        const checkedBoxes = $('input[name=\"selected_species\"]:checked');
        
        if (checkedBoxes.length < 2) {
          // If trying to uncheck and would result in < 2, prevent it
          $(this).prop('checked', true);
          
          // Show warning
          if (typeof Shiny !== 'undefined') {
            Shiny.notifications.show({
              html: '⚠️ At least 2 species must remain selected for MSA analysis.',
              type: 'warning',
              duration: 3000
            });
          }
        }
        
        validateSpeciesSelection();
      });
      
      // Function to attach plot event listeners
      function attachPlotEventListeners() {
        const scorePlotContainer = document.getElementById('score_plot');
        if (scorePlotContainer && !scorePlotContainer.hasAttribute('data-listeners-attached')) {
          let clearHighlightTimeout;
          
          scorePlotContainer.addEventListener('mouseleave', function() {
            // Clear any existing timeout
            if (clearHighlightTimeout) clearTimeout(clearHighlightTimeout);
            
            // Add a small delay to prevent flickering when moving between elements
            clearHighlightTimeout = setTimeout(() => {
              console.log('Plot mouseleave - clearing highlights');
              clearHighlights();
              if (structureViewer) {
                console.log('Calling structure clearHighlight from mouseleave');
                structureViewer.clearHighlight();
              }
            }, 100);
          });
          
          // Cancel clearing if mouse re-enters quickly
          scorePlotContainer.addEventListener('mouseenter', function() {
            if (clearHighlightTimeout) {
              clearTimeout(clearHighlightTimeout);
              clearHighlightTimeout = null;
            }
          });
          
          // Also add event listener to the plotly div inside
          const plotlyDiv = scorePlotContainer.querySelector('.plotly');
          if (plotlyDiv) {
            plotlyDiv.addEventListener('mouseleave', function(e) {
              // Only clear if really leaving the plot area
              const rect = plotlyDiv.getBoundingClientRect();
              const x = e.clientX;
              const y = e.clientY;
              if (x < rect.left || x > rect.right || y < rect.top || y > rect.bottom) {
                console.log('Plotly div mouseleave - clearing highlights');
                clearHighlights();
                if (structureViewer) {
                  structureViewer.clearHighlight();
                }
              }
            });
          }
          
          scorePlotContainer.setAttribute('data-listeners-attached', 'true');
          console.log('Event listeners attached to score plot');
        }
      }
      
      // Shiny message handlers
      Shiny.addCustomMessageHandler('highlight-position', function(message) {
        highlightPosition(message.position);
      });
      
      Shiny.addCustomMessageHandler('clear-highlights', function(message) {
        clearHighlights();
      });
      
      Shiny.addCustomMessageHandler('update-msa', function(message) {
        msaData = message;
        renderMSA(message);
      });
      
      // Function to render MSA
      function renderMSA(data) {
        var container = document.getElementById('msa-viewer');
        if (!container) return;
        
        // Force container to recalculate its width
        container.style.display = 'none';
        container.offsetHeight; // Trigger reflow
        container.style.display = '';
        
        // MSA display settings - make responsive based on container width
        var containerWidth = container.offsetWidth || container.parentElement.offsetWidth;
        var labelWidth = Math.min(100, containerWidth * 0.15); // Dynamic label width, max 15% of container
        var residueWidth = 12; // Reduced width per amino acid position
        var padding = 20; // Account for padding and potential scrollbar
        var availableWidth = containerWidth - labelWidth - padding;
        var maxResidues = Math.floor(availableWidth / residueWidth);
        
        // Set responsive residues per line with better scaling
        var residuesPerLine;
        if (containerWidth < 400) {
          residuesPerLine = Math.max(10, Math.min(maxResidues, 25)); // Very small screens
        } else if (containerWidth < 600) {
          residuesPerLine = Math.max(15, Math.min(maxResidues, 40)); // Small screens
        } else if (containerWidth < 900) {
          residuesPerLine = Math.max(25, Math.min(maxResidues, 60)); // Medium screens
        } else {
          residuesPerLine = Math.max(40, Math.min(maxResidues, 80)); // Large screens
        }
        
        // Adjust label frequency and size based on available space and screen size
        var labelEvery, fontSize, labelFontSize, residueWidthActual;
        if (containerWidth < 480) {
          // Very small screens (mobile)
          labelEvery = 5;
          fontSize = '8px';
          labelFontSize = '7px';
          residueWidthActual = 8;
          residuesPerLine = Math.max(10, Math.min(Math.floor(availableWidth / residueWidthActual), 25));
        } else if (containerWidth < 768) {
          // Small screens (tablet portrait)
          labelEvery = 5;
          fontSize = '9px';
          labelFontSize = '8px';
          residueWidthActual = 10;
          residuesPerLine = Math.max(15, Math.min(Math.floor(availableWidth / residueWidthActual), 40));
        } else if (containerWidth < 900) {
          // Medium screens (tablet landscape)
          labelEvery = 10;
          fontSize = '10px';
          labelFontSize = '9px';
          residueWidthActual = 12;
          residuesPerLine = Math.max(25, Math.min(Math.floor(availableWidth / residueWidthActual), 60));
        } else {
          // Large screens (desktop)
          labelEvery = 20;
          fontSize = '11px';
          labelFontSize = '10px';
          residueWidthActual = 12;
          residuesPerLine = Math.max(40, Math.min(Math.floor(availableWidth / residueWidthActual), 80));
        }
        
        var html = '';
        var maxLength = Math.max(...data.sequences.map(seq => seq.sequence.length));
        
        // Split MSA into blocks
        for (var blockStart = 0; blockStart < maxLength; blockStart += residuesPerLine) {
          var blockEnd = Math.min(blockStart + residuesPerLine, maxLength);
          
          // Add position ruler
          html += '<div class=\"msa-block\" style=\"margin-bottom: 15px;\">';
          html += '<div class=\"msa-ruler\" style=\"margin-bottom: 5px;\">';
          html += '<span class=\"msa-label\" style=\"display: inline-block; width: ' + labelWidth + 'px; font-size: ' + labelFontSize + ';\"></span>';
          
          for (var pos = blockStart; pos < blockEnd; pos++) {
            var displayPos = pos + 1;
            var positionStyle = 'display: inline-block; width: ' + residueWidthActual + 'px; text-align: center; font-size: ' + labelFontSize + '; margin: 0 1px;';
            
            if (displayPos % labelEvery === 0 || displayPos === 1) {
              html += '<span class=\"msa-ruler-label\" style=\"' + positionStyle + ' color: #666; font-weight: bold;\">' + displayPos + '</span>';
            } else {
              html += '<span class=\"msa-ruler-tick\" style=\"' + positionStyle + ' color: #ccc;\">|</span>';
            }
          }
          html += '</div>';
          
          // Add sequences for this block
          data.sequences.forEach(function(seq, index) {
            html += '<div class=\"msa-row\">';
            // Clean sequence name by removing data/CDS/ prefix
            var cleanName = seq.name.replace(/^data\\/CDS\\//, '');
            html += '<span class=\"msa-label\" style=\"display: inline-block; width: ' + labelWidth + 'px; font-weight: bold; color: #2c3e50; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-size: ' + fontSize + ';\">' + cleanName + '</span>';
            html += '<span class=\"msa-sequence\">';
            
            for (var i = blockStart; i < blockEnd && i < seq.sequence.length; i++) {
              var aa = seq.sequence[i];
              var bgColor = getClustalColor(aa, data.sequences, i);
              var margin = containerWidth < 480 ? '0' : '0 1px';
              var padding = containerWidth < 480 ? '0' : '1px 2px';
              html += '<span class=\"msa-position\" data-position=\"' + (i + 1) + '\" style=\"display: inline-block; width: ' + residueWidthActual + 'px; text-align: center; padding: ' + padding + '; margin: ' + margin + '; border-radius: 2px; transition: background-color 0.2s ease; background-color: ' + bgColor + '; color: ' + getTextColor(bgColor) + '; font-size: ' + fontSize + '; font-weight: bold;\">' + aa + '</span>';
            }
            
            html += '</span></div>';
          });
          
          html += '</div>';
        }
        
        container.innerHTML = html;
        
        // Hide loading spinner
        document.getElementById('msa-loading').style.display = 'none';
      }
      
      // Function to get CLUSTAL coloring scheme
      function getClustalColor(aa, sequences, position) {
        // Get all amino acids at this position
        var columnAAs = sequences.map(seq => seq.sequence[position] || '-');
        var uniqueAAs = [...new Set(columnAAs)];
        
        // If all are the same (fully conserved), use strong color
        if (uniqueAAs.length === 1 && uniqueAAs[0] !== '-') {
          return getClustalAAColor(aa, 'strong');
        }
        
        // Calculate conservation level
        var maxCount = 0;
        var aaCount = {};
        columnAAs.forEach(function(residue) {
          aaCount[residue] = (aaCount[residue] || 0) + 1;
          maxCount = Math.max(maxCount, aaCount[residue]);
        });
        
        var conservationLevel = maxCount / sequences.length;
        
        // Strong conservation (>80%)
        if (conservationLevel > 0.8) {
          return getClustalAAColor(aa, 'strong');
        }
        // Weak conservation (>60%)
        else if (conservationLevel > 0.6) {
          return getClustalAAColor(aa, 'weak');
        }
        // Check for similar properties
        else {
          var similarGroups = [
            ['A', 'I', 'L', 'M', 'F', 'W', 'V'], // Hydrophobic
            ['K', 'R'],                            // Positive
            ['E', 'D'],                            // Negative
            ['N', 'Q'],                            // Amide
            ['S', 'T'],                            // Hydroxyl
            ['C'],                                 // Sulfur
            ['G'],                                 // Flexible
            ['P'],                                 // Proline
            ['H', 'Y']                             // Aromatic
          ];
          
          var group = similarGroups.find(g => g.includes(aa));
          if (group) {
            var groupCount = columnAAs.filter(residue => group.includes(residue)).length;
            var groupConservation = groupCount / sequences.length;
            if (groupConservation > 0.6) {
              return getClustalAAColor(aa, 'weak');
            }
          }
          
          return '#FFFFFF'; // No conservation
        }
      }
      
      // CLUSTAL color scheme
      function getClustalAAColor(aa, intensity) {
        var colors = {
          strong: {
            'A': '#80a0f0', 'I': '#80a0f0', 'L': '#80a0f0', 'M': '#80a0f0', 'F': '#80a0f0', 'W': '#80a0f0', 'V': '#80a0f0', // Blue - Hydrophobic
            'K': '#f01505', 'R': '#f01505',         // Red - Positive
            'E': '#c048c0', 'D': '#c048c0',         // Magenta - Negative  
            'N': '#00ff00', 'Q': '#00ff00',         // Green - Amide
            'S': '#00ff00', 'T': '#00ff00',         // Green - Hydroxyl
            'C': '#f08080',                         // Pink - Sulfur
            'G': '#f09048',                         // Orange - Flexible
            'P': '#ffff00',                         // Yellow - Proline
            'H': '#15a4a4', 'Y': '#15a4a4',        // Cyan - Aromatic
            '-': '#FFFFFF'                          // White - Gap
          },
          weak: {
            'A': '#c0d0ff', 'I': '#c0d0ff', 'L': '#c0d0ff', 'M': '#c0d0ff', 'F': '#c0d0ff', 'W': '#c0d0ff', 'V': '#c0d0ff', // Light Blue
            'K': '#ff8a80', 'R': '#ff8a80',         // Light Red
            'E': '#e0a0e0', 'D': '#e0a0e0',         // Light Magenta
            'N': '#80ff80', 'Q': '#80ff80',         // Light Green
            'S': '#80ff80', 'T': '#80ff80',         // Light Green
            'C': '#ffb0b0',                         // Light Pink
            'G': '#ffb080',                         // Light Orange
            'P': '#ffff80',                         // Light Yellow
            'H': '#80d4d4', 'Y': '#80d4d4',        // Light Cyan
            '-': '#FFFFFF'                          // White
          }
        };
        
        return colors[intensity][aa] || '#FFFFFF';
      }
      
      // Function to determine text color based on background
      function getTextColor(bgColor) {
        // Simple algorithm: if background is dark, use white text, otherwise black
        var hex = bgColor.replace('#', '');
        var r = parseInt(hex.substr(0, 2), 16);
        var g = parseInt(hex.substr(2, 2), 16);
        var b = parseInt(hex.substr(4, 2), 16);
        var brightness = (r * 299 + g * 587 + b * 114) / 1000;
        return brightness > 128 ? '#000000' : '#FFFFFF';
      }
      
      // Make MSA responsive to window resize
      window.addEventListener('resize', function() {
        if (msaData) {
          renderMSA(msaData);
        }
      });
      
    ")),
    
    # Mol* (molstar) 3D Structure Viewer Integration
    tags$script(HTML("
      // Mol* Structure Viewer Integration
      class MolstarStructureViewer {
        constructor(containerId) {
          this.container = document.getElementById(containerId);
          this.plugin = null;
          this.currentHighlight = null;
          this.isSpinning = false;
          this.structureData = null;
          this.residueMap = new Map();
          
          this.init();
        }
        
        async init() {
          if (!this.container) {
            console.error('Structure viewer container not found!');
            return;
          }
          
          console.log('Initializing Mol* viewer...');
          
          try {
            // Clear loading message
            this.container.innerHTML = '';
            
            // Check if molstar is loaded
            if (typeof molstar === 'undefined') {
              throw new Error('Mol* library not loaded');
            }
            
            // Create Mol* viewer with proper API
            const spec = {
              target: this.container,
              layoutIsExpanded: false,
              layoutShowControls: false,
              layoutShowRemoteState: false,
              layoutShowSequence: false,
              layoutShowLog: false,
              layoutShowLeftPanel: false,
              viewportShowExpand: false,
              viewportShowSelectionMode: false,
              viewportShowAnimation: false,
              pdbProvider: 'rcsb',
              emdbProvider: 'rcsb'
            };
            
            // Try different molstar initialization methods
            if (molstar.Viewer && molstar.Viewer.create) {
              this.plugin = await molstar.Viewer.create(spec);
            } else if (molstar.createPlugin) {
              this.plugin = await molstar.createPlugin(spec);  
            } else if (window.molstar && window.molstar.Viewer) {
              this.plugin = await window.molstar.Viewer.create(spec);
            } else {
              throw new Error('Mol* API not found. Trying alternative initialization...');
            }
            
            // Setup controls
            this.setupControls();
            
            console.log('Mol* viewer initialized successfully');
            
          } catch (error) {
            console.error('Error initializing Mol* viewer:', error);
            // Fallback to 3Dmol.js viewer
            this.container.innerHTML = '<div style=\"display: flex; align-items: center; justify-content: center; height: 100%; color: #666; flex-direction: column; padding: 20px;\"><i class=\"fa fa-cube\" style=\"font-size: 48px; margin-bottom: 20px; color: #3498db;\"></i><p><strong>3D Structure Viewer</strong></p><p style=\"font-size: 12px; text-align: center; color: #999;\">Using fallback viewer</p><p style=\"font-size: 11px; color: #999; margin-top: 10px;\">Mol* error: ' + error.message + '</p></div>';
            this.initFallbackViewer();
          }
        }
        
        initFallbackViewer() {
          // Initialize 3Dmol.js as fallback
          if (typeof $3Dmol !== 'undefined') {
            console.log('Initializing 3Dmol.js fallback viewer');
            this.container.innerHTML = '';
            this.viewer3d = $3Dmol.createViewer(this.container, {
              defaultcolors: $3Dmol.rasmolElementColors
            });
            this.viewer3d.setBackgroundColor(0xffffff);
            this.useFallback = true;
            this.setupControls();
          } else {
            console.log('No fallback viewer available');
          }
        }
        
        setupControls() {
          const resetBtn = document.getElementById('reset-view');
          const spinBtn = document.getElementById('toggle-spin');
          
          if (resetBtn) {
            resetBtn.addEventListener('click', () => this.resetView());
          }
          
          if (spinBtn) {
            spinBtn.addEventListener('click', () => this.toggleSpin());
          }
        }
        
        async loadStructure(pdbData) {
          console.log('Loading structure...', this.useFallback ? '(using fallback)' : '(using Mol*)');
          
          if (this.useFallback && this.viewer3d) {
            // Use 3Dmol.js fallback
            try {
              this.viewer3d.clear();
              this.viewer3d.addModel(pdbData, 'pdb');
              this.viewer3d.setStyle({}, {cartoon: {color: 'spectrum'}});
              this.viewer3d.zoomTo();
              this.viewer3d.render();
              console.log('Structure loaded successfully with 3Dmol.js');
            } catch (error) {
              console.error('Error loading structure with 3Dmol.js:', error);
            }
            return;
          }
          
          if (!this.plugin) {
            console.error('Mol* plugin not initialized');
            return;
          }
          
          try {
            console.log('Loading structure into Mol*...');
            
            // Clear existing structure
            await this.plugin.clear();
            
            // Parse and load PDB data
            const data = await this.plugin.builders.data.rawData({
              data: pdbData,
              label: 'AlphaFold Structure'
            }, { state: { isGhost: false } });
            
            const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'pdb');
            const model = await this.plugin.builders.structure.createModel(trajectory);
            const structure = await this.plugin.builders.structure.createStructure(model);
            
            // Create cartoon representation with B-factor coloring
            const cartoonRepr = await this.plugin.builders.structure.representation.addRepresentation(structure, {
              type: 'cartoon',
              colorTheme: { name: 'uncertainty', params: {} }, // AlphaFold confidence coloring
              sizeTheme: { name: 'uniform', params: {} }
            });
            
            // Store structure reference
            this.structureData = { structure, cartoonRepr };
            
            // Build residue map for highlighting
            await this.buildResidueMap(structure);
            
            // Focus on structure
            await this.plugin.builders.camera.focusLoci(structure.root);
            
            console.log('Structure loaded successfully');
            
          } catch (error) {
            console.error('Error loading structure:', error);
            // Show error in container
            this.container.innerHTML = '<div style=\"display: flex; align-items: center; justify-content: center; height: 100%; color: red; flex-direction: column;\"><i class=\"fa fa-exclamation-triangle\" style=\"font-size: 48px; margin-bottom: 20px;\"></i><p>Error loading structure</p><p style=\"font-size: 12px;\">' + error.message + '</p></div>';
          }
        }
        
        async buildResidueMap(structure) {
          // Build a map of residue positions for highlighting
          if (!structure || !structure.root) return;
          
          try {
            // For now, we'll implement a simple highlighting approach
            console.log('Building residue map for highlighting');
            
          } catch (error) {
            console.error('Error building residue map:', error);
          }
        }
        
        async highlightResidue(position) {
          console.log('Highlighting residue at position:', position);
          
          if (this.useFallback && this.viewer3d) {
            // Use 3Dmol.js highlighting
            try {
              // Clear previous highlights first
              await this.clearHighlight();
              
              console.log('Applying 3Dmol.js highlight to residue:', position);
              
              // Add prominent highlight with both stick and sphere representation
              this.viewer3d.addStyle({resi: position}, {
                stick: {color: 'red', radius: 1.0},
                sphere: {color: 'red', radius: 2.0, opacity: 0.8}
              });
              
              // Add a label to make it even more prominent
              this.viewer3d.addLabel('Residue ' + position, 
                {resi: position, atom: 'CA'}, 
                {
                  fontSize: 12, 
                  fontColor: 'white', 
                  backgroundColor: 'red',
                  backgroundOpacity: 0.8,
                  borderColor: 'black',
                  borderThickness: 1,
                  borderOpacity: 1.0
                }
              );
              
              this.viewer3d.render();
              this.currentHighlight = position;
              console.log('3Dmol.js highlighting applied successfully');
              
            } catch (error) {
              console.error('Error highlighting with 3Dmol.js:', error);
            }
            return;
          }
          
          if (!this.plugin || !this.structureData) {
            console.log('Cannot highlight: plugin or structure not ready');
            return;
          }
          
          try {
            // Clear previous highlight
            await this.clearHighlight();
            
            console.log('Highlighting residue at position:', position, 'using Mol*');
            
            // Get the structure
            const structure = this.plugin.managers.structure.hierarchy.current.structures[0];
            if (!structure) {
              console.log('No structure available for highlighting');
              return;
            }
            
            // Create selection for the specific residue
            const selection = Script.getStructureSelection(Q => Q.struct.generator.atomGroups({
              'residue-test': Q.core.rel.eq([position, Q.struct.atomProperty.macromolecular.auth_seq_id()])
            }), structure.cell);
            
            // Create highlight representation
            const highlightRepr = await this.plugin.builders.structure.representation.addRepresentation(structure, {
              type: 'ball-and-stick',
              color: 'red',
              size: 'uniform',
              sizeParams: { value: 2.0 }
            }, { tag: 'highlight-residue-' + position });
            
            // Apply selection to the representation
            await this.plugin.managers.structure.component.updateRepresentations(structure, [{
              type: 'ball-and-stick',
              colorTheme: { name: 'uniform', params: { value: 0xff0000 } }, // Red color
              sizeTheme: { name: 'uniform', params: { value: 2.0 } }
            }], selection);
            
            this.currentHighlight = position;
            console.log('Mol* highlighting applied for residue:', position);
            
          } catch (error) {
            console.error('Error highlighting residue with Mol*:', error);
          }
        }
        
        async clearHighlight() {
          console.log('clearHighlight called, useFallback:', this.useFallback);
          
          if (this.useFallback && this.viewer3d) {
            // Clear 3Dmol.js highlights
            try {
              console.log('Clearing 3Dmol.js highlights...');
              
              // Remove all labels first
              this.viewer3d.removeAllLabels();
              
              // Reset all atoms to default cartoon style
              this.viewer3d.setStyle({}, {cartoon: {color: 'spectrum'}});
              
              // If we have a specific residue that was highlighted, reset it explicitly
              if (this.currentHighlight) {
                console.log('Resetting specific residue:', this.currentHighlight);
                const selector = {resi: this.currentHighlight};
                this.viewer3d.setStyle(selector, {cartoon: {color: 'spectrum'}});
              }
              
              // Force a complete re-render
              this.viewer3d.render();
              this.currentHighlight = null;
              console.log('3Dmol.js highlights cleared successfully');
              
            } catch (error) {
              console.error('Error clearing 3Dmol.js highlights:', error);
            }
            return;
          }
          
          if (!this.plugin) return;
          
          try {
            // Clear any existing highlights in Mol*
            console.log('Clearing Mol* highlights');
            
            // Remove any highlight representations by tag
            if (this.currentHighlight) {
              const highlightTag = 'highlight-residue-' + this.currentHighlight;
              const representations = this.plugin.managers.structure.hierarchy.current.representations;
              
              for (const repr of representations) {
                if (repr.cell.tag === highlightTag) {
                  await this.plugin.managers.structure.hierarchy.remove([repr.cell]);
                }
              }
            }
            
            // Alternative approach: reset all representations to default
            const structure = this.plugin.managers.structure.hierarchy.current.structures[0];
            if (structure) {
              // Remove any custom highlight representations and reset to cartoon
              await this.plugin.managers.structure.component.updateRepresentations(structure, [{
                type: 'cartoon',
                colorTheme: { name: 'chain-id' },
                sizeTheme: { name: 'uniform' }
              }]);
            }
            
            this.currentHighlight = null;
            console.log('Mol* highlights cleared');
            
          } catch (error) {
            console.error('Error clearing Mol* highlight:', error);
          }
        }
        
        resetView() {
          if (this.useFallback && this.viewer3d) {
            try {
              this.viewer3d.zoomTo();
              this.viewer3d.render();
              console.log('View reset (3Dmol.js)');
            } catch (error) {
              console.error('Error resetting view (3Dmol.js):', error);
            }
            return;
          }
          
          if (!this.plugin || !this.structureData) return;
          
          try {
            this.plugin.builders.camera.focusLoci(this.structureData.structure.root);
            console.log('View reset');
          } catch (error) {
            console.error('Error resetting view:', error);
          }
        }
        
        toggleSpin() {
          if (this.useFallback && this.viewer3d) {
            try {
              this.isSpinning = !this.isSpinning;
              const btn = document.getElementById('toggle-spin');
              
              if (this.isSpinning) {
                this.viewer3d.spin('y', 1);
                if (btn) btn.textContent = 'Stop Spin';
              } else {
                this.viewer3d.spin(false);
                if (btn) btn.textContent = 'Toggle Spin';
              }
              
              console.log('Spin toggled (3Dmol.js):', this.isSpinning);
            } catch (error) {
              console.error('Error toggling spin (3Dmol.js):', error);
            }
            return;
          }
          
          if (!this.plugin) return;
          
          try {
            this.isSpinning = !this.isSpinning;
            const btn = document.getElementById('toggle-spin');
            
            if (this.isSpinning) {
              // Enable spinning if supported
              if (this.plugin.canvas3d && this.plugin.canvas3d.props) {
                this.plugin.canvas3d.props.camera.spinning = true;
                this.plugin.canvas3d.commit();
              }
              if (btn) btn.textContent = 'Stop Spin';
            } else {
              // Disable spinning
              if (this.plugin.canvas3d && this.plugin.canvas3d.props) {
                this.plugin.canvas3d.props.camera.spinning = false;
                this.plugin.canvas3d.commit();
              }
              if (btn) btn.textContent = 'Toggle Spin';
            }
            
            console.log('Spin toggled:', this.isSpinning);
            
          } catch (error) {
            console.error('Error toggling spin:', error);
          }
        }
      }
      
      // Global structure viewer instance
      let structureViewer = null;
      
      // Initialize structure viewer when DOM is ready
      document.addEventListener('DOMContentLoaded', function() {
        console.log('DOM loaded, checking 3D viewer availability...');
        
        // Wait for scripts to load
        setTimeout(() => {
          // Try to initialize viewer
          console.log('Attempting to initialize 3D structure viewer...');
          try {
            structureViewer = new MolstarStructureViewer('structure-viewer');
            console.log('Structure viewer initialized successfully');
          } catch (error) {
            console.error('Error initializing structure viewer:', error);
            const container = document.getElementById('structure-viewer');
            if (container) {
              container.innerHTML = '<div style=\"display: flex; align-items: center; justify-content: center; height: 100%; color: #666; flex-direction: column; padding: 20px;\"><i class=\"fa fa-cube\" style=\"font-size: 48px; margin-bottom: 20px; color: #3498db;\"></i><p><strong>3D Structure Viewer</strong></p><p style=\"font-size: 12px;\">Initialization error: ' + error.message + '</p><p style=\"font-size: 11px; color: #999; margin-top: 10px;\">Structure will be displayed here after analysis</p></div>';
            }
          }
        }, 3000); // Wait longer for libraries to load
      });
      
      // Fallback initialization for Shiny
      $(document).ready(function() {
        setTimeout(() => {
          if (!structureViewer) {
            console.log('Fallback initialization of structure viewer...');
            try {
              structureViewer = new MolstarStructureViewer('structure-viewer');
              console.log('Structure viewer initialized via fallback');
            } catch (error) {
              console.error('Fallback initialization failed:', error);
            }
          }
          
          // Add mouse leave event listener to score plot to clear highlights
          const scorePlotContainer = document.getElementById('score_plot');
          if (scorePlotContainer) {
            attachPlotEventListeners();
          }
          
          // Set up MutationObserver to detect when plot is rendered/re-rendered
          const observer = new MutationObserver(function(mutations) {
            mutations.forEach(function(mutation) {
              if (mutation.type === 'childList') {
                const scorePlotContainer = document.getElementById('score_plot');
                if (scorePlotContainer && scorePlotContainer.querySelector('.plotly')) {
                  attachPlotEventListeners();
                }
              }
            });
          });
          
          // Start observing changes to the document
          observer.observe(document.body, {
            childList: true,
            subtree: true
          });
        }, 5000);
      });
      
      // Shiny message handlers for structure
      Shiny.addCustomMessageHandler('load-structure', function(message) {
        console.log('Received load-structure message');
        if (structureViewer && message.pdbData) {
          structureViewer.loadStructure(message.pdbData);
        } else {
          console.log('Structure viewer not ready or no PDB data');
        }
      });
      
      Shiny.addCustomMessageHandler('highlight-structure', function(message) {
        if (structureViewer && message.position) {
          structureViewer.highlightResidue(message.position);
        }
      });
      
      Shiny.addCustomMessageHandler('clear-structure-highlight', function(message) {
        console.log('Received clear-structure-highlight message');
        if (structureViewer) {
          console.log('Structure viewer exists, calling clearHighlight');
          structureViewer.clearHighlight();
        } else {
          console.log('No structure viewer available');
        }
      });
    ")),
    
    # Progress indicator overlay
    div(id = "progress-overlay", class = "progress-container",
        div(style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); text-align: center; color: white;",
            div(class = "loading-spinner",
                tags$i(class = "fa fa-spinner fa-spin fa-3x"),
                h3("Running EpicTope...", style = "margin-top: 20px;")
            )
        )
    ),
    
    # Status bar
    fluidRow(
      column(width = 12,
        div(id = "status-container", style = "margin-bottom: 20px;",
            verbatimTextOutput("status")
        )
      )
    ),
    
    # Main content - Score Plot front and center
    fluidRow(
      column(width = 12,
        box(width = NULL, status = "success", solidHeader = TRUE,
            title = tagList(icon("chart-line"), "EpicTope Score Plot"),
            withSpinner(
              plotlyOutput("score_plot", height = "400px"),
              type = 6, color = "#3498db"
            )
        )
      )
    ),
    
    # Structure and MSA visualization row
    fluidRow(
      column(width = 6,
        box(width = NULL, status = "info", solidHeader = TRUE,
            title = tagList(icon("cube"), "3D Structure (Mol*)"),
            div(id = "structure-container", style = "height: 500px; border: 1px solid #ddd; position: relative;",
                div(id = "structure-viewer", style = "width: 100%; height: 100%; background: #ffffff;",
                    div(style = "display: flex; align-items: center; justify-content: center; height: 100%; color: #666; flex-direction: column;",
                        tags$i(class = "fa fa-dna fa-3x", style = "margin-bottom: 15px; color: #3498db;"),
                        "Loading Mol* viewer..."
                    )
                ),
                div(id = "structure-controls", style = "position: absolute; top: 10px; right: 10px; z-index: 2000;",
                    tags$button("Reset View", id = "reset-view", class = "btn btn-sm btn-secondary", 
                               style = "margin-right: 5px;"),
                    tags$button("Toggle Spin", id = "toggle-spin", class = "btn btn-sm btn-secondary")
                ),
                uiOutput("structure_loading_spinner")
            )
        )
      ),
      column(width = 6,
        box(width = NULL, status = "info", solidHeader = TRUE,
            title = tagList(icon("align-left"), "Multiple Sequence Alignment"),
            div(id = "msa-container", style = "height: 500px; overflow-y: auto; overflow-x: hidden; border: 1px solid #ddd; padding: 10px;",
                div(id = "msa-viewer", style = "font-family: 'Courier New', monospace; font-size: 12px;"),
                uiOutput("msa_loading_spinner")
            )
        )
      )
    ),
    
    # Additional tabs for detailed data
    fluidRow(
      column(width = 12,
        tabsetPanel(id = "data_tabs", type = "tabs",
          # Data Table Tab
          tabPanel("📊 Data Table",
            box(width = NULL, status = "primary", solidHeader = TRUE,
                title = tagList(icon("table"), "Detailed Results"),
                withSpinner(
                  DT::dataTableOutput("results_table"),
                  type = 6, color = "#3498db"
                )
            )
          ),
          
          # Export Tab
          tabPanel("💾 Export",
            fluidRow(
              column(width = 12,
                box(width = NULL, status = "success", solidHeader = TRUE,
                    title = tagList(icon("download"), "Export Options"),
                    div(class = "summary-card",
                        h4("Available Downloads"),
                        p("Export your analysis results in various formats:"),
                        div(class = "action-buttons",
                            uiOutput("export_buttons")
                        )
                    )
                )
              )
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  reactive_data <- reactiveValues(
    pfam_domains = NULL,
    uniprot_id = NULL, 
    res_df = NULL, 
    custom_msa = NULL, 
    protein_sequence = NULL,
    uniprot_data = NULL,
    analysis_complete = FALSE,
    analysis_running = FALSE,
    msa_data = NULL,
    structure_file = NULL
  )
  
  # Conditionally render loading spinner for structure
  output$structure_loading_spinner <- renderUI({
    if (!reactive_data$analysis_running) return(NULL)
    div(id = "structure-loading", class = "loading-spinner", 
        style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);",
        tags$i(class = "fa fa-spinner fa-spin fa-2x"),
        p("Loading 3D structure...", style = "margin-top: 10px;")
    )
  })

  # Render protein structure with custom viewer
  observe({
    if(!is.null(reactive_data$structure_file) && file.exists(reactive_data$structure_file)) {
      # Read PDB file
      pdb_content <- paste(readLines(reactive_data$structure_file), collapse = "\n")
      
      # Send PDB data to client-side structure viewer
      session$sendCustomMessage("load-structure", list(pdbData = pdb_content))
    }
  })

  # Conditionally render loading spinner for MSA
  output$msa_loading_spinner <- renderUI({
    if (!reactive_data$analysis_running) return(NULL)
    div(id = "msa-loading", class = "loading-spinner",
        tags$i(class = "fa fa-spinner fa-spin fa-2x"),
        p("Loading MSA...", style = "margin-top: 10px;")
    )
  })
  
  moving_average <- function(x, window_size) {
          if (window_size %% 2 == 0) stop("window_size must be odd")
          n <- length(x)
          ma <- numeric(n)
          half_window <- floor(window_size / 2)
          # Start edge: progressively increase window size from 4 up to (window_size - 1)
          for (i in 1:half_window) {
            ma[i] <- mean(x[1:(i + half_window)])
          }
          # Middle: full window
          for (i in (half_window + 1):(n - half_window)) {
            ma[i] <- mean(x[(i - half_window):(i + half_window)])
          }
          # End edge: progressively decrease window size from (window_size - 1) to 4
          for (i in (n - half_window + 1):n) {
            start_index <- i - (n - i) - half_window
            ma[i] <- mean(x[start_index:n])
          }
          return(ma)
        }

  # Interactive score plot with domain bars
  output$score_plot <- renderPlotly({
    req(reactive_data$res_df)
    
    # Prepare data like the old version
    res_df <- reactive_data$res_df
    
    # Check if required columns exist
    if(!"position" %in% colnames(res_df)) {
        error_plot <- plot_ly(
            x = c(0, 100), 
            y = c(0, 1), 
            type = "scatter", 
            mode = "lines",
            line = list(color = "transparent"),
            showlegend = FALSE,
            source = "score_plot"
        ) %>%
        layout(
            title = "Error: Position column not found in data.",
            xaxis = list(title = "Amino Acid Position", showgrid = FALSE),
            yaxis = list(title = "Score", showgrid = FALSE)
        )
        return(error_plot)
    }
    
    # Handle score column - check for different possible column names
    score_col <- NULL
    if("min" %in% colnames(res_df)) {
        score_col <- "min"
    } else if("min_val" %in% colnames(res_df)) {
        score_col <- "min_val"
    } else if("score" %in% colnames(res_df)) {
        score_col <- "score"
    } else if("norm_score" %in% colnames(res_df)) {
        score_col <- "norm_score"
    } else {
        error_plot <- plot_ly(
            x = c(0, 100), 
            y = c(0, 1), 
            type = "scatter", 
            mode = "lines",
            line = list(color = "transparent"),
            showlegend = FALSE,
            source = "score_plot"
        ) %>%
        layout(
            title = "Error: No score column found in data.",
            xaxis = list(title = "Amino Acid Position", showgrid = FALSE),
            yaxis = list(title = "Score", showgrid = FALSE)
        )
        return(error_plot);
    }
    
    res_df$id <- as.numeric(res_df$position)
    res_df$min_val <- as.numeric(res_df[[score_col]])
    res_df <- res_df[order(res_df$id),]  # Sort by position
    res_df$min_val_ma <- moving_average(res_df$min_val, 7)  # Apply moving average
    plot_data <- res_df[, c("id", "min_val_ma")]
    
    # Reshape data to long format
    plot_data <- reshape(plot_data,
                         varying = list(names(plot_data)[names(plot_data) != "id"]),
                         v.names = "value",
                         idvar = "id",
                         direction = "long")
    
    plot_data$value <- as.numeric(plot_data$value)
    
    # Define axis labels
    y_range <- round(range(as.numeric(plot_data$value), na.rm = TRUE), 2)
    y_breaks <- seq(from = y_range[1], to = y_range[2], by = (y_range[2] - y_range[1]) / 5)
    x_ticks <- seq(0, max(plot_data$id, na.rm = TRUE), by = 100)
    x_labels <- as.character(x_ticks)
    
    # Get amino acid sequence for hover info
    protein_seq <- reactive_data$protein_sequence
    if(!is.null(protein_seq)) {
        aa_vector <- strsplit(protein_seq, "")[[1]]
        plot_data$amino_acid <- aa_vector[plot_data$id]
    } else {
        plot_data$amino_acid <- "N/A"
    }
    
    # Create hover text with amino acid info
    plot_data$hover_text <- paste0(
        "Position: ", plot_data$id, 
        "<br>Amino Acid: ", plot_data$amino_acid,
        "<br>Score: ", round(plot_data$value, 3)
    )
    
    # Calculate y-axis limits (including domain space if needed)
    y_max <- max(plot_data$value, na.rm = TRUE)
    y_min <- min(plot_data$value, na.rm = TRUE)
    y_range <- y_max - y_min
    
    # Initial y-limits (will be updated if domains are present)
    final_y_max <- y_max + (y_range * 0.05)  # Default padding
    
    # Create the ggplot with hover text
    p <- ggplot(plot_data, aes(x = id, y = as.numeric(value), group = 1, text = hover_text)) +
        geom_line(color = "black") +  # Solid line
        scale_x_continuous("Amino Acid Position", breaks = x_ticks, labels = x_labels) +
        scale_y_continuous("Score", breaks = y_breaks, labels = round(y_breaks, 2), limits = c(y_min, final_y_max)) +
        theme_classic() +  # Use theme_classic() for better ggplotly compatibility
        theme(
            panel.grid = element_blank(),  # Remove grid lines
            axis.line = element_line(color = "black")  # Ensures black borders remain in ggplotly
        ) +
        ggtitle(paste("EpicTope Analysis -", reactive_data$uniprot_id))
    
    # Overlay protein domains if available (both Pfam and CATH)
    all_domains <- data.frame()
    
    # Add Pfam domains
    if (!is.null(reactive_data$pfam_domains) && nrow(reactive_data$pfam_domains) > 0) {
      pfam_data <- reactive_data$pfam_domains
      pfam_data$database <- "Pfam"
      pfam_data$domain_type <- "Functional"
      all_domains <- rbind(all_domains, pfam_data[, c("domain_id", "domain_name", "start", "end", "database", "domain_type")])
    }
    
    # Add CATH domains
    if (!is.null(reactive_data$cath_domains) && nrow(reactive_data$cath_domains) > 0) {
      cath_data <- reactive_data$cath_domains
      cath_data$domain_type <- "Structural"
      all_domains <- rbind(all_domains, cath_data[, c("domain_id", "domain_name", "start", "end", "database", "domain_type")])
    }
    
    if (nrow(all_domains) > 0) {
      # Calculate domain bar positioning with proper stacking
      domain_bar_height <- y_range * 0.04  # Slightly smaller bars
      domain_spacing <- y_range * 0.02     # Space between domain rows
      
      # Separate domains by type for stacking
      pfam_domains <- all_domains[all_domains$database == "Pfam", ]
      cath_domains <- all_domains[all_domains$database == "TED", ]
      
      # Create domain data for geom layers with hover text
      domain_data <- data.frame()
      
      # Process Pfam domains (bottom row)
      if (nrow(pfam_domains) > 0) {
        pfam_y_start <- y_max + (y_range * 0.08)
        pfam_colors <- scales::brewer_pal(type = "qual", palette = "Set1")(nrow(pfam_domains))
        
        for (i in seq_len(nrow(pfam_domains))) {
          d <- pfam_domains[i,]
          domain_data <- rbind(domain_data, data.frame(
            xmin = d$start,
            xmax = d$end,
            ymin = pfam_y_start,
            ymax = pfam_y_start + domain_bar_height,
            domain_id = d$domain_id,
            domain_name = d$domain_name,
            start_pos = d$start,
            end_pos = d$end,
            database = d$database,
            domain_type = d$domain_type,
            color = pfam_colors[i],
            hover_text = paste0(
              "Functional Domain: ", d$domain_name, " (", d$domain_id, ")",
              "<br>Database: Pfam",
              "<br>Range: ", d$start, "-", d$end,
              "<br>Length: ", (d$end - d$start + 1), " amino acids"
            ),
            x_center = (d$start + d$end) / 2,
            y_center = pfam_y_start + domain_bar_height/2,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Process CATH domains (top row)
      if (nrow(cath_domains) > 0) {
        # Position CATH domains above Pfam domains if both exist
        cath_y_start <- if (nrow(pfam_domains) > 0) {
          y_max + (y_range * 0.08) + domain_bar_height + domain_spacing
        } else {
          y_max + (y_range * 0.08)
        }
        
        # Create color mapping based on unique domain IDs
        unique_cath_domains <- unique(cath_domains$domain_id)
        cath_colors <- scales::brewer_pal(type = "qual", palette = "Set2")(length(unique_cath_domains))
        names(cath_colors) <- unique_cath_domains
        
        for (i in seq_len(nrow(cath_domains))) {
          d <- cath_domains[i,]
          # Assign color based on domain ID, not row number
          segment_color <- cath_colors[d$domain_id]
          
          domain_data <- rbind(domain_data, data.frame(
            xmin = d$start,
            xmax = d$end,
            ymin = cath_y_start,
            ymax = cath_y_start + domain_bar_height,
            domain_id = d$domain_id,
            domain_name = d$domain_name,
            start_pos = d$start,
            end_pos = d$end,
            database = d$database,
            domain_type = d$domain_type,
            color = segment_color,
            hover_text = paste0(
              "CATH Structural Domain: ", d$domain_name,
              "<br>AlphaFold Model: ", d$domain_id,
              "<br>Database: CATH/TED", 
              "<br>Range: ", d$start, "-", d$end,
              "<br>Length: ", (d$end - d$start + 1), " amino acids"
            ),
            x_center = (d$start + d$end) / 2,
            y_center = cath_y_start + domain_bar_height/2,
            stringsAsFactors = FALSE
          ))
        }
      }
      
      # Update y-axis limits to accommodate stacked domains
      if (nrow(domain_data) > 0) {
        max_domain_y <- max(domain_data$ymax)
        final_y_max <- max_domain_y + (y_range * 0.05)  # Add some padding above domains
        
        # Update the y-axis scale
        p <- p + scale_y_continuous("Score", 
                                   breaks = y_breaks, 
                                   labels = round(y_breaks, 2), 
                                   limits = c(y_min, final_y_max))
      }
      
      # Add domain rectangles with hover text
      p <- p + geom_rect(
        data = domain_data,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, text = hover_text),
        fill = domain_data$color,
        alpha = 0.7,
        color = "black",
        linewidth = 0.5,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
      
      # Add domain labels (abbreviated for space)
      p <- p + geom_text(
        data = domain_data,
        aes(x = x_center, y = y_center, label = ifelse(nchar(domain_name) > 12, 
                                                      paste0(substr(domain_name, 1, 10), "..."), 
                                                      domain_name)),
        size = 2.5,
        fontface = "bold",
        color = "black",
        inherit.aes = FALSE
      )
    }
    
    # Convert to interactive plotly with custom hover text
    fig <- ggplotly(p, tooltip = "text", source = "score_plot")
    
    # Customize hover behavior for domain bars specifically
    if (nrow(all_domains) > 0) {
      # Find domain rectangle traces and enhance their hover info
      for (i in seq_along(fig$x$data)) {
        trace <- fig$x$data[[i]]
        # Check if this trace represents domain rectangles (geom_rect creates fillcolor traces)
        if (!is.null(trace$type) && trace$type == "scatter" && 
            !is.null(trace$fill) && trace$fill == "tonexty") {
          # This is likely a rectangle trace, enhance its hover template
          if (!is.null(trace$text) && length(trace$text) > 0) {
            # Keep the existing text which should contain our hover_text
            trace$hovertemplate <- "%{text}<extra></extra>"
          }
        }
      }
    }
    
    # Register plotly_hover and plotly_unhover events on the plotly object
    fig <- event_register(fig, 'plotly_hover')
    fig <- event_register(fig, 'plotly_unhover')

    return(fig)
  })
  
  # Listen for hover events on the plot
  observe({
    event_data <- event_data("plotly_hover", source = "score_plot")
    if (!is.null(event_data)) {
      # For ggplotly, the position is in the x coordinate
      # Ensure x is numeric before rounding
      if (!is.null(event_data$x) && is.numeric(event_data$x)) {
        position <- round(event_data$x)
        
        # Debug output
        cat("Hovering over position:", position, "\n")
        
        # Send position to client-side JavaScript
        session$sendCustomMessage("highlight-position", list(position = position))
        
        # Highlight position in custom structure viewer
        if (!is.null(reactive_data$structure_file)) {
          cat("Structure file available, highlighting position:", position, "\n")
          
          # Send highlight message to custom structure viewer
          session$sendCustomMessage("highlight-structure", list(position = position))
        } else {
          cat("No structure file available\n")
        }
      } else {
        cat("Invalid or non-numeric position data\n")
      }
    }
  })
  
  # Clear highlights when not hovering
  observe({
    event_data <- event_data("plotly_unhover", source = "score_plot")
    if (!is.null(event_data)) {
      # Debug output
      cat("Plot unhover event triggered\n")
      
      session$sendCustomMessage("clear-highlights", list())
      
      # Clear structure highlights
      if (!is.null(reactive_data$structure_file)) {
        session$sendCustomMessage("clear-structure-highlight", list())
      }
    }
  })
  
  # Results table output
  output$results_table <- DT::renderDataTable({
    if(is.null(reactive_data$res_df) || nrow(reactive_data$res_df) == 0) {
      return(NULL)
    } else {
      df <- reactive_data$res_df

      # Build display_df in requested column order
      display_df <- data.frame(
        Position = df$position,
        Amino_Acid = if("aa" %in% colnames(df)) as.character(df$aa) else NA,
        Min = if("min" %in% colnames(df)) round(as.numeric(df$min), 4) else NA,
        Normalized_Entropy = if("normalized_entropy" %in% colnames(df)) round(as.numeric(df$normalized_entropy), 4) else NA,
        RSA = if("rsa" %in% colnames(df)) round(as.numeric(df$rsa), 4) else NA,
        Inv_Anchor2 = if("inv_anchor2" %in% colnames(df)) round(as.numeric(df$inv_anchor2), 4) else NA,
        SS = if("ss" %in% colnames(df)) as.character(df$ss) else NA
      )

      DT::datatable(display_df, 
                   options = list(
                     scrollX = TRUE,
                     paging = FALSE,
                     searching = TRUE,
                     ordering = TRUE,
                     dom = 'Bfrtip',
                     buttons = c('copy', 'csv', 'excel')
                   ),
                   class = 'cell-border stripe hover',
                   style = 'bootstrap4')
    }
  })
  
  # Export buttons
  output$export_buttons <- renderUI({
    if(reactive_data$analysis_complete) {
      tagList(
        downloadButton("download_csv", "📄 Download CSV", 
                        class = "btn btn-primary btn-block",
                        style = "margin-bottom: 10px;"),
        downloadButton("download_json", "📋 Download JSON", 
                        class = "btn btn-info btn-block",
                        style = "margin-bottom: 10px;"),
        downloadButton("download_plot", "📊 Download Plot", 
                        class = "btn btn-success btn-block",
                        style = "margin-bottom: 10px;")
      )
    } else {
      div(class = "alert alert-info",
          "🔄 Please run analysis first to enable export options.")
    }
  })
  
  # Species count display
  output$species_count_display <- renderText({
    selected <- input$selected_species
    count <- if(is.null(selected)) 0 else length(selected)
    paste(count, "of 7 species selected")
  })
  
  # Handle species selection buttons
  observeEvent(input$select_all_species, {
    updateCheckboxGroupInput(session, "selected_species", 
                           selected = c("bos_taurus", "canis_lupus_familiaris", "gallus_gallus", 
                                        "homo_sapiens", "mus_musculus", "takifugu_rubripes", "xenopus_tropicalis"))
  })
  
  # Prevent deselecting below minimum and validate species selection
  observe({
    if (!is.null(input$selected_species)) {
      num_selected <- length(input$selected_species)
      
      # If user tries to deselect below  2, prevent it
      if (num_selected < 2 && !is.null(isolate(input$selected_species))) {
        # Get the previous valid selection from reactive data or default to human + mouse
        prev_selection <- isolate(input$selected_species)
        if (length(prev_selection) < 2) {
          # Force selection of human + mouse as fallback
          updateCheckboxGroupInput(session, "selected_species", 
                                 selected = c("homo_sapiens", "mus_musculus"))
          showNotification("⚠️ At least 2 species must remain selected for MSA analysis.", 
                           type = "warning", duration = 3)
        }
      }
      
      # Update status based on selection
      if (num_selected < 2) {
        output$status <- renderText("⚠️ Please select at least 2 species for multiple sequence alignment.")
      } else {
        output$status <- renderText(paste("✅", num_selected, "species selected for MSA."))
      }
    }
  })
  
  # Run analysis
  observeEvent(input$run_analysis, {
    # Log initial resource usage
    log_resource_usage("Analysis Start")
    
    # Reset analysis flags
    reactive_data$analysis_complete <- FALSE
    reactive_data$analysis_running <- TRUE
    
    # Show progress overlay
    runjs("$('#progress-overlay').show();")
    
    # Validate inputs
    req(input$uniprot_id)
    if (input$uniprot_id == "") {
      output$status <- renderText("❌ Please enter a UniProt ID.")
      reactive_data$analysis_running <- FALSE
      runjs("$('#progress-overlay').hide();")
      return()
    }
    
    if (is.null(input$selected_species) || length(input$selected_species) < 2) {
      output$status <- renderText("❌ Please select at least 2 species for multiple sequence alignment.")
      reactive_data$analysis_running <- FALSE
      runjs("$('#progress-overlay').hide();")
      return()
    }
    
    # Store validated ID
    reactive_data$uniprot_id <- input$uniprot_id
    query <- reactive_data$uniprot_id
    
    # Log analysis start with UniProt ID
    log_resource_usage("Analysis Start", query)
    
    # Run analysis with progress indicator
    withProgress(message = "🔄 Running EpicTope analysis...", value = 0, {
      tryCatch({
        # Setup and configuration
        incProgress(0.1, detail = "⚙️ Setting up configuration...")
        log_resource_usage("Setup", query)
        setup_files()
        check_config()
        
        # Fetch UniProt data
        incProgress(0.2, detail = "🔍 Fetching UniProt data...")
        log_resource_usage("UniProt Fetch", query)
        uniprot_fields <- c("accession", "id", "gene_names", "xref_alphafolddb", 
                           "sequence", "organism_name", "organism_id")
        uniprot_data <- query_uniProt(query = query, fields = uniprot_fields)
        
        if (is.null(uniprot_data) || nrow(uniprot_data) == 0) {
          stop("No UniProt data found for the provided ID.")
        }
        
        # Store UniProt data and protein sequence
        reactive_data$uniprot_data <- uniprot_data
        reactive_data$protein_sequence <- uniprot_data$Sequence
        
        # Fetch protein domains
        incProgress(0.25, detail = "🧬 Fetching protein domains...")
        log_resource_usage("Domain Fetch", query)
        
        # Fetch Pfam domains
        pfam_domains <- fetch_pfam_domains(query)
        reactive_data$pfam_domains <- pfam_domains
        
        # Fetch CATH domains
        cath_domains <- fetch_cath_domains_ted(query)
        reactive_data$cath_domains <- cath_domains
        
        # Log domain counts
        pfam_count <- if (!is.null(pfam_domains)) nrow(pfam_domains) else 0
        cath_count <- if (!is.null(cath_domains)) length(unique(cath_domains$domain_id)) else 0
        cat("Domain fetching complete: ", pfam_count, " Pfam domains, ", cath_count, " CATH domains\n")
        log_resource_usage("Domain Fetch Complete", query)
          
        # Process AlphaFold structure
        incProgress(0.3, detail = "🏗️ Fetching AlphaFold structure...")
        log_resource_usage("AlphaFold Fetch", query)
        if (is.na(uniprot_data$AlphaFoldDB)) {
          uniprot_data$AlphaFoldDB <- query
        }
            
        alphafold_file <- fetch_alphafold(gsub(";", "", uniprot_data$AlphaFoldDB))
            
            if (is.na(alphafold_file)) {
                stop(paste("No AlphaFold structure found for", query))
            }
            
            # Process structure data
            incProgress(0.4, detail = "⚗️ Processing structure data...")
            log_resource_usage("Structure Processing", query)
            dssp_res <- dssp_command(alphafold_file)
            dssp_df <- parse_dssp(dssp_res)
            
            # Get disorder predictions
            incProgress(0.5, detail = "🔀 Calculating disorder predictions...")
            log_resource_usage("Disorder Prediction", query)
            iupred_df <- iupredAnchor(query)
            
            # Process MSA data
            incProgress(0.6, detail = "🧮 Processing sequence alignment...")
            log_resource_usage("MSA Start", query)
            msa_res <- NULL
            query_for_shannon <- query
            
            # Handle MSA generation with selected species
            # Use default MSA calculation with user-selected species
            # Get sequence from UniProt data
            seq <- Biostrings::AAStringSet(uniprot_data$Sequence, start=NA, end=NA, width=NA, use.names=TRUE)
            
            # Use selected species instead of default species list
            selected_species_list <- input$selected_species
            incProgress(0.7, detail = paste("🔍 Using", length(selected_species_list), "selected species..."))
            
            aa_files <- list.files(cds_folder, pattern = paste0(selected_species_list, ".*\\.all.fa$", collapse = "|"), 
                             ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
            names(aa_files) <- aa_files
              
            # BLAST search
            blast_results <- lapply(aa_files, function(.x) { protein_blast(seq, .x) })
            find_best_match <- function(.x) { head(.x[order(.x$E),], 1) }
            blast_best_match <- lapply(blast_results, find_best_match)
            blast_seqs <- lapply(blast_best_match, fetch_sequences)
            blast_seqs[[query]] <- seq
            blast_stringset  <- Biostrings::AAStringSet(unlist(lapply(blast_seqs, function(.x){.x[[1]]})))
            
            # Multiple sequence alignment
            msa_res <- muscle(blast_stringset)
            msa_res <- Biostrings::AAStringSet(msa_res)
            log_resource_usage("MSA Complete", query)
            
            # Calculate Shannon entropy
            incProgress(0.8, detail = "📊 Calculating Shannon entropy...")
            log_resource_usage("Shannon Entropy", query)
            shannon_df <- shannon_reshape(msa_res, query)
            
            # Join features in dataframe
            features_df <- Reduce(function(x, y) merge(x, y, all=TRUE), 
                                 list(shannon_df, dssp_df, iupred_df), 
                                 accumulate=FALSE)
            
            # Calculate scores
            incProgress(0.9, detail = "🧮 Calculating final scores...")
            log_resource_usage("Final Scoring", query)
            norm_feats_df <- calculate_scores(features_df)
            
            # Merge results
            res_df <- merge(norm_feats_df, features_df)
            
            # Store results in reactive data
            reactive_data$res_df <- res_df
            log_resource_usage("Analysis Complete", query)
            
            # Update UI elements
            incProgress(0.95, detail = "🎨 Updating visualizations...")
            
            # Store structure file for interactive highlighting
            reactive_data$structure_file <- alphafold_file
            
            # Prepare MSA data for custom rendering
            if(!is.null(msa_res)) {
                # Convert MSA to format for client-side rendering
                msa_sequences <- lapply(seq_along(msa_res), function(i) {
                  list(
                    name = as.character(names(msa_res)[i] %||% paste0("Seq", i)),
                    sequence = as.character(msa_res[[i]])
                  )
                })
                reactive_data$msa_data <- msa_sequences
                # Send MSA data to client
                session$sendCustomMessage("update-msa", list(
                  sequences = msa_sequences
                ))
              } else {
                # If MSA is missing, show error in MSA viewer
                runjs("document.getElementById('msa-viewer').innerHTML = '<div style=\'color:red;\'>Error: MSA data missing or malformed.</div>';");
                runjs("document.getElementById('msa-loading').style.display = 'none';");
              }
              
              # Mark analysis as complete
              reactive_data$analysis_complete <- TRUE
              reactive_data$analysis_running <- FALSE
              
              # Success message with domain information
              pfam_count <- if (!is.null(reactive_data$pfam_domains)) nrow(reactive_data$pfam_domains) else 0
              cath_count <- if (!is.null(reactive_data$cath_domains)) length(unique(reactive_data$cath_domains$domain_id)) else 0
              
              domain_info <- if (pfam_count > 0 || cath_count > 0) {
                paste("| Domains:", pfam_count, "Pfam,", cath_count, "CATH")
              } else {
                "| No domains found"
              }
              
              output$status <- renderText(paste("✅ Analysis completed successfully for", query, 
                                               "| Sequence length:", nchar(reactive_data$protein_sequence), "aa", domain_info))
              incProgress(1.0, detail = "🎉 Complete!")
              
          }, error = function(e) {
              # Error handling
              output$status <- renderText(paste("❌ Error during analysis:", e$message))
              reactive_data$analysis_complete <- FALSE
              reactive_data$analysis_running <- FALSE
          }, finally = {
              # Hide progress overlay
              runjs("$('#progress-overlay').hide();")
          })
      })
  })
  
  
  # Clear all inputs and outputs
  observeEvent(input$clear_all, {
    # Reset all reactive values
    reactive_data$uniprot_id <- NULL
    reactive_data$res_df <- NULL
    reactive_data$custom_msa <- NULL
    reactive_data$protein_sequence <- NULL
    reactive_data$uniprot_data <- NULL
    reactive_data$analysis_complete <- FALSE
    reactive_data$analysis_running <- FALSE
    reactive_data$msa_data <- NULL
    reactive_data$structure_file <- NULL
    reactive_data$pfam_domains <- NULL
    reactive_data$cath_domains <- NULL
    
    # Reset UI inputs
    updateTextInput(session, "uniprot_id", value = "")
    updateCheckboxGroupInput(session, "selected_species", 
                           selected = c("bos_taurus", "canis_lupus_familiaris", "gallus_gallus", 
                                        "homo_sapiens", "mus_musculus", "takifugu_rubripes", "xenopus_tropicalis"))
    
    # Clear visualizations
    output$score_plot <- renderPlotly({
      plot_ly(
          x = c(0, 100), 
          y = c(0, 1), 
          type = "scatter", 
          mode = "lines",
          line = list(color = "transparent"),
          showlegend = FALSE,
          source = "score_plot"
      ) %>%
      layout(
          title = "No data available. Please run analysis first.",
          xaxis = list(title = "Amino Acid Position", showgrid = FALSE),
          yaxis = list(title = "Score", showgrid = FALSE)
      )
    })
    output$results_table <- DT::renderDataTable(NULL)
    output$download_ui <- renderUI(NULL)
    output$export_buttons <- renderUI(
      div(class = "alert alert-info", "🔄 Please run analysis first to enable export options.")
    )
    output$status <- renderText("🧹 All inputs and outputs cleared.")
    
    # Clear client-side visualizations
    session$sendCustomMessage("clear-highlights", list())
    session$sendCustomMessage("clear-structure-highlight", list())
    runjs("
        document.getElementById('msa-viewer').innerHTML = '';
        document.getElementById('msa-loading').style.display = 'block';
        // Clear 3D structure viewer
        if (structureViewer) {
            if (structureViewer.useFallback && structureViewer.viewer3d) {
                structureViewer.viewer3d.clear();
            } else if (structureViewer.plugin) {
                structureViewer.plugin.clear();
            }
            structureViewer.currentHighlight = null;
            structureViewer.structureData = null;
        }
    ")
    
    # Hide progress overlay if visible
    runjs("$('#progress-overlay').hide();")
  })
  
  # Enhanced download handlers
  output$download_csv <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0(reactive_data$uniprot_id, "_epitope_analysis_", timestamp, ".csv")
    },
    content = function(file) {
      if(!is.null(reactive_data$res_df)) {
        df <- reactive_data$res_df
        # Convert list columns to character
        for(col in names(df)) {
          if(is.list(df[[col]])) {
            df[[col]] <- vapply(df[[col]], function(x) {
              if(is.null(x)) "" else paste(as.character(x), collapse = ";")
            }, character(1))
          }
        }
        write.csv(df, file = file, row.names = FALSE)
      }
    }
  )
  
  output$download_json <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0(reactive_data$uniprot_id, "_epitope_analysis_", timestamp, ".json")
    },
    content = function(file) {
      if(!is.null(reactive_data$res_df)) {
        analysis_data <- list(
          metadata = list(
            uniprot_id = reactive_data$uniprot_id,
            analysis_date = Sys.time(),
            sequence_length = nchar(reactive_data$protein_sequence),
            protein_info = reactive_data$uniprot_data
          ),
          results = reactive_data$res_df
        )
        jsonlite::write_json(analysis_data, file, pretty = TRUE)
      }
    }
  )
  
  output$download_plot <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0(reactive_data$uniprot_id, "_score_plot_", timestamp, ".png")
    },
    content = function(file) {
      if(!is.null(reactive_data$res_df) && nrow(reactive_data$res_df) > 0) {
        # Moving average function
        moving_average <- function(x, window_size) {
          n <- length(x)
          ma <- numeric(n)
          for (i in 1:min(window_size, n)) {
            mean_value <- mean(x[1:i], na.rm = TRUE)
            if (is.nan(mean_value)) mean_value <- NA
            ma[i] <- mean_value
          }
          for (i in (window_size + 1):n) {
            window_values <- x[(i - window_size + 1):i]
            mean_value <- mean(window_values, na.rm = TRUE)
            if (is.nan(mean_value)) mean_value <- NA
            ma[i] <- mean_value
          }
          return(ma)
        }
        
        res_df <- reactive_data$res_df
        
        # Handle different score column names
        score_col <- NULL
        if("min" %in% colnames(res_df)) {
            score_col <- "min"
        } else if("min_val" %in% colnames(res_df)) {
            score_col <- "min_val"
        } else if("score" %in% colnames(res_df)) {
            score_col <- "score"
        } else if("norm_score" %in% colnames(res_df)) {
            score_col <- "norm_score"
        }
        
        if(!is.null(score_col)) {
            res_df$id <- as.numeric(res_df$position)
            res_df$min_val <- as.numeric(res_df[[score_col]])
            res_df <- res_df[order(res_df$id),]
            res_df$min_val_ma <- moving_average(res_df$min_val, 7)
            
            # Create high-resolution plot for download
            p <- ggplot(res_df, aes(x = id, y = min_val_ma)) +
                geom_line(color = "black", linewidth = 1.2) +
                scale_x_continuous("Amino Acid Position") +
                scale_y_continuous("EpicTope Score") +
                theme_classic() +
                theme(
                    plot.title = element_text(size = 16, hjust = 0.5),
                    axis.title = element_text(size = 14),
                    axis.text = element_text(size = 12),
                    panel.grid = element_blank(),
                    axis.line = element_line(color = "black")
                ) +
                ggtitle(paste("EpicTope Score -", reactive_data$uniprot_id))
            
            ggsave(file, plot = p, width = 12, height = 6, dpi = 300)
        }
      }
    }
  )
  
  # Legacy download handler for compatibility
  output$download_score <- downloadHandler(
    filename = function() {
      paste0(reactive_data$uniprot_id, "_score.csv")
    },
    content = function(file) {
      if(!is.null(reactive_data$res_df)) {
        df <- reactive_data$res_df
        # Convert list columns to character
        for(col in names(df)) {
          if(is.list(df[[col]])) {
            df[[col]] <- vapply(df[[col]], function(x) {
              if(is.null(x)) "" else paste(as.character(x), collapse = ";")
            }, character(1))
          }
        }
        write.csv(df, file = file, row.names = FALSE)
      }
    }
  )
}

shinyApp(ui = ui, server = server)

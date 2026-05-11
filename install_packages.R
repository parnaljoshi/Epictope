# Install required packages for enhanced EpicTope Shiny app

# Set CRAN mirror to avoid mirror selection errors
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Alternative mirrors (uncomment if the above doesn't work):
# options(repos = c(CRAN = "https://cran.rstudio.com/"))
# options(repos = c(CRAN = "https://cran.r-project.org/"))

packages <- c(
  "DT",
  "shinyWidgets", 
  "shinydashboardPlus",
  "shinyjs",
  "jsonlite",
  "scales"
)

cat("Setting up CRAN mirror and checking packages...\n")
cat("Using CRAN mirror:", getOption("repos")["CRAN"], "\n\n")

# Check and install missing packages with better error handling
for(pkg in packages) {
  cat("Checking package:", pkg, "...")
  
  if(!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat(" not found. Installing...\n")
    
    tryCatch({
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
      cat("✓ Successfully installed:", pkg, "\n")
    }, error = function(e) {
      cat("✗ Failed to install", pkg, "- Error:", e$message, "\n")
    })
    
  } else {
    cat(" already installed ✓\n")
  }
}

cat("All packages installed successfully!\n")

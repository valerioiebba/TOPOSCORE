# List all required packages
packages <- c(
  "shiny",
  "DT",
  "ggplot2",
  "dplyr",
  "here",
  "readr",
  "ggnewscale",
  "ggrepel"
)

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))
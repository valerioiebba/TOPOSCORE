library(dplyr)
library(readr)

calculate_toposcore <- function(input_file) {
  # Read the data
  data <- read_csv(input_file)
  
  # Read signature species
  sig1_species <- read.table("data/sig1.txt", skip = 1, stringsAsFactors = FALSE)$V1
  sig2_species <- read.table("data/sig2.txt", skip = 1, stringsAsFactors = FALSE)$V1
  
  # Check which species are present in the data
  available_sig1 <- sig1_species[sig1_species %in% colnames(data)]
  available_sig2 <- sig2_species[sig2_species %in% colnames(data)]
  
  # Warn about missing species
  missing_sig1 <- setdiff(sig1_species, colnames(data))
  missing_sig2 <- setdiff(sig2_species, colnames(data))
  
  if(length(missing_sig1) > 0 || length(missing_sig2) > 0) {
    warning(sprintf("Missing species:\nSIG1: %s\nSIG2: %s",
                   paste(missing_sig1, collapse=", "),
                   paste(missing_sig2, collapse=", ")))
  }
  
  # Function to check if species is present (abundance > 0)
  is_present <- function(x) ifelse(x > 0, 1, 0)
  
  # Calculate scores for each sample using only available species
  results <- data %>%
    rowwise() %>%
    mutate(
      # Count present species using only available species
      SIG1_count = sum(across(all_of(available_sig1), is_present)),
      SIG2_count = sum(across(all_of(available_sig2), is_present)),
      
      # Calculate S score using total counts of available species
      S_score = round((SIG2_count/length(available_sig2) - 
                      SIG1_count/length(available_sig1) + 1)/2, 3),
      
      # Initial classification based on S score
      SIG_class = case_when(
        S_score <= 0.5351 ~ "SIG1",
        S_score >= 0.7911 ~ "SIG2",
        TRUE ~ "Gray"
      ),
      
      # Akkermansia status (if available)
      Akk_status = if("Akkermansia_muciniphila" %in% colnames(data)) {
        case_when(
          Akkermansia_muciniphila >= 4.799 ~ "High",
          Akkermansia_muciniphila == 0 ~ "Zero",
          TRUE ~ "Normal"
        )
      } else {
        "Unknown"
      },
      
      # Updated Toposcore classification
      Toposcore = case_when(
        SIG_class == "SIG1" ~ "SIG1+",
        SIG_class == "SIG2" ~ "SIG2+",
        SIG_class == "Gray" & (Akk_status %in% c("High", "Zero")) ~ "SIG1+",
        SIG_class == "Gray" & (Akk_status == "Normal") ~ "SIG2+",
        TRUE ~ "SIG2+"  # Default to SIG2+ if Akkermansia status is unknown
      )
    ) %>%
    ungroup() %>%
    select(
      Sample_id,
      SIG1_count,
      SIG2_count,
      S_score,
      SIG_class,
      Toposcore,
      # Include OS12 if it exists in the input data
      any_of("OS12")
    )
  
  return(results)
}
library(dplyr)
library(readr)

calculate_toposcore <- function(input_file, taxonomy_database = c("species_Jan21", "SGB_Jan21", "GTDB_r207", "GTDB_r220")) {
  
  # Read the data
  data <- read_tsv(input_file)
  sampleIDs <- data$Sample_id
  data <- data %>% select(-Sample_id)
  
  # Read signature species
  sig1_species <- read.table("data/sig1.txt", stringsAsFactors = FALSE, header = T, sep = "\t")
  sig2_species <- read.table("data/sig2.txt", stringsAsFactors = FALSE, header = T, sep = "\t")
  
  # Check which species are present in the data
  sig1_species$PresentAbsent <- apply(sig1_species, 1, function(x, column = taxonomy_database) {
    species <- unlist(strsplit(as.character(x[column]), ","))
    matched <- intersect(species, colnames(data))
    if(length(matched) != 0) {
      return(1)
    }
    else {
      return(0)
    }
  })
  sig2_species$PresentAbsent <- apply(sig2_species, 1, function(x, column = taxonomy_database) {
    species <- unlist(strsplit(as.character(x[column]), ","))
    matched <- intersect(species, colnames(data))
    if(length(matched) != 0) {
      return(1)
    }
    else {
      return(0)
    }
  })
  
  if(length(which(sig1_species$PresentAbsent == 0)) > 0 || length(which(sig2_species$PresentAbsent == 0)) > 0) {
    warning(sprintf("Missing species:\nSIG1: %s\nSIG2: %s",
                    paste(sig1_species[,taxonomy_database, drop = T][sig1_species$PresentAbsent == 0], collapse=", "),
                    paste(sig2_species[,taxonomy_database, drop = T][sig2_species$PresentAbsent == 0], collapse=", ")))
  }
  
  sig1_species <- sig1_species %>% filter(PresentAbsent == 1)
  sig2_species <- sig2_species %>% filter(PresentAbsent == 1)
  
  count_groups <- function(row, group_df, column) {
    gt_zero <- row > 0
    names(gt_zero) <- names(row)
    
    sum(sapply(group_df[,column, drop = T], function(cols) {
      col_list <- unlist(strsplit(cols, ","))
      col_list <- col_list[col_list %in% names(gt_zero)]
      any(gt_zero[col_list])
    }))
  }
  
  data$SIG1_count <- apply(data[, !sapply(data, is.character)], 1, count_groups, group_df = sig1_species, column = taxonomy_database)
  data$SIG2_count <- apply(data[, !sapply(data, is.character)], 1, count_groups, group_df = sig2_species, column = taxonomy_database)
  
  data$S_score <- round((data$SIG2_count/45 - 
                           data$SIG1_count/37 + 1)/2, 3)
  
  data$SIG_class <- apply(data[, !sapply(data, is.character)], 1, function(x) { 
    if (x['S_score'] <= 0.5351) {
      return('SIG1')
    }
    if (x['S_score'] >= 0.7911) {
      return('SIG2')
    }
    else {
      return("Gray")
    }
  })
  
  if (taxonomy_database != "SGB_Jan21") {
    results <- data %>% rowwise() %>%
      mutate(
        Akk_status = case_when(
          sum(c_across(all_of(grep("muciniphila$", colnames(data), perl = T))), na.rm = TRUE) >= 4.799 ~ "High",
          sum(c_across(all_of(grep("muciniphila$", colnames(data), perl = T))), na.rm = TRUE) == 0 ~ "Zero",
          sum(c_across(all_of(grep("muciniphila$", colnames(data), perl = T))), na.rm = TRUE) > 0 & sum(c_across(all_of(grep("muciniphila$", colnames(data), perl = T))), na.rm = TRUE) < 4.799 ~ "Normal"),
        # Updated Toposcore classification
        Toposcore = case_when(
          SIG_class == "SIG1" ~ "SIG1+",
          SIG_class == "SIG2" ~ "SIG2+",
          SIG_class == "Gray" & (Akk_status %in% c("High", "Zero")) ~ "SIG1+",
          SIG_class == "Gray" & (Akk_status == "Normal") ~ "SIG2+",
          TRUE ~ "SIG2+"))
    results$Sample_id <- sampleIDs
    results <- results %>% select(
      Sample_id, 
      Akk_status,
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
  if (taxonomy_database == "SGB_Jan21") {
    results <- data %>% rowwise() %>%
      mutate(
        Akk_status = case_when(
          sum(c_across(all_of(grep("SGB9226", colnames(data), perl = T))), na.rm = TRUE) >= 4.799 ~ "High",
          sum(c_across(all_of(grep("SGB9226", colnames(data), perl = T))), na.rm = TRUE) == 0 ~ "Zero",
          sum(c_across(all_of(grep("SGB9226", colnames(data), perl = T))), na.rm = TRUE) > 0 & sum(c_across(all_of(grep("SGB9226", colnames(data), perl = T))), na.rm = TRUE) < 4.799 ~ "Normal"),
        # Updated Toposcore classification
        Toposcore = case_when(
          SIG_class == "SIG1" ~ "SIG1+",
          SIG_class == "SIG2" ~ "SIG2+",
          SIG_class == "Gray" & (Akk_status %in% c("High", "Zero")) ~ "SIG1+",
          SIG_class == "Gray" & (Akk_status == "Normal") ~ "SIG2+",
          TRUE ~ "SIG2+"))
    results$Sample_id <- sampleIDs
    results <- results %>% select(
      Sample_id, 
      Akk_status,
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
}

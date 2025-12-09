library(dplyr)
library(readr)

calculate_toposcore <- function(input_file, taxonomy_database = c("species_Jan21", "SGB_Jan21", "GTDB_r207", "GTDB_r220")) {
  
  # Read-in the data for MetaPhlAn format
  if(taxonomy_database == "species_Jan21" | taxonomy_database == "SGB_Jan21") {
    data <- read_tsv(input_file, name_repair = "minimal", comment = "#")
    if(length(grep("\\|", data$clade_name)) == 0) {
      stop("It appears as the taxonomy separator in your table is not in standard MetaPhlAn format, i.e k__|p__")
    }
    else {
      taxa <- data$clade_name
      data <- data %>% select(-clade_name)
      sampleIDs <- colnames(data)
      data <- as_tibble(t(data))
      colnames(data) <- taxa
    }
  }
  
  # Read-in the data for GTDB format
  if(taxonomy_database == "GTDB_r207" | taxonomy_database == "GTDB_r220") {
    data <- read_tsv(input_file, name_repair = "minimal", comment = "#")
    if(length(grep(";", data$clade_name)) == 0) {
      stop("It appears as the taxonomy separator in your table is not in standard GTDB format, i.e d__;p__")
    }
    else {
      SGB9226 <- FALSE
      taxa <- data$clade_name
      data <- data %>% select(-clade_name)
      sampleIDs <- colnames(data)
      data <- as_tibble(t(data))
      colnames(data) <- taxa
      data <- data[,colnames(data)[which(unlist(lapply(colnames(data), function(x) {length(unlist(strsplit(x, ";"))) })) >= 7)]]
      colnames(data) <- lapply(colnames(data), function(x) {unlist(strsplit(x, ";"))[7] } )
      #Sum and remove duplicate columns 
      col_names <- names(data)
      dup_names <- unique(col_names[duplicated(col_names)])
      data_subset <- data[, !col_names %in% dup_names]
      # Add summed columns for each duplicated name
      for (dup_name in dup_names) {
        cols <- which(col_names == dup_name)
        data_subset[[dup_name]] <- rowSums(data[, cols])
      }
      data <- data_subset
    }
  }
  
  if (taxonomy_database == "species_Jan21") {
    SGB9226 <- ifelse(length(grep("SGB9226$", colnames(data))) > 0, TRUE, FALSE)
    if (SGB9226 == TRUE) {
      akk_sgb9226_ab <- data %>% select(grep("SGB9226$", colnames(data)))
    }
    data <- data[,colnames(data)[which(unlist(lapply(colnames(data), function(x) {length(unlist(strsplit(x, "\\|"))) })) >= 7)]]
    colnames(data) <- unlist(lapply(colnames(data), function(x) {unlist(strsplit(x, "\\|"))[7] } ))
    colnames(data) <- gsub("s__", "", colnames(data))
    #Sum and remove duplicate columns 
    col_names <- names(data)
    dup_names <- unique(col_names[duplicated(col_names)])
    data_subset <- data[, !col_names %in% dup_names]
    # Add summed columns for each duplicated name
    for (dup_name in dup_names) {
      cols <- which(col_names == dup_name)
      data_subset[[dup_name]] <- rowSums(data[, cols])
    }
    if (SGB9226 == TRUE) {
      data <- bind_cols(data_subset, akk_sgb9226_ab)
    }
    if (SGB9226 == FALSE) {
      data <- data_subset
    }
  }
  
  if (taxonomy_database == "SGB_Jan21") {
    SGB9226 <- ifelse(length(grep("SGB9226$", colnames(data))) > 0, TRUE, FALSE)
    if(SGB9226 == FALSE) {
      akk <- data %>% select(grep("muciniphila$", colnames(data)))
    }
    data <- data[,colnames(data)[which(unlist(lapply(colnames(data), function(x) {length(unlist(strsplit(x, "\\|"))) })) == 8)]]
    colnames(data) <- unlist(lapply(colnames(data), function(x) {unlist(strsplit(x, "\\|"))[8] } ))
    colnames(data) <- gsub("t__", "", colnames(data))
    col_names <- names(data)
    dup_names <- unique(col_names[duplicated(col_names)])
    if (length(dup_names) > 0) {
      stop(paste0("There are duplicate SGB ids in your table;\n", paste0(dup_names, collapse = ';'),"\n Please check the merging of the profiles prior to uploading the table."))
    }
    if(SGB9226 == FALSE) {
      data <- bind_cols(data, akk)
    }
  }
  
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
  
  if ((taxonomy_database == "species_Jan21" & SGB9226 == FALSE) | taxonomy_database %in% c("GTDB_r207", "GTDB_r220") | (taxonomy_database == "SGB_Jan21" & SGB9226 == FALSE)) {
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
  if ((taxonomy_database == "SGB_Jan21" & SGB9226 == TRUE) | (taxonomy_database == "species_Jan21" & SGB9226 == TRUE)) {
    results <- data %>% rowwise() %>%
      mutate(
        Akk_status = case_when(
          sum(c_across(all_of(grep("SGB9226$", colnames(data), perl = T))), na.rm = TRUE) >= 4.799 ~ "High",
          sum(c_across(all_of(grep("SGB9226$", colnames(data), perl = T))), na.rm = TRUE) == 0 ~ "Zero",
          sum(c_across(all_of(grep("SGB9226$", colnames(data), perl = T))), na.rm = TRUE) > 0 & sum(c_across(all_of(grep("SGB9226$", colnames(data), perl = T))), na.rm = TRUE) < 4.799 ~ "Normal"),
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

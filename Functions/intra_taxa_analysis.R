
intra_taxa_analysis = function(intra_comparisons_list, infos_data, similarity_threshold, 
                               name_col_similarity, name_col_accession, name_col_taxa, verbose = T){
  
  ## Initialisation
  
  # To avoid that the function summarise print any messages
  options(dplyr.summarise.inform = F)
  
  # Creating the columns names with suffixes
  name_col_accession1 = paste0(name_col_accession, "1")
  name_col_accession2 = paste0(name_col_accession, "2")
  name_col_taxa1 = paste0(name_col_taxa, "1")
  name_col_taxa2 = paste0(name_col_taxa, "2")
  
  # Verifying that all necessary columns are in infos_data
  if(length(which(colnames(infos_data) %in% c(name_col_accession, name_col_taxa))) != 2)
    stop(paste0('The table provided in infos_data does not contain columns named ', 
                name_col_accession, ' and ', name_col_taxa, '.'))
  
  # Verifying that all necessary columns are in intra_comparisons_list
  if(any(sapply(intra_comparisons_list, function(x) 
    length(which(colnames(x) %in% c(name_col_taxa1, name_col_taxa2, name_col_accession1, name_col_accession2, name_col_similarity))) != 5)))
    stop(paste0('The list provided in intra_comparisons_list has tables that do not contain columns named ', 
                name_col_taxa1, ', ', name_col_taxa2, ', ', name_col_accession1, ', ', name_col_accession2, ' and ', name_col_similarity, '.'))
  
  
  ## Computing summary statistics about the similarity within each taxa
  
  # Getting the number of sequences per taxa
  number_sequences_compared_per_taxa = setNames(dplyr::summarise(
    .data = dplyr::group_by(.data = infos_data, .data[[name_col_taxa]]), 
    dplyr::n()), c(name_col_taxa, "N"))
  
  # Opposing each sequences of a taxa in a matrix and taking the length of the diagonal to get only single comparisons
  number_comparisons_per_taxa = sapply(split(number_sequences_compared_per_taxa, number_sequences_compared_per_taxa[[name_col_taxa]]), 
                                       function(x) length(matrix(ncol = x$N, nrow = x$N)[lower.tri(matrix(ncol = x$N, nrow = x$N), diag = F)]))
  
  # Converting the reference vector as a dataframe
  reference_taxa = setNames(tibble(names(number_comparisons_per_taxa), number_sequences_compared_per_taxa$N, 
                                   c(number_comparisons_per_taxa)), c(name_col_taxa, "NUMBER_ACCESSION", "TOTAL_NUMBER_COMPARISONS"))
  
  # Printing which operation the function is doing
  if(verbose) cat(paste0("Computing statistics about the similarity per ", name_col_taxa, ": "))
  
  # Computing the summary of the similarity of the sequences within each taxa, for each primers
  intra_taxa_similarity = dplyr::bind_rows(lapply(intra_comparisons_list, function(y) 
    dplyr::bind_rows(lapply(split(y, y[[name_col_taxa1]]), function(x) 
      tibble(MIN_SIMILARITY = summary(x[[name_col_similarity]])[1], Q1_SIMILARITY = summary(x[[name_col_similarity]])[2],
             MEAN_SIMILARITY = summary(x[[name_col_similarity]])[4], Q3_SIMILARITY = summary(x[[name_col_similarity]])[5],
             MAX_SIMILARITY = summary(x[[name_col_similarity]])[6],
             SD_SIMILARITY = sd(x[[name_col_similarity]]))), .id = toupper(name_col_taxa))), .id = "PRIMERS")
  
  # Adding the reference informations
  intra_taxa_similarity = tibble(merge(intra_taxa_similarity, reference_taxa))
  
  # Ordering the final table
  intra_taxa_similarity = subset(intra_taxa_similarity, select = c("PRIMERS", toupper(name_col_taxa), "NUMBER_ACCESSION", "TOTAL_NUMBER_COMPARISONS", 
                                                                   "MIN_SIMILARITY", "Q1_SIMILARITY", "MEAN_SIMILARITY", "Q3_SIMILARITY", 
                                                                   "MAX_SIMILARITY", "SD_SIMILARITY"))
  
  # Printing that this operation is DONE
  if(verbose) cat("DONE\n")
  
  
  ## Computing the intra-taxa analysis for each thresholds
  
  # Creating an empty list to store the results 
  intra_taxa_details = list()
  intra_taxa_summary_list = list()
  
  # Performing the analysis for each similarity threshold provided
  for(similarity in similarity_threshold){
    
    # Printing for which similarity threshold the computation is currently running
    if(verbose) cat(paste0("Computing intra-", name_col_taxa, " analysis for ", similarity, "% similarity: "))
    
    ## Computing the detailed analysis
    
    # Computing multiple statistics for each taxa depending on the similarity, and for each primer
    intra_taxa_details[[100 - similarity]] = dplyr::bind_rows(lapply(intra_comparisons_list, function(y) 
      dplyr::bind_rows(lapply(split(y, y[[name_col_taxa1]]), function(x) 
        
        # Number of cases divided by 2 because vsearch always performs 2 comparisons (i.e. Seq1 - Seq2 & Seq2 - Seq1)
        tibble(NUMBER_ACCESSION = length(unique(c(x[[name_col_accession1]], x[[name_col_accession2]]))),
               NUMBER_INTRA_ERRORS = nrow(x[which(x[[name_col_similarity]] < similarity),])/2, 
               TOTAL_NUMBER_COMPARISONS = nrow(x)/2, 
               PERCENT_INTRA_ERRORS = nrow(x[which(x[[name_col_similarity]] < similarity),]) / nrow(x) * 100,
               MEAN_SIMILARITY_INTRA_ERRORS = mean(x[which(x[[name_col_similarity]] < similarity),][[name_col_similarity]]),
               SD_SIMILARITY_INTRA_ERRORS = sd(x[which(x[[name_col_similarity]] < similarity),][[name_col_similarity]]))
        
      ), .id = toupper(name_col_taxa))), .id = "PRIMERS")
    
    # Naming the first layer of the list with the similarity used
    names(intra_taxa_details)[100 - similarity] = similarity
    
    
    ## Computing the summary analysis
    
    # Splitting the tables obtained above per primer
    intra_taxa_similarity_list = split(intra_taxa_similarity, intra_taxa_similarity$PRIMERS)
    intra_taxa_details_list = split(intra_taxa_details[[100 - similarity]], intra_taxa_details[[100 - similarity]]$PRIMERS)
    
    # Number of wrongly discrimined sequences within a taxa per primer
    number_wrong_discrimined_cases = sapply(intra_taxa_details_list, function(x) sum(x$NUMBER_INTRA_ERRORS))
    
    # Percentage of wrongly discrimined sequences within a taxa over the total number of sequences sharing this taxa per primer
    percent_wrong_discrimined_cases = sapply(intra_taxa_details_list, function(x) 
      sum(x$NUMBER_INTRA_ERRORS) / sum(x$TOTAL_NUMBER_COMPARISONS) * 100)
    
    # Number of taxa having two or more sequences wrongly discrimined per primer
    number_taxa_with_wrong_discrimined_cases = sapply(intra_taxa_similarity_list, function(x) 
      nrow(subset(x, as.numeric(NUMBER_ACCESSION) > 1 & as.numeric(MEAN_SIMILARITY) < similarity)))
    
    # Percentage of taxa having two or more sequences wrongly discrimined over the total number of taxa with 2 sequences or more per primer
    percent_taxa_with_wrong_discrimined_cases = sapply(intra_taxa_similarity_list, function(x) 
      nrow(subset(x, as.numeric(NUMBER_ACCESSION) > 1 & as.numeric(MEAN_SIMILARITY) < similarity)) /
        nrow(subset(x, NUMBER_ACCESSION > 1)) * 100)
    
    # Mean number of wrongly discrimined comparisons per taxa per primer
    mean_number_wrong_discrimined_cases_per_taxa = sapply(intra_taxa_details_list, function(x) 
      mean(x$NUMBER_INTRA_ERRORS, na.rm = T))
    
    # Mean percentage of wrongly discrimined comparisons per taxa over the total number of comparisons within each taxa per primer
    percent_number_wrong_discrimined_cases_per_taxa = sapply(intra_taxa_details_list, function(x) 
      mean(x$PERCENT_INTRA_ERRORS, na.rm = T))
    
    # Mean similarity of wrongly discrimined comparisons per taxa per primer
    mean_similarity_wrong_discrimined_cases_per_taxa = sapply(intra_taxa_details_list, function(x) 
      mean(x$MEAN_SIMILARITY_INTRA_ERRORS, na.rm = T))
    
    # Standard deviation around the mean similarity of wrongly discrimined comparisons per taxa per primer
    sd_similarity_wrong_discrimined_cases_per_taxa = sapply(intra_taxa_details_list, function(x) 
      mean(x$SD_SIMILARITY_INTRA_ERRORS, na.rm = T))
    
    # Overall mean similarity in all comparisons per primer
    overall_mean_similarity = sapply(intra_taxa_similarity_list, function(x) mean(as.numeric(x$MEAN_SIMILARITY), na.rm = T))
    
    # Standard deviation around the mean overall similarity in all comparisons per primer
    overall_sd_similarity = sapply(intra_taxa_similarity_list, function(x) sd(as.numeric(x$MEAN_SIMILARITY)))
    
    # Assembling all these results in a single table
    intra_taxa_summary = tibble::tibble(names(intra_comparisons_list),
                                        number_wrong_discrimined_cases,
                                        round(percent_wrong_discrimined_cases, 3),
                                        number_taxa_with_wrong_discrimined_cases,
                                        round(percent_taxa_with_wrong_discrimined_cases, 3),
                                        round(mean_number_wrong_discrimined_cases_per_taxa, 3),
                                        round(percent_number_wrong_discrimined_cases_per_taxa, 3),
                                        round(mean_similarity_wrong_discrimined_cases_per_taxa, 3),
                                        round(sd_similarity_wrong_discrimined_cases_per_taxa, 3),
                                        round(overall_mean_similarity, 3),
                                        round(overall_sd_similarity, 3))
    
    # Setting the column names
    colnames(intra_taxa_summary) = c("PRIMERS", "NUMBER_INTRA_ERRORS", "PERCENT_INTRA_ERRORS", 
                                     paste0("NUMBER_", name_col_taxa,"_WITH_INTRA_ERRORS"),
                                     paste0("PERCENT_", name_col_taxa,"_WITH_INTRA_ERRORS"), 
                                     paste0("MEAN_NUMBER_INTRA_ERRORS_PER_", name_col_taxa),
                                     paste0("MEAN_PERCENT_INTRA_ERRORS_PER_", name_col_taxa), 
                                     paste0("MEAN_SIMILARITY_INTRA_ERRORS_PER_", name_col_taxa),
                                     paste0("SD_SIMILARITY_INTRA_ERRORS_PER_", name_col_taxa), 
                                     "OVERALL_MEAN_SIMILARITY", "OVERALL_SD_SIMILARITY")
    
    # Replacing any NaN by 0 (NaN created if no data to compute the summary)
    for(i in 2:ncol(intra_taxa_summary)){ intra_taxa_summary[[i]][is.nan(intra_taxa_summary[[i]])] = 0 }
    
    # To avoid losing any primer without wrongly discrimined cases for a certain threshold
    if(nrow(intra_taxa_summary) != length(names(intra_comparisons_list)))
      intra_taxa_summary = tibble::tibble(merge(data.frame(PRIMERS = names(intra_comparisons_list)), intra_taxa_summary, all = T))
    
    # Making the column names upper cases
    colnames(intra_taxa_summary) = toupper(colnames(intra_taxa_summary))
    
    # Putting the result in the list
    intra_taxa_summary_list[[100 - similarity]] = intra_taxa_summary
    
    # Naming the first layer of the list with the similarity used
    names(intra_taxa_summary_list)[100 - similarity] = similarity
    
    # Printing that analysis for this similarity threshold is DONE
    if(verbose) cat("DONE\n")
    
    # Cleaning the memory at each step of the loop
    gc()
    
  }
  
  # To remove NULL elements from the list
  intra_taxa_summary_list = Filter(Negate(is.null), intra_taxa_summary_list)
  intra_taxa_details = Filter(Negate(is.null), intra_taxa_details)
  
  # Returning the results
  return(list(SUMMARY = intra_taxa_summary_list, DETAILS_SIMILARITY = intra_taxa_similarity, 
              DETAILS_PER_TAXA = intra_taxa_details))
  
}

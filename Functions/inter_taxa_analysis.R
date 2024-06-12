
inter_taxa_analysis = function(inter_comparisons_list, infos_data, similarity_threshold, 
                               name_col_similarity, name_col_accession, name_col_taxa, verbose = T){
  
  ## Initialisation
  
  # To avoid that the function summarise print any messages
  options(dplyr.summarise.inform = F)
  
  # Creating the columns names with suffixes
  name_col_taxa1 = paste0(name_col_taxa, "1")
  name_col_taxa2 = paste0(name_col_taxa, "2")
  
  # Verifying that all necessary columns are in infos_data
    if(length(which(colnames(infos_data) %in% c(name_col_accession, name_col_taxa))) != 2)
      stop(paste0('The table provided in infos_data does not contain columns named ', 
                  name_col_accession, ' and ', name_col_taxa, '.'))
  
  # Verifying that the column containing the similarity is in inter_comparisons_list
  if(any(sapply(inter_comparisons_list, function(x) length(which(colnames(x) == name_col_similarity)) != 1)))
    stop(paste0('The list provided in inter_comparisons_list has tables that do not contain a column named ', 
                name_col_similarity, '.'))
  
  # Verifying that all other necessary columns are in inter_comparisons_list
  if(any(sapply(inter_comparisons_list, function(x) length(which(colnames(x) %in% c(name_col_taxa1, name_col_taxa2))) != 2)))
    stop(paste0('The list provided in inter_comparisons_list has tables that do not contain columns named ', 
                name_col_taxa1, ' and ', name_col_taxa2, '.'))
  
  # Warning the user if one or more similarity threshold(s) is/are inferior to the minimum similarity threshold found
  if(min(similarity_threshold) < min(sapply(inter_comparisons_list, function(x) min(x[[name_col_similarity]], na.rm = T))))
    warning('The minimum similarity threshold provided is lower than the minimum similarity in "inter_comparisons_list".')
  
  # Warning the user if one or more similarity threshold(s) is/are superior to the maximum similarity threshold found
  if(max(similarity_threshold) > max(sapply(inter_comparisons_list, function(x) max(x[[name_col_similarity]], na.rm = T))))
    warning('The maximum similarity threshold provided is higher than the maximum similarity in "inter_comparisons_list".')
  
  
  ## Making a dataframe with the number of accession per taxa
  
  if(verbose) cat("\nComputing reference matrix:")
  
  number_accession_per_taxa = sapply(split(infos_data, infos_data[[name_col_taxa]]), 
                                    function(x) length(unique(x[[name_col_accession]])))
  
  number_accession_per_taxa = setNames(tibble::tibble(data.frame(1:length(number_accession_per_taxa),
                                                                 names(number_accession_per_taxa),  
                                                                 number_accession_per_taxa)), 
                                       c(paste0(name_col_taxa, "_INDEX"), name_col_taxa, "NUMBER_ACCESSION"))
  
  
  ## Making a matrix of all combinations
  
  # The sign %o% is calls outer that compute a matrix which corresponds to the pairwise product of two vectors
  multiplied_combinations = number_accession_per_taxa$NUMBER_ACCESSION %o% number_accession_per_taxa$NUMBER_ACCESSION
  
  multiplied_combinations_all = multiplied_combinations 
  
  # To only select the lower triangle of the matrix (unique comparisons without self-comparisons)
  multiplied_combinations[upper.tri(multiplied_combinations, diag = T)] = NA
  
  if(verbose) cat(" DONE\n\n")
  
  
  ## Performing the analysis for each similarity tresholds and each primer
  
  # Creating an empty list to store the results 
  inter_taxa_analysis_summary = list()
  inter_taxa_analysis_details = list()
  
  # Loop to compute the inter-taxa analysis for each similarity threshold
  for(similarity in similarity_threshold){
    
    # Printing for which similarity threshold the computation is currently running
    if(verbose) cat(paste0("Computing inter-taxa analysis for ", similarity, "% similarity\n"))
    if(verbose) cat("---------------------------------------------------------------------\n")
    
    # Printing which operation is currently running
    if(verbose) cat("Computing detailed analysis for each metabarcode: ")
    
    # Initialisation of the vector that will indicate if any undiscrimined cases has been found
    found = 0
    
    # Loop to compute the inter-taxa analysis for each primer
    for(i in 1:length(inter_comparisons_list)){
      
      # Printing which metabarcodes is currently analysed
      if(verbose) cat(paste0(i, "/", length(inter_comparisons_list), " ... "))
      
      
      ## Making a list with the number of cases per taxa
      
      # To get only the comparisons from two different taxa that are higher than the similarity threshold 
      taxa_threshold = inter_comparisons_list[[i]][which(as.numeric(inter_comparisons_list[[i]][[name_col_similarity]]) >= similarity),]
      
      taxa_list = split(taxa_threshold, taxa_threshold[[name_col_taxa1]])
      
      undiscrimined_cases_list = list()
      
      # List with each element corresponding to all cases in which a taxa has been wrongly discrimined from other taxa
      undiscrimined_cases_list = lapply(taxa_list, function(x) sapply(split(x, x[[name_col_taxa2]]), nrow))
      
      if(length(undiscrimined_cases_list) > 0){
        
        ## Making a dataframe with the mean similarity and the standard deviation for each comparison
        
        taxa_similarity = setNames(dplyr::summarise(
          .data = dplyr::group_by(.data = taxa_threshold, .data[[name_col_taxa1]], .data[[name_col_taxa2]]), 
          mean(.data[[name_col_similarity]]), sd(.data[[name_col_similarity]])),
          c(name_col_taxa1, name_col_taxa2, "MEAN_SIMILARITY_INTER_ERRORS", "SD_SIMILARITY_INTER_ERRORS"))
        
        
        ## Making a matrix to remove duplicated cases (i.e. taxa1 - taxa2 & taxa2 - taxa1 = same combination)
        
        matrix_inter_taxa = xtabs(unlist(undiscrimined_cases_list) ~ 
                                   rep(names(undiscrimined_cases_list), times = lengths(undiscrimined_cases_list)) + 
                                   unlist(lapply(undiscrimined_cases_list, names)))
        
        # Removing the duplicated combinations
        matrix_inter_taxa[upper.tri(matrix_inter_taxa, diag = T)] = NA 
        
        inter_taxa = as.data.frame(matrix_inter_taxa)
        
        colnames(inter_taxa) = c(name_col_taxa1, name_col_taxa2, "NUMBER_INTER_ERRORS")
        
        inter_taxa = tibble::tibble(subset(inter_taxa, as.numeric(NUMBER_INTER_ERRORS) > 0))
        
        
        ## Adding taxa index for a fast subset of the reference matrix
        
        inter_taxa_index = merge(inter_taxa, number_accession_per_taxa[,1:2], by.x = name_col_taxa1, 
                                by.y = name_col_taxa, suffixes = c("1","2"))
        
        inter_taxa_index = merge(inter_taxa_index, number_accession_per_taxa[,1:2], by.x = name_col_taxa2, 
                                by.y = name_col_taxa, suffixes = c("1","2"))
        
        inter_taxa_index = dplyr::left_join(inter_taxa, inter_taxa_index, by = c(name_col_taxa1, name_col_taxa2))
        
        
        ## Searching informations for each taxa combination in the reference matrix
        
        total_cases_undiscrimined = list()
        
        for(j in 1:nrow(inter_taxa_index)){
          
          total_cases_undiscrimined[j] = multiplied_combinations_all[inter_taxa_index[j,][[paste0(name_col_taxa, "_INDEX1")]], 
                                                                     inter_taxa_index[j,][[paste0(name_col_taxa, "_INDEX2")]]]
          
        }
        
        
        ## Creating and ordering the final detailed table
        
        inter_taxa_primer = tibble::tibble(cbind(tibble::tibble(PRIMERS = names(inter_comparisons_list)[i]), inter_taxa, 
                                                tibble::tibble(TOTAL_NUMBER_COMPARISONS = unlist(total_cases_undiscrimined),
                                                PERCENT_INTER_ERRORS = inter_taxa$NUMBER_INTER_ERRORS / 
                                                  unlist(total_cases_undiscrimined) * 100)))
        
        inter_taxa_primer = tibble::tibble(merge(inter_taxa_primer, taxa_similarity))
        
        inter_taxa_primer = tibble::tibble(cbind(inter_taxa_primer[,3], inter_taxa_primer[,1:2], inter_taxa_primer[,4:8]))
        
        inter_taxa_primer = inter_taxa_primer[order(-inter_taxa_primer$TOTAL_NUMBER_COMPARISONS, 
                                                  -inter_taxa_primer$PERCENT_INTER_ERRORS),]
        
        if(nrow(inter_taxa_primer) > 0) found = found + 1
        
        else found = 0
        
      }
      
      if(found == 1) inter_taxa_final = inter_taxa_primer
      
      else if(as.numeric(found) > 1) inter_taxa_final = rbind(inter_taxa_final, inter_taxa_primer)
      
      else inter_taxa_final = setNames(tibble::tibble(data.frame(matrix(ncol = 8)))[-1,], 
                                      c("PRIMERS", name_col_taxa1, name_col_taxa2, "NUMBER_INTER_ERRORS", "TOTAL_NUMBER_COMPARISONS", 
                                        "PERCENT_INTER_ERRORS", "MEAN_SIMILARITY_INTER_ERRORS", "SD_SIMILARITY_INTER_ERRORS"))
      
    }
    
    if(verbose) cat("\n\nComputing summary for each metabarcode: ")
    
    
    ## Computing the summary table per primer
    
    # Computing the total number of unique (inferior triangle of the matrix) comparisons possible
    total_number_cases = sum(multiplied_combinations[lower.tri(multiplied_combinations, diag = F)])
    
    # Computing the total number of unique (inferior triangle of the matrix) taxa combinations possible
    total_number_taxa_pair = length(multiplied_combinations[lower.tri(multiplied_combinations, diag = F)])
    
    if(found != 0){
      
      inter_taxa_final_list = split(inter_taxa_final, inter_taxa_final$PRIMERS)
      
      inter_taxa_final_list = lapply(inter_taxa_final_list, function(x) x[!duplicated(x),])
      
      primers = names(inter_taxa_final_list)
      
      # Computing the total number of wrongly undiscrimined cases by each primer
      undiscrimined_cases = sapply(inter_taxa_final_list, function(x) sum(x$NUMBER_INTER_ERRORS))
      
      # Computing the total number of taxa pairs with wrongly undiscrimined cases for each primer
      undiscrimined_taxa_pair = sapply(inter_taxa_final_list, nrow)
      
      # To get the number of comparisons during which the taxa have been successfully discrimined for each primer
      discrimined_taxa_pair = total_number_taxa_pair - undiscrimined_taxa_pair
      
      # Computing the mean number of wrongly undiscrimined cases per taxa pairs, for each primer
      mean_undiscrimined_cases_per_taxa = list()
      for(m in 1:length(inter_taxa_final_list)) { mean_undiscrimined_cases_per_taxa[m] = mean(c(inter_taxa_final_list[[m]]$NUMBER_INTER_ERRORS, 
                                                                                              rep(0, discrimined_taxa_pair[m]))) }
      
      # Computing the mean percentage of wrongly undiscrimined cases per taxa pairs, for each primer
      mean_percent_undiscrimined_cases_per_taxa = list()
      for(n in 1:length(inter_taxa_final_list)) { mean_percent_undiscrimined_cases_per_taxa[n] = mean(c(inter_taxa_final_list[[n]]$PERCENT_INTER_ERRORS, 
                                                                                                      rep(0, discrimined_taxa_pair[n]))) }
      
    }
    
    # In case all sequences belonging to different taxa were successfully discrimined under a certain threshold, for all primers
    else{
      
      undiscrimined_cases = rep(0, length(inter_comparisons_list))
      
      undiscrimined_taxa_pair = rep(0, length(inter_comparisons_list))
      
      mean_undiscrimined_cases_per_taxa = as.list(rep(0, length(inter_comparisons_list)))
      
      mean_percent_undiscrimined_cases_per_taxa = as.list(rep(0, length(inter_comparisons_list)))
      
      primers = names(inter_comparisons_list)
      
    }
    
    # Creating the summary table
    summary_inter_taxa = tibble::tibble(data.frame(primers, undiscrimined_cases,
                                                   round(undiscrimined_cases / total_number_cases * 100, 3),
                                                   undiscrimined_taxa_pair,
                                                   round(undiscrimined_taxa_pair / total_number_taxa_pair * 100, 3),
                                                   round(unlist(mean_undiscrimined_cases_per_taxa), 3),
                                                   round(unlist(mean_percent_undiscrimined_cases_per_taxa), 3)))
    
    # Setting the column names
    colnames(summary_inter_taxa) = c("PRIMERS", "NUMBER_INTER_ERRORS", "PERCENT_INTER_ERRORS", 
                                     paste0("NUMBER_", name_col_taxa, "_PAIRS_WITH_INTER_ERRORS"),
                                     paste0("PERCENT_", name_col_taxa, "_PAIRS_WITH_INTER_ERRORS"), 
                                     paste0("MEAN_NUMBER_INTER_ERRORS_PER_", name_col_taxa, "_PAIRS"),
                                     paste0("MEAN_PERCENT_INTER_ERRORS_PER_", name_col_taxa, "_PAIRS"))
    
    # Ordering the table and removing NA
    summary_inter_taxa = tibble::tibble(merge(tibble::tibble(PRIMERS = names(inter_comparisons_list)), summary_inter_taxa, all = T))
    summary_inter_taxa = summary_inter_taxa[order(summary_inter_taxa$PRIMERS),]
    summary_inter_taxa[is.na(summary_inter_taxa)] = 0
    
    # To avoid losing any primer without wrongly undiscrimined cases for a certain threshold
    if(nrow(summary_inter_taxa) != length(names(inter_comparisons_list)))
      summary_inter_taxa = tibble::tibble(merge(data.frame(PRIMERS = names(inter_comparisons_list)), summary_inter_taxa, all = T))
    
    # Printing that the summary is DONE
    if(verbose) cat("DONE\n---------------------------------------------------------------------\n\n")
    
    
    ## Saving the results in lists
    
    # Creating the final table without duplicated cases
    inter_taxa_final = inter_taxa_final[!duplicated(inter_taxa_final),]
    
    # Making the column names upper cases
    colnames(summary_inter_taxa) = toupper(colnames(summary_inter_taxa))
    colnames(inter_taxa_final) = toupper(colnames(inter_taxa_final))
    
    # Adding the results to the list
    inter_taxa_analysis_summary[[100 - similarity]] = summary_inter_taxa
    inter_taxa_analysis_details[[100 - similarity]] = inter_taxa_final
    
    # Naming the first layer of the lists with the similarity used
    names(inter_taxa_analysis_summary)[100 - similarity] = similarity
    names(inter_taxa_analysis_details)[100 - similarity] = similarity
    
    # Cleaning the memory at each step of the loop
    gc()
    
  }
  
  # To remove NULL elements from the list
  inter_taxa_analysis_summary = Filter(Negate(is.null), inter_taxa_analysis_summary)
  inter_taxa_analysis_details = Filter(Negate(is.null), inter_taxa_analysis_details)
  
  # Returning the results
  return(list(SUMMARY = inter_taxa_analysis_summary, DETAILS_UNDISCRIMINED_PAIRS = inter_taxa_analysis_details))
  
}

multiple_resolution_analyses = function(metabarcodes_data, infos_data, elements_to_subset, similarity_threshold,  
                                        name_col_subset, name_col_taxa, name_col_primer, name_col_accession, 
                                        name_col_fragment, keep_combinations = NULL, program_name, verbose = F){
  
  ## Initialisation
  
  # Converting the vector of elements_to_subset name in a list
  if(class(elements_to_subset) %in% c("character", "numeric", "integer")) subset_list = as.list(elements_to_subset)
  
  # Creating an empty list of the size of the elements_to_subset
  summary_list = vector(mode = "list", length = length(subset_list))
  
  # Creating the columns names with suffixes
  name_col_accession1 = paste0(name_col_accession, "1")
  name_col_accession2 = paste0(name_col_accession, "2")
  name_col_taxa1 = paste0(name_col_taxa, "1")
  name_col_taxa2 = paste0(name_col_taxa, "2")
  
  # Verifying that all necessary columns are in infos_data
  if(length(which(colnames(infos_data) %in% c(name_col_accession, name_col_subset))) != 2)
    stop(paste0('The table provided in infos_data does not contain columns named ', 
                name_col_accession, ' and ', name_col_subset, '.'))
  
  # Verifying that the necessary column is in metabarcodes_data
  if(length(which(colnames(metabarcodes_data) %in% name_col_accession)) != 1)
    stop(paste0('The table provided in metabarcodes_data does not contain a column named ', 
                name_col_accession, '.'))
  
  # Verifying that all elements to subset are present in the corresponding column
  if(!all(unlist(subset_list)[-which(unlist(subset_list) %in% c("All", "Others"))] %in% infos_data[[name_col_subset]]))
    stop(paste0('Some subset given in "elements_to_subset" are not in the column "', name_col_subset, '" of infos_data.'))
  
  
  ## Performing the resolution analysis for each separate elements_to_subset
  
  for(i in 1:length(subset_list)){
    
    if(i > 1) cat("\n\n")
    
    cat(paste0("Computing the resolution analysis for ", paste0(subset_list[[i]], collapse = " + ")))
    cat("\n---------------------------------------------------------------------\n")
    start = Sys.time()
    
    # Getting the names of sequences that belong to the required order
    if(all(subset_list[[i]] == "All")) subset_infos_data = infos_data
    else if(all(subset_list[[i]] == "Others")) subset_infos_data = infos_data[which(!(infos_data[[name_col_subset]] %in% unlist(subset_list))), ]
    else subset_infos_data = infos_data[which(infos_data[[name_col_subset]] %in% subset_list[[i]]), ]
    
    # Getting the corresponding metabarcodes for each primers
    subset_metabarcodes_data = metabarcodes_data[which(metabarcodes_data[[name_col_accession]] %in% subset_infos_data[[name_col_accession]]), ]
    
    # Performing all intra-BIN pairwise comparisons for the given subset
    cat("Performing all intra-BIN pairwise comparisons: ")
    subset_intra_taxa_list = intra_taxa_similarity(metabarcodes_data = subset_metabarcodes_data, infos_data = subset_infos_data,
                                                   name_col_taxa = name_col_taxa, name_col_lower_rank_taxa = NULL, 
                                                   name_col_primer = name_col_primer, name_col_accession = name_col_accession,
                                                   name_col_fragment = name_col_fragment, keep_col_infos_data = name_col_taxa,
                                                   order_col_keeped = NULL, program_name = program_name, verbose = verbose)
    
    # Performing all inter-BIN pairwise comparisons for the given subset
    cat("DONE\nPerforming all inter-BIN pairwise comparisons: ")
    subset_inter_taxa_list = vsearch_pairwise_similarity(metabarcodes_data = subset_metabarcodes_data, 
                                                         infos_data = subset_infos_data, 
                                                         min_similarity_threshold = min(similarity_threshold), 
                                                         max_similarity_threshold = 100, number_rows_considered = c(0,0), 
                                                         name_col_primer = name_col_primer, 
                                                         name_col_accession = name_col_accession, 
                                                         name_col_fragment = name_col_fragment, 
                                                         keep_combinations = NULL, keep_col_infos_data = name_col_taxa, 
                                                         order_col_keeped = NULL, program_name = program_name, verbose = verbose)
    
    # Selecting only the comparisons between sequences from different BIN
    subset_inter_taxa_list = lapply(subset_inter_taxa_list, function(x) x[which(x[[name_col_taxa1]] != x[[name_col_taxa2]]),])
    
    # Computing the intra-BIN analysis summary for all thresholds provided
    cat("DONE\nRunning the intra-BIN analysis for all thresholds provided: ")
    subset_intra_taxa_summary = intra_taxa_analysis(subset_intra_taxa_list, subset_infos_data, similarity_threshold,
                                                    name_col_similarity = "SIMILARITY", name_col_taxa = name_col_taxa, 
                                                    name_col_accession = name_col_accession, verbose = verbose)$SUMMARY
    
    # Computing the inter-BIN analysis summary for all thresholds provided 
    cat("DONE\nRunning the inter-BIN analysis for all thresholds provided: ")
    subset_inter_taxa_summary = inter_taxa_analysis(subset_inter_taxa_list, subset_infos_data, similarity_threshold, 
                                                    name_col_similarity = "SIMILARITY", name_col_taxa = name_col_taxa, 
                                                    name_col_accession = name_col_accession, verbose = verbose)$SUMMARY
    
    # Assembling the tables from both analysis together
    cat("DONE\nAssembling the results of both analysis: ")
    for(j in 1:length(subset_intra_taxa_summary)) { 
      
      subset_summary_temporary = tibble(cbind(tibble(SIMILARITY = names(subset_intra_taxa_summary)[j]),
                                              subset_intra_taxa_summary[[names(subset_intra_taxa_summary)[j]]], 
                                              subset_inter_taxa_summary[[names(subset_intra_taxa_summary)[j]]][,-1]))
      
      subset_summary_temporary$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = 
        subset_summary_temporary$PERCENT_INTRA_ERRORS + subset_summary_temporary$PERCENT_INTER_ERRORS
      
      if(j == 1) subset_summary = subset_summary_temporary
      
      else subset_summary = rbind(subset_summary, subset_summary_temporary)
      
    }
    
    # Ordering the final summary table and adding it to the list
    summary_list[[i]] = subset(subset_summary, select = 
                                 c("SIMILARITY", "PRIMERS", 
                                   "NUMBER_INTRA_ERRORS", "NUMBER_INTER_ERRORS",
                                   "PERCENT_INTRA_ERRORS", "PERCENT_INTER_ERRORS", 
                                   "CUMULATED_PERCENT_ERRORS_PER_ANALYSIS",
                                   paste0("NUMBER_", toupper(name_col_taxa), "_WITH_INTRA_ERRORS"), 
                                   paste0("NUMBER_", toupper(name_col_taxa), "_PAIRS_WITH_INTER_ERRORS"),
                                   paste0("PERCENT_", toupper(name_col_taxa), "_WITH_INTRA_ERRORS"), 
                                   paste0("PERCENT_", toupper(name_col_taxa), "_PAIRS_WITH_INTER_ERRORS"),
                                   paste0("MEAN_NUMBER_INTRA_ERRORS_PER_", toupper(name_col_taxa)), 
                                   paste0("MEAN_NUMBER_INTER_ERRORS_PER_", toupper(name_col_taxa), "_PAIRS"),
                                   paste0("MEAN_PERCENT_INTRA_ERRORS_PER_", toupper(name_col_taxa)), 
                                   paste0("MEAN_PERCENT_INTER_ERRORS_PER_", toupper(name_col_taxa), "_PAIRS")))
    
    # Naming this element of the list with the current order name
    names(summary_list)[i] = paste0(subset_list[[i]], collapse = " + ")
    
    cat("DONE\n---------------------------------------------------------------------\n")
    end = Sys.time() ; duration = difftime(end, start)
    cat(paste("Time taken:", round(duration[[1]], 2), units(duration)))
    
  }
  
  # Binding the summary tables computed for every elements_to_subset together 
  summary_df = dplyr::bind_rows(summary_list, .id = name_col_subset)
  
  # Returning the final summary table containing the result for all elements_to_subset
  return(summary_df)
  
}

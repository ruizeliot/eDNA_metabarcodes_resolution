intra_taxa_similarity = function(metabarcodes_data, infos_data, name_col_taxa, name_col_lower_rank_taxa = NULL, 
                                 name_col_primer, name_col_accession, name_col_fragment, 
                                 keep_col_infos_data = NULL, order_col_keeped = NULL, program_name, verbose = T){
  
  # Setting a general option so that the function summarise does not print any message
  options(dplyr.summarise.inform = F)
  
  # Associating taxonomic infos for each metabarcodes
  metabarcodes_infos = tibble::tibble(merge(metabarcodes_data, infos_data))
  metabarcodes_infos$ROW_ID = 1:nrow(metabarcodes_infos)
  
  # General verification of the validity of the input file
  if(length(unique(metabarcodes_infos[[name_col_primer]])) > 1)
    if(var(sapply(split(metabarcodes_infos, metabarcodes_infos[[name_col_primer]]), nrow)) != 0) 
      stop('The table filled in "metabarcodes_data" must contain equal number of rows accross primers (i.e. same accession).')
  
  # Computing the number of comparisons for each taxa that will be performed to serve as reference
  if(is.null(name_col_lower_rank_taxa)) {
    
    # Getting the number of sequences per taxa (column N)
    number_sequences_compared = setNames(dplyr::summarise(
      .data = dplyr::group_by(.data = metabarcodes_infos[which(metabarcodes_infos[[name_col_primer]] == 
                                                                 unique(metabarcodes_infos[[name_col_primer]])[1]), ],
                              .data[[name_col_taxa]]), dplyr::n()), c(name_col_taxa, "N"))
    
    # Opposing each sequences of a taxa in a matrix and taking the length of the diagonal to get only single comparisons
    number_comparisons = sapply(split(number_sequences_compared, number_sequences_compared[[name_col_taxa]]), 
                                function(x) length(matrix(ncol = x$N, nrow = x$N)[lower.tri(matrix(ncol = x$N, nrow = x$N), diag = F)]))
    
  }
  
  # Computing the number of comparisons for each rank that will be performed to serve as reference
  else {
    
    # Getting the number of sequences per taxa (column N)
    number_sequences_compared = setNames(dplyr::summarise(
      .data = dplyr::group_by(.data = metabarcodes_infos[which(metabarcodes_infos[[name_col_primer]] == 
                                                                 unique(metabarcodes_infos[[name_col_primer]])[1]), ],
                              .data[[name_col_taxa]], .data[[name_col_lower_rank_taxa]]), dplyr::n()), 
      c(name_col_taxa, name_col_lower_rank_taxa, "N"))
    
    # Creating a matrix containing the number of comparisons between each subtaxa (product) and taking the lower triangle
    number_comparisons = sapply(split(number_sequences_compared, number_sequences_compared[[name_col_taxa]]), 
                                function(x) sum(outer(x$N, x$N, "*")[lower.tri(outer(x$N, x$N, "*"), diag = F)]))
    
  }
  
  # Converting the reference vector as a dataframe and removing the taxa without comparisons (1 sequence)
  reference = tibble::tibble(TAXA = names(number_comparisons), TOTAL_NUMBER_CASES = c(number_comparisons))
  reference = subset(reference, TOTAL_NUMBER_CASES != 0)
  
  # Removing the taxa with only one sequences from the metabarcodes_infos file too
  metabarcodes_infos = metabarcodes_infos[which(metabarcodes_infos[[name_col_taxa]] %in% reference$TAXA), ]
  
  # For each remaining taxa
  for(i in 1:nrow(reference)){
    
    # To print the name of the taxa and the number of comparisons to be done by vsearch
    if(verbose) cat(paste0("Retrieving ", reference$TOTAL_NUMBER_CASES[i], " comparisons for ", 
               reference$TAXA[i], " (taxon ", i, "/", length(reference$TAXA), ")\n\n"))
    
    # Printing the details of vsearch if the computation time is a bit long
    if(verbose && reference$TOTAL_NUMBER_CASES[i] > 100000){
      cat("\n----------------------------------------------------------\n")
      cat("Computing the list of relevant comparisons: ")
      internal_verbose = T
    }
    else internal_verbose = F
    
    # Subsetting the metabarcodes to get only those that are in the taxa for which the comparison is done
    taxa_metabarcodes_infos = metabarcodes_infos[which(metabarcodes_infos[[name_col_taxa]] == reference$TAXA[i]),]
    
    # Splitting the table in a list by subtaxa or by row for taxa
    if(is.null(name_col_lower_rank_taxa)) taxa_metabarcodes_infos_list = split(taxa_metabarcodes_infos, taxa_metabarcodes_infos$ROW_ID)
    else taxa_metabarcodes_infos_list = split(taxa_metabarcodes_infos, taxa_metabarcodes_infos[[name_col_lower_rank_taxa]])
    
    # Creating a list containing only accession numbers
    taxa_metabarcodes_infos_list_accession = lapply(taxa_metabarcodes_infos_list, function(x) x$ACCESSION)
    
    # Associating each accession number of a subtaxa with the accession number from all other subtaxa, for all subtaxa
    intra_taxa_combinations = unlist(lapply(1:length(taxa_metabarcodes_infos_list_accession), function(y) 
      unlist(lapply(taxa_metabarcodes_infos_list_accession[[y]], function(x) 
        paste0(x, ":", unlist(taxa_metabarcodes_infos_list_accession[-y]))))))
    
    # Printing the details of vsearch if the computation time is a bit long
    if(verbose && reference$TOTAL_NUMBER_CASES[i] > 100000) cat("DONE\n\n")
    
    # Setting the minimal columns to be keeped in case keep_col_infos_data is NULL
    if(is.null(keep_col_infos_data)) col_kept = c(name_col_taxa, name_col_lower_rank_taxa)
    else col_kept = keep_col_infos_data
    
    # Running vsearch to perform all possible comparisons and then filtering those of interest with "keep_combinations"
    similarity_list_temporary = vsearch_pairwise_similarity(metabarcodes_data = taxa_metabarcodes_infos, 
                                                            infos_data = infos_data, 
                                                            min_similarity_threshold = 0, 
                                                            max_similarity_threshold = 100, 
                                                            number_rows_considered = c(0,0), 
                                                            program_name = program_name,
                                                            name_col_primer = name_col_primer, 
                                                            name_col_accession = name_col_accession, 
                                                            name_col_fragment = name_col_fragment, 
                                                            keep_combinations = intra_taxa_combinations, 
                                                            keep_col_infos_data = col_kept, 
                                                            order_col_keeped = order_col_keeped, 
                                                            verbose = internal_verbose)
    
    # Making sure that there is no comparisons of two sequences from the same subtaxa
    if(!is.null(name_col_lower_rank_taxa)) similarity_list_temporary = lapply(similarity_list_temporary, function(x) 
      x[which(x[[paste0(name_col_lower_rank_taxa, 1)]] != x[[paste0(name_col_lower_rank_taxa, 2)]]),])
    
    # Printing the details of vsearch if the computation time is a bit long
    if(verbose && reference$TOTAL_NUMBER_CASES[i] > 100000)
      cat("\n----------------------------------------------------------\n")
    
    # Binding the table obtained at every iterations below the previous one and naming the list with the primers names 
    if(i == 1) similarity_list = similarity_list_temporary
    else {
      old_names_subset = names(similarity_list_temporary)
      similarity_list = lapply(names(similarity_list), function(x) 
        rbind(similarity_list[[x]], similarity_list_temporary[[x]]))
      names(similarity_list) = old_names_subset
      
    }
    
    # Cleaning the memory at each step of the loop
    gc()
    
  }
  
  # Returning the final table with similarities
  return(similarity_list)
  
}


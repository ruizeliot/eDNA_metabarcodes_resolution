# Requirement: "reshape2", "tibble" & "DECIPHER"

clusters_per_primers = function(infos_data, intra_comparisons_list, inter_comparisons_list, 
                                similarity_threshold, tree_method, name_col_similarity, 
                                name_col_accession, name_col_taxa, nb_processors = NULL, verbose = T){
  
  ## Initialisation
  
  # Creating the columns names with suffixes
  name_col_accession1 = paste0(name_col_accession, "1")
  name_col_accession2 = paste0(name_col_accession, "2")
  name_col_taxa1 = paste0(name_col_taxa, "1")
  name_col_taxa2 = paste0(name_col_taxa, "2")
  
  # Verifying that all necessary columns are in infos_data
  if(length(which(colnames(infos_data) %in% c(name_col_accession, name_col_taxa))) != 2)
    stop(paste0('The table provided in infos_data does not contain columns named ', 
                name_col_taxa, ' and ', name_col_accession, '.'))
  
  # Verifying that all necessary columns are in intra_comparisons_list
  if(any(sapply(intra_comparisons_list, function(x) 
    length(which(colnames(x) %in% c(name_col_taxa1, name_col_taxa2, name_col_accession1, name_col_accession2, name_col_similarity))) != 5)))
    stop(paste0('The list provided in intra_comparisons_list has tables that do not contain columns named ', 
                name_col_taxa1, ', ', name_col_taxa2, ', ', name_col_accession1, ', ', name_col_accession2, ' and ', name_col_similarity, '.'))
  
  # Verifying that all necessary columns are in inter_comparisons_list
  if(any(sapply(inter_comparisons_list, function(x) 
    length(which(colnames(x) %in% c(name_col_taxa1, name_col_taxa2, name_col_accession1, name_col_accession2, name_col_similarity))) != 5)))
    stop(paste0('The list provided in inter_comparisons_list has tables that do not contain columns named ', 
                name_col_taxa1, ', ', name_col_taxa2, ', ', name_col_accession1, ', ', name_col_accession2, ' and ', name_col_similarity, '.'))
  
  
  ## Computing the dataframe of all comparisons
  
  if(verbose) cat(paste0("Creating a table with all possible comparisons: "))
  
  # Creates a matrix with as many rows and columns as accession number, and with each taxa combination in each case
  all_comparisons_matrix = outer(infos_data[[name_col_taxa]], infos_data[[name_col_taxa]], paste, sep = "-") 
  colnames(all_comparisons_matrix) = infos_data[[name_col_accession]]
  rownames(all_comparisons_matrix) = infos_data[[name_col_accession]]
  
  # Converting the matrix into a dataframe
  all_comparisons_df = tibble::tibble(reshape2::melt(all_comparisons_matrix))
  
  # Splitting the taxa column, pasting the ACCESSION columns, and adding a SIMILARITY column
  all_comparisons_df = setNames(tibble::tibble(cbind(tibble::tibble(NA), all_comparisons_df[,1:2], 
                                                     tibble::tibble(TAXA1 = sapply(strsplit(as.character(all_comparisons_df$value), "-"), function(x) x[1]), 
                                                                    TAXA2 = sapply(strsplit(as.character(all_comparisons_df$value), "-"), function(x) x[2]),
                                                                    paste0(all_comparisons_df$Var1, "-", all_comparisons_df$Var2)))), 
                                c(name_col_similarity, name_col_accession1, name_col_accession2, name_col_taxa1, name_col_taxa2, "COMBINATION"))
  
  # Removing the comparison between the same sequences
  all_different_comparisons_df = all_comparisons_df[which(all_comparisons_df[[name_col_accession1]] != 
                                                            all_comparisons_df[[name_col_accession2]]),]
  
  if(verbose) cat("DONE\n\n")
  
  
  ## Computing the number of OTU perceived
  
  # Initilisation of the list that will contain the number of OTU perceived per primer
  clusters_per_similarity = list()
  
  # Looping for each primers (assuming intra_comparisons_list and inter_comparisons_list have the same length)
  for(i in 1:length(intra_comparisons_list)){ 
    
    if(verbose) cat(paste0("Creating the clusters for primer ", names(intra_comparisons_list)[i], ":\n"))
    if(verbose) cat("-------------------------------------------------\n")
    
    # Getting the intra-taxa errors of the primer for the minimum similarity_threshold of similarity (all below will be considered as 0%)
    primer_intra_taxa_errors = intra_comparisons_list[[i]][which(intra_comparisons_list[[i]][[name_col_similarity]] < min(similarity_threshold)),]
    
    # Getting the inter-taxa errors of the primer for the minimum similarity_threshold of similarity (all below will be considered as 0%)
    primer_inter_taxa_errors = inter_comparisons_list[[i]][which(inter_comparisons_list[[i]][[name_col_similarity]] >= min(similarity_threshold)),]
    
    # Binding both results and pasting the ACCESSION columns
    primer_errors = subset(rbind(primer_intra_taxa_errors, primer_inter_taxa_errors), select = 
                             c(name_col_similarity, name_col_accession1, name_col_accession2, name_col_taxa1, name_col_taxa2))
    primer_errors$COMBINATION = paste0(primer_errors[[name_col_accession1]], "-", primer_errors[[name_col_accession2]])
    
    # Replacing the similarities values that don't correspond to intra-taxa errors by 100 (i.e. not wrongly discrimined)
    all_different_comparisons_df[[name_col_similarity]] = replace(all_different_comparisons_df[[name_col_similarity]], 
                                                      which(all_different_comparisons_df[[name_col_taxa1]] == 
                                                              all_different_comparisons_df[[name_col_taxa2]]), 100)
    
    # Replacing the similarities values that don't correspond to inter-taxa errors by 0 (i.e. not wrongly gathered)
    all_different_comparisons_df[[name_col_similarity]] = replace(all_different_comparisons_df[[name_col_similarity]], 
                                                      which(all_different_comparisons_df[[name_col_taxa1]] != 
                                                              all_different_comparisons_df[[name_col_taxa2]]), 0)
    
    # Replacing the similarities values that correspond to errors by the actual value
    similarity_random_sequences_df = rbind(primer_errors, subset(all_different_comparisons_df, 
                                                                 !(COMBINATION %in% primer_errors$COMBINATION)))
    
    # Computing the distance to work with the functions below
    similarity_random_sequences_df$DISTANCE = 100 - similarity_random_sequences_df[[name_col_similarity]]
    
    # Converting into a matrix the distances between each sequences in each case
    similarity_random_sequences_matrix = as.dist(reshape2::dcast(similarity_random_sequences_df, 
                                                                 formula(paste0(name_col_accession1, " ~ ",
                                                                                name_col_accession2)), mean, 
                                                                 value.var = "DISTANCE")[,-1])
    
    # Creating the clusters per levels of similarity
    clusters = DECIPHER::IdClusters(similarity_random_sequences_matrix, method = tree_method,
                          type = "clusters", cutoff = (100 - similarity_threshold),
                          processors = nb_processors, verbose = T)
    
    clusters_per_similarity[[i]] = tibble::tibble(cbind(tibble::tibble(ACCESSION = rownames(clusters)),  
                                                setNames(clusters, paste0("CLUSTERS_", similarity_threshold))))
    
    if(verbose) cat("\n--------------------- DONE ----------------------\n\n")
    
  }
  
  names(clusters_per_similarity) = names(inter_comparisons_list)
  
  return(clusters_per_similarity)
  
}

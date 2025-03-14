# Requirement: "reshape2", "tibble", "Biostrings" & "DECIPHER"

cluster_metabarcodes_decipher = function(unaligned_seqs_per_primer_df = NULL, alignment_per_primer_df = NULL, 
                                         infos_data = NULL, intra_comparisons_list = NULL, inter_comparisons_list = NULL, 
                                         similarity_threshold, tree_methods, max_time_clustering = NULL, name_col_accession, 
                                         name_col_fragment = NULL, name_col_alignment = NULL, name_col_primer = NULL,
                                         name_col_similarity = NULL, name_col_taxa = NULL, nb_threads = NULL, verbose = T){
  
  ## Clustering with DECIPHER on the basis of unaligned fragments with the new method CLUSTERIZE
  
  if(any(tree_methods == "clusterize")){
    
    # Initialisation of the list that will contain the number of clusters per primer and per method
    clusters_per_similarity = list()
    
    # Making sure the Clusterize function exists (depends on DECIPHER versions) and that it is ran separately from other distance based methods
    if(!exists("Clusterize", where = asNamespace("DECIPHER"), inherits = FALSE)) 
      stop("The DECIPHER package needs to be updated as the Clusterize function is not available in this version.")
    if(length(unique(tree_methods)) != 1) warning("The function needs to be run separately with the clusterize method (only this method was applied here).")
    
    # Checking that required arguments are provided if a distance matrix is not directly provided
    if(is.null(unaligned_seqs_per_primer_df)) stop("The argument name_col_similarity is mandatory if not providing a distance matrix.")
    if(is.null(name_col_fragment)) stop("The argument name_col_fragment is mandatory with the clusterize method.")
    if(is.null(name_col_primer)) stop("The argument name_col_primer is mandatory with the clusterize method.")
    
    # Running the Clusterize function for each metabarcode
    primers_list = unique(unaligned_seqs_per_primer_df[[name_col_primer]])
    for(i in 1:length(primers_list)){
      
      if(verbose) cat(paste0("Creating the clusters for primer ", primers_list[i], ":\n"))
      if(verbose) cat("-------------------------------------------------\n")
      
      # Creating the DNAStringSet object with unaligned sequences for the metabarcode considered
      unaligned_seqs_per_primer_subset = unaligned_seqs_per_primer_df[which(unaligned_seqs_per_primer_df[[name_col_primer]] == primers_list[i]),]
      unaligned_seqs_per_primer_subset_dna = Biostrings::DNAStringSet(unaligned_seqs_per_primer_subset[[name_col_fragment]])
      names(unaligned_seqs_per_primer_subset_dna) = unaligned_seqs_per_primer_subset[[name_col_accession]]
      
      # Running the Clusterize function for all cluster provided in ascending order (to get independent clusters)
      clusters = DECIPHER::Clusterize(unaligned_seqs_per_primer_subset_dna, cutoff = sort((100-sort(similarity_threshold))/100))
      
      # Converting to a tibble and renaming columns
      clusters_per_similarity_df = tibble::tibble(cbind(tibble::tibble(ACCESSION = rownames(clusters)),  
                                                        setNames(clusters, paste0("CLUSTERS_", sort(similarity_threshold, decreasing = T)))))
      clusters_per_similarity_df$METHOD = "clusterize"
      
      # Storing tables per metabarcodes in the final list
      clusters_per_similarity[[i]] = clusters_per_similarity_df
      names(clusters_per_similarity)[i] = primers_list[i]
      
      gc()
      if(verbose) cat("\n--------------------- DONE ----------------------\n\n")
      
    }
    
  }
  
  else{
    
    ## Clustering with DECIPHER on the basis of similarity lists computed with VSEARCH
    
    ## Initialisation
    
    # Checking that required arguments are provided if a distance matrix is not directly provided
    if(is.null(infos_data)) stop("The table infos_data is mandatory if not providing a distance matrix.")
    if(is.null(intra_comparisons_list)) stop("The list intra_comparisons_list infos_data is mandatory if not providing a distance matrix.")
    if(is.null(inter_comparisons_list)) stop("The list inter_comparisons_list is mandatory if not providing a distance matrix.")
    if(is.null(name_col_similarity)) stop("The argument name_col_similarity is mandatory if not providing a distance matrix.")
    if(is.null(name_col_taxa)) stop("The argument name_col_taxa is mandatory if not providing a distance matrix.")
    
    # Creating the columns names with suffixes
    name_col_accession1 = paste0(name_col_accession, "1")
    name_col_accession2 = paste0(name_col_accession, "2")
    name_col_taxa1 = paste0(name_col_taxa, "1")
    name_col_taxa2 = paste0(name_col_taxa, "2")
    
    # Checking that names from both list are named and ordered the same
    if(length(setdiff(names(intra_comparisons_list), names(inter_comparisons_list))) != 0 ||
       length(setdiff(names(inter_comparisons_list), names(intra_comparisons_list))) != 0)
      stop("The elements of intra_comparisons_list and inter_comparisons_list are not named and/or ordered the same.")
    
    # Verifying that all necessary columns are in infos_data
    if(length(which(colnames(infos_data) %in% c(name_col_accession, name_col_taxa))) != 2)
      stop(paste0('The table provided in infos_data does not contain columns named ', 
                  name_col_taxa, ' or ', name_col_accession, '.'))
    
    # Verifying that all necessary columns are in intra_comparisons_list
    if(any(sapply(intra_comparisons_list, function(x) 
      length(which(colnames(x) %in% c(name_col_taxa1, name_col_taxa2, name_col_accession1, name_col_accession2, name_col_similarity))) != 5)))
      stop(paste0('The list provided in intra_comparisons_list has tables that do not contain columns named ', 
                  name_col_taxa1, ', ', name_col_taxa2, ', ', name_col_accession1, ', ', name_col_accession2, ' or ', name_col_similarity, '.'))
    
    # Verifying that all necessary columns are in inter_comparisons_list
    if(any(sapply(inter_comparisons_list, function(x) 
      length(which(colnames(x) %in% c(name_col_taxa1, name_col_taxa2, name_col_accession1, name_col_accession2, name_col_similarity))) != 5)))
      stop(paste0('The list provided in inter_comparisons_list has tables that do not contain columns named ', 
                  name_col_taxa1, ', ', name_col_taxa2, ', ', name_col_accession1, ', ', name_col_accession2, ' or ', name_col_similarity, '.'))
    
    
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
    
    gc()
    if(verbose) cat("DONE\n\n")
    
    
    ## Performing the clustering based on VSEARCH comparisons
    
    # Initilisation of the list that will contain the number of clusters per primer
    clusters_per_similarity = list()
    
    # Looping for each primers (assuming intra_comparisons_list and inter_comparisons_list have the same length)
    for(i in 1:length(intra_comparisons_list)){ 
      
      # Initialisation of the list that will contain the number of clusters per method
      clusters_per_similarity_df_list = list()
      
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
      similarity_random_sequences_matrix = reshape2::dcast(similarity_random_sequences_df, 
                                                           formula(paste0(name_col_accession1, " ~ ",
                                                                          name_col_accession2)), mean, 
                                                           value.var = "DISTANCE")
      
      # Renaming the rows with accession numbers
      rownames(similarity_random_sequences_matrix) = similarity_random_sequences_matrix[[1]]  
      similarity_random_sequences_matrix = similarity_random_sequences_matrix[,-1]
      similarity_random_sequences_matrix = as.dist(similarity_random_sequences_matrix)
      
      # Run the function for each method
      for(method in tree_methods){
        
        if(verbose) cat(paste0("Method: ", method, ":\n"))
        
        # Clustering for DNA-based methods using aligned sequences
        if(any(tree_methods %in% c("inexact"))) stop('The method "inexact" from IdCluster is not supported, consider removing it.')
        
        # Clustering for DNA-based methods using aligned sequences
        if(any(tree_methods %in% c("ME", "MP", "ML"))){
          
          # Checking methods chosen are compatible with DECIPHER package version
          if(any(tree_methods %in% c("ME", "MP")) & !exists("Treeline", where = asNamespace("DECIPHER"), inherits = FALSE)) 
            stop('The methods "ME" and "MP" are not implemented in the current version of the DECIPHER package, consider updating it to use them.')
          
          # Checking if the alignment has been provided
          if(is.null(alignment_per_primer_df)) stop('A table with alignments per primer must be provided if using "ML", "ME" and/or "MP" methods.')
          
          # Converting aligned sequences for the current primer to a DNAStringSet
          alignment_per_primer_dna_subset = alignment_per_primer_df[which(alignment_per_primer_df[[name_col_accession]] %in% 
                                                                              labels(similarity_random_sequences_matrix) & 
                                                                              alignment_per_primer_df[[name_col_primer]] %in% names(intra_comparisons_list)[i]),]
          alignment_per_primer_dna_subset_dna = Biostrings::DNAStringSet(alignment_per_primer_dna_subset[[name_col_alignment]])
          names(alignment_per_primer_dna_subset_dna) = alignment_per_primer_dna_subset[[name_col_accession]]
          
          # Filtering the distance matrix so that it contains only sequences on the DNAStringSet
          similarity_random_sequences_matrix_subset = as.matrix(similarity_random_sequences_matrix)
          similarity_random_sequences_matrix_subset = similarity_random_sequences_matrix_subset[
            rownames(similarity_random_sequences_matrix_subset) %in% names(alignment_per_primer_dna_subset_dna),
            colnames(similarity_random_sequences_matrix_subset) %in% names(alignment_per_primer_dna_subset_dna)]
          similarity_random_sequences_matrix_subset = as.dist(similarity_random_sequences_matrix_subset)
          
          # Verifying that their are the same number of sequences in the DNAStringSet and the distace matrix
          if(length(unique(names(alignment_per_primer_dna_subset_dna))) != length(unique(labels(similarity_random_sequences_matrix_subset))))
            stop("You must provide the same aligned sequences than those provided in pairwise similarity files.")
          
          # Verifying that all sequences are aligned except if using the function with methods set to swarm only
          if(length(unique(Biostrings::width(alignment_per_primer_dna_subset_dna))) != 1)
            stop('Sequences must be aligned (not the same length here).')
          
          # Arranging names of the DNAStringSet in the same order than the names of the distance matrix
          alignment_per_primer_dna_subset_dna = alignment_per_primer_dna_subset_dna[
            match(labels(similarity_random_sequences_matrix_subset), names(alignment_per_primer_dna_subset_dna))]
          
          # Creating the clusters per levels of similarity
          if(exists("Treeline", where = asNamespace("DECIPHER"), inherits = FALSE)){
            
            if(!is.null(max_time_clustering)) 
              clusters = DECIPHER::Treeline(myDistMatrix = similarity_random_sequences_matrix_subset, 
                                            myXStringSet = alignment_per_primer_dna_subset_dna,
                                            method = method, maxTime = max_time_clustering,
                                            type = "clusters", cutoff = (100 - similarity_threshold),
                                            processors = nb_threads, verbose = verbose)
            
            else clusters = DECIPHER::Treeline(myDistMatrix = similarity_random_sequences_matrix_subset, 
                                               myXStringSet = alignment_per_primer_dna_subset_dna,
                                               method = method, type = "clusters", cutoff = (100 - similarity_threshold),
                                               processors = nb_threads, verbose = verbose)
            
          }
          
          
          else{ 
            
            if(!is.null(max_time_clustering)) 
              warning('The argument "max_time_clustering" cannot be used with older version of the DECIPHER package, consider updating it.')
            
            clusters = DECIPHER::IdClusters(myDistMatrix = similarity_random_sequences_matrix_subset, 
                                            myXStringSet = alignment_per_primer_dna_subset_dna,
                                            method = method, type = "clusters", cutoff = (100 - similarity_threshold),
                                            processors = nb_threads, verbose = verbose)
            
          }
          
        }
        
        # Clustering for methods only based on a distance matrix
        else{ 
          
          if(exists("Treeline", where = asNamespace("DECIPHER"), inherits = FALSE))
            clusters = DECIPHER::Treeline(myDistMatrix = similarity_random_sequences_matrix, method = method,
                                          type = "clusters", cutoff = (100 - similarity_threshold),
                                          processors = nb_threads, verbose = verbose)
          
          else clusters = DECIPHER::IdClusters(myDistMatrix = similarity_random_sequences_matrix, method = method,
                                               type = "clusters", cutoff = (100 - similarity_threshold),
                                               processors = nb_threads, verbose = verbose)
          
        }
        
        # Converting to a tibble and renaming columns
        clusters_per_similarity_df = tibble::tibble(cbind(tibble::tibble(ACCESSION = rownames(clusters)),  
                                                          setNames(clusters, paste0("CLUSTERS_", similarity_threshold))))
        clusters_per_similarity_df$METHOD = method
        
        clusters_per_similarity_df_list[[method]] = clusters_per_similarity_df
        
      }
      
      clusters_per_similarity[[i]] = do.call(rbind, clusters_per_similarity_df_list)
      names(clusters_per_similarity)[i] = names(intra_comparisons_list)[i]
      
      gc()
      if(verbose) cat("\n--------------------- DONE ----------------------\n\n")
      
    }
    
  }
  
  return(clusters_per_similarity)
  
}
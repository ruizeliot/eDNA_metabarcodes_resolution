
vsearch_pairwise_similarity = function(metabarcodes_data, infos_data = NULL, min_similarity_threshold, 
                                       max_similarity_threshold, number_rows_considered = c(0,0), 
                                       program_name, name_col_order, name_col_primer, name_col_accession, 
                                       name_col_fragment, keep_combinations = NULL, keep_col_infos_data = NULL, 
                                       order_col_keep = NULL, verbose = T){
  
  # Verifying that all necessary columns are in data_metabarcodes
  if(length(which(colnames(metabarcodes_data) %in% c(name_col_accession, name_col_primer, name_col_fragment))) != 
     length(c(name_col_accession, name_col_primer, name_col_fragment)))
    stop('One or more columns specified in "name_col_accession" and/or "name_col_primer" and/or "name_col_fragment" are not in "metabarcodes_data".')
  
  # Verifying that a table with infos about each sequences has been provided if necessary
  if(is.null(infos_data) && !is.null(keep_col_infos_data)) 
    stop('The argument "keep_col_infos_data" is not NULL but no table has been provided in "infos_data" (NULL).')
  
  # Verifying that a table with infos about each sequences has been provided if necessary
  if(!is.null(order_col_keep) && is.null(keep_col_infos_data)) 
    stop('The argument "order_col_keep" is not NULL but no columns have been provided in "keep_col_infos_data" (NULL).')
  
  # Verifying that all necessary columns are in infos_data
  if(is.null(keep_col_infos_data)) 
    if(length(which(colnames(infos_data) %in% c(name_col_accession, keep_col_infos_data))) != 
       length(c(name_col_accession, keep_col_infos_data)))
      stop('One or more columns specified in "name_col_accession" and/or "keep_col_infos_data" are not in "infos_data".')
  else 
    if(length(which(colnames(infos_data) %in% c(name_col_accession))) != 1)
      stop('One or more columns specified in "name_col_accession" are not in "infos_data".')
  
  # Verifying that the order of columns is in the right format
  if(!is.null(order_col_keep))
    if(any(!(names(table(gsub("[[:alpha:]]|[_]", "", order_col_keep))) %in% c(1,2))))
      stop('All column names in "order_col_keep" must be pasted to the suffixes 1 or 2 depending on if it depends on the query (1) or the target (2).')
  
  # Verifying that the order of columns is in the right format
  if(!is.null(order_col_keep))
    if(length(which(name_col_accession == gsub("[0-9]+", "", order_col_keep))) != 2)
      warning('The columns name_col_accession must be present twice in "order_col_keep" (e.g. "AAA1", "AAA2").')
  
  # Verifying that the order of columns is in the right format
  if(!is.null(order_col_keep) && !is.null(keep_col_infos_data))
    if(length(which(keep_col_infos_data %in% unique(gsub("[0-9]+", "", order_col_keep)))) != length(keep_col_infos_data)) {
      if(length(which(keep_col_infos_data %in% unique(gsub("[0-9]+", "", order_col_keep)))) == 0)
      { stop('No columns specified in "keep_col_infos_data" were specified in "order_col_keep".') }
      else { warning('Not the same number of columns specified in "keep_col_infos_data" were specified in "order_col_keep".') }}
  
  # Splitting the dataframe containing the metabarcodes in a list with each elements corresponding to a set of primers
  data_split = split(tibble::tibble(metabarcodes_data), metabarcodes_data[[name_col_primer]])
  
  # Duplicating the list to get its structure but replace its elements by the results
  data_list = data_split
  
  # Loop to compute all pairwise similarities for each primers iteratively
  for(i in 1:length(data_split)){
    
    # Printing the name of the primers for which the calculation is made at each steps
    start = Sys.time()
    if(verbose){
      cat(paste0("ANALYSIS OF FRAGMENTS AMPLIFIED BY PRIMER: ", names(data_split[i]), "\n"))
      cat("-------------------------------------------------------------------\n")
      cat("\n")
    }
    
    # Creating a temporary FASTA file with the accession number as name and the metabarcodes as sequences
    fasta_file = c(rbind(paste0(">", data_split[[i]][[name_col_accession]]),
                         FRAGMENT = data_split[[i]][[name_col_fragment]]))
    write.table(fasta_file, file = "temporary_fasta_file.fasta", row.names = F, col.names = F, quote = F)
    
    # Running the vsearch program using a string in the system function
    if(!verbose) system(paste0(program_name, ' --usearch_global ', "temporary_fasta_file.fasta", 
                               ' --db ', "temporary_fasta_file.fasta", ' --threads ', 0, ' --maxid ', max_similarity_threshold/100,
                               ' --self --id ', min_similarity_threshold/100, ' --minseqlength ', 0, 
                               ' --maxaccepts ', number_rows_considered[1], ' --maxrejects ', number_rows_considered[2],
                               ' --iddef 2 --minwordmatches 0 --userfields query+target+id --maxaccepts 0 --userout ', 
                               "temporary_pairwise_similarity.txt"), intern = F, show.output.on.console = F)
    
    # Same but reprinting the messages of vsearch on the R console 
    else system(paste0(program_name, ' --usearch_global ', "temporary_fasta_file.fasta", 
                       ' --db ', "temporary_fasta_file.fasta", ' --threads ', 0, ' --maxid ', max_similarity_threshold/100,
                       ' --self --id ', min_similarity_threshold/100, ' --minseqlength ', 0, 
                       ' --maxaccepts ', number_rows_considered[1], ' --maxrejects ', number_rows_considered[2],
                       ' --iddef 2 --minwordmatches 0 --userfields query+target+id --userout ', 
                       "temporary_pairwise_similarity.txt"), intern = F)
    
    # Reading the text file containing the results of all comparisons for a primer
    similarity_file = suppressWarnings(data.table::fread("temporary_pairwise_similarity.txt"))
    
    # Placing the temporary FASTA file in the temporary directory (deleted when R session will be terminated)
    file.rename(from = paste0(getwd(), "/temporary_fasta_file.fasta"), to = paste0(tempdir(), "/temporary_fasta_file.fasta"))
    
    # Placing the file containing the results in the temporary directory (deleted when R session will be terminated)
    file.rename(from = paste0(getwd(), "/temporary_pairwise_similarity.txt"), to = paste0(tempdir(), "/temporary_pairwise_similarity.txt"))
    
    # If vsearch returned some results
    if(nrow(similarity_file) > 0){
      
      # Renaming the first column containing the similarity values
      name_col_V1 = paste0(name_col_accession, "1")
      colnames(similarity_file)[which(colnames(similarity_file) == "V1")] = name_col_V1
      
      # Renaming the second column containing the similarity values
      name_col_V2 = paste0(name_col_accession, "2")
      colnames(similarity_file)[which(colnames(similarity_file) == "V2")] = name_col_V2
      
      # Renaming the last column containing the similarity values
      colnames(similarity_file)[which(colnames(similarity_file) == "V3")] = "SIMILARITY"
      
      # Keeping only certain combinations
      if(!is.null(keep_combinations)) {
        
        similarity_file$ID = paste0(similarity_file[[name_col_V1]], ":", similarity_file[[name_col_V2]])
        
        similarity_file = subset(similarity_file, ID %in% keep_combinations)
        
      }
      
      if(!is.null(keep_col_infos_data) || length(keep_col_infos_data) != 0){
        
        # Associating the taxonomy infos with all the accession in the first column
        similarity_file = tibble::tibble(merge(subset(infos_data, select = unique(c(name_col_accession, keep_col_infos_data))), 
                                               similarity_file, by.x = name_col_accession, by.y = name_col_V1, suffixes = c("1","2")))
        colnames(similarity_file)[which(colnames(similarity_file) %in% c(name_col_accession, keep_col_infos_data))] = 
          paste0(colnames(similarity_file)[which(colnames(similarity_file) %in% c(name_col_accession, keep_col_infos_data))], "1")
        
        # Associating the taxonomy infos with all the accession in the previously second column (now its the first)
        similarity_file = tibble::tibble(merge(subset(infos_data, select = unique(c(name_col_accession, keep_col_infos_data))), 
                                               similarity_file, by.x = name_col_accession, by.y = name_col_V2, suffixes = c("1","2")))
        colnames(similarity_file)[which(colnames(similarity_file) %in% c(name_col_accession, keep_col_infos_data))] = 
          paste0(colnames(similarity_file)[which(colnames(similarity_file) %in% c(name_col_accession, keep_col_infos_data))], "2")
        
      }
      
      # Ordering the columns by alphabetical order (but keeping SIMILARITY first)
      if(is.null(order_col_keep)) similarity_file = tibble::tibble(cbind(similarity_file[,ncol(similarity_file)], 
                                                                           similarity_file[,-ncol(similarity_file)]
                                                                           [,order(colnames(similarity_file)[-ncol(similarity_file)])]))
      
      # Ordering the columns by custom order (but keeping SIMILARITY first)
      else similarity_file = subset(similarity_file, select = c("SIMILARITY", order_col_keep))
      
      # Adding the result for the primer in the list
      data_list[[i]] = similarity_file
      
    }
    
    # If vsearch doest not return any result, create a table with NA everywhere so it does not stop the function
    else data_list[[i]] = tibble::tibble(setNames(data.frame(matrix(ncol = length(c("SIMILARITY", order_col_keep)))),
                                                  c("SIMILARITY", order_col_keep)))
    
    # Printing the time making all the comparisons took for each primers (increase with the length of the primer)
    if(verbose){
      
      end = Sys.time()
      
      duration = difftime(end, start)
      
      cat("\n")
      
      cat("-------------------------------------------------------------------\n")
      
      cat(paste("Time taken:", round(duration[[1]], 2), units(duration), "\n"))
      
      cat("\n")
      
      cat("\n")
      
    }
    
  }
  
  # Returning a list with each element corresponding to one primer
  return(data_list)
  
}

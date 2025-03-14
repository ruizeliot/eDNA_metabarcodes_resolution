
gene_extraction = function(position_data, dna_data, genes, name_col_accession){
  
  # Defining the names of the columns
  all_name_start_col = paste0("START_", genes)
  all_name_end_col = paste0("END_", genes)
  all_name_orientation_col = paste0("ORIENTATION_", genes)
  all_necessary_columns = c(name_col_accession, all_name_start_col, all_name_end_col, all_name_orientation_col)
  
  # Verifying all necessary columns are in position_data
  if(!all(all_necessary_columns %in% colnames(position_data)))
    stop(paste0('The table provided in "position_data" must contain columns named: ', 
                paste0(all_necessary_columns, collapse = ", ")))
  
  # Verifying that all positions have been provided in position_data
  position_data = subset(position_data, select = all_necessary_columns)
  if(nrow(position_data[!complete.cases(position_data),]) != 0)
    stop('The table provided in "position_data" must not contain any missing positions.')
  
  # Verifying that position_data and dna_data are corresponding in term of number of sequences
  if(nrow(position_data) != length(dna_data)) 
    stop("Both position data and dna data must contain exactly the same accession numbers.")
  
  # Ordering the table to get the accession numbers in the alphabetical numbers
  position_data = position_data[order(position_data[[name_col_accession]]),]
  
  # Same for the DNA sequences to get them in the same order than the table above
  dna_data = dna_data[order(names(dna_data))]
  
  # Defining empty list to store each fragment
  fragment_list = vector(mode = "list", length = length(genes))
  
  # Looping for each genes
  for(j in 1:length(genes)){
    
    # Clean the memory at each iterations
    gc()
    
    # Getting the name of the column for the current gene
    name_start_col = all_name_start_col[j]
    name_end_col = all_name_end_col[j]
    name_orientation_col = all_name_orientation_col[j]
    
    # For each sequences
    for(i in 1:length(dna_data)){
      
      # As the mitogenome is circular, and we only have linear sequences in NCBI, 
      # a portion of a gene can be at the end of a fragment and the other at the beginning
      if(position_data[[name_start_col]][i] > position_data[[name_end_col]][i]) {
        
        # In this case, paste the fragment from the start position (START_GENE) to the end (width), 
        # and from the begining (1) to the end (END_GENE) together to get the full gene
        fragment_list_temporary = paste0(substr(dna_data[i], position_data[[name_start_col]][i], Biostrings::width(dna_data[i])), 
                                         substr(dna_data[i], 1, position_data[[name_end_col]][i]))
        
        # If it is specified that it is the complementary strand (3' to 5') by NCBI, 
        # conversion in DNAStringSet first, extrapolation of the reverse completement fragment, 
        # and then conversion back in character string using paste
        if(position_data[[name_orientation_col]][i] == "complement") fragment_list_temporary = 
            BiocGenerics::paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(fragment_list_temporary)))
        
        # Storing in the list 
        fragment_list[[j]][[i]] = fragment_list_temporary
        
        # Naming of the fragment with the accession number
        names(fragment_list[[j]][[i]]) = names(dna_data[i])
        
      }
      
      # Same but just taking the fragment between the START_GENE and END_GENE if START_GENE is lower than END_GENE
      else {
        
        fragment_list_temporary = substr(dna_data[i], position_data[[name_start_col]][i], position_data[[name_end_col]][i])
        
        if(position_data[[name_orientation_col]][i] == "complement") fragment_list_temporary = 
            BiocGenerics::paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(fragment_list_temporary)))
        
        fragment_list[[j]][[i]] = fragment_list_temporary
        
        names(fragment_list[[j]][[i]]) = names(dna_data[i])
        
      }
      
    }
    
  }
  
  # Naming the list by gene name
  names(fragment_list) = genes
  
  # Putting all sequences into columns
  extracted_genes = setNames(tibble::tibble(data.frame(do.call(cbind, lapply(fragment_list, function(x) 
    t(as.data.frame(unname(x))))))), paste0("FRAGMENT_", genes))
  
  # Adding a column with the names of the sequences
  extracted_genes = tibble::tibble(cbind(setNames(tibble::tibble(unlist(lapply(fragment_list[[1]], names))), 
                                                  name_col_accession), extracted_genes))
  
  # Returning the final table
  return(extracted_genes)
  
}

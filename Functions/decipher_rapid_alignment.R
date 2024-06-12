
decipher_rapid_alignment = function(data, name_col_group, name_col_accession, name_col_fragment,  
                                    progressive_saving = F, saving_name = NULL, 
                                    nb_processors = NULL, verbose = T, ...){
  
  # Creating a list from the table
  data = split(data, data[[name_col_group]])
  
  # Verifying that a name has been provided to save the alignment if needed
  if(is.null(saving_name) && progressive_saving)
    stop('A name must be provided in "saving_name" to save the alignment.')
  
  # Checking if the column names exist in the input
  if(!all(sapply(data, function(x) all(c(name_col_fragment, name_col_group, name_col_group) %in% colnames(x)))))
    stop('The dataframes must contain columns provided in the "name_col_fragment", "name_col_accession" and "name_col_group" arguments.')
  
  # For each element of name_col_group
  for(i in 1:length(data)){
    
    # Printing infos about each step
    if(verbose) cat("---------------------------------------------------------\n\n")
    if(verbose) cat(paste0("Performing the alignment for ", names(data)[i], ":\n\n"))

    # Converting the DNA sequences (class character) in the Biostrings::DNAStringSet format
    dna_sequence = Biostrings::DNAStringSet(data[[i]][[name_col_fragment]])
    
    # Naming by accession numbers
    names(dna_sequence) = data[[i]][[name_col_accession]]
    
    # If there is only one sequence, it is not possible to perform an alignment with other sequences so it returns the same sequence.
    if(length(dna_sequence) == 1) final = setNames(tibble::tibble(data.frame(names(dna_sequence), names(data)[[i]], paste(dna_sequence))), 
                                                   c(toupper(name_col_accession), toupper(name_col_group), "ALIGNED_SEQUENCES"))
    
    # If there is multiple sequence
    else{
      
      # Construct a chained guide tree with one leag per sequence, ordered with width of sequences
      guide_tree = lapply(order(Biostrings::width(dna_sequence), decreasing = TRUE), function(x) {
        attr(x, "height") = 0
        attr(x, "label") = names(dna_sequence)[x]
        attr(x, "members") = 1L
        attr(x, "leaf") = TRUE
        x})
      attr(guide_tree, "height") = 0.5
      attr(guide_tree, "members") = length(dna_sequence)
      class(guide_tree) = "dendrogram"
      
      # Align sequences using this chained guide tree, which does not requires more iterations or refinements processes
      aligned_sequences = DECIPHER::AlignSeqs(dna_sequence, guideTree = guide_tree, iterations = 0, 
                                    refinements = 0, verbose = verbose, processors = nb_processors, ...)
      
      # Assembling all infos in a dataframe
      final = setNames(tibble::tibble(data.frame(names(dna_sequence), names(data)[[i]], paste(aligned_sequences))), 
                       c(toupper(name_col_accession), toupper(name_col_group), "ALIGNED_SEQUENCES"))
      
    }
    
    # Saves the result per split
    if(progressive_saving) write.csv(final, paste0(saving_name, "_", names(data)[i], "_aligned.csv"), row.names = F)
    
    # Add progresively each table on top of each other
    if(i == 1) final_output = final
    else final_output = rbind(final_output, final)
    
  }
  
  # Return the final table under the tibble::tibble format
  return(tibble::tibble(final_output))
  
}

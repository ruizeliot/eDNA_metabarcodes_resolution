
metabarcode_extraction = function(data_primer, data_aligned, remove_primers = T, max_mismatch_whole_primer, indels, 
                                   name_col_alignment, name_col_accession, name_col_primer_name, 
                                   name_col_orientation, name_col_forward_primer, name_col_reverse_primer,
                                   length1_col_name, length2_col_name = NULL, verbose = T){
  
  # Verifying that all necessary columns are in data_aligned
  if(length(which(colnames(data_aligned) %in% c(name_col_accession, name_col_alignment))) != 2)
    stop(paste0('The table provided in data_aligned does not contain columns named ', 
                name_col_accession, ' and ', name_col_alignment, '.'))
  
  # Verifying that all necessary columns are in data_primer
  data_primer_columns = c(name_col_primer_name, name_col_orientation, name_col_forward_primer, 
                          name_col_reverse_primer, length1_col_name, length2_col_name)
  if(is.null(length2_col_name)) data_primer_columns = data_primer_columns[-length(data_primer_columns)]
  data_primer_columns_message = paste0(paste0(data_primer_columns[-length(data_primer_columns)], collapse = ", "),
                                       " and ", data_primer_columns[length(data_primer_columns)])
  if(length(which(colnames(data_primer) %in% data_primer_columns)) != length(data_primer_columns))
    stop(paste0('The table provided in data_primer does not contain columns named ', data_primer_columns_message, '.'))
  
  if(!any(unique(data_primer[[name_col_orientation]]) %in% c("1 x 2", "2 x 1")) || 
     length(unique(data_primer[[name_col_orientation]])) > 2)
    stop(paste0('The column of data_primer "', name_col_orientation, '" can only contain "1 x 2" and/or "2 x 1".'))
  
  # Small internal function to attribute an index to each nucleotide in an alignment without counting the gaps ("-")
  # splited_dna = an aligned sequence formatted as a vector with each character (i.e. nucleotide or gap) stored in a new element of the vector
  index_letter = function(splited_dna){
    
    # Initialisation of the number of letters already found (i.e. 0 nucleotide before starting reading the alignment)
    letter = 0
    
    # Initialisation of an empty list to store all positions
    pos_list = list()
    
    # For each nucleotide
    for(r in 1:length(splited_dna)){
      
      # Assign a position 0 if the character is a gap
      if(splited_dna[r] == "-") pos_list[r] = 0
      
      else {
        
        if(letter == 0){
          
          letter = 1
          
          # Or assign the position 1 if it is the first nucleotide in the alignment
          pos_list[r] = letter
          
        }
        
        else{
          
          letter = letter + 1
          
          # Or assign the position corresponding to the number of nucleotide already found +1
          pos_list[r] = letter
          
        }
        
      }
      
    }
    
    # Convert the positions stored in a list in a vector (unlist)
    letter_positions = unlist(pos_list)
    
    # Name each position with the corresponding characted in the alignment
    names(letter_positions) = splited_dna
    
    # Return the named positions of each character
    return(letter_positions)
    
  }
  
  # Creating an empty data.frame that will contain all extracted metabarcodes in the end
  metabarcode_data = data.frame()
  
  # For every primer given in data_primer (1 per row)
  for(i in 1:nrow(data_primer)){
    
    # Reinitialisation of the infos concerning the search after each new primer set
    position_start_aligned_list = list()
    position_end_aligned_list = list()
    primer_found = list()
    
    # Getting the name of the primer
    name_primer = data_primer[i,][[name_col_primer_name]]
    
    # Reinitialisation of the number of loop for each new primer set
    number_loop = 0
    
    # Printing infos about what the function is doing and starting of the timer
    if(verbose) start = Sys.time()
    if(verbose) cat(paste0("EXTRACTING METABARCODES USING: ", name_primer, "\n"))
    if(verbose) cat("-------------------------------------------------------------------\n")
    
    # Getting the forward and reverse primer depending on the orientation of the original DNA strand 
    # (5' to 3' = 1 x 2 // 3' to 5' = 2 x 1)
    if(data_primer[i,][[name_col_orientation]] == "2 x 1") {
      forward_primer = data_primer[[name_col_reverse_primer]][i]
      reverse_primer = BiocGenerics::paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(
        data_primer[[name_col_forward_primer]][i]))) }
    else {
      forward_primer = data_primer[[name_col_forward_primer]][i]
      reverse_primer = BiocGenerics::paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(
        data_primer[[name_col_reverse_primer]][i]))) }
    
    # Calculating the total number of nucleotides in both forward and reverse primers
    primer_length = nchar(forward_primer) + nchar(reverse_primer)
    
    # Calculating the known length of the metabarcode which will be expected
    if(is.null(length2_col_name)) known_length_fragment = mean(c(data_primer[i,][[length1_col_name]]), na.rm = T) 
    else known_length_fragment = mean(c(data_primer[i,][[length1_col_name]], data_primer[i,][[length2_col_name]]), na.rm = T) 
    
    # Checking that the length has been successfully computed
    if(length(known_length_fragment) == 0 || any(is.na(known_length_fragment)))
      stop('The argument "length1_col_name" (and "length2_col_name" if not NULL) return no length or NA.')
    
    # For each sequence in the alignment
    for(j in 1:nrow(data_aligned)){
      
      # Counting the number of loops to print the progression as a percentage
      if(number_loop == 0) number_loop = 1
      else number_loop = number_loop + 1
      if(verbose && number_loop > ceiling(nrow(data_aligned) / 100) && number_loop %% ceiling(nrow(data_aligned) / 100) == 0) 
        cat(paste0(round(number_loop / nrow(data_aligned) * 100), "% ... "))
      
      # Attributing a index to each nucleotides using the internal function index_letter
      nucleotide_index = index_letter(strsplit(data_aligned[[name_col_alignment]][j], "")[[1]])
      
      # After storing the alignment infos in nucleotide_index, removing of the gaps to get back to normal sequences and ease the search of each primer
      unaligned_sequence = gsub("-", "", data_aligned[j, ][[name_col_alignment]])
      
      # Using the function matchPattern to find all potential hydribised forward primer, allowing some mismatches and 1 indel or not
      forward_primer_search = Biostrings::matchPattern(forward_primer, unaligned_sequence, 
                                                       max.mismatch = max_mismatch_whole_primer, with.indels = indels)
      
      # Using the function matchPattern to find all potential hydribised reverse primer, allowing some mismatches and 1 indel or not
      reverse_primer_search = Biostrings::matchPattern(reverse_primer, unaligned_sequence, 
                                                       max.mismatch = max_mismatch_whole_primer, with.indels = indels)
      
      # Spotting infered matches outside the given unaligned sequence
      forward_primer_out_of_unaligned_sequence = suppressWarnings(as.numeric(sapply(1:length(forward_primer_search), function(x) 
        tryCatch(BiocGenerics::paste(forward_primer_search[x]), error = function(e) x))))
      reverse_primer_out_of_unaligned_sequence = suppressWarnings(as.numeric(sapply(1:length(reverse_primer_search), function(x) 
        tryCatch(BiocGenerics::paste(reverse_primer_search[x]), error = function(e) x))))
      
      # Getting the position of matches outside the given unaligned sequence
      forward_primer_out_of_unaligned_sequence = forward_primer_out_of_unaligned_sequence[!is.na(forward_primer_out_of_unaligned_sequence)]
      reverse_primer_out_of_unaligned_sequence = reverse_primer_out_of_unaligned_sequence[!is.na(reverse_primer_out_of_unaligned_sequence)]
      
      # Removing the matches outside the given unaligned sequence
      if(length(forward_primer_out_of_unaligned_sequence) != 0) 
        forward_primer_search = forward_primer_search[-forward_primer_out_of_unaligned_sequence]
      if(length(reverse_primer_out_of_unaligned_sequence) != 0) 
        reverse_primer_search = reverse_primer_search[-reverse_primer_out_of_unaligned_sequence]
      
      # If more than one potential hydribised forward primer has been found
      if(length(forward_primer_search) >= 2){
        
        # Get the sequence of all potential hydribised forward primer found
        matched_forward_patterns = sapply(1:length(forward_primer_search), function(x) 
          BiocGenerics::paste(forward_primer_search[x]))
        
        # Count the number of matched nucleotides for each of them
        nb_matched_nucleotides_forward = apply(Biostrings::hasLetterAt(matched_forward_patterns, forward_primer, 1:nchar(forward_primer)), 
                                               1, function(x) length(which(x)))
        
        # Saves the initial position of hydribised forward primer in the names of scores
        names(nb_matched_nucleotides_forward) = 1:length(nb_matched_nucleotides_forward)
        
      }
      
      # If more than one potential hydribised reverse primer has been found
      if(length(reverse_primer_search) >= 2){
        
        # Get the sequence of all potential hydribised reverse primer found
        matched_reverse_patterns = sapply(1:length(reverse_primer_search), function(x) 
          BiocGenerics::paste(reverse_primer_search[x]))
        
        # Count the number of matched nucleotides for each of them
        nb_matched_nucleotides_reverse = apply(Biostrings::hasLetterAt(matched_reverse_patterns, reverse_primer, 1:nchar(reverse_primer)), 
                                               1, function(x) length(which(x)))
        
        # Saves the initial position of hydribised reverse primer in the names of scores
        names(nb_matched_nucleotides_reverse) = 1:length(nb_matched_nucleotides_reverse)
        
      }
      
      # If primers need to be conserved
      if(!remove_primers){
        
        # Searching the starting position of potential hydribised forward primer
        position_start = Biostrings::start(forward_primer_search)
        
        # Searching the ending position of potential hydribised reverse primer
        position_end = Biostrings::end(reverse_primer_search)
        
      }
      
      # If primers do not need to be conserved
      else{
        
        # Same but taking the ending position of the hydribised forward primer + 1 to get the next nucleotide which is the first of the metabarcode
        position_start = Biostrings::end(forward_primer_search) + 1
        
        # Same but taking the starting position of the hydribised reverse primer - 1 to get the previous nucleotide which is the last of the metabarcode
        position_end = Biostrings::start(reverse_primer_search) - 1
        
      }
      
      # Removing the remaining matches of forward primers outside the given unaligned sequence
      if(length(which(position_start > max(nucleotide_index) | position_start <= 0)) != 0){
        if(length(forward_primer_search) >= 2) nb_matched_nucleotides_forward = 
            nb_matched_nucleotides_forward[-which(position_start > max(nucleotide_index) | position_start <= 0)]
        position_start = position_start[-which(position_start > max(nucleotide_index) | position_start <= 0)]
      }
      
      # Removing the remaining matches of reverse primers outside the given unaligned sequence
      if(length(which(position_end > max(nucleotide_index) | position_end <= 0)) != 0){
        if(length(reverse_primer_search) >= 2) nb_matched_nucleotides_reverse = 
            nb_matched_nucleotides_reverse[-which(position_end > max(nucleotide_index) | position_end <= 0)]
        position_end = position_end[-which(position_end > max(nucleotide_index) | position_end <= 0)]
      }
      
      # Indicating the number of matched nucleotides with each hydribised forward primer starting position
      if(length(position_start) >= 2) position_start = setNames(position_start, nb_matched_nucleotides_forward)
      else if(length(position_start) == 1) position_start = setNames(position_start, "1")
      
      # Indicating the number of matched nucleotides with each hydribised reverse primer ending position
      if(length(position_end) >= 2) position_end = setNames(position_end, nb_matched_nucleotides_reverse)
      else if(length(position_end) == 1) position_end = setNames(position_end, "1")
      
      # Depending on if the position found are NA or numbers, we can known which primer has been found in the unaligned sequence
      if(!is.na(position_start[1]) && !is.na(position_end[1])) primer_found[j] = "all"
      else if(is.na(position_start[1]) && !is.na(position_end[1])) primer_found[j] = "no forward"
      else if(!is.na(position_start[1]) && is.na(position_end[1])) primer_found[j] = "no reverse"
      else if(is.na(position_start[1]) && is.na(position_end[1])) primer_found[j] = "none"
      
      # If the primers where found on multiple places (multiple combination possible) on the unaligned sequences
      if(sum(c(length(position_start), length(position_end))) > 2 && length(position_start) != 0
         && length(position_end) != 0){

        # Creates a vector will all possible combinations between starting and ending positions separated by ":"
        combination_position = suppressWarnings(levels(interaction(position_start, position_end, sep = ":")))

        # Creates combinations in the same order but with the number of matched nucleotides
        score_combination_position = suppressWarnings(levels(interaction(names(position_start), names(position_end), sep = ":")))
        
        # Sum the number of matched nucleotides by each hydribised forward and reverse primers
        sum_score_combination_position = sapply(1:length(score_combination_position), function(x)
          sum(as.numeric(sub("\\:.*", "", score_combination_position[[x]])),
              as.numeric(sub(".*:", "", score_combination_position[[x]]))))

        # Saves the initial position of each combination in the names of the sum of scores
        names(sum_score_combination_position) = 1:length(sum_score_combination_position)

        # Order this vector by best sum of scores so that its names gives the position of the best combinations first
        order_combination_nucleotides_matched = as.numeric(names(sort(sum_score_combination_position, decreasing = T)))

        # Reorder the combination by best sum of scores
        combination_position = combination_position[order_combination_nucleotides_matched]

        # For every possible combinations
        for(k in 1:length(combination_position)){

          # Extracts the starting position
          start_position = as.numeric(sub("\\:.*", "", combination_position[k]))

          # Extracts the ending position
          end_position = as.numeric(sub(".*:", "", combination_position[k]))

          # Extract the starting position in the alignment based on the position found in the unaligned sequence using the nucleotide index
          start_position_alignment = which(nucleotide_index == start_position)

          # Extract the ending position in the alignment based on the position found in the unaligned sequence using the nucleotide index
          end_position_alignment = which(nucleotide_index == end_position)
          
          # Cutting the alignment between the starting and ending positions to get the isolated fragment
          if(end_position_alignment < start_position_alignment)
            aligned_fragment = substr(data_aligned[[name_col_alignment]][j], end_position_alignment, start_position_alignment)
          else aligned_fragment = substr(data_aligned[[name_col_alignment]][j], start_position_alignment, end_position_alignment)
          
          # Removing the gaps from the alignment
          fragment = gsub("-", "", aligned_fragment)

          # Storing the length of each fragment extracted with each combinations of positions
          if(k == 1) length_fragment = nchar(fragment)

          else length_fragment = c(length_fragment, nchar(fragment))

        }

        # Choosing the best combinations as the one that gives the fragment with the closest length from the known length
        # and the highest number of matched nucleotides (which.min choose the first element if equal and they are sorted)
        best_combination = combination_position[which.min(abs(length_fragment - known_length_fragment))]

        # Storing the best starting position in the alignment
        position_start_aligned_list[j] = as.numeric(sub("\\:.*", "", best_combination))

        # Storing the best ending position in the alignment
        position_end_aligned_list[j] = as.numeric(sub(".*:", "", best_combination))

      }

      # If the primer was not found twice or more
      else {

        # If the starting position has been found once, store the position in the alignment
        if(!is.na(position_start[1])) position_start_aligned_list[j] = which(nucleotide_index == position_start[1])

        # Or attribute NA if the primer has not been found in the unaligned sequence
        else position_start_aligned_list[j] = NA

        # Same for the ending position
        if(!is.na(position_end[1])) position_end_aligned_list[j] = which(nucleotide_index == position_end[1])

        else position_end_aligned_list[j] = NA

      }

    }

    # Convert in a vector starting positions for all sequences stored in a list
    position_start_aligned = unlist(position_start_aligned_list)

    # Convert in a vector ending positions for all sequences stored in a list
    position_end_aligned = unlist(position_end_aligned_list)

    # Remove the starting positions that are NA
    position_start_aligned = position_start_aligned[!is.na(position_start_aligned)]

    # Remove the ending positions that are NA
    position_end_aligned = position_end_aligned[!is.na(position_end_aligned)]

    # Choose the most frequent starting position as the most probable one by creating a frequency table
    best_position_start_aligned = as.numeric(names(sort(table(position_start_aligned), decreasing = T))[1])

    # Choose the most frequent ending position as the most probable one by creating a frequency table
    best_position_end_aligned = as.numeric(names(sort(table(position_end_aligned), decreasing = T))[1])

    # Stop the function if no primer has been found in any of the sequences
    if(length(best_position_start_aligned) == 0 || length(best_position_end_aligned) == 0)
      stop("No primers found in alignment with such parameters.")

    # Otherwise, continue running the function
    else{
      
      # Cut the alignment using the best starting and ending positions determined above
      if(best_position_end_aligned < best_position_start_aligned){
        warning("Best ending position in the alignment is behind the best starting position.")
        all_metabarcodes = substr(data_aligned[[name_col_alignment]], best_position_end_aligned, best_position_start_aligned)
      }
      else all_metabarcodes = substr(data_aligned[[name_col_alignment]], best_position_start_aligned, best_position_end_aligned)

      # Remove the gaps from the alignment
      all_metabarcodes = gsub("-", "", all_metabarcodes)

      # Store the best positions
      best_positions = paste0(best_position_start_aligned, ":", best_position_end_aligned)

      # Store the state of the search of primers in the unaligned sequences
      all_primer_found = unlist(primer_found)

      # If metabarcodes could be extracted with the best positions
      if(length(all_metabarcodes) != 0){

        # Calculate the length of the fragment with primer by adding the length of primers or not
        if(remove_primers) length_with_primers = nchar(all_metabarcodes) + primer_length
        else length_with_primers = nchar(all_metabarcodes)

        # Store all infos in a table
        if(nrow(metabarcode_data) == 0) metabarcode_data = tibble::tibble(data.frame(ACCESSION = data_aligned[[name_col_accession]],
                                                                                     PRIMERS = name_primer,
                                                                                     PRIMER_FOUND = all_primer_found ,
                                                                                     POSITION_IN_ALIGNMENT = best_positions,
                                                                                     LENGTH_WITH_PRIMERS = length_with_primers,
                                                                                     LENGTH = nchar(all_metabarcodes),
                                                                                     FRAGMENT = all_metabarcodes))

        # Bind all tables established for each primer set (each i loop)
        else metabarcode_data = rbind(metabarcode_data, tibble::tibble(data.frame(ACCESSION = data_aligned[[name_col_accession]],
                                                                                  PRIMERS = name_primer,
                                                                                  PRIMER_FOUND = all_primer_found,
                                                                                  POSITION_IN_ALIGNMENT = best_positions,
                                                                                  LENGTH_WITH_PRIMERS = length_with_primers,
                                                                                  LENGTH = nchar(all_metabarcodes),
                                                                                  FRAGMENT = all_metabarcodes)))
      }

      # Or stop the function if no metabarcodes could be cut with the best positions
      else stop("No metabarcodes could be extracted with the best positions found.")

    }

    # Print the duration of the whole process per primer
    if(verbose) {
      end = Sys.time()
      duration = difftime(end, start)
      cat("\n")
      cat("-------------------------------------------------------------------\n")
      cat(paste("Time taken:", round(duration[[1]], 2), units(duration), "\n"))
      cat("\n")
      cat("\n")
    }

    # Cleaning the memory after restarting the operation for another primer
    gc()
    
  }
  
  # Returns the final dataframe with all results
  return(metabarcode_data)
  
}

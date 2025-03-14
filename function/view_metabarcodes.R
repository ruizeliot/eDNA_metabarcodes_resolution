
view_metabarcodes = function(data, manual_data = NULL, data_primer, fragment_data, nb_sequences = 1, 
                             show_all = F, size_margin = 15, name_col_fragment, name_col_accession, 
                             name_col_manual_primer_name = NULL, name_col_manual_sequences = NULL,
                             name_col_primer_name, name_col_orientation, name_col_forward_primer, 
                             name_col_reverse_primer, length1_col_name, length2_col_name = NULL, 
                             name_col_data_metabarcodes = "FRAGMENT", name_col_data_primers = "PRIMERS", 
                             name_col_data_found = "PRIMER_FOUND", name_col_data_length = "LENGTH", 
                             name_col_data_length_with_primers = "LENGTH_WITH_PRIMERS"){
  
  # Verifying that all necessary columns are in data
  if(length(which(colnames(data) %in% c(name_col_accession, name_col_data_primers, name_col_data_metabarcodes,
                                        name_col_data_found, name_col_data_length, name_col_data_length_with_primers))) != 6)
    stop(paste0('The table provided in data does not contain columns named ', 
                name_col_accession, ', ', name_col_data_primers, ', ', name_col_data_metabarcodes, ', ', 
                name_col_data_found, ', ', name_col_data_length, ' and ', name_col_data_length_with_primers, '.'))
  
  # Verifying that all necessary columns are in fragment_data
  if(length(which(colnames(fragment_data) %in% c(name_col_accession, name_col_fragment))) != 2)
    stop(paste0('The table provided in fragment_data does not contain columns named ', 
                name_col_accession, ' and ', name_col_fragment, '.'))
  
  # Verifying that all necessary columns are in manual_data if provided
  if(!is.null(manual_data))
    if(length(which(colnames(manual_data) %in% c(name_col_accession, name_col_manual_primer_name, 
                                                 name_col_manual_sequences))) != 3)
      stop(paste0('The table provided in manual_data does not contain columns named ', 
                  name_col_accession, ', ', name_col_manual_primer_name, ' and ', name_col_manual_sequences, '.'))
  
  # Verifying that all necessary columns are in data_primer
  data_primer_columns = c(name_col_primer_name, name_col_orientation, name_col_forward_primer, 
                          name_col_reverse_primer, length1_col_name, length2_col_name)
  if(is.null(length2_col_name)) data_primer_columns = data_primer_columns[-length(data_primer_columns)]
  data_primer_columns_message = paste0(paste0(data_primer_columns[-length(data_primer_columns)], collapse = ", "),
                                       " and ", data_primer_columns[length(data_primer_columns)])
  if(length(which(colnames(data_primer) %in% data_primer_columns)) != length(data_primer_columns))
    stop(paste0('The table provided in data_primer does not contain columns named ', data_primer_columns_message, '.'))
  
  # Verifying that there is the same number of metabarcodes accross every primer set
  if(!is.null(manual_data) && var(sapply(split(manual_data, manual_data[[name_col_manual_primer_name]]), 
                                               function(x) nrow(x))) != 0)
    stop("The number of metabarcodes is not the same accross every primer set.")
  
  # Verifying that data contains the same accession numbers than manual_data
  if(!is.null(manual_data) && length(setdiff(sort(unique(data[[name_col_accession]])), 
                                             sort(unique(manual_data[[name_col_accession]])))) != 0){
    
    warning(paste0('Some accession numbers (', length(setdiff(sort(unique(data[[name_col_accession]])), 
                                                              sort(unique(manual_data[[name_col_accession]])))),
                   ') of data are not present in manual_data. The table data was subsetted accordingly'))
    
    data = data[which(data[[name_col_accession]] %in% unique(manual_data[[name_col_accession]])), ]
    
  }
  
  # Running a loop for every primer
  primers = unique(data[[name_col_data_primers]])
  for(j in 1:length(primers)){
    
    # Extracting the automatically extracted metabarcodes relative to the current primer of interest in the loop
    subset_data_primer = data[which(data[[name_col_data_primers]] == primers[j]),]
    
    # Extracting the manually extracted metabarcodes relative to the current primer of interest in the loop
    if(!is.null(manual_data)) subset_manual_data = manual_data[which(manual_data[[name_col_manual_primer_name]] == primers[j]), ]
    
    # Extracting the information relative to this primer
    subset_infos_primer = data_primer[which(data_primer[[name_col_primer_name]] == primers[j]), ]
    
    # Extracting all situations that happened when searching primers (e.g. both found, only one or none)
    cases = table(subset_data_primer[[name_col_data_found]])
    
    # Subsetting the metabarcodes per types of cases in a list
    subset_data_primer_cases = lapply(names(cases), function(name_case) 
      subset_data_primer[which(subset_data_primer[[name_col_data_found]] == name_case), ][[name_col_accession]])
    
    # Verifying that there is enough sequences in each cases
    if(any(sapply(subset_data_primer_cases, length) < nb_sequences)) 
      stop(paste0('The argument "nb_sequences" must be inferior than ', 
                  min(sapply(subset_data_primer_cases, length)), ' (less than ', 
                      nb_sequences, ' sequences in some cases).'))
    
    # For each of these cases
    for(i in 1:length(cases)){
      
      # Subsetting randomly a given number of sequences that were in this case
      random_accession = sample(subset_data_primer_cases[[i]], nb_sequences)
      
      # Getting the forward and reverse primer depending on the orientation of the original DNA strand 
      # (5' to 3' = 1 x 2 // 3' to 5' = 2 x 1)
      if(subset_infos_primer[[name_col_orientation]] == "2 x 1"){
        forward_primer = BiocGenerics::paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(
          subset_infos_primer[[name_col_forward_primer]])))
        reverse_primer = subset_infos_primer[[name_col_reverse_primer]] }
      else{
        forward_primer = subset_infos_primer[[name_col_forward_primer]]
        reverse_primer = BiocGenerics::paste(Biostrings::reverseComplement(Biostrings::DNAStringSet(
          subset_infos_primer[[name_col_reverse_primer]]))) }
      
      # Getting the complete gene corresponding to the random accession numbers chosen above
      whole_fragment = fragment_data[which(fragment_data[[name_col_accession]] %in% random_accession),][[name_col_fragment]]
      
      # Getting the manually extracted fragments corresponding to the random accession numbers chosen above
      if(!is.null(manual_data)) manual_fragment = subset_manual_data[which(subset_manual_data[[name_col_accession]] %in% random_accession),
                                                              ][[name_col_manual_sequences]]
      
      # Getting the automatically extracted fragments corresponding to the random accession numbers chosen above
      automatic_fragment = subset_data_primer[which(subset_data_primer[[name_col_accession]] %in% 
                                                      random_accession), ][[name_col_data_metabarcodes]]
      
      # Aligning all sequences together, with the manually extracted fragment
      if(!is.null(manual_data)) aligned_sequences = DECIPHER::AlignSeqs(Biostrings::DNAStringSet(c(whole_fragment, manual_fragment, 
                                                                                                     automatic_fragment, forward_primer,
                                                                                                     reverse_primer)), verbose = F)
      
      # Aligning all sequences together, without the manually extracted fragment
      else aligned_sequences = DECIPHER::AlignSeqs(Biostrings::DNAStringSet(c(whole_fragment, automatic_fragment, forward_primer,
                                                                                reverse_primer)), verbose = F)
      
      # Getting the positions of the start and the end of the primers to show a certain number of nucleotides around the primers
      if(subset_infos_primer[[name_col_orientation]] == "2 x 1") {
        
        first_pos = which(strsplit(BiocGenerics::paste(aligned_sequences)[length(aligned_sequences)], "")[[1]] != "-")[1]
        
        last_pos = tail(which(strsplit(BiocGenerics::paste(aligned_sequences)[length(aligned_sequences) - 1], "")[[1]] != "-"), 1)
        
        subset_aligned_sequences = Biostrings::DNAStringSet(substr(BiocGenerics::paste(aligned_sequences), 
                                                                   first_pos - size_margin, last_pos + size_margin))
        
      }
      
      # Same with reversed primers
      else {
        
        first_pos = which(strsplit(BiocGenerics::paste(aligned_sequences)[length(aligned_sequences) - 1], "")[[1]] != "-")[1]
        
        last_pos = tail(which(strsplit(BiocGenerics::paste(aligned_sequences)[length(aligned_sequences)], "")[[1]] != "-"), 1)
        
        subset_aligned_sequences = Biostrings::DNAStringSet(substr(BiocGenerics::paste(aligned_sequences), 
                                                                   first_pos - size_margin, last_pos + size_margin))
        
      }
      
      # Naming the different aligned sequences
      if(!is.null(manual_data)) names(aligned_sequences) = c(paste0(toupper(names(cases)[i]), " - Whole - ", random_accession), 
                                                             paste0(toupper(names(cases)[i]), " - Manual - ", random_accession),
                                                             paste0(toupper(names(cases)[i]), " - Automatic - ", random_accession), 
                                                             "Forward primer", "Reverse primer")
      else names(aligned_sequences) = c(paste0(toupper(names(cases)[i]), " - Whole - ", random_accession), 
                                        paste0(toupper(names(cases)[i]), " - Automatic - ", random_accession), 
                                        "Forward primer", "Reverse primer")
      
      # Naming the different regions of interests in the aligned sequences
      if(!is.null(manual_data)) names(subset_aligned_sequences) = c(paste0("Whole - ", random_accession), paste0("Manual - ", random_accession),
                                                                    paste0("Automatic - ", random_accession), "Forward primer", "Reverse primer")
      else names(subset_aligned_sequences) = c(paste0("Whole - ", random_accession), paste0("Automatic - ", random_accession), 
                                               "Forward primer", "Reverse primer")
      
      # Showing the alignments for each case 
      cat("---------------------------------------------------------------------------------------------------------------------")
      
      cat("\n")
      
      cat(paste("Sequences with", names(cases)[i], "primer(s) for:", primers[j]))
      
      cat("\n")
      
      cat("\n")
      
      # Showing the alignments for each case in the console
      print(subset_aligned_sequences)
      
      # Showing the alignments for each case in an HTML viewer
      if(show_all) {
        
        DECIPHER::BrowseSeqs(aligned_sequences)
        
        Sys.sleep(1)
        
      }
      
      cat("\n")
      
      cat("\n")
      
    }
    
    # Showing the summary of the search done by the function metabarcodes_extraction
    cat("---------------------------------------------------------------------------------------------------------------------")
    
    cat("\n")
    
    cat("\n")
    
    cat(paste0("PERCENTAGE SEQUENCES WITH BOTH PRIMERS FOUND: ", round(nrow(subset_data_primer[which(
      subset_data_primer[[name_col_data_found]] == "all"), ]) / 
        nrow(fragment_data) * 100, 2), "%"))
    
    cat("\n")
    
    cat("\n")
    
    cat(paste(toupper(length1_col_name), " FOR", primers[j], ":", data_primer[which(data_primer[[name_col_primer_name]] == 
                                                                          primers[j]),][[length1_col_name]]))
    
    if(!is.null(length2_col_name)){
      
      cat("\n")
      
      cat("\n")
      
      cat(paste(toupper(length2_col_name), " FOR", primers[j], ":", data_primer[which(data_primer[[name_col_primer_name]] == 
                                                                          primers[j]),][[length2_col_name]]))
      
    }
    
    cat("\n")
    
    cat("\n")
    
    if(is.null(manual_data)) cat("LENGTH WITHOUT PRIMERS FOR", primers[j], ":")
    else cat("LENGTH WITHOUT PRIMERS FOR AUTOMATICALLY EXTRACTED", primers[j], ":")
       
    cat("\n")
    
    if(is.null(manual_data)) print(summary(subset_data_primer[[name_col_data_length]]))
    else print(summary(subset_data_primer[[name_col_data_length]]))
    
    cat("\n")
    
    cat("\n")
    
    if(is.null(manual_data)) cat("LENGTH WITH PRIMERS FOR", primers[j], ":")
    else cat("LENGTH WITH PRIMERS FOR MANUALLY EXTRACTED", primers[j], ":")
    
    cat("\n")
    
    if(is.null(manual_data)) print(summary(subset_data_primer[[name_col_data_length_with_primers]]))
    else print(summary(nchar(subset_manual_data[[name_col_manual_sequences]])))
    
    cat("\n")
    
    cat("\n")
    
  }
  
}

# Requirement: "tibble::tibble"

cluster_metabarcodes_mothur_swarm = function(alignment_per_primer_df, output_path, dist_threshold, cutoff_dist_matrix,
                                       methods = "all", name_col_accession, name_col_alignment, name_col_primers, 
                                       mothur_path, vsearch_path, swarm_path, nb_threads = 1, verbose = T){
  
  ## Initialisation
  
  # Function to fill missing columns (OTU) when binding SWARM outputs at different thresholds together
  fill_missing_columns = function(df, all_columns) {
    missing_cols = setdiff(all_columns, colnames(df))
    df[missing_cols] = NA
    df[all_columns]
  }
  
  # Function to replace duplicated ID number in MOTHUR output to prevent it from crashing if the name is too long
  replace_values = function(x, replacement_table, name_col_accession) {
    for (i in seq_along(replacement_table$ID_SHORT)) {
      x = gsub(replacement_table$ID_SHORT[i], replacement_table[[name_col_accession]][i], x, fixed = T)
    }
    return(x)
  }
  
  # Defining abbreviations to be used for MOTHUR methods
  abbreviations_methods = setNames(c("opti", "furthest", "nearest", "average", "agc", "dgc", "swarm"), 
                                   c("opti_mcc", "fn", "nn", "an", "agc", "dgc", "swarm"))
  if(any(methods == "all")) methods = abbreviations_methods
  
  # Defining output paths and empty list to store outputs
  output_path = gsub("/$", "", output_path)
  mothur_path = trimws(mothur_path)
  vsearch_path = trimws(vsearch_path)
  swarm_path = trimws(swarm_path)
  mothur_cluster_list = list()
  
  # Verifying that all necessary columns are in alignment_per_primer_df
  if(length(which(colnames(alignment_per_primer_df) %in% c(name_col_accession, name_col_primers, name_col_alignment))) != 3)
    stop(paste0('The table provided in alignment_per_primer_df does not contain columns named ', 
                name_col_accession, ", ", name_col_primers, ' or ', name_col_alignment, '.'))
  
  # Verifying accession numbers do not contain separators used later in MOTHUR outputs
  if(any(grepl("[,]|[|]", alignment_per_primer_df[[name_col_accession]])))
    stop('The column provided in name_col_accession must not contain "," and "|".')
  
  
  ## Preparing aligned sequences per metabarcodes under a format accepted by MOTHUR
  
  for(primer in unique(alignment_per_primer_df[[name_col_primers]])){
    
    if(verbose) cat(paste0("Creating the clusters for primer ", primer, ":\n"))
    if(verbose) cat("-------------------------------------------------\n")

    # Subset of the initial dataset with aligned sequences for the current primer considered
    alignment_per_primer_df_subset = alignment_per_primer_df[which(alignment_per_primer_df[[name_col_primers]] == primer),]
    initial_accessions = unique(alignment_per_primer_df_subset[[name_col_accession]])
    
    # Verifying that all sequences are aligned except if using the function with methods set to swarm only
    if(length(methods) == 1){ if(methods != "swarm" && length(unique(nchar(alignment_per_primer_df_subset[[name_col_alignment]]))) != 1){
      stop('Sequences must be aligned (not the same length here) for use with another value of methods than methods == "swarm".')}}
    else if(length(methods) != 1) if(length(unique(nchar(alignment_per_primer_df_subset[[name_col_alignment]]))) != 1)
      stop('Sequences must be aligned (not the same length here) for use with another value of methods than methods == "swarm".')
    
    # Indicating all duplicated aligned sequences in the name of a single representative one
    alignment_dedup_metabarcode = tibble::tibble(aggregate(as.formula(paste0(name_col_accession, " ~ ", name_col_alignment)), 
                                                   data = alignment_per_primer_df_subset, FUN = function(x) paste0(x, collapse = "|")))
    
    if(any(methods == "swarm")){
      
      # Unaligning sequences as this automatically done by SWARM (but it needs only sequences without ambiguous base)
      alignment_per_primer_df_subset$UNALIGNED_FRAGMENT = gsub("-", "", alignment_per_primer_df_subset[[name_col_alignment]])
      swarm_current_metabarcode = alignment_per_primer_df_subset[!grepl("[^ATGC]", alignment_per_primer_df_subset$UNALIGNED_FRAGMENT),]
      max_length_metabarcode_swarm = max(nchar(swarm_current_metabarcode$UNALIGNED_FRAGMENT), na.rm = T)
      initial_accessions = unique(swarm_current_metabarcode[[name_col_accession]])
      
      # Indicating all duplicated unaligned sequences in the name of a single representative one
      swarm_dedup_metabarcode = tibble::tibble(aggregate(as.formula(paste0(name_col_accession, " ~ UNALIGNED_FRAGMENT")), 
                                                             data = swarm_current_metabarcode, FUN = function(x) paste0(x, collapse = "|")))
      
      # Creating a shorter accession number when it is too long because SWARM fails when names are too long
      swarm_dedup_metabarcode$ID_SHORT = sub("[|].*", "", swarm_dedup_metabarcode[[name_col_accession]])
      swarm_dedup_metabarcode$ID = sapply(1:nrow(swarm_dedup_metabarcode), function(i) 
        ifelse(nchar(swarm_dedup_metabarcode[[name_col_accession]][i]) <= 1500, 
               swarm_dedup_metabarcode[[name_col_accession]][i], swarm_dedup_metabarcode$ID_SHORT[i]))
      replacement_table_swarm = swarm_dedup_metabarcode[which(nchar(swarm_dedup_metabarcode[[name_col_accession]]) > 1500),]
      
      # Saving the fasta file with sequences for SWARM with the argument size set to 1 as no information on amplicon abundance in PCR
      swarm_dedup_metabarcode_fasta = c(rbind(paste0(">", swarm_dedup_metabarcode$ID, ";size=1"),
                                                  FRAGMENT = swarm_dedup_metabarcode$UNALIGNED_FRAGMENT))
      write.table(swarm_dedup_metabarcode_fasta, 
                  file = paste0(output_path, "/", primer, "_swarm.fasta"), 
                  row.names = F, col.names = F, quote = F)
      
    }
    
    if(length(methods) != 1 || length(methods[which(methods == "swarm")]) == 0){
      
        # Saving the fasta file containing all unique alignments
        alignment_dedup_metabarcode_fasta = c(rbind(paste0(">", alignment_dedup_metabarcode[[name_col_accession]]),
                                                    FRAGMENT = alignment_dedup_metabarcode[[name_col_alignment]]))
        write.table(alignment_dedup_metabarcode_fasta, 
                    file = paste0(output_path, "/", primer, " alignment.fasta"), 
                    row.names = F, col.names = F, quote = F)
        
        # Running MOTHUR functions to compute the count and distance tables
        system(sprintf(paste0(mothur_path, " ", '"#unique.seqs(fasta=%s/%s alignment.fasta)"'), output_path, primer), ignore.stdout = T)
        system(sprintf(paste0(mothur_path, " ", '"#dist.seqs(fasta=%s/%s alignment.unique.fasta, cutoff=%s, processors=%s)"'), 
                       output_path, primer, cutoff_dist_matrix, nb_threads), ignore.stdout = T)
        
    }
    
    
    ## Creating the clusters for each primer, each method, and each distance threshold
    
    otu_expanded_dist_threshold_methods = list()
    
    for(method in methods){
      
      if(verbose) cat(paste0("Creating clusters using the method: ", method))
      
      # Running the clustering with SWARM but using the same output format than for MOTHUR
      if(method == "swarm"){
        
        otu_table_swarm_list = list()
        
        # Running SWARM for each distance threshold as it is not possible to run it for multiple ones at the same time
        for(i in 1:length(dist_threshold)){

          system(paste0(swarm_path, " -r -z -t ", nb_threads, " -d ", floor(max_length_metabarcode_swarm*dist_threshold[i]),
                        " -o ", output_path, "/", primer, ".swarm_", dist_threshold[i], ".list ", output_path, "/", primer, "_swarm.fasta"),
                 show.output.on.console = F)
          
          # Formatting the output of SWARM so that is is exactly similar to a MOTHUR output
          otu_table_swarm = tibble::tibble(read.delim(paste0(output_path, "/", primer, ".swarm_", dist_threshold[i], ".list"), header = F))
          otu_table_swarm$V1 = dist_threshold[i]
          colnames(otu_table_swarm)[-c(1:2)] = paste0("Otu", 1:(ncol(otu_table_swarm)-2))
          colnames(otu_table_swarm)[1:2] = c("label", "numASVs")
          
          # Replacing short accession numbers with the ID of all those that are duplicated per group if the name was too long for SWARM
          otu_table_swarm[] = lapply(otu_table_swarm, function(col) {
            if (is.character(col)) replace_values(col, replacement_table_swarm, name_col_accession) else col})
          otu_table_swarm_list[[i]] = otu_table_swarm
          
        }
        
      }
      
      # Running the clustering with VSEARCH functions included in MOTHUR
      else if(method %in% c("agc", "dgc"))
        system(sprintf(paste0(mothur_path, " ", '"#cluster(fasta=%s/%s alignment.unique.fasta, count=%s/%s alignment.count_table, method=%s, cutoff=%s, vsearch=%s)"'), 
                       output_path, primer, output_path, primer, method, paste0(dist_threshold, collapse = "-"), vsearch_path), ignore.stdout = T)
      
      # Running the clustering with the OPTICLUST algorithm included in MOTHUR
      else if(method == "opti")
        system(sprintf(paste0(mothur_path, " ", '"#cluster(column=%s/%s alignment.unique.dist, count=%s/%s alignment.count_table, method=%s, cutoff=%s)"'),
                       output_path, primer, output_path, primer, method, paste0(dist_threshold, collapse = "-")), ignore.stdout = T)
      
      # Running the clustering with others algorithms included in MOTHUR except "unique"
      else system(sprintf(paste0(mothur_path, " ", '"#cluster(column=%s/%s alignment.unique.dist, count=%s/%s alignment.count_table, method=%s, cutoff=%s)"'),
                          output_path, primer, output_path, primer, method, max(dist_threshold, na.rm = T)), ignore.stdout = T)

      # Binding all tables under MOTHUR format returned by SWARM by adding NA for non-attributed OTU
      if(method == "swarm") otu_table = do.call(rbind, lapply(otu_table_swarm_list, fill_missing_columns, 
                                                                       all_columns = unique(unlist(lapply(otu_table_swarm_list, colnames)))))
      
      # Reading the table directly returned by MOTHUR containing clusters for other methods
      else otu_table = tibble::tibble(read.delim(paste0(output_path, "/", primer, " alignment.unique.", 
                                                        names(abbreviations_methods[which(abbreviations_methods == method)]), ".list")))
      
      # Selecting rows containing the current similarity threshold
      if(any(colnames(otu_table) == "numASVs")) otu_table_clean = subset(otu_table, label != "unique", select = -c(numASVs))
      else otu_table_clean = subset(otu_table, label != "unique", select = -c(numOtus))
      
      # Determine if clusters other than unique sequences could be found
      if(nrow(otu_table_clean) != 0){
        
        # Expand the table into separate rows for each accession number
        otu_table_long = do.call(rbind, lapply(2:ncol(otu_table_clean), function(i) {
          tibble::tibble(DIST_THRESHOLD = otu_table_clean$label,
                 CLUSTER = names(otu_table_clean)[i], ACCESSION = otu_table_clean[[i]])}))
        
        # Converting the cluster column  in character format to a numeric column
        otu_table_long$CLUSTER = as.numeric(gsub("ASV|Otu", "", otu_table_long$CLUSTER))
        
        # Expand the table into separate rows for each accession number
        otu_table_long[[name_col_accession]] = strsplit(otu_table_long[[name_col_accession]], "[,|]")
        otu_expanded = do.call(rbind, lapply(1:nrow(otu_table_long), function(i) {
          tibble::tibble(CLUSTER = otu_table_long$CLUSTER[i], 
                         DIST_THRESHOLD = paste0("CLUSTERS_", 100 - (as.numeric(otu_table_long$DIST_THRESHOLD[i]) * 100)),
                         ACCESSION = gsub(";size=1", "", otu_table_long[[name_col_accession]][[i]]))}))
        otu_expanded = otu_expanded[!is.na(otu_expanded[[name_col_accession]]),]
        
        # Retransform to a wide format with all clusters created per column
        otu_expanded_dist_threshold = tibble::tibble(reshape(data.frame(otu_expanded), idvar = name_col_accession, 
                                                     timevar = "DIST_THRESHOLD", direction = "wide"))
        colnames(otu_expanded_dist_threshold) = gsub("CLUSTER\\.", "", colnames(otu_expanded_dist_threshold))
        
        # Adding missing cluster columns as NA columns if they could not be delineated at some distance thresholds
        missing_cols = setdiff(c(name_col_accession, paste0("CLUSTERS_", 100 - (dist_threshold) * 100)), colnames(otu_expanded_dist_threshold))
        if(length(missing_cols) != 0){
          for(col in missing_cols) { otu_expanded_dist_threshold[[col]] = NA }
          otu_expanded_dist_threshold = otu_expanded_dist_threshold[, c(name_col_accession, paste0("CLUSTERS_", 100 - (dist_threshold) * 100))]}
        
      }
      
      # If clusters other than unique sequences could not be found, create and empty dataframe
      else otu_expanded_dist_threshold = setNames(tibble::tibble(cbind(tibble::tibble(ID = unique(unlist(otu_table[-c(1,2)]))), 
                                                                data.frame(matrix(NA, nrow = 1, ncol = length(dist_threshold))))),
                                                   c(name_col_accession, paste0("CLUSTERS_", 100 - (dist_threshold) * 100)))
      
      # List of missing accession number since their distance is greater than the maximum limit
      missing_accessions = unique(c(setdiff(unique(otu_expanded_dist_threshold[[name_col_accession]]), initial_accessions),
                                    setdiff(initial_accessions, unique(otu_expanded_dist_threshold[[name_col_accession]]))))
      
      # Adding missing accession numbers as new clusters for each threshold (inferior to the maximum threshold checked)
      if(length(missing_accessions) != 0){
        
        # Compute the maximum index of cluster per column
        max_values = sapply(otu_expanded_dist_threshold[2:ncol(otu_expanded_dist_threshold)], function(x) 
          if (is.numeric(x)) max(x, na.rm = TRUE) else NA)
        
        # Create one row for each new accession number, by attributing the maximum for each column + the new index number (increment)
        otu_expanded_dist_new = do.call(rbind, lapply(seq_along(missing_accessions), function(i) {
          accession = missing_accessions[i]
          increment = i  
          new_row = max_values + increment 
          data.frame(ACCESSION = accession, as.list(new_row), check.names = FALSE)
        }))
        
        # Ensure columns match before merging
        otu_expanded_dist_new[names(otu_expanded_dist_threshold)] = otu_expanded_dist_new
        otu_expanded_dist_threshold = rbind(otu_expanded_dist_threshold, otu_expanded_dist_new)
        
      }
      
      # Creating a list per method with clusters for each similarity threshold
      otu_expanded_dist_threshold$METHOD = method
      otu_expanded_dist_threshold_methods[[method]] = otu_expanded_dist_threshold
      
      if(verbose) cat(" - DONE\n")
      
    }
    
    # Binding the table with all methods together, and putting them in a list
    mothur_cluster_list[[primer]] = do.call("rbind", otu_expanded_dist_threshold_methods)
    
    if(verbose) cat("\n--------------------- DONE ----------------------\n\n")
    
  }
  
  # Removing the MOTHUR logfile created each time MOTHUR function is run and returning the final list
  file.remove(grep("logfile", list.files(pattern = "mothur"), value = T))
  return(mothur_cluster_list)
  
}

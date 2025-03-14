perception_mock_communities = function(clusters_per_primers, nb_random_community, 
                                       nb_seq_sampled, nb_seq_per_taxa = NULL, clusters_cols_prefix,
                                       name_col_taxa, name_col_accession, output_path){
  
  # Creates an output directory if it does not already exists
  if(!dir.exists(output_path)) dir.create(output_path)
  
  number_loop = 0
  
  # Checking that all elements of the list contains the same same sequences affiliated to the same taxa
  if(stats::var(sapply(clusters_per_primers, function(x) length(unique(x[[name_col_taxa]])))) != 0 ||
     stats::var(sapply(clusters_per_primers, function(x) length(unique(x[[name_col_accession]])))) != 0)
    stop('Each element of the list "clusters_per_primers" must contain the same taxa and sequence names.')
  
  # Creates the reference taxonomic tablee and get the thresholds from column names on the basis of the first element
  clusters_per_primers = clusters_per_primers[sort(names(clusters_per_primers))]
  taxo_data = tibble::tibble(clusters_per_primers[[1]])
  thresholds = as.numeric(sub(clusters_cols_prefix, "", grep(paste0("^", clusters_cols_prefix), colnames(
    clusters_per_primers[[1]]), value = T)))
  
  for(i in 1:length(nb_seq_sampled)){
    
    nb_perceived_OTU_random_community_temp = list()
    
    for(j in 1:nb_random_community){
      
      # Computing the total number of taxa of a certain community
      if(is.null(nb_seq_per_taxa)) nb_seq_per_taxa = c(seq(0.5, 1, by = 0.1), seq(1, 2, by = 0.2)[-1])
      nb_taxa_community = nb_seq_sampled[[i]] * sample(nb_seq_per_taxa, 1)
      
      # Subsetting a certain amount of taxa to represent all taxa present in the community
      random_taxa = sample(unique(taxo_data[[name_col_taxa]]), nb_taxa_community)
      random_taxa_df = taxo_data[taxo_data[[name_col_taxa]] %in% random_taxa,]
      
      # Selecting a random number of sequences retrieved and sequenced per taxa
      random_seq_per_taxa_df = do.call(rbind, lapply(split(random_taxa_df, random_taxa_df[[name_col_taxa]]), function(x)
        x[sample(nrow(x), sample(1:nrow(x), 1)), ]))
      
      # Selecting a determined number of sequences that were kept for analysis (e.g. good sequencing quality)
      random_sequences_df = random_seq_per_taxa_df[sample(nrow(random_seq_per_taxa_df), nb_seq_sampled[[i]], replace = T), ]

      # Subsetting the list with the clusters per primers and thresholds
      clusters_per_primers_random = lapply(clusters_per_primers, function(x) 
        x[x[[name_col_accession]] %in% random_sequences_df[[name_col_accession]],])

      # Computing the number of perceived OTU per primers
      nb_perceived_OTU_per_primers = tibble(setNames(do.call(cbind, lapply(clusters_per_primers_random, function(y) 
        data.frame(apply(y[, grep(paste0("^", clusters_cols_prefix), colnames(y))], 2, function(x) length(unique(x)))))), names(clusters_per_primers_random)))

      # Saving the characteristics of the community and the perception depending on the primer set in a list
      nb_perceived_OTU_per_primers = tibble(cbind(tibble(SEQ = nrow(random_sequences_df)), 
                                                  tibble(TAXA = length(unique(random_sequences_df[[name_col_taxa]]))),
                                                  tibble(SIMILARITY = thresholds),
                                                  nb_perceived_OTU_per_primers))

      # Saving the result in a dataframe
      if(j == 1) nb_perceived_OTU_per_primers_all_sim = nb_perceived_OTU_per_primers
      else nb_perceived_OTU_per_primers_all_sim = rbind(nb_perceived_OTU_per_primers_all_sim, 
                                                        nb_perceived_OTU_per_primers)
      
    }
    
    # Saving the result per sequence number interval to avoid getting a too large file in the end
    write.csv(nb_perceived_OTU_per_primers_all_sim, paste0(output_path, "/random_com_", nb_seq_sampled[i], "_seq.csv"), row.names = F)
    
    # Printing the progression of the simulation as a percentage
    if(number_loop == 0) number_loop = 1
    else number_loop = number_loop + 1
    if(number_loop > ceiling(length(nb_seq_sampled) / 100) && number_loop %% 
       ceiling(length(nb_seq_sampled) / 100) == 0) 
      cat(paste0(round(number_loop / length(nb_seq_sampled) * 100), "% ... "))
    
  }
  
  # Opening each file per nb_seq_sampled to summarise the mean number of taxa, OTU and diversity bias in mock communities
  files_communities = list.files(paste0(output_path, "/"), pattern = "_seq[.]csv")
  for(i in 1:length(files_communities)){
    
    # Performing the analysis independently for each similarity threshold so that a list per cutoff is returned
    community_result = read.csv(paste0(paste0(output_path, "/"), files_communities[i]))
    community_result_list = split(community_result, community_result$SIMILARITY)
    
    if(i == 1) summary_per_threshold = setNames(lapply(community_result_list, function(x) 
      tibble(NB_SEQ = x[1,1][[1]],
             MEAN_NB_TAXA = mean(x$TAXA),
             SD_NB_TAXA = sd(x$TAXA),
             PRIMERS = primers,
             MEAN_NB_OTU = apply(x[,-c(1:3)], 2, mean),
             SD_NB_OTU = apply(x[,-c(1:3)], 2, sd),
             MEAN_PERCENT_PERCEIVED = apply(apply(x[,-c(1:3)], 2, function(y) 
               y / x$TAXA * 100), 2, mean),
             SD_PERCENT_PERCEIVED = apply(apply(x[,-c(1:3)], 2, function(y) 
               y / x$TAXA * 100), 2, sd))), names(community_result_list))
    
    else summary_per_threshold = setNames(lapply(1:length(community_result_list), function(x) rbind(summary_per_threshold[[x]], tibble(
      NB_SEQ = community_result_list[[x]][1,1][[1]],
      MEAN_NB_TAXA = mean(community_result_list[[x]]$TAXA),
      SD_NB_TAXA = sd(community_result_list[[x]]$TAXA),
      PRIMERS = primers,
      MEAN_NB_OTU = apply(community_result_list[[x]][,-c(1:3)], 2, mean),
      SD_NB_OTU = apply(community_result_list[[x]][,-c(1:3)], 2, sd),
      MEAN_PERCENT_PERCEIVED = apply(apply(community_result_list[[x]][,-c(1:3)], 2, function(y) 
        y / community_result_list[[x]]$TAXA * 100), 2, mean),
      SD_PERCENT_PERCEIVED = apply(apply(community_result_list[[x]][,-c(1:3)], 2, function(y) 
        y / community_result_list[[x]]$TAXA * 100), 2, sd)))), names(community_result_list))
    
  }
  
  # Converting the list to a single table by grouping list elements together and adding the threshold in the SIMILARITY column
  summary_per_threshold_combined = do.call(rbind, Map(function(df, name) {
    df$SIMILARITY = name
    df}, summary_per_threshold, names(summary_per_threshold)))
  summary_per_threshold_combined = cbind(summary_per_threshold_combined[,9], summary_per_threshold_combined[,1:8])
  
  return(tibble::tibble(summary_per_threshold_combined))
  
}

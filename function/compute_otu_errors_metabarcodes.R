compute_otu_errors_metabarcodes = function(clusters_taxa_list, name_col_method, name_col_taxa, clusters_cols_prefix, verbose = T){
  
  # Initialize lists to store results
  intra_errors_summary = list()
  inter_errors_summary = list()
  
  # Compute error rates independently for each metabarcode
  for(primer_name in names(clusters_taxa_list)){
    
    if(verbose) cat(paste0("Computing errors rates for primer ", primer_name, ":\n"))
    if(verbose) cat("-------------------------------------------------\n")
    
    cluster_cols = grep(paste0("^", clusters_cols_prefix), names(clusters_taxa_list[[primer_name]]), value = TRUE)
    
    # Iterate over each cluster column corresponding to different similarity thresholds
    for(cluster_col in cluster_cols){
      
      similarity_threshold = sub(clusters_cols_prefix, "", cluster_col)
      intra_list = list()
      inter_list = list()
      
      # Iterate over each clustering method to compute independently false-positive and false-negative error rates
      for(method in unique(clusters_taxa_list[[primer_name]][[name_col_method]])){
        
        # Creates a specific table for each similarity threshold and clustering method
        df_method = clusters_taxa_list[[primer_name]][clusters_taxa_list[[primer_name]][[name_col_method]] == method, ]
        df_method = df_method[!is.na(df_method[[cluster_col]]), ]
        
        if(nrow(df_method) > 0){
          
          # Compute the number of intra-BIN errors (false-positives) as those belonging to two different clusters
          bin_summary = aggregate(df_method[[cluster_col]], by = list(df_method[[name_col_taxa]]), function(x) as.integer(length(unique(x)) != 1))
          num_intra = aggregate(df_method[[cluster_col]], by = list(df_method[[name_col_taxa]]), function(x) length(unique(x)))
          
          # Compute counts and percentages about intra-BIN errors (false-positives)
          num_bins_with_intra = sum(unlist(bin_summary$x))  
          percent_bins_with_intra = (num_bins_with_intra * 100) / nrow(bin_summary)
          mean_num_intra_per_bin = mean(num_intra$x)
          mean_percent_intra_per_bin = mean(num_intra$x * 100 / nrow(num_intra))
          
          # Regroup all summary statistics about intra-BIN errors (false-positives) into a table
          intra_list[[method]] = data.frame(PRIMERS = primer_name, METHOD = method, NUMBER_INTRA_ERRORS = num_bins_with_intra,
                                            PERCENT_INTRA_ERRORS = percent_bins_with_intra, MEAN_NUMBER_INTRA_ERRORS_PER_BIN = mean_num_intra_per_bin,
                                            MEAN_PERCENT_INTRA_ERRORS_PER_BIN = mean_percent_intra_per_bin)
          
          # Compute counts and percentages about inter-BIN errors (false-negatives) 
          cluster_summary = aggregate(df_method[[name_col_taxa]], by = list(df_method[[cluster_col]]), function(x) as.integer(length(unique(x)) != 1))
          num_inter = aggregate(df_method[[name_col_taxa]], by = list(df_method[[cluster_col]]), function(x) length(unique(x)))
          
          # Regroup all summary statistics about inter-BIN errors (false-negatives) into a table
          num_clusters_with_inter = sum(unlist(cluster_summary$x))  
          percent_clusters_with_inter = (num_clusters_with_inter * 100) / nrow(cluster_summary)
          mean_num_inter_per_cluster = mean(num_inter$x)
          mean_percent_inter_per_cluster = mean(num_inter$x * 100 / nrow(num_inter))
          
          # Regroup all summary statistics about inter-BIN errors (false-negatives) into a table
          inter_list[[method]] = data.frame(PRIMERS = primer_name, METHOD = method, NUMBER_INTER_ERRORS = num_clusters_with_inter,
                                            PERCENT_INTER_ERRORS = percent_clusters_with_inter, MEAN_NUMBER_INTER_ERRORS_PER_CLUSTER = mean_num_inter_per_cluster,
                                            MEAN_PERCENT_INTER_ERRORS_PER_CLUSTER = mean_percent_inter_per_cluster)
        } 
        
        # Creates empty tables if clusters could not be created for certain similarity thresholds and methods
        else{
          
          intra_list[[method]] = data.frame(PRIMERS = primer_name, METHOD = method, NUMBER_INTRA_ERRORS = NA, PERCENT_INTRA_ERRORS = NA, 
                                            MEAN_NUMBER_INTRA_ERRORS_PER_BIN = NA, MEAN_PERCENT_INTRA_ERRORS_PER_BIN = NA)
          
          inter_list[[method]] = data.frame(PRIMERS = primer_name, METHOD = method, NUMBER_INTER_ERRORS = NA, PERCENT_INTER_ERRORS = NA,
            MEAN_NUMBER_INTER_ERRORS_PER_CLUSTER = NA, MEAN_PERCENT_INTER_ERRORS_PER_CLUSTER = NA)
          
        }
        
      }
      
      # Append results for the current metabarcode and similarity threshold to the intra-BIN errors (false-positives) summary table
      intra_errors_summary[[similarity_threshold]] = tibble::tibble(rbind(
        intra_errors_summary[[similarity_threshold]], do.call(rbind, intra_list)))
      
      # Append results for the current metabarcode and similarity threshold to the inter-BIN errors (false-negatives) summary table
      inter_errors_summary[[similarity_threshold]] = tibble::tibble(rbind(
        inter_errors_summary[[similarity_threshold]], do.call(rbind, inter_list)))
      
    }
    
    if(verbose) cat("\n--------------------- DONE ----------------------\n\n")
    
  }
  
  return(list(INTRA_ERRORS = intra_errors_summary, INTER_ERRORS = inter_errors_summary))
  
}

## Replacement: 'tibble' + 'stringr'

# Info: Automatically determines which new informations can replace missing or incorrect infos, for differently sized tables.

# Note: The ID column can be anything (e.g. TAXID, SPECIES) but its values neeed to be unique.
# Note: A deduplification is made on the basis of the id if necessary with the function consensus deduplification.

add_infos = function(old_data, new_data, id){
  
  consensus = list()
  
  if(length(setdiff(colnames(old_data), colnames(new_data))) != 0){
    
    new_data[setdiff(colnames(old_data), colnames(new_data))] = NA
    
    new_data = new_data[match(colnames(old_data), colnames(new_data))]
    
  }
  
  replace_partial_df_rows = function(data_old, data_new, common_reference, common_column){
    
    for(i in 1:length(unique(common_reference))){
      
      data_to_replace = data_old[ , which(colnames(data_old) == 
                                            common_column):ncol(data_old)][data_old[ , which(colnames(data_old) == 
                                                                                               common_column):ncol(data_old)][[common_column]] %in% 
                                                                             common_reference[i], ]
      
      if(i == 1){
        
        pos_replaced = row.names(data_to_replace)
        
        replaced_data = cbind(data_old[as.numeric(row.names(data_to_replace)),][ , 1:which(colnames(data_old) == common_column) - 1], data_new[i, ])
        
      }
      
      else{
        
        pos_replaced = c(pos_replaced, row.names(data_to_replace))
        
        replaced_data = rbind(replaced_data, 
                              cbind(data_old[as.numeric(row.names(data_to_replace)),][ , 1:which(colnames(data_old) == common_column) - 1],
                                    data_new[i, ]))
        
      } 
      
    }
    
    data_old[pos_replaced, ] = replaced_data
    
    data_old
    
  }
  
  if(length(unique(new_data[[id]])) != length(new_data[[id]])) new_data = suppressWarnings(consensus_deduplification(new_data, id))
  
  old_data_deduplicated = suppressWarnings(consensus_deduplification(old_data, id))
  
  new_data = new_data[ , which(colnames(new_data) == id):ncol(old_data)] 
  
  new_data = new_data[complete.cases(new_data[[id]]), ]
  
  data_merged = merge(new_data, old_data_deduplicated, all = T)
  
  duplicates = data_merged[which(duplicated(data_merged[[id]]) | duplicated(data_merged[[id]], fromLast = TRUE)), ]
  
  deduplicated_df = duplicates[!duplicated(duplicates[[id]]),]
  
  for(i in 1:length(unique(duplicates[[id]]))){
    
    data = duplicates[duplicates[[id]] == unique(duplicates[[id]])[i],]
    
    for(j in 1:ncol(data)){ 
      
      col = data[[j]][!is.na(data[[j]])]
      
      col = col[col != "NA"]
      
      {if(length(unique(col)) == 1) consensus_info = unique(col)
        
        else if(length(unique(col)) > 1) consensus_info = paste(unique(col), collapse = " or ")
        
        else if(length(unique(col)) == 0) consensus_info = NA}
      
      if(length(grep("or", consensus_info)) != 0){
        
        different_values = strsplit(consensus_info, split = " ")[[1]]
        
        different_values = different_values[!different_values == "or"]
        
        if(any(duplicated(different_values))) {
          
          unique_values = different_values[!duplicated(different_values)]
          
          consensus_info = paste(unique_values, collapse = " or ")
          
        }
        
      }
      
      consensus[j] = consensus_info
      
    }
    
    if(anyNA(deduplicated_df)) deduplicated_df[is.na(deduplicated_df)] = "NAN"
    
    deduplicated_df[deduplicated_df[[id]] == deduplicated_df[[id]][i], ][1,] = unlist(consensus)
    
  }
  
  if(all(apply(deduplicated_df, 1, function(x) all(is.na(x))))) final_output = old_data
  
  else{
    
    deduplicated_df[deduplicated_df == "NAN"] = NA
    
    deduplicated_df = deduplicated_df[which(!apply(deduplicated_df, 1, function(x) all(is.na(x)))),]
    
    final_output = suppressWarnings(tibble::tibble(replace_partial_df_rows(data_old = data.frame(old_data), data_new = deduplicated_df, 
                                                                   common_reference = deduplicated_df[!is.na(deduplicated_df[[id]]),][[id]], common_column = id)))
    
  }
  
  final_output
  
}
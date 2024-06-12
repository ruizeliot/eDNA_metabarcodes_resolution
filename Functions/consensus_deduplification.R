## Replacement : 'tibble'

# Info: Gather element of the ID column and check if different values are present for the same ID. 
# If so, it removes the NA and combine the different values by separating (e.g. "x or y or z").

consensus_deduplification = function(data, id){
  
  consensus_infos = list()
  
  if(anyNA(data[[id]])) data_no_na = data[!is.na(data[[id]]),]
  
  else data_no_na = data
  
  duplicates = data_no_na[which(duplicated(data_no_na[[id]]) | duplicated(data_no_na[[id]], fromLast = TRUE)), ]
  
  if(nrow(duplicates) != 0){
    
    deduplicated_df = duplicates[!duplicated(duplicates[[id]]),]
    
    not_duplicates = data_no_na[!(data_no_na[[id]] %in% deduplicated_df[[id]]), ]
    
    deduplicated_df[] = mapply(FUN = as, deduplicated_df, sapply(data_no_na, class), SIMPLIFY = FALSE)
    
    for(i in 1:length(unique(duplicates[[id]]))){
      
      duplicates_partitioned = duplicates[duplicates[[id]] == deduplicated_df[[id]][i], ]
      
      differences = apply(duplicates_partitioned, 2, function(x) length(unique(x)))
      
      if(all(differences == 1) || nrow(duplicates_partitioned) == 1) 
        deduplicated_df[deduplicated_df[[id]] == deduplicated_df[[id]][i], ] = duplicates_partitioned[1,]
      
      else {
        
        for(j in 1:ncol(data_no_na)){ 
          
          col = duplicates_partitioned[,j][[1]]
          
          col_differences = apply(duplicates_partitioned, 2, function(x) length(unique(x))) 
          
          if(all(col_differences == 1)) consensus = unique(col)
          
          else {
            
            col_no_na = col[!is.na(col)]
            
            if(length(unique(col_no_na)) == 1) consensus = unique(col)
            
            else consensus = paste(unique(col), collapse = " or ")
            
          }
          
          consensus_infos[j] = consensus
          
        }
        
        unique_row = t(data.frame(unlist(consensus_infos)))
        
        colnames(unique_row) = colnames(data_no_na)
        
        deduplicated_df = data.frame(deduplicated_df)
        
        deduplicated_df[deduplicated_df[[id]] == deduplicated_df[[id]][i], ] = unique_row
        
      }
      
    }
    
    deduplicated_data = rbind(deduplicated_df, not_duplicates)
    
    if(anyNA(data[[id]])){
      
      data_na = data[is.na(data[[id]]),]
      
      data_na = data_na[!duplicated(data_na),]
      
      deduplicated_data = rbind(deduplicated_data, data_na)
      
    }
    
    tibble::tibble(deduplicated_data)
    
  }
  
  else{
    
    warning("All elements of the column are unique.")
    
    tibble::tibble(data)
    
  }
  
}
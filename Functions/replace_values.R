## Replacement: 'tibble'

# Info: Small function to replace the values of a big dataframe by the corrected values obtained indepedently, for multiple variables.

# Note: Use total_replacement = "NA" if the initial data must be keeped for all columns in case of conflict with the new data.
# Note: Otherwise, provide the column names of the columns needing total replacement. Provide "all" to select all columns.

replace_values = function(data_original, data_model, variables_original, variables_model = NULL, 
                          id_original, id_model = NULL, total_replacement){
  
  if(is.null(variables_model)) variables_model = variables_original
  
  if(is.null(id_model)) id_model = id_original
  
  if(total_replacement == "all") total_replacement = variables_original
  
  if(!all(variables_original %in% colnames(data_original)) || !all(variables_model %in% colnames(data_model)))
    stop("Some columns are not in the original dataframes.")
  
  if(length(variables_original) != length(variables_model)) stop("Different number of model variables and variables to replace.")
  
  for(i in 1:length(variables_original)){
    
    if(nrow(data_model) > 0){
      
      for(m in 1:nrow(data_model)){
        
        pos = list()
        
        position = which(data_model[m, ][[id_model]] == data_original[[id_original]])
        
        if(variables_original[i] %in% total_replacement) position_corrected = position
        
        else{
          
          for(n in 1:length(position)){ 
            
            if(is.na(data_original[[variables_original[i]]][position[n]])) {
              
              pos[n] = position[n]
              
            }
            
            else pos[n] = NA
            
          }
          
          position_corrected = unlist(pos)
          
          position_corrected = position_corrected[!is.na(position_corrected)]
          
        }
        
        data_original[[variables_original[i]]] = replace(data_original[[variables_original[i]]], position_corrected, 
                                                         data_model[m, ][[variables_model[i]]])
        
      }
      
    }
    
  }
  
  return(tibble::tibble(data_original))
  
}
# Requirement : 'taxizedb' + 'tibble::tibble'

# Info: Get accepted species name, species rank & higher ranks names from taxid (rapidity made possible by constructing a local database).

id_to_classification = function(id, require_id = NULL, db = "ncbi"){
  
  print_time = function(start_time, end_time, message){
    
    time = round(end_time - start_time, 2)
    
    cat(paste0(message, " - Time taken : ", time, " seconds", "\n", ""))
    
    cat("\n")
    
  }
  
  na_rm = function(x) {unlist(x)[!is.na(unlist(x))]}
  
  superkingdom = list()
  kingdom = list()
  phylum = list()
  class = list()
  order = list()
  family = list()
  genus = list()
  required_id = list()
  
  start_time_retrieving = Sys.time()
  
  output = taxizedb::classification(id, db)
  
  end_time_retrieving = Sys.time()
  
  cat(paste0("Data retrieved in ", round(end_time_retrieving - start_time_retrieving, 2), " seconds\n"))
  
  start_time_processing = Sys.time()
  
  for(i in 1:length(id)){
    
    table = data.frame(ranks = unlist(output[[as.character(id[i])]]$rank), 
                       names = unlist(output[[as.character(id[i])]]$name),
                       ids = unlist(output[[as.character(id[i])]]$id))
    
    table = table[!duplicated(table$names), ]
    
    if(length(table) == 0) {
      
      superkingdom[i] = NA 
      kingdom[i] = NA
      phylum[i] = NA
      class[i] = NA
      order[i] = NA 
      family[i] = NA
      genus[i] = NA
      
    }
    
    else{
      
      if(length(subset(table, ranks == "superkingdom")$names) == 0) superkingdom[i] = NA 
      
      else superkingdom[i] = subset(table, ranks == "superkingdom")$names
      
      if(length(subset(table, ranks == "kingdom")$names) == 0) kingdom[i] = NA
      
      else kingdom[i] = subset(table, ranks == "kingdom")$names
      
      if(length(subset(table, ranks == "phylum")$names) == 0) phylum[i] = NA
      
      else phylum[i] = subset(table, ranks == "phylum")$names
      
      if(length(subset(table, ranks == "class")$names) == 0) class[i] = NA
      
      else class[i] = subset(table, ranks == "class")$names
      
      if(length(subset(table, ranks == "order")$names) == 0) order[i] = NA 
      
      else order[i] = subset(table, ranks == "order")$names
      
      if(length(subset(table, ranks == "family")$names) == 0) family[i] = NA
      
      else family[i] = subset(table, ranks == "family")$names
      
      if(length(subset(table, ranks == "genus")$names) == 0) genus[i] = NA
      
      else genus[i] = subset(table, ranks == "genus")$names
      
      if(!is.null(require_id)) {
        
        if(length(subset(table, ranks == require_id)$names) == 0) required_id[i] = NA
        
        else required_id[i] = subset(table, ranks == require_id)$id
        
      }
      
    }
    
  }
  
  
  species = taxizedb::taxid2name(id, db)
  
  rank = taxizedb::taxid2rank(id, db)
  
  end_time_processing = Sys.time()
  
  cat(paste0("Data processed in ", round(end_time_processing - start_time_processing, 2), " seconds\n"))
  
  output = tibble::tibble(data.frame(ID = id, SPECIES = species, RANK = rank, GENUS = unlist(genus), FAMILY = unlist(family), ORDER = unlist(order),
                             CLASS = unlist(class), PHYLUM = unlist(phylum), KINGDOM = unlist(kingdom), 
                             SUPERKINGDOM = unlist(superkingdom)))
  
  
  if(!is.null(require_id)) tibble::tibble(cbind(output, REQUIRED_ID = unlist(required_id))) 
  
  else output
  
  output
  
}
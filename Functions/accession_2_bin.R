## Requirement: 'bold' + 'tibble'

accession_2_bin = function(accession, type = "accession", division_number = NULL, csv_name = NULL, progressive_saving = F, final_saving = T,
                           remove_ncbi_version = F, unique = F, multiple_accession_list = F){
  
  if(is.null(division_number) && !multiple_accession_list)
    stop('The parameter "division_number" must be provided if "multiple_accession_list" is set to FALSE.')
  
  if(is.null(csv_name) && progressive_saving || is.null(csv_name) && final_saving)
    stop('The parameter "csv_name" must be provided if "progressive_saving" and/or "final_saving" is set to TRUE.')
  
  if(!is.null(division_number) && division_number > length(accession)) 
    stop(paste0('The parameter "division_number" must be inferior the length of accession (', length(accession), ' here).'))
  
  if(!(type %in% c("accession", "bold_id"))) stop('The parameter "type" must be either "accession" or "bold_id".')
  
  if(remove_ncbi_version) accession = stringr::word(accession, 1, sep = stringr::fixed("."))
  
  if(unique) accession = unique(accession)
  
  if(multiple_accession_list) accession_string_splitted = accession
  
  else{
    
    accession_groups = rep(seq(ceiling(length(accession)/division_number)), each = division_number)[1:length(accession)]
    
    accession_string_splitted = split(accession, accession_groups)
    
  }
  
  for(i in 1:length(accession_string_splitted)){
    
    start = Sys.time()
    
    cat(paste0("Searching for subset ", i, "/", length(accession_string_splitted), "\n"))
    
    found_bin = tryCatch(bold::bold_specimens(ids = accession_string_splitted[[i]]), 
                         error = function(e) data.frame(ACCESSION = accession_string_splitted[[i]], 
                                                        BIN = NA, BOLD_TAXID = NA, BOLD_SPECIES = NA, 
                                                        BOLD_GENUS = NA, BOLD_FAMILY = NA, BOLD_ORDER = NA, SIMILARITY = NA))
    
    if(type == "bold_id") found_bin = found_bin[match(accession_string_splitted[[i]], found_bin[,1]),]
    
    else found_bin = found_bin[match(accession_string_splitted[[i]], found_bin[,2]),]
    
    if(multiple_accession_list) found_bin = found_bin[which(found_bin$bin_uri != "")[1], ]
    
    if(ncol(found_bin) == 7) final = found_bin
    
    else{
      
      if(type == "accession") found_data = data.frame(ACCESSION = found_bin[,2], BIN = found_bin[,8], 
                              BOLD_TAXID = found_bin[,21], BOLD_SPECIES = found_bin[,22], 
                              BOLD_GENUS = found_bin[,20], BOLD_FAMILY = found_bin[,16], 
                              BOLD_ORDER = found_bin[,14], SIMILARITY = 1)
      
      else found_data = data.frame(ACCESSION = found_bin[,1], BIN = found_bin[,8], 
                                   BOLD_TAXID = found_bin[,21], BOLD_SPECIES = found_bin[,22], 
                                   BOLD_GENUS = found_bin[,20], BOLD_FAMILY = found_bin[,16], 
                                   BOLD_ORDER = found_bin[,14], SIMILARITY = 1)
      
      if(ncol(data.frame(apply(found_data, 2, function(x) gsub("^$|^ $", NA, x)))) == ncol(found_data)){
        
        found_data = data.frame(apply(found_data, 2, function(x) gsub("^$|^ $", NA, x)))
        
      }
      
      found_data$BOLD_TAXID = as.numeric(found_data$BOLD_TAXID)
      
      if(!multiple_accession_list && length(found_data$ACCESSION) != length(accession_string_splitted[[i]])) {
        
        all_accession = data.frame(ACCESSION = accession_string_splitted[[i]], 
                                   BIN = NA, BOLD_TAXID = NA, BOLD_SPECIES = NA, BOLD_GENUS = NA, BOLD_FAMILY = NA, 
                                   BOLD_ORDER = NA, SIMILARITY = NA)
        
        final = merge(all_accession, found_data, by = "ACCESSION", all.y = F)
        
        final = final[,-which(grepl(".x", colnames(final)))]
        
        colnames(final) = gsub(".y", "", colnames(final))
        
      }
      
      else final = found_data
      
    }
    
    end = Sys.time()
    
    duration = difftime(end, start)
    
    cat(paste("Time taken:", round(duration[[1]], 2), units(duration), "\n"))
    
    cat("------------------------------------\n")
    
    cat("\n")
    
    if(progressive_saving) write.csv(final, paste0(csv_name, "_", i, ".csv"), row.names = F)
    
    if(i == 1) final_output = final
    
    else final_output = rbind(final_output, final)
    
  }
  
  if(final_saving) write.csv(final_output, paste0(csv_name, ".csv"), row.names = F)
  
  if(type == "bold_id") {
    
    colnames(final_output)[1] = "BOLD_ID"
    
    final_output = subset(final_output, select = -SIMILARITY)
    
  }
  
  return(tibble::tibble(final_output))
  
}
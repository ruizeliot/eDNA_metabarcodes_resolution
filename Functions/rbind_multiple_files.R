## Requirement: 'tibble' + 'DECIPHER'

rbind_multiple_files = function(fixed_name, path = NULL, verbose = T, file_extension = "csv"){
  
  if(is.null(path)) found_names = list.files()[grep(fixed_name, list.files())]
  
  else found_names = list.files(path)[grep(fixed_name, list.files(path))]
  
  if(is.null(path)) found_names = found_names[grep(paste0(".", file_extension), found_names)]
  
  else found_names = paste0(path, "/", found_names[grep(paste0(".", file_extension), found_names)])
  
  for(i in 1:length(found_names)){
    
    start = Sys.time()
    
    if(file_extension == "csv") data = read.csv(found_names[i], header = T, sep = ",")
    
    else {
      
      data = Biostrings::readDNAStringSet(found_names[i]) 
      
      data = tibble::tibble(data.frame(ACCESSION = names(data), PRIMER_SET = gsub(paste0(path, "/|", fixed_name, "|.", file_extension), "", found_names[i]),
                               LENGTH = Biostrings::width(data), SEQUENCES = BiocGenerics::paste(data)))
      
    }
    
    end = Sys.time()
    
    duration = difftime(end, start)
    
    if(verbose){
      
      cat(paste0("Reading the csv for ", found_names[i], " DONE\n"))
      
      cat(paste("Time taken:", round(duration[[1]], 2), units(duration), "\n"))
      
      cat("------------------------------------\n")
      
      cat("\n")
      
    }
    
    if(i == 1) final_output = data
    
    else final_output = rbind(final_output, data)
    
  }
  
  return(tibble::tibble(final_output))
  
}

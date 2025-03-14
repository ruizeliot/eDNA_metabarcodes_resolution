## Replacement: 'tibble'

# Info: Gather element of the ID column and check if different values are present for the same ID. 
# If so, it removes the NA and combine the different values by separating (e.g. "x or y or z").

consensus_deduplication = function(data, id){
  
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
  
  if(length(unique(new_data[[id]])) != length(new_data[[id]])) new_data = suppressWarnings(consensus_deduplication(new_data, id))
  
  old_data_deduplicated = suppressWarnings(consensus_deduplication(old_data, id))
  
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



## Replacement: 'worms' + 'tibble' + 'taxizedb' + 'stringr'

# Info: Automatically detect the columns with missing info. 
# Info: Based on references columns, it searches for all superior taxa in GBIF, ITIS, Catalog of Life and World Flora Online local databases.
# Info: Conflictual taxa are both conserved and separated with an "or", to avoid splitting a taxa in multiple synonyms.

# Note: BASE_TAXA can be set to automatic for dataframes obtained after treatment with the other functions (i.e. same column names).
# Note: Otherwise, it must be a vector with the reference the name of the reference column for each column with NA.

complete_taxonomy = function(data, BASE_TAXA = "automatic"){
  
  wormsbynames_modified = function(taxon_names, ids = FALSE, match = FALSE, verbose = TRUE, 
            chunksize = 50, like = "false", marine_only = "true", 
            sleep_btw_chunks_in_sec = 0.1) 
  {
    stopifnot(inherits(taxon_names, "character"))
    if (match) {
      ids = TRUE
    }
    search_options = paste0("like=", like, "&marine_only=", 
                             marine_only)
    my_worms = list()
    request = "http://www.marinespecies.org/rest/AphiaRecordsByNames"
    wrapname = gsub(" ", "%20", taxon_names)
    chunk = split(wrapname, ceiling(seq_along(taxon_names)/chunksize))
    chunkid = split(1:length(taxon_names), ceiling(seq_along(taxon_names)/chunksize))
    if(verbose) cat("REQUESTING ", length(taxon_names), " ITEMS BY NAME from World Register of Marine Species (www.marinespecies.org), ", 
        format(Sys.time(), "%d/%m/%Y %X"), " (CC-BY)\n", 
        sep = "")
    for (round in 1:length(chunk)) {
      if (verbose) {
        cat(sprintf("%62s", paste0("chunk ", 
                                   round, "/", length(chunk))), "\n")
      }
      m = paste(paste(paste0("scientificnames[]=", chunk[[round]]), 
                       collapse = "&"), search_options, sep = "&")
      r = httr::GET(paste(request, m, sep = "?"))
      Sys.sleep(sleep_btw_chunks_in_sec)
      stopifnot(r$status_code == 200)
      r_parsed = httr::content(r, as = "parsed")
      for (i in 1:length(r_parsed)) {
        w_index = unlist(chunkid[round])[i]
        if (length(r_parsed[[i]]) == 0) {
          my_worms[[w_index]] = NA
          if (verbose) cat(sprintf("%-46s       %-40s", taxon_names[w_index], 
                      "no match"), "\n")
        }
        else if (length(r_parsed[[i]]) == 1) {
          my_worms[[w_index]] = r_parsed[[i]][[1]]
        }
        else {
          l = length(r_parsed[[i]])
          my_worms[[w_index]] = r_parsed[[i]][[l]]
          for (j in 1:l) {
            if (r_parsed[[i]][[j]]$status == "accepted") {
              my_worms[[w_index]] = r_parsed[[i]][[j]]
            }
          }
        }
      }
    }
    non.null.list = lapply(my_worms, lapply, function(x) ifelse(is.null(x), 
                                                                 NA, x))
    worms = plyr::rbind.fill(lapply(non.null.list, as.data.frame, 
                               stringsAsFactors = F))
    worms$NA. = NULL
    if (ids) {
      worms = cbind(data.frame(id = 1:nrow(worms), name = taxon_names, 
                                stringsAsFactors = F), worms)
    }
    if (verbose) {
      cat("by names ........................................... DONE\n")
    }
    if (match) {
      nonefound = is.na(worms[, "AphiaID"])
      failed_species = taxon_names[nonefound]
      if (length(failed_species) > 0) {
        failed_worms = wormsbymatchnames(failed_species, 
                                          verbose = verbose, ids = FALSE)
        worms[nonefound, c(F, F, rep(T, ncol(failed_worms)))] = failed_worms
      }
      else {
        if (verbose) cat("  Nothing to match.\n")
      }
    }
    return(worms)
  }
  
  fct_env = new.env()
  
  na_columns = colSums(is.na(data)) 
  
  names_na_columns = names(na_columns[na_columns > 0])
  
  if(anyNA(names_na_columns)) stop("One of the column is named NA.")
  
  if(BASE_TAXA != "automatic") warning("You must provide a vector with the name of the reference column for each column with NA.")
  
  else{
    
    if(ncol(random_infos_data) == 10)
      if(!all(colnames(random_infos_data) == c("ID", "SPECIES", "RANK", "GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM")))
        stop('You must provide a dataframe with the following column names: "ID", "SPECIES", "RANK" (facultative column), "GENUS" (facultative column), "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM" (in capital) and they must be in this order.')
    
    if(ncol(random_infos_data) == 9)
      if(!all(colnames(random_infos_data) == c("ID", "SPECIES", "GENUS", "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM")) ||
         !all(colnames(random_infos_data) == c("ID", "SPECIES", "RANK", "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM")))
         stop('You must provide a dataframe with the following column names: "ID", "SPECIES", "RANK" (facultative column), "GENUS" (facultative column), "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM" (in capital) and they must be in this order.')
    
    else if(ncol(random_infos_data) == 8) 
      if(!all(colnames(random_infos_data) == c("ID", "SPECIES", "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM")))
      stop('You must provide a dataframe with the following column names: "ID", "SPECIES", "RANK" (facultative column), "GENUS" (facultative column), "FAMILY", "ORDER", "CLASS", "PHYLUM", "KINGDOM", "SUPERKINGDOM" (in capital) and they must be in this order.')
                                                   
  }
  
  database_taxonomy = function(initial_data, database, column_to_correct, BASE_TAXA = "automatic") {
    
    if(BASE_TAXA == "automatic") {
      
      if(column_to_correct == "FAMILY") BASE_TAXA = "SPECIES"
      
      else if(column_to_correct == "ORDER") BASE_TAXA = "FAMILY"
      
      else if(column_to_correct == "CLASS") BASE_TAXA = "ORDER"
      
      else if(column_to_correct == "PHYLUM") BASE_TAXA = "CLASS"
      
      else if(column_to_correct == "KINGDOM") BASE_TAXA = "PHYLUM"
      
      else if(column_to_correct == "SUPERKINGDOM") BASE_TAXA = "KINGDOM" 
      
    }
    
    missing_data = initial_data[which(!complete.cases(initial_data[column_to_correct])),]
    
    if(nrow(missing_data) != 0){
      
      if(database == "worms"){
        
        worms_valid_id = tryCatch(wormsbynames_modified(as.character(missing_data[[BASE_TAXA]]), verbose = F), 
                                  error = function(e) NULL)
        
        worms_valid_id = worms_valid_id$valid_AphiaID[!is.na(as.numeric(worms_valid_id$valid_AphiaID))]
        
        if(length(worms_valid_id) == 0) data_corrected = initial_data
        
        else{
          
          worms_new_infos = worms::wormsbyid(worms_valid_id, verbose = F)
          
          worms_new_infos = data.frame(SPECIES = worms_new_infos$scientificname, GENUS = worms_new_infos$genus, 
                                       FAMILY = worms_new_infos$family, ORDER = worms_new_infos$order, CLASS = worms_new_infos$class, 
                                       PHYLUM = worms_new_infos$phylum, KINGDOM = worms_new_infos$kingdom)
          
          data_corrected = add_infos(initial_data, worms_new_infos, id = BASE_TAXA)
          
        }
        
      }
      
      else {
        
        id_data = taxizedb::name2taxid(missing_data[[BASE_TAXA]][which(!is.na(missing_data[[BASE_TAXA]]))], db = database, out_type = "summary")
        
        id_data = id_data[!duplicated(id_data$name),]
        
        if(length(id_data$id) != 0) {
          
          new_classification = id_to_classification(unique(id_data$id), db = database)
          
          if(BASE_TAXA == "SPECIES" && database == "gbif") { new_classification["GENUS"] = stringr::word(new_classification$GENUS, 1) }
          
          data_corrected = add_infos(initial_data, new_classification, id = BASE_TAXA)
          
          
        }
        
        else { data_corrected = initial_data }
        
      }
      
    }
    
    else { data_corrected = initial_data }
    
    data_corrected
    
  }
  
  results = function(data_corrected, message, ok = F){
    
    na_col = colSums(is.na(data_corrected)) 
    
    na_col = na_col[na_col > 0]
    
    cat("\n")
    
    cat(paste(message, "refinement:\n"))
    
    cat("\n")
    
    if(ok) {
      
      end_vector = rep(0, length(names_na_columns))
      
      names(end_vector) = names_na_columns
      
      print(end_vector)
      
    }
    
    else print(na_col)
    
    cat("\n")
    
  }
  
  show_differences = function(data_original, data_new, data_before, database_name){
    
    before_na = sum(sapply(X = data_before, FUN = function(x) sum(is.na(x))))
    
    after_na = sum(sapply(X = data_new, FUN = function(x) sum(is.na(x))))
    
    if(after_na == 0){
      
      cat("\n")
      
      results(data_new, paste(database_name, names_na_columns[i]), ok = T)
      
      cat("All taxonomic informations were retrieved!\n")
      
      cat("\n")
      
    }
    
    else if(before_na == after_na) {
      
      cat("\n")
      
      cat(paste(database_name, names_na_columns[i], "refinement: No information added\n"))
      
      cat("\n")
      
    }
    
    else {
      
      results(data_new, paste(database_name, names_na_columns[i]))
      
    }
    
  }
  
  results(data, "Before")
  
  for(i in 1:length(names_na_columns)){
    
    if(i == 1){
      
      assign(paste("gbif", 1, sep = "_"), data, envir = fct_env)
      
      assign(paste("gbif", (i + 1), sep = "_"), database_taxonomy(get(paste("gbif", (i), sep = "_"), envir = fct_env), 
                                                                  "gbif", column_to_correct = names_na_columns[i]), envir = fct_env)
      
      show_differences(data, get(paste("gbif", (i + 1), sep="_"), envir = fct_env), 
                       get(paste("gbif", (i), sep = "_"), envir = fct_env), database_name = "'GBIF'")
      
    }
    
    else if(i > 1) { 
      
      assign(paste("gbif", (i + 1), sep = "_"), database_taxonomy(get(paste("Tour", (i), sep = "_"), envir = fct_env), 
                                                                  "gbif", column_to_correct = names_na_columns[i]), envir = fct_env) 
      
      show_differences(data, get(paste("gbif", (i + 1), sep="_"), envir = fct_env), 
                       get(paste("Tour", (i), sep = "_"), envir = fct_env), database_name = "'GBIF'")
      
    }
    
    if(anyNA(get(paste("gbif", (i + 1), sep = "_"), envir = fct_env)[[names_na_columns[i]]])) {
      
      assign(paste("itis", (i + 1), sep = "_"), database_taxonomy(get(paste("gbif", (i + 1), sep = "_"), envir = fct_env), 
                                                                  "itis", column_to_correct = names_na_columns[i]), envir = fct_env)
      
      show_differences(data, get(paste("itis", (i + 1), sep="_"), envir = fct_env), 
                       get(paste("gbif", (i + 1), sep = "_"), envir = fct_env), database_name = "'ITIS'")
      
      if(anyNA(get(paste("itis", (i + 1), sep = "_"), envir = fct_env)[[names_na_columns[i]]])) {
        
        assign(paste("col", (i + 1), sep = "_"), database_taxonomy(get(paste("itis", (i + 1), sep = "_"), envir = fct_env), 
                                                                   "col", column_to_correct = names_na_columns[i]), envir = fct_env)
        
        show_differences(data, get(paste("col", (i + 1), sep="_"), envir = fct_env), 
                         get(paste("itis", (i + 1), sep = "_"), envir = fct_env), database_name = "'Catalog of Life'")
        
        if(anyNA(get(paste("col", (i + 1), sep = "_"), envir = fct_env)[[names_na_columns[i]]])) {
          
          assign(paste("wfo", (i + 1), sep = "_"), database_taxonomy(get(paste("col", (i + 1), sep = "_"), envir = fct_env), 
                                                                     "wfo", column_to_correct = names_na_columns[i]), envir = fct_env)
          
          show_differences(data, get(paste("wfo", (i + 1), sep="_"), envir = fct_env), 
                           get(paste("col", (i + 1), sep = "_"), envir = fct_env), database_name = "'World Flora Online'")
          
          if(anyNA(get(paste("wfo", (i + 1), sep = "_"), envir = fct_env)[[names_na_columns[i]]])) {
            
            assign(paste("worms", (i + 1), sep = "_"), database_taxonomy(get(paste("wfo", (i + 1), sep = "_"), envir = fct_env), 
                                                                         "worms", column_to_correct = names_na_columns[i]), envir = fct_env)
            
            show_differences(data, get(paste("worms", (i + 1), sep="_"), envir = fct_env), 
                             get(paste("wfo", (i + 1), sep = "_"), envir = fct_env), database_name = "'WoRMS'")
            
            {if(anyNA(get(paste("worms", (i + 1), sep = "_"), envir = fct_env)[[names_na_columns[i]]])) {
              
              assign(paste("Tour", (i + 1), sep = "_"), get(paste("worms", (i + 1), sep = "_"), envir = fct_env), envir = fct_env)
              
              cat("\n")
              
              cat(paste("Not all", names_na_columns[i], "informations were retrieved despite searching in 5 databases.\n"))
              
              cat("\n")
              
            }
              
              else assign(paste("Tour", (i + 1), sep = "_"), get(paste("worms", (i + 1), sep = "_"), envir = fct_env), envir = fct_env)}
            
          }
          
          else assign(paste("Tour", (i + 1), sep = "_"), get(paste("wfo", (i + 1), sep = "_"), envir = fct_env), envir = fct_env)
          
        }
        
        else assign(paste("Tour", (i + 1), sep = "_"), get(paste("col", (i + 1), sep = "_"), envir = fct_env), envir = fct_env)
        
      }
      
      else assign(paste("Tour", (i + 1), sep = "_"), get(paste("itis", (i + 1), sep = "_"), envir = fct_env), envir = fct_env)
      
    }
    
    else assign(paste("Tour", (i + 1), sep = "_"), get(paste("gbif", (i + 1), sep = "_"), envir = fct_env), envir = fct_env)
    
  }
  
  final_data = suppressWarnings(get(paste("Tour", (i + 1), sep = "_"), envir = fct_env))
  
  final_data[final_data == "Not assigned"] = NA
  
  final_data[final_data == "incertae sedis"] = NA
  
  if(sum(sapply(X = final_data, FUN = function(x) sum(is.na(x)))) != 0){
    
    results(data, "Before")
    
    results(final_data, "After last")
    
  }
  
  tibble::tibble(final_data)
  
}

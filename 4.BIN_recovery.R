##### INITIALISATION #####

## Loading required packages and custom functions

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(bold)
library(tidyverse)


## Loading fish mitogenomes and associated taxonomic informations

act_mitogenomes = tibble(read.csv("Data/mitogenomes_taxonomy/Actinopterygii - Mitogenomes taxonomy.csv"))
act_mitogenomes 


## Loading extracted metabarcodes

filtered_metabarcodes_corrected = tibble(read.csv("Data/metabarcodes/Metabarcodes - N correction - Length filtering.csv"))
filtered_metabarcodes_corrected


## Loading the primer data informations

primer_data_no_ps1 = tibble(read.csv("Data/primers_infos/Primer data - All infos - No PS1.csv"))
primer_data_no_ps1





##### FIRST SEARCH OF BIN #####

## BLAST of the Genbank accession numbers in BOLD (which automatically retrieves Genbank barcodes)

actinopterygii_accession = unique(filtered_metabarcodes_corrected$ACCESSION)

start1 = Sys.time()
bin_first_search = accession_2_bin(actinopterygii_accession, division_number = 500, csv_name = "Data/BIN/BIN - First search")
end1 = Sys.time()
difftime(end1, start1) # Duration : 20mn

bin_first_search = tibble(read.csv("Data/BIN/BIN - First search.csv"))
bin_first_search


## Checking the number of BIN retrieved with this rapid method

missing_bin_first_search = subset(bin_first_search, is.na(BIN))
missing_bin_first_search 

nrow(missing_bin_first_search) / length(actinopterygii_accession) * 100 # 26% still missing





##### SECOND SEARCH OF BIN (BEST MATCH ID) #####

## BLAST of the COI barcodes in BOLD -> checking the best match ID

missing_sequences_first_search = subset(filtered_metabarcodes_corrected, ACCESSION %in% missing_bin_first_search$ACCESSION & 
                                          PRIMERS == "FishF1-FishR1")
missing_sequences_first_search

start2 = Sys.time()
bold_identification = bold_identify(paste(missing_sequences_first_search$FRAGMENT), db = "COX1_SPECIES")
end2 = Sys.time()
difftime(end2, start2)    # Duration: 7h30
names(bold_identification) = missing_sequences_first_search$ACCESSION # Naming list elements with Genbank accession numbers


## Getting position of missing and not missing elements in list

missing_pos = which(unlist(lapply(bold_identification, is.null)))
{if(length(missing_pos) != 0) not_missing_pos = setdiff(seq(1:length(bold_identification)), missing_pos)}
{if(length(missing_pos) == 0) not_missing_pos = seq(1:length(bold_identification))}


## Extracting found elements

found_table = bold_identification[not_missing_pos]
found_id = sapply(found_table, function(x) x[1,]$ID)


## Searching in BOLD for informations on the best match ID 

start3 = Sys.time()
second_search = accession_2_bin(accession = found_id, type = "bold_id",
                                division_number = ifelse(length(found_id) < 501, length(found_id), 500), final_saving = F)
end3 = Sys.time()
difftime(end3, start3)    # Duration: 50s
second_search 


## Extracting similarity and adding info the second_search table

similarity = unlist(lapply(found_table, function(x) x[1,6]))
similarity_table = tibble(data.frame(ACCESSION = names(similarity), BOLD_ID = found_id, SIMILARITY = as.numeric(similarity)))
similarity_table

second_search = tibble(merge(second_search, similarity_table))
second_search = second_search[!duplicated(second_search),]
second_search


## Constructing a table containing all accession initially searched, found or not

second_search_final = tibble(cbind(second_search[,2], second_search[,1], second_search[,4:9], second_search[,3]))
second_search_final

if(length(missing_pos) != 0){
  
  missing_df = data.frame(ACCESSION = names(bold_identification)[missing_pos], BOLD_ID = NA, BIN = NA, 
                          BOLD_TAXID = NA, BOLD_SPECIES = NA, BOLD_GENUS = NA, BOLD_FAMILY = NA, 
                          BOLD_ORDER = NA, SIMILARITY = NA)
  
  second_search_final = tibble(rbind(second_search_final, missing_df))
  
}

second_search_final


## Saving the result

write.csv(second_search_final, "Data/BIN/BIN - Second search.csv", row.names = F)
second_search_final = tibble(read.csv("Data/BIN/BIN - Second search.csv"))
second_search_final





##### THIRD SEARCH OF BIN (BEST MATCHED BIN) #####

## Getting sequences still missing a BIN and building a list of all found bold identifier for each one of them

not_first_pos = which(is.na(second_search_final$BIN) & !is.na(second_search_final$BOLD_ID))

all_found_id = lapply(bold_identification[setdiff(c(1:length(bold_identification)), missing_pos)], function(x) x[,1])

still_missing_bin = all_found_id[not_first_pos]


## Searching informations for the first BIN found, and not the first row (best match ID does not always have a BIN)

start4 = Sys.time()
third_search = accession_2_bin(accession = still_missing_bin, type = "bold_id", final_saving = F, multiple_accession_list = T)
end4 = Sys.time()
difftime(end4, start4) # Duration: 12.5mn
third_search


## Retrieving all similarities found for other match than the first

all_similarity = lapply(bold_identification[not_first_pos], function(x) x[,6])

no_id_found = which(unlist(lapply(all_similarity, is.null)))
no_id_found

{if(length(no_id_found) != 0) id_found = setdiff(seq(1:length(all_similarity)), no_id_found)}
{if(length(no_id_found) == 0) id_found = seq(1:length(all_similarity))}

all_similarity_found = all_similarity[id_found] 
all_similarity_found 

corresponding_id = all_found_id[which(names(all_found_id) %in% names(all_similarity_found))]
corresponding_id


## Retrieving the best matched BIN

position_first_bin = sapply(corresponding_id, function(x) which(x %in% third_search$BOLD_ID)[1]) # 1st best matched BIN
position_first_bin

accession_instead_id = which(is.na(position_first_bin))
accession_instead_id 

all_similarity_found_corrected = all_similarity_found[-accession_instead_id]
corresponding_id_corrected = corresponding_id[-accession_instead_id]
position_first_bin_corrected = position_first_bin[-accession_instead_id]

similarity_list_corrected = list()
corresponding_id_list = list()

for(i in 1:length(all_similarity_found_corrected)) { 
  
  similarity_list_corrected[i] = all_similarity_found_corrected[[i]][position_first_bin_corrected[i]] 
  
  corresponding_id_list[i] = corresponding_id_corrected[[i]][position_first_bin_corrected[i]]
  
}


## Compiling results in a dataframe

similarity_third_search = tibble(data.frame(ACCESSION = names(all_similarity_found_corrected), BOLD_ID = unlist(corresponding_id_list),
                                            SIMILARITY = unlist(similarity_list_corrected)))
similarity_third_search

third_search_final = tibble(merge(third_search, similarity_third_search))
third_search_final = tibble(cbind(third_search_final[,8], third_search_final[,c(1:7,9)]))
third_search_final = third_search_final[!duplicated(third_search_final),]
third_search_final

third_search_bin_found = subset(bin_second_search_added, ACCESSION %in% third_search_final$ACCESSION & is.na(BIN))
third_search_final = subset(third_search_final, ACCESSION %in% third_search_bin_found$ACCESSION)


## Saving the results

write.csv(third_search_final, "Data/BIN/BIN - Third search.csv", row.names = F)
third_search_final = tibble(read.csv("Data/BIN/BIN - Third search.csv"))
third_search_final





##### CREATING DATABASES DEPENDING ON BIN #####

## Adding progressively new informations found during each searches

bin_second_search_added = replace_values(data_original = bin_first_search, data_model = second_search_final[,-2], 
                                         variables_original = colnames(bin_first_search)[2:8], 
                                         id_original = "ACCESSION", total_replacement = "all")

bin_database = replace_values(data_original = bin_second_search_added, data_model = third_search_final[,-2], 
                              variables_original = colnames(bin_second_search_added)[2:8], 
                              id_original = "ACCESSION", total_replacement = "all")


## Saving the full dataframe

write.csv(bin_database, "Data/BIN/BIN database - Full.csv", row.names = F)
bin_database = tibble(read.csv("Data/BIN/BIN database - Full.csv"))
bin_database


## Removing all unsure/unknown identifications and saving

bin_database_corrected = subset(bin_database, !is.na(BIN) & SIMILARITY == 1)

write.csv(bin_database_corrected, "Data/BIN/BIN database - Corrected.csv", row.names = F)
bin_database_corrected = tibble(read.csv("Data/BIN/BIN database - Corrected.csv"))
bin_database_corrected


## Some stats on the data

length(unique(bin_database_corrected$BIN)) ; length(unique(bin_database_corrected$ACCESSION)) # 2669 verified BIN & 5438 sequences
length(unique(bin_database_corrected$BIN)) / length(unique(bin_database$BIN)) * 100 # 95% of all BIN found
length(unique(bin_database_corrected$ACCESSION)) / length(unique(bin_database$ACCESSION)) * 100 # 89% of all fish sequences retained


## Getting only metabarcodes associated with a verified BIN

filtered_metabarcodes_final = tibble(merge(bin_database_corrected[,1:2], filtered_metabarcodes_corrected))
write.csv(filtered_metabarcodes_final, "Data/metabarcodes/Filtered metabarcodes final.csv", row.names = F)
filtered_metabarcodes_final = tibble(read.csv("Data/metabarcodes/Filtered metabarcodes final.csv"))
filtered_metabarcodes_final


## Associated taxonomic informations and BIN for these sequences

act_mitogenomes_bin_corrected = tibble(merge(bin_database_corrected[,1:2], act_mitogenomes[,-13]))
write.csv(act_mitogenomes_bin_corrected, "Data/taxonomy_BIN_infos/Actinopterygii - Taxonomy - BIN.csv", row.names = F)
act_mitogenomes_bin_corrected = tibble(read.csv("Data/taxonomy_BIN_infos/Actinopterygii - Taxonomy - BIN.csv"))
act_mitogenomes_bin_corrected


## Constructing a reference database for the mean length of each metabarcodes

mean_length_metabarcodes = filtered_metabarcodes_final %>% group_by(PRIMERS) %>% 
  summarize(MEAN_LENGTH = mean(LENGTH), MEAN_LENGTH_PRIMERS = mean(LENGTH_WITH_PRIMERS))

mean_length_metabarcodes = tibble(merge(subset(primer_data_no_ps1, select = c("SHORT_NAME", "GENE")), 
                                        mean_length_metabarcodes, by.x = "SHORT_NAME", by.y = "PRIMERS"))
colnames(mean_length_metabarcodes)[1] = "PRIMERS"
mean_length_metabarcodes

write.csv(mean_length_metabarcodes, "Data/primers_infos/Mean length data.csv", row.names = F)
mean_length_metabarcodes = tibble(read.csv("Data/primers_infos/Mean length data.csv"))
mean_length_metabarcodes %>% print(n = 22)



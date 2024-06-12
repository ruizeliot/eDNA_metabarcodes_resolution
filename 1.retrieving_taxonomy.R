##### INITIALISATION #####

## Installing the functions designed for this study directly from Github 

#library(remotes) ; install_github("Eliot-RUIZ/eDNAevaluation", upgrade = "never")
#library(eDNAevaluation) -> NOT POSSIBLE NOW BECAUSE PACKAGE IS PRIVATE


## In case of error, even after updating R (see package "installr"), please install the function locally and load dependencies manually

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(taxizedb)


## Download of the databases -> may take 10mn each

#db_download_ncbi()
#db_download_itis()
#db_download_gbif()
#db_download_wfo()
#db_download_col() 





##### DATA RETRIEVING #####

## Data donwload

# Note: Data downloaded in FASTA and sorted by TAXID from the NCBI nucleotide database by typing 
# "actinopterygii complete genome" and clicking on "mitochondrion" in the left panel.


## Data reading

actinopterygii = readDNAStringSet("Data/NCBI/Actinopterygii - Complete mitochondrial sequences.fasta") 
actinopterygii


## Accession number extraction

actinopterygii_accession = word(names(actinopterygii), 1)
head(actinopterygii_accession)


## Extraction of the names and sequences in the fasta file

actinopterygii_sequences_names = names(actinopterygii)                                            # Extraction of sequences names
actinopterygii_sequences = paste(actinopterygii)                                                  # Extraction of sequences
actinopterygii_sequences_table = tibble(data.frame(ACCESSION = actinopterygii_accession, 
                                                   COMPLETE_NAMES = actinopterygii_sequences_names, 
                                                   SEQUENCES = actinopterygii_sequences))
actinopterygii_sequences_table





##### TAXON NAME FILTERING #####

## Counting and removing of the mention UNVERIFIED

length(actinopterygii_sequences_names[grep("unverified|UNVERIFIED", actinopterygii_sequences_names)]) # 76 Unverified IDs

actinopterygii_species_not_unverified = gsub(paste0(c("UNVERIFIED: ", "UNVERIFIED ", "unverified: ", "unverified "), 
                                                    collapse = "|"),"", actinopterygii_sequences_names)                               

actinopterygii_species = tibble(data.frame(SPECIES = word(actinopterygii_species_not_unverified, 2, 3), 
                                           COMPLETE_NAMES = actinopterygii_sequences_names))
actinopterygii_species


## Removing of hydrids

nrow(actinopterygii_species[grep(" x ", actinopterygii_species$COMPLETE_NAMES), ])
actinopterygii_species[grep(" x ", actinopterygii_species$COMPLETE_NAMES), ] = NA
nrow(actinopterygii_species[grep(" hybrid", actinopterygii_species$COMPLETE_NAMES), ])
actinopterygii_species[grep(" hybrid", actinopterygii_species$COMPLETE_NAMES), ] = NA                 # 188 hybrids


## Removing uncertain identifications

nrow(actinopterygii_species[grep(" cf.", actinopterygii_species$COMPLETE_NAMES), ])                   # 17 uncertain IDs
actinopterygii_species[grep(" cf.", actinopterygii_species$COMPLETE_NAMES), ] = NA


## Removing subspecies and varieties

nrow(actinopterygii_species[grep(" var.", actinopterygii_species$COMPLETE_NAMES), ])                  # 30 varieties
actinopterygii_species[grep(" var.", actinopterygii_species$COMPLETE_NAMES), ] = NA


## Removing species with affinities to a genus

nrow(actinopterygii_species[grep(" aff.", actinopterygii_species$COMPLETE_NAMES), ])                  # 16 species with affinities
actinopterygii_species[grep(" aff.", actinopterygii_species$COMPLETE_NAMES), ] = NA


## Removing of unidentified species

nrow(actinopterygii_species[grep(" sp.", actinopterygii_species$COMPLETE_NAMES), ])                   # 336 unidentified sequences
actinopterygii_species[grep(" sp.", actinopterygii_species$COMPLETE_NAMES), ] = NA


## Removing lines with NA only

actinopterygii_species_corrected = actinopterygii_species[!apply(actinopterygii_species, 1, function(X) all(is.na(X))), ]
actinopterygii_species_corrected                                                                      # 10,934 sequences from specimens


## Taxid number retrievement from NCBI local database (might not be the same depending on NCBI database updating)

actinopterygii_taxid_table = name2taxid(actinopterygii_species_corrected$SPECIES, out_type = "summary")
actinopterygii_taxid_table

actinopterygii_taxid_merged = tibble(merge(actinopterygii_taxid_table, actinopterygii_species_corrected, 
                                              by.x = "name", by.y = "SPECIES", all = T))
colnames(actinopterygii_taxid_merged)[1:2] = c("SPECIES", "TAXID")
actinopterygii_taxid_merged = actinopterygii_taxid_merged[,c(2,1,3)]
actinopterygii_taxid_merged$TAXID = as.numeric(actinopterygii_taxid_merged$TAXID)

nrow(actinopterygii_taxid_merged) == nrow(actinopterygii_species_corrected) # Not the same number of rows

setdiff(actinopterygii_species_corrected$COMPLETE_NAMES, actinopterygii_taxid_corrected$COMPLETE_NAMES)
setdiff(actinopterygii_taxid_corrected$COMPLETE_NAMES, actinopterygii_species_corrected$COMPLETE_NAMES)

actinopterygii_taxid_merged_duplicated =
  actinopterygii_taxid_merged[which(duplicated(actinopterygii_taxid_merged$COMPLETE_NAMES) | 
                                      duplicated(actinopterygii_taxid_merged$COMPLETE_NAMES, fromLast = T)),]
actinopterygii_taxid_merged_duplicated # Ariomma indicum with two TAXID duplicated
id_to_classification(as.numeric(actinopterygii_taxid_merged_duplicated$TAXID)) # Only Ariomma indica is accepted

actinopterygii_taxid = actinopterygii_taxid_merged[-which(duplicated(
  actinopterygii_taxid_merged$COMPLETE_NAMES)), ] # Removing the first row containing Ariomna indicum

nrow(actinopterygii_taxid) == nrow(actinopterygii_species_corrected) # Same number of rows


## Removing supplementary characters from names (might not be the same depending on NCBI database updating)

actinopterygii_taxid[!complete.cases(actinopterygii_taxid), ] # Taxid not found for 29 specimens
actinopterygii_taxid$SPECIES = replace(actinopterygii_taxid$SPECIES, grep("TPA_asm", actinopterygii_taxid$SPECIES), 
                                       c("Astyanax mexicanus", "Astyanax aeneus", "Psalidodon fasciatus"))
actinopterygii_taxid$SPECIES = gsub(",", "", actinopterygii_taxid$SPECIES)

actinopterygii_taxid_table2 = name2taxid(actinopterygii_taxid$SPECIES, out_type = "summary")
actinopterygii_taxid_table2

actinopterygii_taxid_corrected = tibble(merge(actinopterygii_taxid_table2, actinopterygii_taxid [,-1], 
                                           by.x = "name", by.y = "SPECIES", all = T))
colnames(actinopterygii_taxid_corrected)[1:2] = c("SPECIES", "TAXID")
actinopterygii_taxid_corrected = actinopterygii_taxid_corrected[,c(2,1,3)]
actinopterygii_taxid_corrected$TAXID = as.numeric(actinopterygii_taxid_corrected$TAXID)
actinopterygii_taxid_corrected

nrow(actinopterygii_taxid_corrected) == nrow(actinopterygii_taxid) # Not the same number of rows

actinopterygii_taxid_corrected[which(duplicated(
  actinopterygii_taxid_corrected$COMPLETE_NAMES)), ] # Verifying that it is still Ariomma indicum causing problems

actinopterygii_taxid_corrected = actinopterygii_taxid_corrected[-which(duplicated(
  actinopterygii_taxid_corrected$COMPLETE_NAMES)), ]
actinopterygii_taxid_corrected

nrow(actinopterygii_taxid_corrected) == nrow(actinopterygii_taxid) # Same number of rows


##  Correcting manually names mispelling from the NCBI browser online (might not be the same depending on NCBI database updating)

actinopterygii_taxid_corrected[!complete.cases(actinopterygii_taxid_corrected), ]$SPECIES 
# Searching manually each names on NCBI to get name without mispelling and TAXID

accepted_species_names_NCBI = c(rep("Candidia barbata", 3), "Chaetodon lineolatus", rep("Epinephelus quoyans", 2), 
                                "Ecsenius bicolor", rep("Vanmanenia hainanensis", 3), rep("Labracoglossa argenteiventris", 2),
                                rep("Pseudexostoma yunnanense", 3), "Sebastes schlegelii", rep("Bangana rendahli", 2), 
                                rep("Sinocrossocheilus labiatus", 2), rep("Sturisomatichthys panamensis", 2), 
                                rep("Torquigener brevipinnis", 2), "Trachipterus trachypterus")

accepted_species_names_NCBI_taxid = name2taxid(accepted_species_names_NCBI)

actinopterygii_corrected_names = data.frame(INITIAL_SPECIES = actinopterygii_taxid_corrected
                                            [!complete.cases(actinopterygii_taxid_corrected), ]$SPECIES,
                                            SPECIES = accepted_species_names_NCBI,
                                            TAXID = accepted_species_names_NCBI_taxid)
actinopterygii_corrected_names

actinopterygii_taxid_corrected = replace_values(cbind(actinopterygii_taxid_corrected, INITIAL_SPECIES = actinopterygii_taxid_corrected$SPECIES), 
                                                actinopterygii_corrected_names, variables_original = c("TAXID", "SPECIES"), 
                                                id_original = "INITIAL_SPECIES", total_replacement = "all")
actinopterygii_taxid_corrected[!complete.cases(actinopterygii_taxid_corrected), ]             # All species with TAXID
actinopterygii_taxid_corrected                                                                # 10,598 sequences belonging to species





##### RETRIEVING AND CORRECTING TAXONOMY #####

## Species name and classification retrievement from the local NCBI database using the taxid number

actinopterygii_taxid_unique = actinopterygii_taxid_corrected[!duplicated(actinopterygii_taxid_corrected$TAXID), ]
actinopterygii_taxid_unique                                                                   # 3,384 different species

actinopterygii_taxonomy = id_to_classification(actinopterygii_taxid_unique$TAXID)             # Duration : 40 seconds
actinopterygii_taxonomy 


## Correcting NCBI species names

# Removing subspecies
subset(actinopterygii_taxonomy, RANK != "species")
actinopterygii_taxonomy[["SPECIES"]] = word(actinopterygii_taxonomy[["SPECIES"]], 1, 2) # Taking only the first two words

# Last checkings (all must be empty)
actinopterygii_taxonomy[which(duplicated(actinopterygii_taxonomy[["SPECIES"]]) | 
                                duplicated(actinopterygii_taxonomy[["SPECIES"]], fromLast = TRUE)), ]       # No ambiguous TAXID
actinopterygii_taxonomy$SPECIES[which(sapply(strsplit(actinopterygii_taxonomy$SPECIES, " "), length) != 2)] # Only 2 first words
actinopterygii_taxonomy[which(!grepl("[^A-Za-z]", actinopterygii_taxonomy[["SPECIES"]])),]                  # Only letters
actinopterygii_taxonomy[which(grepl("hybrid", actinopterygii_taxonomy[["SPECIES"]])),]                      # No hybrids
actinopterygii_taxonomy[which(grepl(" sp.", actinopterygii_taxonomy[["SPECIES"]])),]                        # No unidentified species


## Automatic completion of taxonomic informations

# Stats on missing informations
na_column = colSums(is.na(actinopterygii_taxonomy)) ; na_column[na_column > 0]                      # Number of NA per column
na_row = sum(!complete.cases(actinopterygii_taxonomy)) ; na_row                                     # Number of lines with NA
na_row / nrow(actinopterygii_taxonomy) * 100                                                        # 5% of sequences with missing taxonomical infos

# Automatic completion of taxonomic informations from multiple databases 
actinopterygii_taxonomy_completed = complete_taxonomy(actinopterygii_taxonomy)                      # Duration: 2s for 3500 species
actinopterygii_taxonomy_completed

# Correcting conflicts manually (might not be the same depending on taxonomy database updating)
lengths(regmatches(actinopterygii_taxonomy_completed, gregexpr(" or ", actinopterygii_taxonomy_completed)))   # 3383 conflicts
actinopterygii_taxonomy_completed[with(actinopterygii_taxonomy_completed, grepl(" or ", paste(CLASS, KINGDOM))), ]     

actinopterygii_taxonomy_completed$CLASS = gsub("Actinopteri or Actinopterygii|Actinopterygii or Actinopteri|Actinopteri|Actinopterygii or Teleostei|Teleostei or Actinopterygii", 
                                               "Actinopterygii", actinopterygii_taxonomy_completed$CLASS)

actinopterygii_taxonomy_completed$KINGDOM = gsub("Metazoa or Animalia|Animalia or Metazoa", "Metazoa", 
                                                 actinopterygii_taxonomy_completed$KINGDOM)

lengths(regmatches(actinopterygii_taxonomy_completed, gregexpr(" or ", actinopterygii_taxonomy_completed)))  # 0 conflicts

# Adding manually remaining infos (might not be the same depending on taxonomy database updating)
colSums(is.na(actinopterygii_taxonomy_completed))

still_missing_families = actinopterygii_taxonomy_completed[!complete.cases(actinopterygii_taxonomy_completed),]
still_missing_families

missing_position_families = which(!complete.cases(actinopterygii_taxonomy_completed))
missing_position_families

actinopterygii_taxonomy_completed$FAMILY = replace(actinopterygii_taxonomy_completed$FAMILY, 
                                                   missing_position_families[1], "Percalatidae")
actinopterygii_taxonomy_completed$FAMILY = replace(actinopterygii_taxonomy_completed$FAMILY, 
                                                   missing_position_families[-1], "Paedocyprididae")

# Checking if all taxonomic infos were successfully completed
colSums(is.na(actinopterygii_taxonomy_completed))                   



##### CREATING THE FULL DATABASE #####

## Generalization of species informations to all sequences

colnames(actinopterygii_taxonomy_completed)[1] = "TAXID"                                              # Changing the column name for merging
actinopterygii_merged = tibble(merge(actinopterygii_taxonomy_completed, actinopterygii_taxid_corrected, by = "TAXID"))
length(unique(actinopterygii_merged[,2][[1]])) == length(unique(actinopterygii_taxonomy_completed$SPECIES)) # Same number of species
actinopterygii_merged = actinopterygii_merged[,-11]                                                       # Removing duplicated column
colnames(actinopterygii_merged)[2] = "SPECIES"
actinopterygii_merged


## Assembly with the sequences

act_mitogenomes = tibble(merge(actinopterygii_merged, actinopterygii_sequences_table))
act_mitogenomes = act_mitogenomes[,-c(1,4,12)]
act_mitogenomes


## Removoving replicated (after reviewing) accession numbers (NC_xxxx)

act_mitogenomes = act_mitogenomes[!grepl("NC_", act_mitogenomes$ACCESSION),]


## Removing useless accession version

subset(act_mitogenomes, duplicated(word(ACCESSION, 1, sep = fixed(".")))) # No sequences with in 2 different versions
act_mitogenomes$ACCESSION = word(act_mitogenomes$ACCESSION, 1, sep = fixed("."))
act_mitogenomes


## Removing some sequences to match the dataset used in following articles (filtering on ecological data)

removed_sequences_accession = c("KF386025", "KM257863", "KU674798", "MH463445", "MH463450", "MH463451", "MH479388", 
                                "MH479389", "MH479394", "MH479396", "MH479398", "MH479399", "MH479400")
act_mitogenomes = subset(act_mitogenomes, !(ACCESSION %in% removed_sequences_accession))
act_mitogenomes


## Saving/reading the data

write.csv(act_mitogenomes, "Data/mitogenomes_taxonomy/Actinopterygii - Mitogenomes taxonomy.csv", row.names = F)
act_mitogenomes = tibble(read.csv("Data/mitogenomes_taxonomy/Actinopterygii - Mitogenomes taxonomy.csv"))
act_mitogenomes 

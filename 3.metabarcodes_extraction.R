##### INITIALISATION #####

## Installing the functions designed for this study directly from Github 

#library(remotes) ; install_github("Eliot-RUIZ/eDNAevaluation", upgrade = "never")
#library(eDNAevaluation) -> NOT POSSIBLE NOW BECAUSE PACKAGE IS PRIVATE


## In case of error, even after updating R (see package "installr"), please install the function locally and load dependencies manually

miceadds::source.all("Functions not commented yet")
library(tibble)
library(DECIPHER)
library(stringr)


## Loading the primer data informations

primer_data = tibble(read.csv("Data/primers_infos/Primer data - All infos.csv"))
primer_data


## Loading the unaligned genes

fragment_tab = tibble(read.csv("Data/genes/Isolated genes.csv"))
fragment_tab


## Loading the aligned genes

alignment_12S = tibble(read.csv("Data/genes/complete_gene_12S_aligned.csv"))
alignment_12S 

alignment_16S = tibble(read.csv("Data/genes/complete_gene_16S_aligned.csv"))
alignment_16S 

alignment_COI = tibble(read.csv("Data/genes/complete_gene_COI_aligned.csv"))
alignment_COI 

alignment_CytB = tibble(read.csv("Data/genes/complete_gene_CytB_aligned.csv"))
alignment_CytB 


## Loading manually extracted metabarcodes for reference against automatic extraction

metabarcodes_12S_JDD = rbind_multiple_files("Actino_12S-", "Data/manual_extraction", file_extension = "fas")
metabarcodes_12S_JDD # Teleo1 not done + 3 missing accession numbers

metabarcodes_16S_JDD = rbind_multiple_files("Actino_16S-", "Data/manual_extraction", file_extension = "fas")
metabarcodes_16S_JDD 

metabarcodes_CytB_JDD = rbind_multiple_files("Actino_CytB-", "Data/manual_extraction", file_extension = "fas")
metabarcodes_CytB_JDD # Fish2b & Fish2deg assembled, as well as L14735c & L14735c2, FishCB not done

metabarcodes_CytB_JDD = rbind(tibble(cbind(tibble(PRIMER_SET = "Fish2b"), subset(metabarcodes_CytB_JDD, PRIMER_SET == "Fish2b&deg")[,c(1,3:4)])),
                              tibble(cbind(tibble(PRIMER_SET = "Fish2deg"), subset(metabarcodes_CytB_JDD, PRIMER_SET == "Fish2b&deg")[,c(1,3:4)])),
                              tibble(cbind(tibble(PRIMER_SET = "L14735c"), subset(metabarcodes_CytB_JDD, PRIMER_SET == "Fish2b&deg")[,c(1,3:4)])),
                              tibble(cbind(tibble(PRIMER_SET = "L14735c2"), subset(metabarcodes_CytB_JDD, PRIMER_SET == "L14735c&c2")[,c(1,3:4)])),
                              subset(metabarcodes_CytB_JDD, PRIMER_SET %in% c("L14841", "L14912"))[,c(2,1,3,4)])
metabarcodes_CytB_JDD

# COI metabarcodes manual extraction not done





##### METABARCODE EXTRACTION #####

## Extracting 12S metabarcodes ##

start1 = Sys.time()
metabarcodes_12S = metabarcode_extraction(subset(primer_data, GENE == "12S"), alignment_12S, 
                                          max_mismatch_whole_primer = 5, indels = F, # FALSE because wrong fragment from Teleo1 otherwise
                                          name_col_alignment = "ALIGNED_SEQUENCES", name_col_accession = "ACCESSION", 
                                          name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                                          name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                                          length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH")
end1 = Sys.time()
difftime(end1, start1) # Duration: 7mn
metabarcodes_12S

view_metabarcodes(subset(metabarcodes_12S, PRIMERS != "Teleo1"), metabarcodes_12S_JDD, 
                  data_primer = primer_data, fragment_data = fragment_tab, size_margin = 10, 
                  name_col_fragment = "FRAGMENT_12S", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH")
# Error in the manual alignment of 12SV5 and Teleo2

view_metabarcodes(subset(metabarcodes_12S, PRIMERS == "Teleo1"), 
                  data_primer = primer_data, fragment_data = fragment_tab, nb_sequences = 3, show_all = T, 
                  name_col_fragment = "FRAGMENT_CytB", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH") 
# Checking Teleo1 without manual confirmation -> difficult to align but position in alignment is right and length too


## Extracting 16S metabarcodes ##

start2 = Sys.time()
metabarcodes_16S = metabarcode_extraction(subset(primer_data, GENE == "16S"), alignment_16S, 
                                          max_mismatch_whole_primer = 5, indels = F, # FALSE To be consistent with 12S
                                          name_col_alignment = "ALIGNED_SEQUENCES", name_col_accession = "ACCESSION", 
                                          name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                                          name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                                          length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH")
end2 = Sys.time()
difftime(end2, start2) # Duration: 8mn

view_metabarcodes(metabarcodes_16S, metabarcodes_16S_JDD, data_primer = primer_data, 
                  fragment_data = fragment_tab, size_margin = 10, 
                  name_col_fragment = "FRAGMENT_16S", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH")



## Extracting CytB metabarcodes ##

start3 = Sys.time()
metabarcodes_CytB = metabarcode_extraction(subset(primer_data, GENE == "CytB"), alignment_CytB, 
                                           max_mismatch_whole_primer = 5, indels = T, 
                                           name_col_alignment = "ALIGNED_SEQUENCES", name_col_accession = "ACCESSION", 
                                           name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                                           name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                                           length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH")
end3 = Sys.time()
difftime(end3, start3) # Duration: 5mn
metabarcodes_CytB

view_metabarcodes(subset(metabarcodes_CytB, PRIMERS != "FishCB"), metabarcodes_CytB_JDD, 
                  data_primer = primer_data, fragment_data = fragment_tab, size_margin = 10, 
                  name_col_fragment = "FRAGMENT_CytB", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH") 
# Problem of L14841 with the alignment of primers

view_metabarcodes(subset(metabarcodes_CytB, PRIMERS == "L14841"), metabarcodes_CytB_JDD, 
                  data_primer = primer_data, fragment_data = fragment_tab, nb_sequences = 8, show_all = T, 
                  name_col_fragment = "FRAGMENT_CytB", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH") 
# The problem is due to the alignment of the 5' tip of the forward primer
# But the fragment seems okay (right size and similar to the fragment manually extracted)

view_metabarcodes(subset(metabarcodes_CytB, PRIMERS == "FishCB"), 
                  data_primer = primer_data, fragment_data = fragment_tab, nb_sequences = 4, show_all = T, 
                  name_col_fragment = "FRAGMENT_CytB", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH") 
# Checking FishCB without manual confirmation


## Extracting COI metabarcodes ##

metabarcodes_PS1 = metabarcode_extraction(subset(primer_data, SHORT_NAME == "PS1"), alignment_COI, # ERROR: Primer not found for PS1 !
                                          max_mismatch_whole_primer = 5, indels = T,
                                          name_col_alignment = "ALIGNED_SEQUENCES", name_col_accession = "ACCESSION", 
                                          name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                                          name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                                          length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH") 

start4 = Sys.time()
metabarcodes_COI = metabarcode_extraction(subset(primer_data, SHORT_NAME != "PS1" & GENE == "COI"), alignment_COI, 
                                          max_mismatch_whole_primer = 5, indels = T, 
                                          name_col_alignment = "ALIGNED_SEQUENCES", name_col_accession = "ACCESSION", 
                                          name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                                          name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                                          length1_col_name = "ARTICLE_LENGTH", length2_col_name = "ZHANG_LENGTH")
end4 = Sys.time()
difftime(end4, start4) # Duration: 2mn
metabarcodes_COI 

view_metabarcodes(metabarcodes_COI, data_primer = subset(primer_data, SHORT_NAME != "PS1"), 
                  fragment_data = fragment_tab, size_margin = 10, 
                  name_col_fragment = "FRAGMENT_COI", name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMER_SET", name_col_manual_sequences = "SEQUENCES",
                  name_col_primer_name = "SHORT_NAME", name_col_orientation = "PRIMER_USED",
                  name_col_forward_primer = "FORWARD", name_col_reverse_primer = "REVERSE",
                  length1_col_name = "ARTICLE_LENGTH")


## Removing the PS1 primer from the infos file

primer_data_no_ps1 = subset(primer_data, SHORT_NAME != "PS1")
write.csv(primer_data_no_ps1, "Data/primers_infos/Primer data - All infos - No PS1.csv", row.names = F)
primer_data_no_ps1 = tibble(read.csv("Data/primers_infos/Primer data - All infos - No PS1.csv"))
primer_data_no_ps1





##### METABARCODE FILTERING #####

## Gross Filtering of extreme (1/2 inferior to Q1 or 2 * superior to Q3) metabarcodes ##

all_metabarcodes = rbind(metabarcodes_12S, metabarcodes_16S, metabarcodes_CytB, metabarcodes_COI)
all_metabarcodes

filtered_metabarcodes_list = lapply(split(all_metabarcodes, all_metabarcodes$PRIMERS), function(x) 
  x[which(x$LENGTH >= 1/2 * summary(x$LENGTH)[2][[1]] & x$LENGTH <= 2 * summary(x$LENGTH)[5][[1]]), ])
sapply(filtered_metabarcodes_list, nrow)

filtered_accession = lapply(filtered_metabarcodes_list, function(x) x$ACCESSION)

common_accession = Reduce(intersect, filtered_accession) 

filtered_metabarcodes = subset(all_metabarcodes, ACCESSION %in% common_accession)

write.csv(filtered_metabarcodes, "Data/metabarcodes/Metabarcodes - Length filtering.csv", row.names = F)
filtered_metabarcodes = tibble(read.csv("Data/metabarcodes/Metabarcodes - Length filtering.csv"))
filtered_metabarcodes

length(unique(filtered_metabarcodes$ACCESSION)) # 6140 sequences retained after length correction


## Second filtering of sequences with anormal number of N (>10% of mean length)

filtered_metabarcodes$NUMBER_N = str_count(filtered_metabarcodes$FRAGMENT, "N")
table(filtered_metabarcodes$NUMBER_N)

metabarcodes_number_n_list = split(filtered_metabarcodes, filtered_metabarcodes$PRIMERS)

max_number_n = sapply(metabarcodes_number_n_list, function(x) round(0.1 * mean(x$LENGTH))) # Max number of N = 10% of mean length
for(i in 1:length(max_number_n)){ 
  
  if(i == 1) metabarcodes_number_n_corrected = subset(metabarcodes_number_n_list[[i]], NUMBER_N < max_number_n[i]) 
  
  else metabarcodes_number_n_corrected = rbind(metabarcodes_number_n_corrected, 
                                               subset(metabarcodes_number_n_list[[i]], NUMBER_N < max_number_n[i]))
  
}
metabarcodes_number_n_corrected

filtered_accession_number_n = lapply(split(metabarcodes_number_n_corrected, metabarcodes_number_n_corrected$PRIMERS), function(x) x$ACCESSION)
common_accession_number_n = Reduce(intersect, filtered_accession_number_n) 
filtered_metabarcodes_corrected = subset(filtered_metabarcodes, ACCESSION %in% common_accession_number_n, select = -NUMBER_N)

filtered_metabarcodes_corrected 

write.csv(filtered_metabarcodes_corrected, "Data/metabarcodes/Metabarcodes - N correction - Length filtering.csv", row.names = F)
filtered_metabarcodes_corrected = tibble(read.csv("Data/metabarcodes/Metabarcodes - N correction - Length filtering.csv"))
filtered_metabarcodes_corrected

length(unique(filtered_metabarcodes_corrected$ACCESSION))        # 6089 sequences retained after the N correction





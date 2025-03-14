##### INITIALISATION #####

## Loading required packages and custom functions

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(data.table)


## Loading fish mitogenomes and associated taxonomic informations

act_mitogenomes = tibble(read.csv("Data/mitogenomes_taxonomy/Actinopterygii - Mitogenomes taxonomy.csv"))
act_mitogenomes 

act_mitogenomes_dna = setNames(DNAStringSet(act_mitogenomes$SEQUENCES), act_mitogenomes$ACCESSION)
act_mitogenomes_dna





##### EXTRACTING GENES POSITIONS #####

## Retrieving gene positions in the genbank file

position_gb = fread("Data/NCBI/Infos Actinopterygii mitogenome.gb", header = F, sep = "") 
position_gb


## Preparing the differents inputs of the function gene_position

genes = c("12S", "16S", "COI", "CytB") 

# List of the different names given to the gene in "product" or "gene" on NCBI (separated by | which means or)
# This list has been obtained by trial and error until all cases in which no position could be retrieved was explained
pattern_list = setNames(list(paste0('12S ribosomal RNA"|12S ribosomal RNA subunit"|12 ribosomal RNA"|12S rivbosomal RNA"|',
                                    'small subunit ribosomal RNA"|s-RNA"|s-rRNA"|rrnS"|12S rRNA"|12SrRNA"|12S"|',
                                    'small ribosomal RNA|12 rRNA"'), 
                             paste0('16S ribosomal RNA"|16S ribosomal RNA subunit"|16S rivbosomal RNA"|16 ribosomal RNA"|',
                                    '16S ribosamal RNA"|large subunit ribosomal RNA"|l-RNA"|l-rRNA"|rrnL"|16S rRNA"|16S"|',
                                    'large ribosomal RNA|l6S ribosomal RNA"|16 rRNA"'), 
                             paste0('oxidase subunit I"|ocidase subunit I"|oxidase subunit 1"|oxidase subunit-I"|',
                                    'oxidase subunit-1"|cox1"|CoxI"|oxidase I"|oxidase 1"|oxidase I;|COX1"|COI"|',
                                    'oxidase subunit idase subunit I"|oxidase subunits I"|oxidase subunits 1"|',
                                    'oxydase subunit 1"|oxydase subunit I"'), 
                             paste0('ytochrome b"|ytochorome b"|ytochrome-b"|ytochrome-B"|ytohrome b"|ytchorome b"|CYTB"|',
                                    'ytochome b"|cytB"|cytochrome b;|Cyt b"|cytb"|Cytb"|Cyt b"')), genes)

# Rowname(s) (left column on NCBI) indicating the rows in which the gene position are stored
type_gene_list = setNames(list("rRNA |gene |mRNA ", "rRNA |gene |mRNA ", "CDS |gene ", "CDS |gene "), genes)


## Running the function

start = Sys.time()
position_tab = gene_position(position_gb, genes, pattern_list, type_gene_list)
end = Sys.time()
difftime(end, start)  # Duration: 1.5mn
position_tab


## Checking length of sequences missing some fragments 

summary(nchar(subset(act_mitogenomes, ACCESSION %in% position_tab$COMPLETE_POSITIONS$ACCESSION)$SEQUENCES)) # Min = 15000
summary(nchar(subset(act_mitogenomes, ACCESSION %in% position_tab$MISSING_POSITIONS$ACCESSION)$SEQUENCES)) # Mostly uncomplete mitogenome


## Manual verification of plausible complete sequences with missing genes -> Example code for NC_030485 (same for 9 other)

missing_positions_sequences = act_mitogenomes_dna[names(act_mitogenomes_dna) %in% position_tab$MISSING_POSITIONS$ACCESSION]
missing_positions_sequences

normal_length_missing_accession = names(missing_positions_sequences[width(missing_positions_sequences) > 15000]) # > 15000 = plausible
normal_length_missing_accession

normal_length_missing_positions = subset(position_tab$MISSING_POSITIONS, ACCESSION %in% normal_length_missing_accession)
normal_length_missing_positions

manual_split_gb = split(position_gb, findInterval(1:nrow(position_gb), which(grepl("LOCUS", position_gb$V1)) + 1))
example_check_gb = lapply(manual_split_gb[-1], function(x) grep("NC_030485", x$V1))
example_check_gb = example_check_gb[lapply(example_check_gb, length) > 0]
example_check_gb = names(example_check_gb[lapply(example_check_gb, function(x) which(3 %in% x)) > 0])
View(manual_split_gb[-1][[as.numeric(na.omit(example_check_gb))]])

# Reasons of missing positions:
# Same name 12S & 16S: NC_024573 // KJ643927
# Uncomplete mitogenome (missing CytB): MG599474 // NC_033859 // KU291530
# Uncomplete mitogenome (missing 12S & 16S): KT867088 // KT867089 // NC_029216 
# Uncomplete mitogenome (missing 12S): U62532 // AM919428 // U62532


## Removing sequences with abscence of some positions explained

verified_missing_gene = c("MG599474", "NC_033859", "KU291530", "NC_024573", "KJ643927", "KT867088", 
                          "KT867089", "NC_029216", "U62532", "AM919428", "U62532")

normal_length_missing_positions = subset(normal_length_missing_positions, !(ACCESSION %in% verified_missing_gene))
normal_length_missing_positions # All plausible complete mitogenome with missing positions were explained above 

# If no reason is found, just modify the gene_position() function by adding the unpredicted case to the list of accepted variations/errors


## Removing replicated (after reviewing) accession numbers (NC_)

complete_positions = position_tab$COMPLETE_POSITIONS[!grepl("_", position_tab$COMPLETE_POSITIONS$ACCESSION),]
complete_positions


## Saving the positions

write.csv(complete_positions, "Data/genes/Gene positions.csv", row.names = F)
complete_positions = tibble(read.csv("Data/genes/Gene positions.csv"))
complete_positions





##### EXTRACTING GENES #####

## Preparing the position file

world_dna_verified = act_mitogenomes_dna[names(act_mitogenomes_dna) %in% complete_positions$ACCESSION]
position_verified = subset(complete_positions, ACCESSION %in% names(world_dna_verified))
position_verified = position_verified[order(position_verified$ACCESSION),]
position_verified

# Adding 100nt before CytB to include L14735c & L14735c2 primers which are in tRNA Glu, the gene before CytB
position_verified = tibble(cbind(position_verified[,1:10], START_CytB = position_verified$START_CytB - 100, position_verified[,12:13])) 


## Preparing the sequence database

correct_dna_data = DNAStringSet(act_mitogenomes$SEQUENCES)
names(correct_dna_data) = act_mitogenomes$ACCESSION
correct_dna_data = correct_dna_data[names(correct_dna_data) %in% names(world_dna_verified)]
correct_dna_data


## Extracting metabarcodes

start = Sys.time()
fragment_tab = gene_extraction(position_verified, correct_dna_data, genes = genes, name_col_accession = "ACCESSION")
fragment_tab
end = Sys.time() ; difftime(end, start)  # Duration: 50s


## Filtering metabarcodes with aberrant length

length_fragment = tibble(cbind(fragment_tab[,1], apply(fragment_tab[,-1], c(1,2), nchar)))
length_fragment

apply(length_fragment[,-1], 2, summary)

fragment_tab = fragment_tab[-which(length_fragment$FRAGMENT_12S < 500 | length_fragment$FRAGMENT_12S > 2000 | 
                                     length_fragment$FRAGMENT_16S < 1000 | length_fragment$FRAGMENT_16S > 3000), ]
fragment_tab 


## Removing the complement sequences that will cause differences between non-specific amplicons and metabarcodes

complement_accession = subset(complete_positions, ORIENTATION_12S == "complement" | ORIENTATION_16S == "complement" | 
                                ORIENTATION_COI == "complement" | ORIENTATION_CytB == "complement")$ACCESSION

fragment_tab = subset(fragment_tab, !(ACCESSION %in% complement_accession))


## Saving the genes extracted

write.csv(fragment_tab, "Data/genes/Isolated genes.csv", row.names = F)
fragment_tab = tibble(read.csv("Data/genes/Isolated genes.csv"))
fragment_tab




##### ALIGNMENT OF EACH COMPLETE GENE #####

## Converting the table with fragments to a long format

fragment_tab_long = fragment_tab %>% pivot_longer(!ACCESSION, values_to = "FRAGMENT", names_to = "GENE")
fragment_tab_long$GENE = word(fragment_tab_long$GENE, 2, sep = fixed("_"))
fragment_tab_long


## Alignment of each genes

start = Sys.time()
alignment_all_genes = decipher_rapid_alignment(fragment_tab, name_col_group = "GENE", 
                                               name_col_accession = "ACCESSION",
                                               name_col_fragment = "FRAGMENT",  
                                               progressive_saving = T, final_saving = F, 
                                               global_name = "Data/genes/complete_gene")  
end = Sys.time()
difftime(end, start)  # Duration: 9h


## Visualising of some sequences for each alignments

BrowseSeqs(DNAStringSet(subset(alignment_all_genes, GENE == "12S")$ALIGNED_SEQUENCES[1:20]))

BrowseSeqs(DNAStringSet(subset(alignment_all_genes, GENE == "16S")$ALIGNED_SEQUENCES[1:20]))

BrowseSeqs(DNAStringSet(subset(alignment_all_genes, GENE == "COI")$ALIGNED_SEQUENCES[1:20]))

BrowseSeqs(DNAStringSet(subset(alignment_all_genes, GENE == "CytB")$ALIGNED_SEQUENCES[1:20]))


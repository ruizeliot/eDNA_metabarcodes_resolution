##### INITIALISATION #####

## Loading required packages and custom functions

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(ggplot2)
library(extrafont)


## Loading required packages here

library(ape)
library(ggtree)
library(adegenet)
library(randomcoloR)


## Loading all needed files

act_mitogenomes_bin_corrected = tibble(read.csv("Data/taxonomy_BIN_infos/Actinopterygii - Taxonomy - BIN.csv"))
act_mitogenomes_bin_corrected

filtered_metabarcodes_final = tibble(read.csv("Data/metabarcodes/Filtered metabarcodes final.csv"))
filtered_metabarcodes_final

primers = sort(unique(filtered_metabarcodes_final$PRIMERS))
primers





##### CREATING A TREE FOR AN EXAMPLE FAMILY #####

## Creating the dataset files

salmonidae_dataset = subset(act_mitogenomes_bin_corrected, FAMILY == "Salmonidae")
salmonidae_dataset

salmonidae_dna = subset(filtered_metabarcodes_final, ACCESSION %in% salmonidae_dataset$ACCESSION)
salmonidae_dna

salmonidae_dna_taxo = tibble(merge(salmonidae_dataset, salmonidae_dna))
salmonidae_dna_taxo


## Creating the NJ tree

salmonidae_dna_taxo_barcode = subset(salmonidae_dna_taxo, PRIMERS == "FishF1-FishR1")
salmonidae_dna_taxo_barcode

salmonidae_fasta_file = c(rbind(paste0(">", salmonidae_dna_taxo_barcode$BIN, " - ", salmonidae_dna_taxo_barcode$ACCESSION),
                          FRAGMENT = salmonidae_dna_taxo_barcode$FRAGMENT))
write.table(salmonidae_fasta_file, file = "Output/example_resolution_analysis/Salmonidae barcode.fasta", 
            row.names = F, col.names = F, quote = F)

dna = fasta2DNAbin(file = "Output/example_resolution_analysis/Salmonidae barcode.fasta", quiet = T)
distance_matrix_barcode = dist.dna(dna)
phylo_barcode = nj(distance_matrix_barcode)


## Settings different discernable colours per BIN 

group_bins = lapply(split(salmonidae_dna_taxo_barcode, salmonidae_dna_taxo_barcode$BIN), function(x) paste0(x$BIN, " - ", x$ACCESSION))
group_bins

phylo_barcode = groupOTU(phylo_barcode, group_bins)

set.seed(2) # Change seed to change random colour
colours = as.list(distinctColorPalette(length(group_bins)))

for(i in 1:length(colours)){
  
  if(i == 1) barcode_colours = setNames(rep(colours[[i]], length(group_bins[[i]])), 
                                        rep(names(group_bins)[i], length(group_bins[[i]])))
  
  else barcode_colours = c(barcode_colours, setNames(rep(colours[[i]], length(group_bins[[i]])), 
                                                     rep(names(group_bins)[i], length(group_bins[[i]]))))
  
}

options(ignore.negative.edge = T)
tree_barcode = ggtree(phylo_barcode, right = F) +
  geom_tiplab(aes(label = "▬", colour = group), align = T, size = 4, linesize = NA) + # UNICODE OF LONG DASH: &#x25AC
  geom_treescale(x = 0, y = -40, offset = 2) +
  theme(legend.position = 'none') +
  scale_color_manual(values = barcode_colours)
tree_barcode + hexpand(0.5)





##### RESOLUTION ANALYSIS #####

## Inter-BIN errors heatmap

similarity_comparisons_list_salmonidae = 
  vsearch_pairwise_similarity(metabarcodes_data = salmonidae_dna, 
                              infos_data = salmonidae_dataset, 
                              min_similarity_threshold = 99, 
                              max_similarity_threshold = 100,
                              name_col_primer = "PRIMERS", 
                              name_col_accession = "ACCESSION", 
                              name_col_fragment = "FRAGMENT", 
                              keep_col_infos_data = c("BIN", "SPECIES", "GENUS", "FAMILY", "ORDER"), 
                              order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                                   "BIN1", "BIN2", 
                                                   "SPECIES1", "SPECIES2",  
                                                   "GENUS1", "GENUS2",  
                                                   "FAMILY1", "FAMILY2",  
                                                   "ORDER1", "ORDER2"),
                              # Change it with the path + name of your VSEARCH program
                              program_name = "vsearch-2.17.1-win-x86_64.exe")

inter_bin_salmonidae_list = lapply(similarity_comparisons_list_salmonidae, function(x) x[which(x$BIN1 != x$BIN2), ])

inter_bin_errors = data.frame(matrix(ncol = length(primers), nrow = length(labels(dna))))
colnames(inter_bin_errors) = primers
row.names(inter_bin_errors) = labels(dna)
inter_bin_errors[is.na(inter_bin_errors)] = "No error"
inter_bin_errors

for(i in 1:length(primers)){
  
  if(nrow(inter_bin_salmonidae_list[[i]]) != 0){
    
    inter_bin_errors_names = c(paste0(inter_bin_salmonidae_list[[i]]$BIN1, " - ", inter_bin_salmonidae_list[[i]]$ACCESSION1),
                               paste0(inter_bin_salmonidae_list[[i]]$BIN2, " - ", inter_bin_salmonidae_list[[i]]$ACCESSION2))
    
    inter_bin_errors_names = inter_bin_errors_names[!duplicated(inter_bin_errors_names)]
    
    inter_bin_errors[[primers[i]]] = replace(inter_bin_errors[[primers[i]]], 
                                             which(rownames(inter_bin_errors) %in% inter_bin_errors_names), "Inter-BIN error")
    
  }
  
}

inter_bin_errors = inter_bin_errors[, names(sort(sapply(inter_bin_salmonidae_list, nrow)))] # To order by number of error
inter_bin_errors

gheatmap(tree_barcode, inter_bin_errors, offset = 0.01, width = 2.5, 
         colnames_angle = -90, colnames_offset_y = -3, hjust = 0,
         font.size = 4, color = NA, family = "Segoe UI Semilight") +
  scale_fill_manual(values = c("#FF1300", "#23D33F")) +
  theme(legend.position = 'none')

ggsave("Graph/example_resolution_analysis/Inter-BIN errors - Salmonidae - 99 - Automatic.png", dpi = 900, type = "cairo", 
       width = 21.40479, height = 22.38375, units = "cm")


## Intra-BIN errors heatmap

intra_bin_comparisons_salmonidae = intra_taxa_similarity(metabarcodes_data = salmonidae_dna, 
                                                         infos_data = salmonidae_dataset, 
                                                         name_col_taxa = "BIN", 
                                                         name_col_primer = "PRIMERS", 
                                                         name_col_accession = "ACCESSION", 
                                                         name_col_fragment = "FRAGMENT", 
                                                         keep_col_infos_data = c("BIN", "SPECIES", "GENUS", "FAMILY", "ORDER"), 
                                                         order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                                                              "BIN1", "BIN2", 
                                                                              "SPECIES1", "SPECIES2",  
                                                                              "GENUS1", "GENUS2",  
                                                                              "FAMILY1", "FAMILY2",  
                                                                              "ORDER1", "ORDER2"),
                                                         # Change it with the path + name of your VSEARCH program
                                                         program_name = "vsearch-2.17.1-win-x86_64.exe")

var(sapply(intra_bin_comparisons_salmonidae, nrow)) == 0 # Equal number of rows for all primers as expected (so variance of 0)

intra_bin_salmonidae_list = lapply(intra_bin_comparisons_salmonidae, function(x) 
  x[which(x$BIN1 == x$BIN2 & x$SIMILARITY < 99),])

intra_bin_errors = data.frame(matrix(ncol = length(primers), nrow = length(labels(dna))))
colnames(intra_bin_errors) = primers
row.names(intra_bin_errors) = labels(dna)
intra_bin_errors[is.na(intra_bin_errors)] = "No error"

for(i in 1:length(primers)){
  
  if(nrow(intra_bin_salmonidae_list[[i]]) != 0){
    
    intra_bin_errors_names = c(paste0(intra_bin_salmonidae_list[[i]]$BIN1, " - ", intra_bin_salmonidae_list[[i]]$ACCESSION1),
                               paste0(intra_bin_salmonidae_list[[i]]$BIN2, " - ", intra_bin_salmonidae_list[[i]]$ACCESSION2))
    
    intra_bin_errors_names = intra_bin_errors_names[!duplicated(intra_bin_errors_names)]
    
    intra_bin_errors[[primers[i]]] = replace(intra_bin_errors[[primers[i]]], 
                                             which(rownames(intra_bin_errors) %in% intra_bin_errors_names), "Intra-BIN error")
    
  }
  
}

bin_name = word(rownames(intra_bin_errors), 1)
bin_name_unique = rownames(intra_bin_errors)[-which(duplicated(bin_name) | duplicated(bin_name, fromLast = T))]
intra_bin_errors[which(rownames(intra_bin_errors) %in% bin_name_unique),] = "NA" # Just NA to get it white automatically
intra_bin_errors = intra_bin_errors[, names(sort(sapply(intra_bin_salmonidae_list, nrow)))] # To order by number of error

gheatmap(tree_barcode, intra_bin_errors, offset = 0.01, width = 2.5, 
         colnames_angle = -90, colnames_offset_y = -3, hjust = 0,
         font.size = 4, color = NA, family = "Segoe UI Semilight") +
  scale_fill_manual(values = c("#FF1300", "blue",  "#23D33F")) +
  theme(legend.position = 'none')

ggsave("Graph/example_resolution_analysis/Intra-BIN errors - Salmonidae - 99 - Automatic.png", dpi = 900, type = "cairo", 
       width = 21.40479, height = 22.38375, units = "cm")





##### OTHER FIGURES TO HELP CREATE THE FINAL ONE ON PHOTOSHOP #####

## Coloured branches scaled

tree_barcode_color = ggtree(phylo_barcode, aes(color = group), right = F) +
  geom_tiplab(aes(label = "▬", colour = group), align = T, size = 4, linesize = NA) + # UNICODE: &#x25AC
  geom_treescale(x = 0, y = -40, offset = 2) +
  theme(legend.position = 'none') +
  scale_color_manual(values = barcode_colours)
tree_barcode_color + hexpand(0.5)

gheatmap(tree_barcode_color, intra_bin_errors, offset = 0.01, width = 2.5, 
         colnames_angle = -90, colnames_offset_y = -3, hjust = 0,
         font.size = 4, color = NA, family = "Segoe UI Semilight") +
  scale_fill_manual(values = c("#FF1300", "blue",  "#23D33F")) +
  theme(legend.position = 'none')

ggsave("Graph/example_resolution_analysis/Salmonidae - 99 - Colored tree.png", dpi = 900, type = "cairo", 
       width = 21.40479, height = 22.38375, units = "cm")


## Coloured branches not scaled

tree_barcode_color = ggtree(phylo_barcode, aes(color = group), right = F, branch.length = "none") +
  geom_tiplab(aes(label = "▬", colour = group), align = T, size = 4, linesize = NA) + # UNICODE: &#x25AC
  geom_treescale(x = 0, y = -40, offset = 2) +
  theme(legend.position = 'none') +
  scale_color_manual(values = barcode_colours)
tree_barcode_color + hexpand(0.5)

gheatmap(tree_barcode_color, intra_bin_errors, offset = 3, width = 2.5, 
         colnames_angle = -90, colnames_offset_y = -3, hjust = 0,
         font.size = 4, color = NA, family = "Segoe UI Semilight") +
  scale_fill_manual(values = c("#FF1300", "blue",  "#23D33F")) +
  theme(legend.position = 'none')

ggsave("Graph/example_resolution_analysis/Salmonidae - 99 - Colored tree - Not scaled.png", dpi = 900, type = "cairo", 
       width = 21.40479, height = 22.38375, units = "cm")


## Not coloured branches not scaled

tree_barcode_color = ggtree(phylo_barcode, right = F, branch.length = "none") +
  geom_tiplab(aes(label = "▬", colour = group), align = T, size = 4, linesize = NA) + # UNICODE: &#x25AC
  geom_treescale(x = 0, y = -40, offset = 2) +
  theme(legend.position = 'none') +
  scale_color_manual(values = barcode_colours)
tree_barcode_color + hexpand(0.5)

gheatmap(tree_barcode_color, intra_bin_errors, offset = 3, width = 2.5, 
         colnames_angle = -90, colnames_offset_y = -3, hjust = 0,
         font.size = 4, color = NA, family = "Segoe UI Semilight") +
  scale_fill_manual(values = c("#FF1300", "blue",  "#23D33F")) +
  theme(legend.position = 'none')

ggsave("Graph/example_resolution_analysis/Salmonidae - 99 - Not colored tree - Not scaled.png", dpi = 900, type = "cairo", 
       width = 21.40479, height = 22.38375, units = "cm")

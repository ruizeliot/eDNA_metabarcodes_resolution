##### INITIALISATION #####

## Loading required packages and custom functions

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(ggplot2)
library(extrafont)
library(tidyverse)
library(patchwork)
library(ggtext)


## Default theme for the following ggplot graph

theme_ = function(base_family = "Segoe UI Semilight", ...){
  theme_bw(base_family = base_family, ...) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust=0.5, vjust=0.5, face='bold',size=18),
      axis.title.y = element_text(margin=unit(c(0,0.2,0,0), "cm")),
      axis.text.y = element_text(margin=unit(c(0.1,0.1,0.1,0.1), "cm")),
      axis.ticks.length = unit(0.05, "in"))
}


## Preparing the fonts for good PDF rendering

library(showtext)

font_add("Segoe UI Semilight", regular = "C:/Windows/Fonts/segoeuil.ttf")
font_add("Segoe UI Semibold", regular = "C:/Windows/Fonts/seguisb.ttf")
font_add("Segoe UI", regular = "C:/Windows/Fonts/segoeui.ttf", bold = "C:/Windows/Fonts/segoeuib.ttf",
         italic = "C:/Windows/Fonts/segoeuii.ttf", bolditalic = "C:/Windows/Fonts/segoeuiz.ttf")
font_add("Segoe UI Semibold", regular = "C:/Windows/Fonts/seguisb.ttf", italic = "C:/Windows/Fonts/seguisbi.ttf")
font_add("Segoe UI Semilight", regular = "C:/Windows/Fonts/segoeuil.ttf", italic = "C:/Windows/Fonts/seguisli.ttf")

showtext_opts(dpi = 600) # Must match the DPI used when saving graphs
showtext_auto(enable = TRUE)



## Loading all needed files

act_mitogenomes_bin_corrected = tibble(read.csv("Data/taxonomy_BIN_infos/Actinopterygii - Taxonomy - BIN.csv"))
act_mitogenomes_bin_corrected

filtered_metabarcodes_final = tibble(read.csv("Data/metabarcodes/Filtered metabarcodes final.csv"))
filtered_metabarcodes_final

mean_length_metabarcodes = tibble(read.csv("Data/primers_infos/Mean length data.csv"))
mean_length_metabarcodes %>% print(n = 22)

primers = sort(unique(filtered_metabarcodes_final$PRIMERS))
primers

filtered_metabarcodes_final_infos = tibble(merge(filtered_metabarcodes_final, mean_length_metabarcodes))
filtered_metabarcodes_final_infos

intra_bin_comparisons = readRDS("ADLIFISH_eruiz/clustering_ml_decipher/Intra-BIN - Raw comparisons - All primers.rds") # Duration: 15s
intra_bin_comparisons[1:2]

inter_bin_comparisons = readRDS("ADLIFISH_eruiz/clustering_ml_decipher/Inter-BIN - Raw comparisons - All primers.rds") # Duration: 1mn
inter_bin_comparisons[1:2]

intra_bin_summary = readRDS("Output/intra_BIN/Intra-BIN analysis - Summary - All primers.rds")
intra_bin_summary[1:2]

inter_bin_summary = readRDS("Output/inter_BIN/Inter-BIN analysis - Summary - All primers.rds")
inter_bin_summary[1:2]




##### CLUSTERING METABARCODES USING VARIOUS METHODS BASED ON PAIRWISE OR GLOBAL SIMILARITY #####

## Aligning each metabarcode with the default number of iterations and refinements for MOTHUR

start0 = Sys.time()
for(i in 1:length(primers)){
  
  cat(paste0("Aligning metabarcodes for ", primers[i], ":\n"))
  
  alignment_per_primer_data_subset = subset(filtered_metabarcodes_final_infos, PRIMERS == primers[i])
  alignment_per_primer_data_subset$ALIGNMENT_PER_PRIMER = as.character(AlignSeqs(DNAStringSet(
    alignment_per_primer_data_subset$FRAGMENT), verbose = F, processors = 14))
  
  if(i == 1) alignment_per_primer_data = alignment_per_primer_data_subset
  else alignment_per_primer_data = rbind(alignment_per_primer_data, alignment_per_primer_data_subset)
  
}
end0 = Sys.time()
difftime(end0, start0) # Duration: 1.7h on a standard laptop

write.csv(alignment_per_primer_data, "Output/comparison_clustering/Filtered metabarcodes final with alignment per primer.csv", row.names = F)
alignment_per_primer_data = tibble(read.csv("Output/comparison_clustering/Filtered metabarcodes final with alignment per primer.csv"))
alignment_per_primer_data


## Clustering based on with hierarchical methods implemented in DECIPHER using on pairwise similarity

start1 = Sys.time()
clusters_decipher_pairwise_similarity = cluster_metabarcodes_decipher(infos_data = act_mitogenomes_bin_corrected, 
                                                                      intra_comparisons_list = intra_bin_comparisons, 
                                                                      inter_comparisons_list = inter_bin_comparisons, 
                                                                      similarity_threshold = 90:99, 
                                                                      tree_methods = c("single", "complete", "UPGMA", "WPGMA", "NJ"), 
                                                                      name_col_similarity = "SIMILARITY",
                                                                      name_col_accession = "ACCESSION", 
                                                                      name_col_taxa = "BIN",
                                                                      nb_threads = 32)
end1 = Sys.time()
difftime(end1, start1) # Duration: 1.4h on a standard laptop

saveRDS(clusters_decipher_pairwise_similarity, "Output/comparison_clustering/Clusters per primers list - Simple clustering methods.rds")
clusters_decipher_pairwise_similarity = readRDS("Output/comparison_clustering/Clusters per primers list - Simple clustering methods.rds")
clusters_decipher_pairwise_similarity[1:2]


## Clustering based on with the new Clusterize method implemented in DECIPHER using on pairwise similarity

start2 = Sys.time()
clusters_accession_list_clusterize_methods = cluster_metabarcodes_decipher(unaligned_seqs_per_primer_df = filtered_metabarcodes_final_infos,
                                                                           similarity_threshold = 90:99, 
                                                                           tree_methods = "clusterize", 
                                                                           name_col_fragment = "FRAGMENT",
                                                                           name_col_accession = "ACCESSION", 
                                                                           name_col_primer = "PRIMERS",
                                                                           nb_threads = 14)
end2 = Sys.time()
difftime(end2, start2) # Duration: 50mn on a standard laptop

saveRDS(clusters_accession_list_clusterize_methods, "Output/comparison_clustering/Clusters per primers list - Clusterize method.rds")
clusters_accession_list_clusterize_methods = readRDS("Output/comparison_clustering/Clusters per primers list - Clusterize method.rds")
clusters_accession_list_clusterize_methods[1:2]


## Clustering based on heuristic methods implemented in SWARM and MOTHUR using global and pairwise similarity, respectively

start3 = Sys.time()
mothur_swarm_heuristic_clustering = 
  cluster_metabarcodes_mothur_swarm(alignment_per_primer_df = alignment_per_primer_data,
                                    dist_threshold = seq(0.01, 0.10, by = 0.01), cutoff_dist_matrix = 0.2,
                                    output_path = "Output/mothur_output", 
                                    method = c("opti", "agc", "dgc", "swarm"),
                                    name_col_accession = "ACCESSION", name_col_alignment = "ALIGNMENT_PER_PRIMER", name_col_primers = "PRIMERS",
                                    mothur_path = "mothur.exe", vsearch_path = "vsearch-2.17.1-win-x86_64.exe",
                                    swarm_path = "swarm.exe", nb_threads = 14)
end3 = Sys.time()
difftime(end3, start3) # Duration: 2.5h on a standard laptop
mothur_swarm_heuristic_clustering[1:2]


## Clustering based on hierarchical methods implemented in MOTHUR using the full distance matrix to avoid getting problems when making trees

start4 = Sys.time()
mothur_hierarchical_clustering = 
  cluster_metabarcodes_mothur_swarm(alignment_per_primer_df = alignment_per_primer_data,
                                    dist_threshold = seq(0.01, 0.10, by = 0.01), cutoff_dist_matrix = 1,
                                    output_path = "Output/mothur_output_simple_methods", method = c("furthest", "nearest", "average"),
                                    name_col_accession = "ACCESSION", name_col_alignment = "ALIGNMENT_PER_PRIMER", name_col_primers = "PRIMERS",
                                    mothur_path = "mothur.exe", vsearch_path = "vsearch-2.17.1-win-x86_64.exe",
                                    swarm_path = "swarm.exe", nb_threads = 14)
end4 = Sys.time()
difftime(end4, start4) # Duration: 2.5h on a standard laptop
mothur_hierarchical_clustering[1:2]


## Binding the clusters obtained with all methods together

all_clusters_per_methods = list()

for(i in 1:length(primers)){
  
  all_primers_names = c(names(mothur_swarm_heuristic_clustering)[i], names(mothur_hierarchical_clustering)[i],
                        names(clusters_decipher_pairwise_similarity)[i], names(clusters_decipher_clusterize)[i])
  
  if(length(unique(all_primers_names)) != 1) warning("List elements names are not the same.")
  
  all_clusters_per_methods[[i]] = rbind(mothur_swarm_heuristic_clustering[[i]], mothur_hierarchical_clustering[[i]],
                                        clusters_decipher_pairwise_similarity[[i]], clusters_decipher_clusterize[[i]])
  
  names(all_clusters_per_methods)[i] = unique(all_primers_names)
  
}

saveRDS(all_clusters_per_methods, "Output/comparison_clustering/All clusters per thresholds - MOTHUR & SWARM & DECIPHER.rds")
all_clusters_per_methods = readRDS("Output/comparison_clustering/All clusters per thresholds - MOTHUR & SWARM & DECIPHER.rds")
all_clusters_per_methods[1:2]




##### COMPUTING THE NUMBER OF FALSE-POSITIVES AND FALSE-NEGATIVES FROM OTU OBTAINED WITH VARIOUS CLUSTERING METHODS #####

## Adding the BIN to the clusters table

all_clusters_per_methods_bins = lapply(1:length(all_clusters_per_methods), function(i)
  tibble(merge(subset(filtered_metabarcodes_final[,1:3], PRIMERS == names(all_clusters_per_methods)[i]),
               all_clusters_per_methods[[i]], by = "ACCESSION")))
names(all_clusters_per_methods_bins) = names(all_clusters_per_methods)
all_clusters_per_methods_bins[1:2]


## Computing the number of intra-BIN (false-positive) and inter-BIN (false-negative) errors

start5 = Sys.time()
clustering_errors_mothur = compute_otu_errors_metabarcodes(all_clusters_per_methods_bins, name_col_method = "METHOD", 
                                                     name_col_taxa = "BIN", clusters_cols_prefix = "CLUSTERS_")
end5 = Sys.time()
difftime(end5, start5) # Duration: 8mn on a standard laptop
clustering_errors_mothur[[1]][[2]]


## Assembling summary tables together

clustering_intra_bin_summary_df = bind_rows(clustering_errors_mothur$INTRA_ERRORS, .id = "SIMILARITY")
clustering_intra_bin_summary_df

clustering_inter_bin_summary_df = bind_rows(clustering_errors_mothur$INTER_ERRORS, .id = "SIMILARITY")
clustering_inter_bin_summary_df

clustering_intra_bin_summary_df = clustering_intra_bin_summary_df[order(
  clustering_intra_bin_summary_df$SIMILARITY, clustering_intra_bin_summary_df$PRIMERS, clustering_intra_bin_summary_df$METHOD),]
clustering_intra_bin_summary_df

clustering_inter_bin_summary_df = clustering_inter_bin_summary_df[order(
  clustering_inter_bin_summary_df$SIMILARITY, clustering_inter_bin_summary_df$PRIMERS, clustering_inter_bin_summary_df$METHOD),]
clustering_inter_bin_summary_df

clustering_intra_inter_bin_summary = tibble(cbind(clustering_intra_bin_summary_df, clustering_inter_bin_summary_df[,-c(1:3)]))

clustering_intra_inter_bin_summary = subset(clustering_intra_inter_bin_summary, select = 
                                   c("SIMILARITY", "PRIMERS", "METHOD", "NUMBER_INTRA_ERRORS", "PERCENT_INTRA_ERRORS", 
                                     "NUMBER_INTER_ERRORS", "PERCENT_INTER_ERRORS"))
clustering_intra_inter_bin_summary

clustering_intra_inter_bin_summary$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = 
  clustering_intra_inter_bin_summary$PERCENT_INTRA_ERRORS +
  clustering_intra_inter_bin_summary$PERCENT_INTER_ERRORS

clustering_intra_inter_bin_summary = tibble(merge(mean_length_metabarcodes[,1:3], clustering_intra_inter_bin_summary))


## Changing metabarcodes names

clustering_intra_inter_bin_summary_subset = subset(clustering_intra_inter_bin_summary, !(PRIMERS %in% c("Fish2deg", "L14735c2")))

clustering_intra_inter_bin_summary_subset$PRIMERS = replace(clustering_intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(clustering_intra_inter_bin_summary_subset$PRIMERS %in% "Fish2b"), "Fish2b/deg")

clustering_intra_inter_bin_summary_subset$PRIMERS = replace(clustering_intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(clustering_intra_inter_bin_summary_subset$PRIMERS %in% "L14735c"), "L14735c/c2")

clustering_intra_inter_bin_summary_subset$PRIMERS = replace(clustering_intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(clustering_intra_inter_bin_summary_subset$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")


## Merging the results of clusterings to the table summarizing results from the intra/inter-BIN analyses

intra_bin_summary_df = bind_rows(intra_bin_summary, .id = "SIMILARITY")
intra_bin_summary_df = intra_bin_summary_df[order(intra_bin_summary_df$SIMILARITY),]
intra_bin_summary_df

inter_bin_summary_df = bind_rows(inter_bin_summary, .id = "SIMILARITY")
inter_bin_summary_df = inter_bin_summary_df[order(inter_bin_summary_df$SIMILARITY),]
inter_bin_summary_df

intra_inter_bin_summary = tibble(cbind(intra_bin_summary_df, inter_bin_summary_df[,-c(1:2)]))
intra_inter_bin_summary = subset(intra_inter_bin_summary, select = 
                                   c("SIMILARITY", "PRIMERS", "NUMBER_INTRA_ERRORS", "PERCENT_INTRA_ERRORS", 
                                     "NUMBER_INTER_ERRORS", "PERCENT_INTER_ERRORS"))
intra_inter_bin_summary$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = 
  intra_inter_bin_summary$PERCENT_INTRA_ERRORS + intra_inter_bin_summary$PERCENT_INTER_ERRORS

intra_inter_bin_summary = tibble(merge(mean_length_metabarcodes[,1:3], intra_inter_bin_summary))

intra_inter_bin_summary_subset = subset(intra_inter_bin_summary, !(PRIMERS %in% c("Fish2deg", "L14735c2")))

intra_inter_bin_summary_subset$PRIMERS = replace(intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(intra_inter_bin_summary_subset$PRIMERS %in% "Fish2b"), "Fish2b/deg")

intra_inter_bin_summary_subset$PRIMERS = replace(intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(intra_inter_bin_summary_subset$PRIMERS %in% "L14735c"), "L14735c/c2")

intra_inter_bin_summary_subset$PRIMERS = replace(intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(intra_inter_bin_summary_subset$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")
intra_inter_bin_summary_subset

all_methods_grouped = tibble(rbind(clustering_intra_inter_bin_summary_subset,
                                   cbind(intra_inter_bin_summary_subset[,1:4], tibble(METHOD = "raw sim"), intra_inter_bin_summary_subset[,5:9])))
all_methods_grouped


## Preparing the table for plotting

all_methods_grouped_filtered = all_methods_grouped %>% filter(METHOD != "agc") %>% group_by(PRIMERS, METHOD) %>%
  mutate(SIMILARITY = as.numeric(SIMILARITY),
         BEST_THRESHOLD = ifelse(min(CUMULATED_PERCENT_ERRORS_PER_ANALYSIS, na.rm = TRUE) == CUMULATED_PERCENT_ERRORS_PER_ANALYSIS, TRUE, FALSE)) %>%
  mutate(COUNT_TRUE = sum(BEST_THRESHOLD, na.rm = TRUE),
         MIN_SIMILARITY_TRUE = ifelse(COUNT_TRUE > 1, min(SIMILARITY[BEST_THRESHOLD], na.rm = TRUE), NA_real_),
         BEST_THRESHOLD = ifelse(COUNT_TRUE > 1 & BEST_THRESHOLD & SIMILARITY != MIN_SIMILARITY_TRUE, FALSE, BEST_THRESHOLD)) %>%
  select(-COUNT_TRUE, -MIN_SIMILARITY_TRUE) %>% ungroup() %>%
  mutate(PRIMERS_GENE = factor(paste0(PRIMERS, " (", GENE, ")"), levels = 
                                 c("FishF1-R1 (COI)", "Minibar (COI)", "L14912 (CytB)", "L14841 (CytB)", "L14735c/c2 (CytB)", "FishCB (CytB)", "Teleo2 (12S)", 
                                   "MiFish (12S)", "Fish16S (16S)", "Ac12S (12S)", "AcMDB (12S)", "16SFD (16S)", "Ve16S (16S)", "Vert16S (16S)", "12SV5 (12S)", 
                                   "12SF1R1 (12S)", "Ac16S (16S)", "L2513 (16S)", "Teleo1 (12S)", "Fish2b/deg (CytB)")))

scale_color_distinct = c('black', 'yellow', '#b15928', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a', '#1f78b4', '#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6')

methods_full_names = setNames(c("MOTHUR (global sim.)\nUPGMA", "DECIPHER (pairwise sim.)\nFurthest neighbor", "MOTHUR (global sim.)\nVSEARCH", "MOTHUR (global sim.)\nFurthest neighbor", 
           "MOTHUR (global sim.)\nNearest neighbor", "DECIPHER (pairwise sim.)\nNeighbor joining", "MOTHUR (global sim.)\nOPTICLUST", "Intra/inter-BIN analyses\n(pairwise sim.)",
           "DECIPHER (pairwise sim.)\nNearest neighbor", "SWARM (pairwise sim.)", "DECIPHER (pairwise sim.)\nUPGMA", "DECIPHER (pairwise sim.)\nWPGMA", "DECIPHER (pairwise sim.)\nClusterize"),
         c("average", "complete", "dgc", "furthest", "nearest", "NJ", "opti", "raw sim", "single", "swarm", "UPGMA", "WPGMA", "clusterize"))

ordered_methods_full_names = c("Intra/inter-BIN analyses\n(pairwise sim.)", "SWARM (pairwise sim.)", "DECIPHER (pairwise sim.)\nClusterize", "DECIPHER (pairwise sim.)\nNeighbor joining", 
  "DECIPHER (pairwise sim.)\nNearest neighbor", "DECIPHER (pairwise sim.)\nFurthest neighbor", "DECIPHER (pairwise sim.)\nUPGMA", 
  "DECIPHER (pairwise sim.)\nWPGMA", "MOTHUR (global sim.)\nOPTICLUST", "MOTHUR (global sim.)\nVSEARCH", "MOTHUR (global sim.)\nNearest neighbor", 
  "MOTHUR (global sim.)\nFurthest neighbor", "MOTHUR (global sim.)\nUPGMA")

all_methods_grouped_filtered$METHOD = methods_full_names[all_methods_grouped_filtered$METHOD]
all_methods_grouped_filtered$METHOD = factor(all_methods_grouped_filtered$METHOD, levels = ordered_methods_full_names)
all_methods_grouped_filtered


## Saving a summary of the minimal global error rate and its associated similarity threshold

all_methods_grouped_filtered_summary_clustering = all_methods_grouped_filtered %>% filter(METHOD != "Raw pairwise similarity") %>% group_by(PRIMERS) %>% 
  filter(CUMULATED_PERCENT_ERRORS_PER_ANALYSIS == min(CUMULATED_PERCENT_ERRORS_PER_ANALYSIS, na.rm = T)) %>%
  mutate(PROGRAM = word(METHOD, 1), SIMILARITY_TYPE = ifelse(grepl("pairwise", METHOD), "Pairwise"), METHOD = ifelse(grepl("Raw", METHOD), "Direct comparison", word(METHOD, 2, sep = fixed("\n")))) %>%
  select(c(PRIMERS, GENE, PROGRAM, SIMILARITY_TYPE, METHOD, SIMILARITY, PERCENT_INTRA_ERRORS, PERCENT_INTER_ERRORS, CUMULATED_PERCENT_ERRORS_PER_ANALYSIS)) %>%
  arrange(CUMULATED_PERCENT_ERRORS_PER_ANALYSIS) %>% mutate(across(where(is.numeric), ~ round(., 2)))
all_methods_grouped_filtered_summary_clustering

write.csv(all_methods_grouped_filtered_summary_clustering, "Output/comparison_clustering/Best clustering methods per metabarcode.csv", row.names = F)


## Plotting the number of false-positives per clustering methods and similarity thresholds

ggplot(all_methods_grouped_filtered, aes(x = as.numeric(SIMILARITY), y = as.numeric(PERCENT_INTRA_ERRORS),
                                                      color = METHOD, fill = METHOD, group = METHOD)) + 
  facet_wrap(~PRIMERS_GENE, scales = "free_y") + geom_point(shape = 21) + geom_line() + theme_() +
  scale_color_manual(values = scale_color_distinct, guide = "legend") +
  scale_fill_manual(values = scale_color_distinct) + 
  scale_x_continuous(breaks = c(90, 93, 96, 99), labels = function(x) paste0(x, '%')) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))), 
                     minor_breaks = seq(1, 50, 1), labels = function(x) paste0(x, '%')) +
  theme(legend.title = element_text(family = "Segoe UI Semibold"), 
        strip.text = element_text(family = "Segoe UI Semibold"), legend.key.spacing.y = unit(0.3, "cm"),
        panel.grid.major.y = element_line(size = 0.1, color = "darkgrey"),
        panel.grid.minor.y = element_line(size = 0.1, color = "darkgrey")) +
  labs(x = "Clustering threshold", y = "Percentage of over-splitting errors")
  
ggsave("Graph/comparison_clustering/Comparison clustering oversplitting - All methods.png", 
       device = png, dpi = 600, width = 25.53229, height = 14.39333, units = "cm")


## Plotting the number of false-negatives per clustering methods and similarity thresholds

ggplot(all_methods_grouped_filtered, aes(x = as.numeric(SIMILARITY), y = as.numeric(PERCENT_INTER_ERRORS),
                                color = METHOD, fill = METHOD, group = METHOD)) + 
  facet_wrap(~PRIMERS_GENE, scales = "free_y") + geom_point(shape = 21) + geom_line() + theme_() +
  scale_color_manual(values = scale_color_distinct, guide = "legend") +
  scale_fill_manual(values = scale_color_distinct) + 
  scale_x_continuous(breaks = c(90, 93, 96, 99), labels = function(x) paste0(x, '%')) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 5) * 1.1)))), 
                     minor_breaks = seq(5, 85, 5), labels = function(x) paste0(x, '%')) +
  theme(legend.title = element_text(family = "Segoe UI Semibold"), 
        strip.text = element_text(family = "Segoe UI Semibold"), legend.key.spacing.y = unit(0.3, "cm"),
        panel.grid.major.y = element_line(size = 0.1, color = "darkgrey"),
        panel.grid.minor.y = element_line(size = 0.1, color = "darkgrey")) +
  labs(x = "Clustering threshold", y = "Percentage of over-merging errors")

ggsave("Graph/comparison_clustering/Comparison clustering overmerging - All methods.png", 
       device = png, dpi = 600, width = 25.53229, height = 14.39333, units = "cm")


## Plotting the sum of false-positives and false-negatives per clustering methods and similarity thresholds

ggplot(all_methods_grouped_filtered, aes(x = as.numeric(SIMILARITY), y = as.numeric(CUMULATED_PERCENT_ERRORS_PER_ANALYSIS),
                                         color = METHOD, group = METHOD)) + 
  facet_wrap(~PRIMERS_GENE, scales = "free_y") + geom_point() + geom_line() + theme_() +
  geom_point(data = subset(all_methods_grouped_filtered, BEST_THRESHOLD), shape = 21, fill = "white", stroke = 1) +
  scale_color_manual(values = scale_color_distinct, guide = "legend") +
  scale_fill_manual(values = scale_color_distinct) + 
  scale_x_continuous(breaks = c(90, 93, 96, 99), labels = function(x) paste0(x, '%')) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 5) * 1.1)))), 
                     minor_breaks = seq(5, 85, 5), labels = function(x) paste0(x, '%')) +
  theme(legend.title = element_text(family = "Segoe UI Semibold"), 
        strip.text = element_text(family = "Segoe UI Semibold"), legend.key.spacing.y = unit(0.3, "cm"),
        panel.grid.major.y = element_line(size = 0.1, color = "darkgrey"),
        panel.grid.minor.y = element_line(size = 0.1, color = "darkgrey")) +
  labs(x = "Clustering threshold", y = "Sum of over-splitting and over-merging percentages")

ggsave("Graph/comparison_clustering/Comparison clustering cumulated overmerging and oversplitting errors - All methods.png", 
       device = png, dpi = 600, width = 25.53229, height = 14.39333, units = "cm")

# ggsave("Graph/comparison_clustering/Comparison clustering cumulated overmerging and oversplitting errors - All methods.pdf", 
#        device = cairo_pdf, width = 25.53229, height = 14.39333, units = "cm")

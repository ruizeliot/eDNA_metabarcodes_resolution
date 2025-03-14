##### INITIALISATION #####

## Loading required packages and custom functions

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(ggplot2)
library(extrafont)
library(tidyverse)
library(data.table)
library(patchwork)


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


## Function for adding a strip appearance to a text in ggplot

element_textbox = function(...) {
  el = element_text(...)
  class(el) = c("element_textbox", class(el))
  el
}

element_grob.element_textbox = function(element, ...) {
  text_grob = NextMethod()
  rect_grob = element_grob(calc_element("strip.background", theme_bw()))
  ggplot2:::absoluteGrob(grid::gList(element_grob(calc_element("strip.background", theme_bw())), text_grob),
                         height = grid::grobHeight(text_grob), width = grid::unit(1, "npc"))
}


## Loading all needed files

act_mitogenomes_bin_corrected = tibble(read.csv("Data/taxonomy_BIN_infos/Actinopterygii - Taxonomy - BIN.csv"))
act_mitogenomes_bin_corrected

filtered_metabarcodes_final = tibble(read.csv("Data/metabarcodes/Filtered metabarcodes final.csv"))
filtered_metabarcodes_final

bin_database_corrected = tibble(read.csv("Data/BIN/BIN database - Corrected.csv"))
bin_database_corrected

act_mitogenomes_bin_full_infos = tibble(merge(act_mitogenomes_bin_corrected, bin_database_corrected))
act_mitogenomes_bin_full_infos

mean_length_metabarcodes = tibble(read.csv("Data/primers_infos/Mean length data.csv"))
mean_length_metabarcodes %>% print(n = 22)

primers = sort(unique(filtered_metabarcodes_final$PRIMERS))
primers






##### INTRA-BIN ANALYSIS #####

## Comparing all sequences sharing same BIN

start1 = Sys.time()
intra_bin_comparisons = intra_taxa_similarity(metabarcodes_data = filtered_metabarcodes_final, 
                                              infos_data = act_mitogenomes_bin_full_infos, 
                                              name_col_taxa = "BIN", 
                                              name_col_primer = "PRIMERS", 
                                              name_col_accession = "ACCESSION", 
                                              name_col_fragment = "FRAGMENT", 
                                              keep_col_infos_data = c("BIN", 
                                                                     "SPECIES", "BOLD_SPECIES", 
                                                                     "GENUS", "BOLD_GENUS", 
                                                                     "FAMILY", "BOLD_FAMILY", 
                                                                     "ORDER", "BOLD_ORDER"), 
                                              order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                                                   "BIN1", "BIN2", 
                                                                   "SPECIES1", "BOLD_SPECIES1", 
                                                                   "SPECIES2", "BOLD_SPECIES2", 
                                                                   "GENUS1", "BOLD_GENUS1", 
                                                                   "GENUS2", "BOLD_GENUS2", 
                                                                   "FAMILY1", "BOLD_FAMILY1", 
                                                                   "FAMILY2", "BOLD_FAMILY2", 
                                                                   "ORDER1", "BOLD_ORDER1", 
                                                                   "ORDER2", "BOLD_ORDER2"),
                                              # Change it with the path + name of your VSEARCH program
                                              program_name = "vsearch-2.17.1-win-x86_64.exe")
end1 = Sys.time()
difftime(end1, start1) # Duration: 29mn

saveRDS(intra_bin_comparisons, "Output/intra_BIN/Intra-BIN - Raw comparisons - All primers.rds") # Duration: 30s
intra_bin_comparisons = readRDS("Output/intra_BIN/Intra-BIN - Raw comparisons - All primers.rds") # Duration: 15s
intra_bin_comparisons[1:2]


## Running all analysis for all thresholds of similarity

start2 = Sys.time()
intra_bin_analysis_all_similarities = intra_taxa_analysis(intra_bin_comparisons, act_mitogenomes_bin_full_infos, 
                                                          similarity_threshold = 90:99, name_col_similarity = "SIMILARITY", 
                                                          name_col_accession = "ACCESSION", name_col_taxa = "BIN")
end2 = Sys.time()
difftime(end2, start2) # Duration: 5mn


## Saving and reading the summary of the similarity per BIN

intra_bin_similarity = intra_bin_analysis_all_similarities$DETAILS_SIMILARITY

write.csv(intra_bin_similarity, "Output/intra_BIN/Intra-BIN analysis - Similarity.csv", row.names = F)
intra_bin_similarity = tibble(read.csv("Output/intra_BIN/Intra-BIN analysis - Similarity.csv"))
intra_bin_similarity


## Saving and reading all the results of the summary of the intra-BIN analysis for all thresholds

intra_bin_summary = intra_bin_analysis_all_similarities$SUMMARY

saveRDS(intra_bin_summary, "Output/intra_BIN/Intra-BIN analysis - Summary - All primers.rds")
intra_bin_summary = readRDS("Output/intra_BIN/Intra-BIN analysis - Summary - All primers.rds")
intra_bin_summary[1:2]


## Saving and reading all the results of the details of the intra-BIN analysis for all thresholds

intra_bin_details = intra_bin_analysis_all_similarities$DETAILS_PER_TAXA

saveRDS(intra_bin_details, "Output/intra_BIN/Intra-BIN analysis - Details - All primers.rds")
intra_bin_details = readRDS("Output/intra_BIN/Intra-BIN analysis - Details - All primers.rds")
intra_bin_details[1:2]





##### INTER-BIN ANALYSIS #####

## Runing comparisons between all sequences in the database 

start1 = Sys.time()
similarity_comparisons_list = vsearch_pairwise_similarity(metabarcodes_data = filtered_metabarcodes_final, 
                                                          infos_data = act_mitogenomes_bin_full_infos, 
                                                          min_similarity_threshold = 90, 
                                                          max_similarity_threshold = 100,
                                                          name_col_primer = "PRIMERS", 
                                                          name_col_accession = "ACCESSION", 
                                                          name_col_fragment = "FRAGMENT", 
                                                          keep_col_infos_data = c("BIN", 
                                                                                 "SPECIES", "BOLD_SPECIES", 
                                                                                 "GENUS", "BOLD_GENUS", 
                                                                                 "FAMILY", "BOLD_FAMILY", 
                                                                                 "ORDER", "BOLD_ORDER"), 
                                                          order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                                                               "BIN1", "BIN2", 
                                                                               "SPECIES1", "BOLD_SPECIES1", 
                                                                               "SPECIES2", "BOLD_SPECIES2", 
                                                                               "GENUS1", "BOLD_GENUS1", 
                                                                               "GENUS2", "BOLD_GENUS2", 
                                                                               "FAMILY1", "BOLD_FAMILY1", 
                                                                               "FAMILY2", "BOLD_FAMILY2", 
                                                                               "ORDER1", "BOLD_ORDER1", 
                                                                               "ORDER2", "BOLD_ORDER2"),
                                                          # Change it with the path + name of your VSEARCH program
                                                          program_name = "vsearch-2.17.1-win-x86_64.exe")
end1 = Sys.time()
difftime(end1, start1) # Duration: 59mn


## Filtering comparisons from different BIN to get the inter-BIN

inter_bin_comparisons = lapply(similarity_comparisons_list, function(x) x[which(x$BIN1 != x$BIN2),])


## Saving and reading all the comparisons retained

saveRDS(inter_bin_comparisons, "Output/inter_BIN/Inter-BIN - Raw comparisons - All primers.rds") # Duration: 1.5mn
inter_bin_comparisons = readRDS("Output/inter_BIN/Inter-BIN - Raw comparisons - All primers.rds") # Duration: 1mn
inter_bin_comparisons[1:2]


## Detailed and summarised inter-BIN analysis at all thresholds

start4 = Sys.time()
inter_bin_analysis_all_similarities = inter_taxa_analysis(inter_bin_comparisons, act_mitogenomes_bin_full_infos, 
                                                          similarity_threshold = 90:99, name_col_similarity = "SIMILARITY", 
                                                          name_col_accession = "ACCESSION", name_col_taxa = "BIN")
end4 = Sys.time()
difftime(end4, start4) # Duration: 36mn


## Saving and reading all the summary of the inter-BIN analysis for all thresholds

inter_bin_summary = inter_bin_analysis_all_similarities$SUMMARY

saveRDS(inter_bin_summary, "Output/inter_BIN/Inter-BIN analysis - Summary - All primers.rds")
inter_bin_summary = readRDS("Output/inter_BIN/Inter-BIN analysis - Summary - All primers.rds")
inter_bin_summary[1:2]


## Saving and reading all the details of the inter-BIN analysis for all thresholds

inter_bin_details = inter_bin_analysis_all_similarities$DETAILS_UNDISCRIMINED_PAIRS

saveRDS(inter_bin_details, "Output/inter_BIN/Inter-BIN analysis - Details - All primers.rds")
inter_bin_details = readRDS("Output/inter_BIN/Inter-BIN analysis - Details - All primers.rds")
inter_bin_details[1:2]





##### PLOTTING DETAILS ON INTRA-BIN RESOLUTION #####

## Plotting of details for intra-BIN resolution at 97%

intra_bin_graph_97 = intra_bin_details$`97`

wrongly_discrimined_97 = sapply(split(intra_bin_graph_97, intra_bin_graph_97$PRIMERS), function(x)
  sum(x$NUMBER_INTRA_ERRORS, na.rm = T) / sum(x$TOTAL_NUMBER_COMPARISONS, na.rm = T) * 100)
wrongly_discrimined_97 = tibble(data.frame(PRIMERS = names(wrongly_discrimined_97),
                                           PERCENT_INTRA_BIN_ERRORS = wrongly_discrimined_97))
wrongly_discrimined_97 = wrongly_discrimined_97[order(wrongly_discrimined_97$PERCENT_INTRA_BIN_ERRORS), ]

wrongly_discrimined_97$PERCENT_INTRA_BIN_ERRORS_FACTOR = paste0(round(wrongly_discrimined_97$PERCENT_INTRA_BIN_ERRORS, 2), "%")

intra_bin_graph_97 = tibble(merge(wrongly_discrimined_97, intra_bin_graph_97))
intra_bin_graph_97$PERCENT_INTRA_BIN_ERRORS_FACTOR = factor(intra_bin_graph_97$PERCENT_INTRA_BIN_ERRORS_FACTOR,
                                                               levels = unique(wrongly_discrimined_97$PERCENT_INTRA_BIN_ERRORS_FACTOR))
intra_bin_graph_97

intra_bin_graph_97 = subset(intra_bin_graph_97, !(PRIMERS %in% c("L14735c2", "Fish2deg")))

intra_bin_graph_97$PRIMERS = replace(intra_bin_graph_97$PRIMERS, 
                                          which(intra_bin_graph_97$PRIMERS %in% "Fish2b"), "Fish2b/deg")

intra_bin_graph_97$PRIMERS = replace(intra_bin_graph_97$PRIMERS, 
                                          which(intra_bin_graph_97$PRIMERS %in% "L14735c"), "L14735c/c2")

intra_bin_graph_97$PRIMERS = replace(intra_bin_graph_97$PRIMERS, 
                                          which(intra_bin_graph_97$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")

ggplot(subset(intra_bin_graph_97, PERCENT_INTRA_ERRORS != 0 & PRIMERS != "L14735c2"), 
       aes(x = PRIMERS, y = PERCENT_INTRA_ERRORS, fill = PERCENT_INTRA_BIN_ERRORS)) +
  facet_grid(. ~ PERCENT_INTRA_BIN_ERRORS_FACTOR, drop = T, scales = "free") +
  scale_y_continuous(labels = function(x) paste0(x, '%')) + geom_violin(adjust = 1/3, color = "red", size = 0.2) +
  scale_fill_gradient2(mid = "#fff5f0", high = "#cb181d") +
  geom_hline(aes(yintercept = 0), color = "#6BD400", size = 1) +
  theme_() + theme(axis.title.x = element_blank(), legend.position = "none", strip.text = element_text(size = 7.5, family = "Segoe UI Semibold"),
                   plot.subtitle = element_textbox(hjust = 0.5, family = "Segoe UI Semibold", margin = margin(t = 5, b = 3)),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(y = "Percentage of sequences wrongly discrimined per BIN",
       subtitle = "Total percentage of intra-BIN errors at 97%")

ggsave("Graph/intra_BIN/Intra-BIN 97 percent details violin.png", device = png, dpi = 600, type = "cairo",
       width = 25.53229, height = 14.39333, units = "cm")


## Plotting of details for intra-BIN resolution at 99%

intra_bin_graph_99 = intra_bin_details$`99`

wrongly_discrimined_99 = sapply(split(intra_bin_graph_99, intra_bin_graph_99$PRIMERS), function(x)
  sum(x$NUMBER_INTRA_ERRORS, na.rm = T) / sum(x$TOTAL_NUMBER_COMPARISONS, na.rm = T) * 100)
wrongly_discrimined_99 = tibble(data.frame(PRIMERS = names(wrongly_discrimined_99),
                                           PERCENT_INTRA_BIN_ERRORS = wrongly_discrimined_99))
wrongly_discrimined_99 = wrongly_discrimined_99[order(wrongly_discrimined_99$PERCENT_INTRA_BIN_ERRORS), ]

wrongly_discrimined_99$PERCENT_INTRA_BIN_ERRORS_FACTOR = paste0(round(wrongly_discrimined_99$PERCENT_INTRA_BIN_ERRORS, 2), "%")

intra_bin_graph_99 = tibble(merge(wrongly_discrimined_99, intra_bin_graph_99))
intra_bin_graph_99$PERCENT_INTRA_BIN_ERRORS_FACTOR = factor(intra_bin_graph_99$PERCENT_INTRA_BIN_ERRORS_FACTOR,
                                                            levels = unique(wrongly_discrimined_99$PERCENT_INTRA_BIN_ERRORS_FACTOR))
intra_bin_graph_99

intra_bin_graph_99 = subset(intra_bin_graph_99, !(PRIMERS %in% c("L14735c2", "Fish2deg")))

intra_bin_graph_99$PRIMERS = replace(intra_bin_graph_99$PRIMERS, 
                                     which(intra_bin_graph_99$PRIMERS %in% "Fish2b"), "Fish2b/deg")

intra_bin_graph_99$PRIMERS = replace(intra_bin_graph_99$PRIMERS, 
                                     which(intra_bin_graph_99$PRIMERS %in% "L14735c"), "L14735c/c2")

intra_bin_graph_99$PRIMERS = replace(intra_bin_graph_99$PRIMERS, 
                                     which(intra_bin_graph_99$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")

ggplot(subset(intra_bin_graph_99, PERCENT_INTRA_ERRORS != 0 & PRIMERS != "L14735c2"), 
       aes(x = PRIMERS, y = PERCENT_INTRA_ERRORS, fill = PERCENT_INTRA_BIN_ERRORS)) +
  facet_grid(. ~ PERCENT_INTRA_BIN_ERRORS_FACTOR, drop = T, scales = "free") +
  scale_y_continuous(labels = function(x) paste0(x, '%')) + geom_violin(adjust = 1/3, color = "red", size = 0.2) +
  scale_fill_gradient2(mid = "#fff5f0", high = "#cb181d") +
  geom_hline(aes(yintercept = 0), color = "#6BD400", size = 1) +
  theme_() + theme(axis.title.x = element_blank(), legend.position = "none", strip.text = element_text(size = 7.5, family = "Segoe UI Semibold"),
                   plot.subtitle = element_textbox(hjust = 0.5, family = "Segoe UI Semibold", margin = margin(t = 5, b = 3)),
                   axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(y = "Percentage of sequences wrongly discrimined per BIN",
       subtitle = "Total percentage of intra-BIN errors at 99%")

ggsave("Graph/intra_BIN/Intra-BIN 99 percent details violin.png", device = png, dpi = 600, type = "cairo",
       width = 25.53229, height = 14.39333, units = "cm")




##### PLOTTING THE RELATIVE IMPORTANCE OF INTER-BIN AND INTRA-BIN ERRORS DEPENDING ON THE PRIMER #####

## Taxonomic resolution analysis summary at 97%

resolution_analysis_97 = tibble(cbind(intra_bin_summary$`97`, inter_bin_summary$`97`[,-1]))

taxonomic_resolution_summary_97 = subset(resolution_analysis_97, select = c("PRIMERS", "NUMBER_INTRA_ERRORS", "PERCENT_INTRA_ERRORS", 
                                                                            "NUMBER_INTER_ERRORS", "PERCENT_INTER_ERRORS"))

taxonomic_resolution_summary_97$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = 
  taxonomic_resolution_summary_97$PERCENT_INTRA_ERRORS +
  taxonomic_resolution_summary_97$PERCENT_INTER_ERRORS

taxonomic_resolution_summary_97 = tibble(merge(mean_length_metabarcodes[,1:3], taxonomic_resolution_summary_97))

taxonomic_resolution_summary_97 = taxonomic_resolution_summary_97[order(taxonomic_resolution_summary_97$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS),]

taxonomic_resolution_summary_97 = subset(taxonomic_resolution_summary_97, !(PRIMERS %in% c("L14735c2", "Fish2deg")))

taxonomic_resolution_summary_97$PRIMERS = replace(taxonomic_resolution_summary_97$PRIMERS, 
                                                  which(taxonomic_resolution_summary_97$PRIMERS %in% "Fish2b"), "Fish2b/deg")

taxonomic_resolution_summary_97$PRIMERS = replace(taxonomic_resolution_summary_97$PRIMERS, 
                                                  which(taxonomic_resolution_summary_97$PRIMERS %in% "L14735c"), "L14735c/c2")

taxonomic_resolution_summary_97$PRIMERS = replace(taxonomic_resolution_summary_97$PRIMERS, 
                                                  which(taxonomic_resolution_summary_97$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")

taxonomic_resolution_summary_97$PRIMERS = factor(taxonomic_resolution_summary_97$PRIMERS, levels = taxonomic_resolution_summary_97$PRIMERS)
taxonomic_resolution_summary_97 %>% print(n = 22)

relative_misidentified_97 = ggplot(taxonomic_resolution_summary_97) + 
  geom_bar(aes(x = PRIMERS, y = CUMULATED_PERCENT_ERRORS_PER_ANALYSIS, fill = "Inter-BIN"), stat = "identity") + 
  geom_bar(aes(x = PRIMERS, y = PERCENT_INTRA_ERRORS, fill = "Intra-BIN"), stat = "identity") + 
  theme_() + labs(x = "", y = "Cumulated percentage of misidentified cases in each analysis") +
  scale_y_continuous(labels = function(x) paste0(x, '%')) +
  scale_fill_manual(values = c("#41CDFF", "#FF4141")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = unit(c(0.05,0.05,-0.3,0.05), "cm"), size = 12), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 14),
        legend.title = element_blank(), legend.position = c(0.12,0.915), legend.text = element_text(size = 12)) 

relative_info_taxo_97 = ggplot(taxonomic_resolution_summary_97, aes(x = PRIMERS, y = MEAN_LENGTH, group = 1)) + 
  geom_line(size = 0.7, color = "black") +
  geom_point(color = "#41AAFF", size = 4) + geom_point(color = "black") + 
  theme_minimal(base_family = "Segoe UI Semilight") + ylab("Mean length") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(margin = unit(c(0,0.3,0,0), "cm"), size = 14), 
        axis.text.y = element_text(size = 12), legend.position = "top", legend.text = element_text(size = 12),
        legend.title = element_text(family = "Segoe UI Semibold", size = 14)) +
  coord_fixed(ratio = 1 / 100) +
  geom_col(aes(x = PRIMERS, y = -110), fill = "white", width = 1) +
  geom_hline(aes(yintercept = 0), color = "#EBEBEB", size = 1.1) +
  geom_col(aes(x = PRIMERS, y = -100, fill = GENE), width = 0.9) +
  geom_col(aes(x = PRIMERS, y = -20), fill = "white", width = 1) +
  scale_fill_manual(values = c("#009E73", "#F5CF00", "#D55E00", "#B600FF")) + labs(fill = "REGION")

(relative_info_taxo_97 / relative_misidentified_97) 


## Taxonomic resolution analysis summary at 99%

resolution_analysis_99 = tibble(cbind(intra_bin_summary$`99`, inter_bin_summary$`99`[,-1]))

taxonomic_resolution_summary_99 = subset(resolution_analysis_99, select = c("PRIMERS", "NUMBER_INTRA_ERRORS", "PERCENT_INTRA_ERRORS", 
                                                                            "NUMBER_INTER_ERRORS", "PERCENT_INTER_ERRORS"))

taxonomic_resolution_summary_99$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = taxonomic_resolution_summary_99$PERCENT_INTRA_ERRORS +
  taxonomic_resolution_summary_99$PERCENT_INTER_ERRORS

taxonomic_resolution_summary_99 = tibble(merge(mean_length_metabarcodes[,1:3], taxonomic_resolution_summary_99))

taxonomic_resolution_summary_99 = taxonomic_resolution_summary_99[order(taxonomic_resolution_summary_99$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS),]

taxonomic_resolution_summary_99 = subset(taxonomic_resolution_summary_99, !(PRIMERS %in% c("L14735c2", "Fish2deg")))

taxonomic_resolution_summary_99$PRIMERS = replace(taxonomic_resolution_summary_99$PRIMERS, 
                                                  which(taxonomic_resolution_summary_99$PRIMERS %in% "Fish2b"), "Fish2b/deg")

taxonomic_resolution_summary_99$PRIMERS = replace(taxonomic_resolution_summary_99$PRIMERS, 
                                                  which(taxonomic_resolution_summary_99$PRIMERS %in% "L14735c"), "L14735c/c2")

taxonomic_resolution_summary_99$PRIMERS = replace(taxonomic_resolution_summary_99$PRIMERS, 
                                                  which(taxonomic_resolution_summary_99$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")

taxonomic_resolution_summary_99$PRIMERS = factor(taxonomic_resolution_summary_99$PRIMERS, levels = taxonomic_resolution_summary_99$PRIMERS)
taxonomic_resolution_summary_99 %>% print(n = 22)

relative_misidentified_99 = ggplot(taxonomic_resolution_summary_99) + 
  geom_bar(aes(x = PRIMERS, y = CUMULATED_PERCENT_ERRORS_PER_ANALYSIS, fill = "Inter-BIN"), stat = "identity") + 
  geom_bar(aes(x = PRIMERS, y = PERCENT_INTRA_ERRORS, fill = "Intra-BIN"), stat = "identity") + 
  theme_() + labs(x = "", y = "Cumulated percentage of misidentified cases in each analysis") +
  scale_y_continuous(labels = function(x) paste0(x, '%')) +
  scale_fill_manual(values = c("#41CDFF", "#FF4141")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = unit(c(0.05,0.05,-0.3,0.05), "cm"), size = 12), 
        axis.text.y = element_text(size = 12), axis.title.y = element_blank(), legend.text = element_text(size = 12),
        legend.title = element_blank(), legend.position = c(0.12,0.915)) 

relative_info_taxo_99 = ggplot(taxonomic_resolution_summary_99, aes(x = PRIMERS, y = MEAN_LENGTH, group = 1)) + 
  geom_line(size = 0.7, color = "black") +
  geom_point(color = "#41AAFF", size = 4) + geom_point(color = "black") + 
  theme_minimal(base_family = "Segoe UI Semilight") + ylab("Mean length") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12), legend.position = "top", legend.text = element_text(size = 12),
        legend.title = element_text(family = "Segoe UI Semibold", size = 14)) +
  coord_fixed(ratio = 1 / 100) +
  geom_col(aes(x = PRIMERS, y = -110), fill = "white", width = 1) +
  geom_hline(aes(yintercept = 0), color = "#EBEBEB", size = 1.1) +
  geom_col(aes(x = PRIMERS, y = -100, fill = GENE), width = 0.9) +
  geom_col(aes(x = PRIMERS, y = -20), fill = "white", width = 1) +
  scale_fill_manual(values = c("#009E73", "#F5CF00", "#D55E00", "#B600FF")) + labs(fill = "REGION")

(relative_info_taxo_99 / relative_misidentified_99) 


## Saving the last two graphs together

(relative_info_taxo_97 / relative_misidentified_97) | (relative_info_taxo_99 / relative_misidentified_99) 

ggsave("Graph/summary_resolution_analysis/Summary taxonomic resolution 97 & 99.png", dpi = 600, device = png,
       width = 30, height = 22.38375, units = "cm")





##### PLOTTING ALL RESULTS OF THE INTRA-BIN AND INTER-BIN ANALYSES #####

## Assembling summary tables together

intra_bin_summary_df = bind_rows(intra_bin_summary, .id = "SIMILARITY")
intra_bin_summary_df

inter_bin_summary_df = bind_rows(inter_bin_summary, .id = "SIMILARITY")
inter_bin_summary_df

intra_bin_summary_df = intra_bin_summary_df[order(intra_bin_summary_df$SIMILARITY),]
intra_bin_summary_df

inter_bin_summary_df = inter_bin_summary_df[order(inter_bin_summary_df$SIMILARITY),]
inter_bin_summary_df

intra_inter_bin_summary = tibble(cbind(intra_bin_summary_df, inter_bin_summary_df[,-c(1:2)]))

intra_inter_bin_summary = subset(intra_inter_bin_summary, select = 
                                   c("SIMILARITY", "PRIMERS", "NUMBER_INTRA_ERRORS", "PERCENT_INTRA_ERRORS", 
                                     "NUMBER_INTER_ERRORS", "PERCENT_INTER_ERRORS"))
intra_inter_bin_summary

intra_inter_bin_summary$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = 
  intra_inter_bin_summary$PERCENT_INTRA_ERRORS +
  intra_inter_bin_summary$PERCENT_INTER_ERRORS

intra_inter_bin_summary = tibble(merge(mean_length_metabarcodes[,1:3], intra_inter_bin_summary))


## Assembling the same metabarcodes with different primers together

intra_inter_bin_summary_subset = subset(intra_inter_bin_summary, !(PRIMERS %in% c("Fish2deg", "L14735c2")))

intra_inter_bin_summary_subset$PRIMERS = replace(intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(intra_inter_bin_summary_subset$PRIMERS %in% "Fish2b"), "Fish2b/deg")

intra_inter_bin_summary_subset$PRIMERS = replace(intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(intra_inter_bin_summary_subset$PRIMERS %in% "L14735c"), "L14735c/c2")

intra_inter_bin_summary_subset$PRIMERS = replace(intra_inter_bin_summary_subset$PRIMERS, 
                                                 which(intra_inter_bin_summary_subset$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")


## Computing the minimum error threshold

best_errors_per_primer = sapply(split(intra_inter_bin_summary_subset, intra_inter_bin_summary_subset$PRIMERS), function(x) 
  min(x$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS))
best_errors_per_primer = sort(best_errors_per_primer)
best_errors_per_primer

best_gene_per_primer = sapply(split(intra_inter_bin_summary_subset, intra_inter_bin_summary_subset$PRIMERS), function(x) 
  x[which.min(x$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS),]$GENE)
best_gene_per_primer = best_gene_per_primer[match(names(best_errors_per_primer), names(best_gene_per_primer))]
best_gene_per_primer

best_similarity_per_primer = sapply(split(intra_inter_bin_summary_subset, intra_inter_bin_summary_subset$PRIMERS), function(x) 
  x[which.min(x$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS),]$SIMILARITY)
best_similarity_per_primer = best_similarity_per_primer[match(names(best_errors_per_primer), names(best_similarity_per_primer))]
best_similarity_per_primer


## Creating the table with the star indicating the best similarity threshold

best_similarity_per_primer_graph = tibble(cbind(subset(intra_inter_bin_summary_subset, select = c("PRIMERS", "SIMILARITY")), 
                                                tibble(BEST_THRESHOLD = "")))
best_similarity_per_primer_graph

for(i in 1:length(best_similarity_per_primer)){
  
  best_similarity_per_primer_graph$BEST_THRESHOLD = 
    replace(best_similarity_per_primer_graph$BEST_THRESHOLD,
            which(best_similarity_per_primer_graph$PRIMERS == names(best_similarity_per_primer)[i] & 
                    best_similarity_per_primer_graph$SIMILARITY == best_similarity_per_primer[i]), "*")
  
}

best_similarity_per_primer_graph


## Changing the facet names to indicate the minimum cumulated errors

best_errors_per_primer_name = setNames(paste0(names(best_errors_per_primer), " (", best_gene_per_primer, "): ", 
                                              round(best_errors_per_primer, 2), "%"), names(best_errors_per_primer))
best_errors_per_primer_name


## Creating a table with all features for plotting

intra_inter_bin_summary_graph = tibble(merge(best_similarity_per_primer_graph, intra_inter_bin_summary_subset))

intra_inter_bin_summary_graph$SIMILARITY = as.numeric(intra_inter_bin_summary_graph$SIMILARITY)

intra_inter_bin_summary_graph$PRIMERS = factor(intra_inter_bin_summary_graph$PRIMERS, levels = names(best_errors_per_primer))
intra_inter_bin_summary_graph


## Creating and saving the overall plot

ggplot(intra_inter_bin_summary_graph) +
  facet_wrap(~PRIMERS, scales = "free_y", nrow = 4, ncol = 5, labeller = as_labeller(best_errors_per_primer_name)) +
  geom_bar(aes(x = SIMILARITY, y = CUMULATED_PERCENT_ERRORS_PER_ANALYSIS, fill = "Inter-BIN"), stat = "identity") + 
  geom_bar(aes(x = SIMILARITY, y = PERCENT_INTRA_ERRORS, fill = "Intra-BIN"), stat = "identity") + 
  labs(x = "", y = "Stacked percentage of misidentified cases in each type of analysis") +
  scale_x_continuous(breaks = 90:99, labels = function(x) paste0(x, '%')) +
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(0, (max(x) + 1) * 1.1)))), 
                     minor_breaks = seq(1, 50, 1), labels = function(x) paste0(x, '%')) +
  scale_fill_manual(values = c("#41CDFF", "#FF4141")) +
  theme_() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), legend.title = element_blank(), 
                   legend.position = "none", strip.text = element_text(family = "Segoe UI Semibold"),
                   panel.grid.major.y = element_line(size = 0.1, color = "darkgrey"),
                   panel.grid.minor.y = element_line(size = 0.1, color = "darkgrey")) +
  geom_text(aes(label = BEST_THRESHOLD, x = SIMILARITY, y = CUMULATED_PERCENT_ERRORS_PER_ANALYSIS), 
            position = position_dodge(width = 0.9), vjust = 0.2, size = 7, family = "Segoe UI Semibold")

ggsave("Graph/summary_resolution_analysis/Summary taxonomic resolution - All thresholds.png", 
       device = png, dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, units = "cm")


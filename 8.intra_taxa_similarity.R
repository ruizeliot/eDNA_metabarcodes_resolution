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





##### INTRA-BIN SUMMARY #####

## Loading intra-BIN comparisons

intra_bin_comparisons = readRDS("Output/intra_BIN/Intra-BIN - Raw comparisons - All primers.rds") # Duration: 15s
intra_bin_comparisons[1:2]


## Computing statistics

median_intra_bin = sapply(intra_bin_comparisons, function(x) sapply(split(x, x$BIN1), function(y) median(y$SIMILARITY)))
median_intra_bin_df = pivot_longer(tibble(cbind(tibble(BIN = rownames(median_intra_bin)), tibble(as.data.frame(median_intra_bin)))), 
                                   !BIN, "PRIMERS", values_to = "MEDIAN")

mad_intra_bin = sapply(intra_bin_comparisons, function(x) sapply(split(x, x$BIN1), function(y) mad(y$SIMILARITY)))
mad_intra_bin_df = pivot_longer(tibble(cbind(tibble(BIN = rownames(mad_intra_bin)), tibble(as.data.frame(mad_intra_bin)))), 
                                !BIN, "PRIMERS", values_to = "MAD")

mean_intra_bin = sapply(intra_bin_comparisons, function(x) sapply(split(x, x$BIN1), function(y) mean(y$SIMILARITY)))
mean_intra_bin_df = pivot_longer(tibble(cbind(tibble(BIN = rownames(mean_intra_bin)), tibble(as.data.frame(mean_intra_bin)))), 
                                 !BIN, "PRIMERS", values_to = "MEAN")

sd_intra_bin = sapply(intra_bin_comparisons, function(x) sapply(split(x, x$BIN1), function(y) sd(y$SIMILARITY)))
sd_intra_bin_df = pivot_longer(tibble(cbind(tibble(BIN = rownames(sd_intra_bin)), tibble(as.data.frame(sd_intra_bin)))), 
                               !BIN, "PRIMERS", values_to = "SD")

similarity_per_bin = tibble(cbind(median_intra_bin_df, mad_intra_bin_df[,3], mean_intra_bin_df[,3], sd_intra_bin_df[,3]))

write.csv(similarity_per_bin, "Output/intra_BIN/Similarity per BIN.csv", row.names = F)
similarity_per_bin = tibble(read.csv("Output/intra_BIN/Similarity per BIN.csv"))
similarity_per_bin





##### INTRA-GENUS ANALYSIS AND SUMMARY #####

## Comparing all sequences sharing same genus

start4 = Sys.time()
intra_genus_comparisons = intra_taxa_similarity(metabarcodes_data = filtered_metabarcodes_final, 
                                                infos_data = act_mitogenomes_bin_full_infos, 
                                                name_col_taxa = "GENUS", 
                                                name_col_lower_rank_taxa = "BIN",
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
end4 = Sys.time()
difftime(end4, start4) # Duration: 11mn


## Saving and the reading the output

saveRDS(intra_genus_comparisons, "Output/intra_genus/Intra-genus - Raw comparisons - All primers.rds")
intra_genus_comparisons = readRDS("Output/intra_genus/Intra-genus - Raw comparisons - All primers.rds")
intra_genus_comparisons[1:2]


## Computing statistics

median_intra_genus = sapply(intra_genus_comparisons, function(x) sapply(split(x, x$GENUS1), function(y) median(y$SIMILARITY)))
median_intra_genus_df = pivot_longer(tibble(cbind(tibble(GENUS = rownames(median_intra_genus)), tibble(as.data.frame(median_intra_genus)))), 
                                     !GENUS, "PRIMERS", values_to = "MEDIAN")

mad_intra_genus = sapply(intra_genus_comparisons, function(x) sapply(split(x, x$GENUS1), function(y) mad(y$SIMILARITY)))
mad_intra_genus_df = pivot_longer(tibble(cbind(tibble(GENUS = rownames(mad_intra_genus)), tibble(as.data.frame(mad_intra_genus)))), 
                                  !GENUS, "PRIMERS", values_to = "MAD")

mean_intra_genus = sapply(intra_genus_comparisons, function(x) sapply(split(x, x$GENUS1), function(y) mean(y$SIMILARITY)))
mean_intra_genus_df = pivot_longer(tibble(cbind(tibble(GENUS = rownames(mean_intra_genus)), tibble(as.data.frame(mean_intra_genus)))), 
                                   !GENUS, "PRIMERS", values_to = "MEAN")

sd_intra_genus = sapply(intra_genus_comparisons, function(x) sapply(split(x, x$GENUS1), function(y) sd(y$SIMILARITY)))
sd_intra_genus_df = pivot_longer(tibble(cbind(tibble(GENUS = rownames(sd_intra_genus)), tibble(as.data.frame(sd_intra_genus)))), 
                                 !GENUS, "PRIMERS", values_to = "SD")

similarity_per_genus = tibble(cbind(median_intra_genus_df, mad_intra_genus_df[,3], mean_intra_genus_df[,3], sd_intra_genus_df[,3]))

write.csv(similarity_per_genus, "Output/intra_genus/Similarity per genus.csv", row.names = F)
similarity_per_genus = tibble(read.csv("Output/intra_genus/Similarity per genus.csv"))
similarity_per_genus





##### INTRA-FAMILY ANALYSIS AND SUMMARY #####

## Comparing all sequences sharing same family

start7 = Sys.time()
intra_family_comparisons = intra_taxa_similarity(metabarcodes_data = filtered_metabarcodes_final, 
                                                 infos_data = act_mitogenomes_bin_full_infos, 
                                                 name_col_taxa = "FAMILY", 
                                                 name_col_lower_rank_taxa = "GENUS",
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
end7 = Sys.time()
difftime(end7, start7) # Duration: 17mn


## Saving and reading the output

saveRDS(intra_family_comparisons, "Output/intra_family/Intra-family - Raw comparisons - All primers.rds") # Duration: 1.1mn
intra_family_comparisons = readRDS("Output/intra_family/Intra-family - Raw comparisons - All primers.rds") # Duration: 1.1mn
intra_family_comparisons[1:2]


## Computing statistics

median_intra_family = sapply(intra_family_comparisons, function(x) sapply(split(x, x$FAMILY1), function(y) median(y$SIMILARITY)))
median_intra_family_df = pivot_longer(tibble(cbind(tibble(FAMILY = rownames(median_intra_family)), tibble(as.data.frame(median_intra_family)))), 
                                      !FAMILY, "PRIMERS", values_to = "MEDIAN")

mad_intra_family = sapply(intra_family_comparisons, function(x) sapply(split(x, x$FAMILY1), function(y) mad(y$SIMILARITY)))
mad_intra_family_df = pivot_longer(tibble(cbind(tibble(FAMILY = rownames(mad_intra_family)), tibble(as.data.frame(mad_intra_family)))), 
                                   !FAMILY, "PRIMERS", values_to = "MAD")

mean_intra_family = sapply(intra_family_comparisons, function(x) sapply(split(x, x$FAMILY1), function(y) mean(y$SIMILARITY)))
mean_intra_family_df = pivot_longer(tibble(cbind(tibble(FAMILY = rownames(mean_intra_family)), tibble(as.data.frame(mean_intra_family)))), 
                                    !FAMILY, "PRIMERS", values_to = "MEAN")

sd_intra_family = sapply(intra_family_comparisons, function(x) sapply(split(x, x$FAMILY1), function(y) sd(y$SIMILARITY)))
sd_intra_family_df = pivot_longer(tibble(cbind(tibble(FAMILY = rownames(sd_intra_family)), tibble(as.data.frame(sd_intra_family)))), 
                                  !FAMILY, "PRIMERS", values_to = "SD")

similarity_per_family = tibble(cbind(median_intra_family_df, mad_intra_family_df[,3], mean_intra_family_df[,3], sd_intra_family_df[,3]))

write.csv(similarity_per_family, "Output/intra_family/Similarity per family.csv", row.names = F)
similarity_per_family = tibble(read.csv("Output/intra_family/Similarity per family.csv"))
similarity_per_family





##### INTRA-ORDER ANALYSIS AND SUMMARY #####

## Performing the intra-order analysis primer by primer to free the memory at each step (saturation otherwise)

start9 = Sys.time()
for(i in 1:length(primers)){
  
  cat(paste0("Perfoming intra-order analysis for ", primers[i]))
  cat("\n-----------------------------------------------------\n\n")
  
  intra_order_comparisons = intra_taxa_similarity(subset(filtered_metabarcodes_final, PRIMERS == primers[i]),
                                                  infos_data = act_mitogenomes_bin_full_infos, 
                                                  name_col_taxa = "ORDER", 
                                                  name_col_lower_rank_taxa = "FAMILY",
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
  
  cat("\nSaving the file: ")
  write.csv(intra_order_comparisons[[primers[i]]], paste0("Output/intra_order/Intra-order - ", 
                                                          primers[i], ".csv"), row.names = F)
  cat("DONE")
  
  cat("\nCleaning the memory: ") ; gc() ; cat("DONE")
  cat("\n\n-----------------------------------------------------\n\n")
  
}
end9 = Sys.time()
difftime(end9, start9) # Duration: 45mn


## Reading the output

start12 = Sys.time()
intra_order_comparisons = list()
for(i in 1:length(primers)){ cat(paste0(i, "/", length(primers), ": "))
  intra_order_comparisons[[i]] = tibble(read.csv(paste0("Output/intra_order/Intra-order - ",
                                                        primers[i], ".csv")))
  names(intra_order_comparisons)[i] = primers[i] ; cat("DONE\n")}
end12 = Sys.time()
difftime(end12, start12) # Duration: 3.3mn
intra_order_comparisons[1:2]


## Computing statistics

median_intra_order = sapply(intra_order_comparisons, function(x) sapply(split(x, x$ORDER1), function(y) median(y$SIMILARITY)))
median_intra_order_df = pivot_longer(tibble(cbind(tibble(ORDER = rownames(median_intra_order)), tibble(as.data.frame(median_intra_order)))), 
                                     !ORDER, "PRIMERS", values_to = "MEDIAN")

mad_intra_order = sapply(intra_order_comparisons, function(x) sapply(split(x, x$ORDER1), function(y) mad(y$SIMILARITY)))
mad_intra_order_df = pivot_longer(tibble(cbind(tibble(ORDER = rownames(mad_intra_order)), tibble(as.data.frame(mad_intra_order)))), 
                                  !ORDER, "PRIMERS", values_to = "MAD")

mean_intra_order = sapply(intra_order_comparisons, function(x) sapply(split(x, x$ORDER1), function(y) mean(y$SIMILARITY)))
mean_intra_order_df = pivot_longer(tibble(cbind(tibble(ORDER = rownames(mean_intra_order)), tibble(as.data.frame(mean_intra_order)))), 
                                   !ORDER, "PRIMERS", values_to = "MEAN")

sd_intra_order = sapply(intra_order_comparisons, function(x) sapply(split(x, x$ORDER1), function(y) sd(y$SIMILARITY)))
sd_intra_order_df = pivot_longer(tibble(cbind(tibble(ORDER = rownames(sd_intra_order)), tibble(as.data.frame(sd_intra_order)))), 
                                 !ORDER, "PRIMERS", values_to = "SD")

similarity_per_order = tibble(cbind(median_intra_order_df, mad_intra_order_df[,3], mean_intra_order_df[,3], sd_intra_order_df[,3]))

write.csv(similarity_per_order, "Output/intra_order/Similarity per order.csv", row.names = F)
similarity_per_order = tibble(read.csv("Output/intra_order/Similarity per order.csv"))
similarity_per_order





##### INTRA-TAXA VISUALISATIONS #####

## Computing each boxplot statistics separately, because ggplot crashes with 60 millions rows to treat

start13 = Sys.time()
all_intra_taxa_similarities = rbind(tibble(cbind(TAXA = paste0("BIN (", nrow(intra_bin_comparisons[[1]]) / 2, " comparisons)"),
                                                 bind_rows(lapply(intra_bin_comparisons, function(x) tibble(SIMILARITY = x$SIMILARITY)), .id = "PRIMERS"))),
                                    tibble(cbind(TAXA = paste0("GENUS (", nrow(intra_genus_comparisons[[1]]) / 2, " comparisons)"),
                                                 bind_rows(lapply(intra_genus_comparisons, function(x) tibble(SIMILARITY = x$SIMILARITY)), .id = "PRIMERS"))),
                                    tibble(cbind(TAXA = paste0("FAMILY (", nrow(intra_family_comparisons[[1]]) / 2, " comparisons)"),
                                                 bind_rows(lapply(intra_family_comparisons, function(x) tibble(SIMILARITY = x$SIMILARITY)), .id = "PRIMERS"))),
                                    tibble(cbind(TAXA = paste0("ORDER (", nrow(intra_order_comparisons[[1]]) / 2, " comparisons)"),
                                                 bind_rows(lapply(intra_order_comparisons, function(x) tibble(SIMILARITY = x$SIMILARITY)), .id = "PRIMERS"))))
end13 = Sys.time()
difftime(end13, start13) # Duration: 1.5mn
all_intra_taxa_similarities

start14 = Sys.time()
all_intra_taxa_boxplot = bind_rows(lapply(split(all_intra_taxa_similarities, all_intra_taxa_similarities$TAXA), 
                                          function(x) tibble(bind_rows(lapply(split(x, x$PRIMERS), function(y)
                                            setNames(data.frame(t(boxplot(as.numeric(y$SIMILARITY), plot = F)$stats)), 
                                                     c("MIN", "LOW", "MID", "TOP", "MAX"))), .id = "PRIMERS"))), .id = "TAXA")
end14 = Sys.time()
difftime(end14, start14) # Duration: 3.8mn
all_intra_taxa_boxplot

all_intra_taxa_boxplot$TAXA = factor(all_intra_taxa_boxplot$TAXA, levels = c(paste0("ORDER (", nrow(intra_order_comparisons[[1]]) / 2, " comparisons)"), 
                                                                             paste0("FAMILY (", nrow(intra_family_comparisons[[1]]) / 2, " comparisons)"), 
                                                                             paste0("GENUS (", nrow(intra_genus_comparisons[[1]]) / 2, " comparisons)"), 
                                                                             paste0("BIN (", nrow(intra_bin_comparisons[[1]]) / 2, " comparisons)")))

all_intra_taxa_boxplot$PRIMERS = replace(all_intra_taxa_boxplot$PRIMERS, which(all_intra_taxa_boxplot$PRIMERS == "FishF1-FishR1"), "FishF1/R1")
all_intra_taxa_boxplot

start15 = Sys.time()
all_intra_taxa_boxplot_outliers = bind_rows(lapply(split(all_intra_taxa_similarities, all_intra_taxa_similarities$TAXA), 
                                                   function(x) bind_rows(lapply(split(x, x$PRIMERS), function(y)
                                                     tibble(OUTLIER = boxplot(as.numeric(y$SIMILARITY), plot = F)$out)
                                                   ), .id = "PRIMERS")), .id = "TAXA")
end15 = Sys.time()
difftime(end15, start15) # Duration: 1mn
all_intra_taxa_boxplot_outliers

all_intra_taxa_boxplot_outliers = all_intra_taxa_boxplot_outliers[!duplicated(all_intra_taxa_boxplot_outliers),]
all_intra_taxa_boxplot_outliers

intra_taxa_boxplot_full_stats = left_join(all_intra_taxa_boxplot, all_intra_taxa_boxplot_outliers)

write.csv(intra_taxa_boxplot_full_stats, "Output/summary_intra_taxa_similarity/All similarities boxplot stats.csv", row.names = F)
intra_taxa_boxplot_full_stats = tibble(read.csv("Output/summary_intra_taxa_similarity/All similarities boxplot stats.csv"))
intra_taxa_boxplot_full_stats


## Intra-taxa boxplot of mean similarity per taxa

# Note: This kind of representation resolves the problem of inequal number of comparisons per taxa, but taxa with only few sequences
# Note: are more likely to contain chimeric sequences.

similarity_stat_all_taxa = tibble(rbind(cbind(TAXA = paste0("BIN (", length(unique(similarity_per_bin$BIN)), ")"), similarity_per_bin[,-1]),
                                        cbind(TAXA = paste0("GENUS (", length(unique(similarity_per_genus$GENUS)), ")"), similarity_per_genus[,-1]),
                                        cbind(TAXA = paste0("FAMILY (", length(unique(similarity_per_family$FAMILY)), ")"), similarity_per_family[,-1]),
                                        cbind(TAXA = paste0("ORDER (", length(unique(similarity_per_order$ORDER)), ")"), similarity_per_order[,-1])))

similarity_stat_all_taxa$TAXA = factor(similarity_stat_all_taxa$TAXA, levels = c(paste0("ORDER (", length(unique(similarity_per_order$ORDER)), ")"), 
                                                                                 paste0("FAMILY (", length(unique(similarity_per_family$FAMILY)), ")"), 
                                                                                 paste0("GENUS (", length(unique(similarity_per_genus$GENUS)), ")"), 
                                                                                 paste0("BIN (", length(unique(similarity_per_bin$BIN)), ")")))

similarity_stat_all_taxa = tibble(merge(mean_length_metabarcodes, similarity_stat_all_taxa))
similarity_stat_all_taxa

similarity_stat_all_taxa_subset = subset(similarity_stat_all_taxa, !(PRIMERS %in% c("Fish2deg", "L14735c2")))

similarity_stat_all_taxa_subset$PRIMERS = replace(similarity_stat_all_taxa_subset$PRIMERS, 
                                                  which(similarity_stat_all_taxa_subset$PRIMERS %in% "Fish2b"), "Fish2b/deg")

similarity_stat_all_taxa_subset$PRIMERS = replace(similarity_stat_all_taxa_subset$PRIMERS, 
                                                  which(similarity_stat_all_taxa_subset$PRIMERS %in% "L14735c"), "L14735c/c2")

similarity_stat_all_taxa_subset$PRIMERS = replace(similarity_stat_all_taxa_subset$PRIMERS, 
                                                  which(similarity_stat_all_taxa_subset$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")

length_ordered_primers = sapply(split(similarity_stat_all_taxa_subset, similarity_stat_all_taxa_subset$PRIMERS), function(x) 
  round(mean(x$MEAN_LENGTH)))
length_ordered_primers = sort(length_ordered_primers)
length_ordered_primers

length_ordered_primers_graph_name = setNames(paste0(names(length_ordered_primers), " (", 
                                                    length_ordered_primers, " bp)"), names(length_ordered_primers))
length_ordered_primers_graph_name

order_taxonomic_resolution = c("FishF1-R1", "Minibar", "L14912", "L14841", "L14735c/c2", "FishCB", "Teleo2", "MiFish", 
                               "Fish16S", "Ac12S", "AcMDB", "16SFD", "Ve16S", "Vert16S", "12SV5", "12SF1R1", "Ac16S", 
                               "L2513", "Teleo1", "Fish2b/deg")

similarity_stat_all_taxa_subset$PRIMERS = factor(similarity_stat_all_taxa_subset$PRIMERS, levels = order_taxonomic_resolution)

similarity_stat_all_taxa_subset = similarity_stat_all_taxa_subset %>% group_by(PRIMERS, TAXA) %>% 
  mutate(FIRST_QUANTILE_PRIMERS = quantile(MEAN, 0.25, na.rm = T),
         THIRD_QUANTILE_PRIMERS = quantile(MEAN, 0.75, na.rm = T))
quantile_per_primers = similarity_stat_all_taxa_subset[,c(1,5,10:11)]
quantile_per_primers = quantile_per_primers[!duplicated(quantile_per_primers),]
quantile_per_primers

best_cutoff_taxa = tibble(as.data.frame(sapply(split(quantile_per_primers, quantile_per_primers$PRIMERS), function(x)
  round(sort(setNames(x$THIRD_QUANTILE_PRIMERS, x$TAXA), decreasing = T)[-1] +
            (sort(setNames(x$FIRST_QUANTILE_PRIMERS, x$TAXA), decreasing = T)[-4] -
               sort(setNames(x$THIRD_QUANTILE_PRIMERS, x$TAXA), decreasing = T)[-1])/2, 1))))
best_cutoff_taxa$LIMIT = c("BIN - GENUS", "GENUS - FAMILY", "FAMILY - ORDER")

best_cutoff_taxa_long = best_cutoff_taxa %>% pivot_longer(!LIMIT, names_to = "PRIMERS", values_to = "THRESHOLD")
best_cutoff_taxa_long$PRIMERS = factor(best_cutoff_taxa_long$PRIMERS, levels = order_taxonomic_resolution)

best_cutoff_taxa_graph = bind_rows(lapply(split(best_cutoff_taxa_long, best_cutoff_taxa_long$PRIMERS), function(x) 
  tibble(GRAPH = paste0(paste0(sort(x$THRESHOLD, decreasing = T), collapse = "% > "), "%"))), .id = "PRIMERS")
best_cutoff_taxa_graph$PRIMERS = factor(best_cutoff_taxa_graph$PRIMERS, levels = order_taxonomic_resolution)

ggplot(similarity_stat_all_taxa_subset, aes(x = TAXA, y = MEAN)) +
  facet_wrap(~PRIMERS, nrow = 2, ncol = 10) +
  geom_hline(data = best_cutoff_taxa_long, aes(yintercept = THRESHOLD), color = "red") +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(aes(colour = TAXA, fill = TAXA), outlier.size = 0.7) +
  geom_boxplot(aes(fill = TAXA), outlier.colour = NA) +
  scale_y_continuous(labels = function(x) paste0(x, '%'), limits = c(40, 103), minor_breaks = seq(40, 100, 5),
                     breaks = seq(40, 100, 10)) + labs(y = "Mean intra-taxa similarity") +
  theme_() + theme(plot.caption = element_text(hjust = 0.5), legend.position = "bottom", 
                   axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                   strip.text = element_text(family = "Segoe UI Semibold"),
                   legend.title = element_text(family = "Segoe UI Semibold"),
                   panel.grid.major.y = element_line(size = 0.3, color = "darkgrey"),
                   panel.grid.minor.y = element_line(size = 0.3, color = "darkgrey")) +
  geom_text(data = best_cutoff_taxa_graph, aes(label = GRAPH, x = 2.5, y = 103), 
            size = 1.9, color = "red", family = "Segoe UI Semibold")

ggsave("Graph/intra_taxa_similarity/Intra-taxa - Mean intra-taxa similarity.png", dpi = 600, 
       device = png, type = "cairo", width = 25.53229, height = 14.39333, units = "cm")


## Plotting of all boxplots together, representing the distribution of similarities for all comparisons

# Note: This graph is biased because taxa with a lot of comparisons (e.g. large number of infra-taxa and/or sequences per taxa) 
# Note: will overwhelm the other taxa with less comparisons. 

intra_taxa_boxplot_full_stats$TAXA = factor(intra_taxa_boxplot_full_stats$TAXA, levels = c(unique(intra_taxa_boxplot_full_stats$TAXA)[4], 
                                                                                           unique(intra_taxa_boxplot_full_stats$TAXA)[2],
                                                                                           unique(intra_taxa_boxplot_full_stats$TAXA)[3], 
                                                                                           unique(intra_taxa_boxplot_full_stats$TAXA)[1]))
intra_taxa_boxplot_full_stats

intra_taxa_boxplot_full_stats_subset = subset(intra_taxa_boxplot_full_stats, !(PRIMERS %in% c("Fish2deg", "L14735c2")))

intra_taxa_boxplot_full_stats_subset$PRIMERS = replace(intra_taxa_boxplot_full_stats_subset$PRIMERS, 
                                                  which(intra_taxa_boxplot_full_stats_subset$PRIMERS %in% "Fish2b"), "Fish2b/deg")

intra_taxa_boxplot_full_stats_subset$PRIMERS = replace(intra_taxa_boxplot_full_stats_subset$PRIMERS, 
                                                  which(intra_taxa_boxplot_full_stats_subset$PRIMERS %in% "L14735c"), "L14735c/c2")

intra_taxa_boxplot_full_stats_subset$PRIMERS = replace(intra_taxa_boxplot_full_stats_subset$PRIMERS, 
                                                  which(intra_taxa_boxplot_full_stats_subset$PRIMERS %in% "FishF1/R1"), "FishF1-R1")

intra_taxa_boxplot_full_stats_subset$PRIMERS = factor(intra_taxa_boxplot_full_stats_subset$PRIMERS, levels = order_taxonomic_resolution)

intra_taxa_boxplot_full_stats_subset = intra_taxa_boxplot_full_stats_subset %>% group_by(PRIMERS, TAXA) %>% 
  mutate(FIRST_QUANTILE_PRIMERS = unique(LOW), THIRD_QUANTILE_PRIMERS = unique(TOP))
quantile_per_primers_full = intra_taxa_boxplot_full_stats_subset[,c(1,2,9:10)]
quantile_per_primers_full = quantile_per_primers_full[!duplicated(quantile_per_primers_full),]
quantile_per_primers_full

best_cutoff_taxa_full = tibble(as.data.frame(sapply(split(quantile_per_primers_full, quantile_per_primers_full$PRIMERS), function(x)
  round(sort(setNames(x$THIRD_QUANTILE_PRIMERS, x$TAXA), decreasing = T)[-1] +
          (sort(setNames(x$FIRST_QUANTILE_PRIMERS, x$TAXA), decreasing = T)[-4] -
             sort(setNames(x$THIRD_QUANTILE_PRIMERS, x$TAXA), decreasing = T)[-1])/2, 1))))
best_cutoff_taxa_full$LIMIT = c("BIN - GENUS", "GENUS - FAMILY", "FAMILY - ORDER")

best_cutoff_taxa_full_long = best_cutoff_taxa_full %>% pivot_longer(!LIMIT, names_to = "PRIMERS", values_to = "THRESHOLD")
best_cutoff_taxa_full_long$PRIMERS = factor(best_cutoff_taxa_full_long$PRIMERS, levels = order_taxonomic_resolution)

best_cutoff_taxa_full_graph = bind_rows(lapply(split(best_cutoff_taxa_full_long, best_cutoff_taxa_full_long$PRIMERS), function(x) 
  tibble(GRAPH = paste0(paste0(sort(x$THRESHOLD, decreasing = T), collapse = "% > "), "%"))), .id = "PRIMERS")
best_cutoff_taxa_full_graph$PRIMERS = factor(best_cutoff_taxa_full_graph$PRIMERS, levels = order_taxonomic_resolution)

ggplot(intra_taxa_boxplot_full_stats_subset) +
  facet_wrap(~PRIMERS, nrow = 2, ncol = 10) +
  geom_hline(data = best_cutoff_taxa_full_long, aes(yintercept = THRESHOLD), color = "red") +
  geom_point(aes(x = TAXA, y = OUTLIER, colour = TAXA), size = 0.7) +
  geom_errorbar(aes(x = TAXA, ymin = MIN, ymax = MAX), width = 0.5) +
  geom_boxplot(data = intra_taxa_boxplot_full_stats_subset[,1:7][!duplicated(intra_taxa_boxplot_full_stats_subset[,1:7]), ], 
               aes(x = TAXA, ymin = MIN, lower = LOW, middle = MID, upper = TOP, 
                   ymax = MAX, fill = TAXA), stat = "identity", width = 0.8) +
  scale_y_continuous(labels = function(x) paste0(x, '%'), limits = c(40, 103), minor_breaks = seq(40, 100, 5),
                     breaks = seq(40, 100, 10)) + labs(y = "Mean intra-taxa similarity") +
  theme_() + theme(plot.caption = element_text(hjust = 0.5), legend.position = "bottom", 
                   axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), 
                   strip.text = element_text(family = "Segoe UI Semibold"),
                   legend.title = element_text(family = "Segoe UI Semibold"),
                   panel.grid.major.y = element_line(size = 0.3, color = "darkgrey"),
                   panel.grid.minor.y = element_line(size = 0.3, color = "darkgrey")) +
  geom_text(data = best_cutoff_taxa_full_graph, aes(label = GRAPH, x = 2.5, y = 103), 
            size = 1.9, color = "red", family = "Segoe UI Semibold")

ggsave("Graph/intra_taxa_similarity/Intra-taxa - All similarities.png", dpi = 600, 
       device = png, type = "cairo", width = 25.53229, height = 14.39333, units = "cm")



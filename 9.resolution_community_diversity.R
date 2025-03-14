##### INITIALISATION #####

## Loading required packages and custom functions

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(ggplot2)
library(extrafont)
library(tidyverse)


## Loading other packages required for this script

library(reshape2)
library(elementalist)
library(gridExtra)
library(grid)


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

mean_length_metabarcodes = tibble(read.csv("Data/primers_infos/Mean length data.csv"))
mean_length_metabarcodes %>% print(n = 22)

primers = sort(unique(filtered_metabarcodes_final$PRIMERS))
primers


## Loading intra-BIN comparisons

intra_bin_comparisons = readRDS("Output/intra_BIN/Intra-BIN - Raw comparisons - All primers.rds") # Duration: 15s
intra_bin_comparisons[1:2]


## Loading inter-BIN comparisons

inter_bin_comparisons = readRDS("Output/inter_BIN/Inter-BIN - Raw comparisons - All primers.rds") # Duration: 1mn
inter_bin_comparisons[1:2]





##### SIMULATING A CLUSTERING PER PRIMER #####

## Assigning each sequences to a cluster according to the 10 thresholds

start1 = Sys.time()
clusters_accession_list = clusters_per_primers_decipher(infos_data = act_mitogenomes_bin_corrected, 
                                                        intra_comparisons_list = intra_bin_comparisons, 
                                                        inter_comparisons_list = inter_bin_comparisons, 
                                                        similarity_threshold = 90:99, 
                                                        tree_method = "NJ", 
                                                        name_col_similarity = "SIMILARITY",
                                                        name_col_accession = "ACCESSION", 
                                                        name_col_taxa = "BIN",
                                                        nb_processors = 35)
end1 = Sys.time()
difftime(end1, start1) # Duration: 2h on a computation cluster


## Saving/reading the result

saveRDS(clusters_accession_list, "Output/community_simulations/NJ clusters per primers list.rds")
clusters_primers_thresholds_list = readRDS("Output/community_simulations/NJ clusters per primers list.rds")
clusters_primers_thresholds_list[1:2]




##### SIMULATING RANDOM COMMUNITIES AND CHECKING THE NUMBER OF OTU PERCEIVED #####

## Adding BINs to the output with NJ OTUs per sequences and thresholds

clusters_primers_thresholds_list_bin = lapply(1:length(clusters_primers_thresholds_list), function(i)
  tibble(merge(act_mitogenomes_bin_corrected[,1:2], clusters_primers_thresholds_list[[i]], by = "ACCESSION")))
names(clusters_primers_thresholds_list_bin) = names(clusters_primers_thresholds_list)
clusters_primers_thresholds_list_bin[1:2]


## Performing mock community simulations

start2 = Sys.time()
perception_mock_communities(clusters_taxa_list = clusters_primers_thresholds_list_bin, 
                              nb_random_community = 1000, nb_seq_sampled = seq(10, 1000, by = 10),
                            output_path = "Output/community_simulations/random_communities_results")
end2 = Sys.time()
difftime(end2, start2) # Duration: 4.5h 


## Preparing the summary file for plotting

perception_per_primers_final = tibble(merge(mean_length_metabarcodes[,1:2], perception_per_primers))
perception_per_primers_final

perception_per_primers_final = subset(perception_per_primers_final, PRIMERS != "L14735c2")
perception_per_primers_final$PRIMERS = replace(perception_per_primers_final$PRIMERS, 
                                               which(perception_per_primers_final$PRIMERS == "L14735c"), 
                                               "L14735c & L14735c2")
perception_per_primers_final


## Saving the summary file

write.csv(perception_per_primers_final, "Output/community_simulations/Summary perception analysis.csv", row.names = F)
perception_per_primers_final = tibble(read.csv("Output/community_simulations/Summary perception analysis.csv"))
perception_per_primers_final





##### PLOTTING THE RESULTS FOR EACH THRESHOLDS #####

## Preparing the data for plotting

perception_per_primers_final = subset(perception_per_primers_final, PRIMERS != "Fish2deg")

perception_per_primers_final$PRIMERS = replace(perception_per_primers_final$PRIMERS, 
                                                 which(perception_per_primers_final$PRIMERS %in% "Fish2b"), "Fish2b/deg")

perception_per_primers_final$PRIMERS = replace(perception_per_primers_final$PRIMERS, 
                                                 which(perception_per_primers_final$PRIMERS %in% "L14735c & L14735c2"), "L14735c/c2")

perception_per_primers_final$PRIMERS = replace(perception_per_primers_final$PRIMERS, 
                                                 which(perception_per_primers_final$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")


## Defining a graphic function easily modifiable

plot_perception_analysis = function(community_data, colours_per_gene_list, upper_padding_y_axis_list, threshold){
  
  community_data_threshold = subset(community_data, SIMILARITY == threshold)
  
  nb_seq = unique(community_data_threshold$NB_SEQ)
  
  community_list_gene = split(community_data_threshold, community_data_threshold$GENE)
  
  gene_names = names(community_list_gene)
  
  graph_gene_list = list()
  
  order_misidentifications_list = list()
  
  best_level_misidentifications_list = list()
  
  for(i in 1:length(community_list_gene)){
    
    ### Optional (to print the best level)
    
    # order_misidentifications = round(sort(sapply(split(community_list_gene[[i]], community_list_gene[[i]]$PRIMERS), function(x) 
    #   abs(c(100 - x$MEAN_PERCENT_PERCEIVED)[which.min(abs(c(100 - x$MEAN_PERCENT_PERCEIVED) - 0))]))), 2)
    # 
    # sign_misidentifications_unsorted = sapply(split(community_list_gene[[i]], community_list_gene[[i]]$PRIMERS), function(x) 
    #   ifelse(sign(min(x$MEAN_PERCENT_PERCEIVED - 100)) == -1, "-", ""))
    # 
    # sign_misidentifications_sorted = sign_misidentifications_unsorted[order(match(names(sign_misidentifications_unsorted), 
    #                                                                               names(order_misidentifications)))]
    # 
    # sign_order_misidentifications = setNames(paste0(sign_misidentifications_sorted, order_misidentifications, "%"),
    #                                          names(order_misidentifications))
    # 
    # sign_order_misidentifications = replace(sign_order_misidentifications, grep("-0%", sign_order_misidentifications), "0%")
    # 
    # order_misidentifications_list[[i]] = order_misidentifications
    # 
    # best_level_misidentifications = sort(sapply(split(community_list_gene[[i]], community_list_gene[[i]]$PRIMERS), function(x)
    #   x[which.min(abs(x$MEAN_PERCENT_PERCEIVED - 100)), ]$NB_SEQ))
    # 
    # best_level_misidentifications_list[[i]] = best_level_misidentifications
    # 
    # names_primers = list()
    # 
    # for(j in 1:length(best_level_misidentifications)) {
    #   
    #   names_primers[j] = paste0(names(best_level_misidentifications)[j], "\n \n", best_level_misidentifications[j], ": ")
    #   
    # }
    # 
    # community_list_gene[[i]]$PRIMERS = rep(paste0(unlist(names_primers)[order(unlist(names_primers))],
    #                                               sign_order_misidentifications[order(names(sign_order_misidentifications))]),
    #                                        each = length(nb_seq))
    # 
    # best_level_misidentifications = best_level_misidentifications[order(match(names(best_level_misidentifications), 
    #                                                                           names(order_misidentifications)))]
    # 
    # names(order_misidentifications) = paste0(names(order_misidentifications), "\n \n", 
    #                                          best_level_misidentifications, ": ", sign_order_misidentifications)
    # 
    # community_list_gene[[i]]$PRIMERS = factor(community_list_gene[[i]]$PRIMERS, levels = names(order_misidentifications))
    
    ### End of the optional part
    
    colours_ratio = colours_per_gene_list[[which(names(colours_per_gene_list) == gene_names[i])]]
    
    graph_gene_list[[i]] = ggplot(community_list_gene[[i]], aes(x = NB_SEQ, y = MEAN_PERCENT_PERCEIVED, 
                                                                group = PRIMERS, color = PRIMERS, fill = PRIMERS)) +
      geom_hline(aes(yintercept = 100), color = "black", size = 1) +  
      geom_ribbon(aes(ymin = MEAN_PERCENT_PERCEIVED - SD_PERCENT_PERCEIVED, 
                      ymax = MEAN_PERCENT_PERCEIVED + SD_PERCENT_PERCEIVED), alpha = 0.1, color = NA) +
      geom_line(size = 0.5) +
      scale_y_continuous(labels = function(x) paste0(x, "%"),
                         limits = c(NA, max(subset(perception_per_primers_final, SIMILARITY == threshold)$MEAN_PERCENT_PERCEIVED) +
                                      upper_padding_y_axis_list[[which(names(upper_padding_y_axis_list) == gene_names[i])]])) +
      scale_color_manual(values = colours_ratio) + scale_fill_manual(values = colours_ratio) +
      guides(color = guide_legend(label.position = "top", nrow = 1, title = paste0(gene_names[i], "  "), title.vjust = 0.5, order = 1),
             fill = guide_legend(label.position = "top", nrow = 1, title = paste0(gene_names[i], "  "), title.vjust = 0.5, order = 1)) +
      theme_() + theme(plot.margin = unit(c(1.2,0.2,0.2,0.2), "cm"),
                       axis.title.x = element_blank(), axis.title.y = element_blank(),
                       legend.title = element_text(size = 15, family = "Segoe UI Semibold"),
                       legend.box.background = element_rect_round(color = "black", fill = "white",
                                                                  size = 0.5, radius = unit(0.1, "snpc")),
                       legend.text = element_text(size = 9.5, family = "Segoe UI Semibold"), # size = 7 if printing the best level
                       legend.box.margin = margin(2, 1, 4, 1), legend.spacing.y = unit(0.3, 'cm'),
                       legend.direction = "horizontal", legend.position = c(0.5, 1),
                       legend.key.height = unit(0.2, "cm"))
    
  }
  
  return(graph_gene_list)
  
}


## Defining the colours for each genes

colour_per_genes = list(c("#BC00B5", "#E7E300", "#757575", "#00CF15", "#00B6D8", "#004FE7", "#F38F00"),
                        c("#BC00B5", "#F38F00", "#E7E300", "#004FE7", "#00CF15", "#757575"),
                        c("#FF0000", "#00B6D8"),
                        c("#FF0000", "#BC00B5", "#00CF15", "#757575", "#00B6D8", "#F38F00"))
names(colour_per_genes) = c("12S", "16S", "COI", "CytB")
colour_per_genes


## Plotting and saving the result of the analysis for the 90% threshold

padding_per_genes = list(8,6,1,5)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_90 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 90)

perception_ratio_90_graph = grid.arrange(perception_ratio_90[[1]], perception_ratio_90[[2]],
                                         perception_ratio_90[[3]], perception_ratio_90[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (90% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 90.png", perception_ratio_90_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 91% threshold

padding_per_genes = list(8,6,1,5)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_91 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 91)

perception_ratio_91_graph = grid.arrange(perception_ratio_91[[1]], perception_ratio_91[[2]],
                                         perception_ratio_91[[3]], perception_ratio_91[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (91% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 91.png", perception_ratio_91_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 92% threshold

padding_per_genes = list(8,6,1,5)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_92 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 92)

perception_ratio_92_graph = grid.arrange(perception_ratio_92[[1]], perception_ratio_92[[2]],
                                         perception_ratio_92[[3]], perception_ratio_92[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (92% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 92.png", perception_ratio_92_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 93% threshold

padding_per_genes = list(8,6,1,5)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_93 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 93)

perception_ratio_93_graph = grid.arrange(perception_ratio_93[[1]], perception_ratio_93[[2]],
                                         perception_ratio_93[[3]], perception_ratio_93[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (93% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 93.png", perception_ratio_93_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 94% threshold

padding_per_genes = list(8,6,1,6)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_94 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 94)

perception_ratio_94_graph = grid.arrange(perception_ratio_94[[1]], perception_ratio_94[[2]],
                                         perception_ratio_94[[3]], perception_ratio_94[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (94% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 94.png", perception_ratio_94_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 95% threshold

padding_per_genes = list(8,6,1,5)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_95 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 95)

perception_ratio_95_graph = grid.arrange(perception_ratio_95[[1]], perception_ratio_95[[2]],
                                         perception_ratio_95[[3]], perception_ratio_95[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (95% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 95.png", perception_ratio_95_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 96% threshold

padding_per_genes = list(8,6,1,6)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_96 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 96)

perception_ratio_96_graph = grid.arrange(perception_ratio_96[[1]], perception_ratio_96[[2]],
                                         perception_ratio_96[[3]], perception_ratio_96[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (96% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 96.png", perception_ratio_96_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 97% threshold

padding_per_genes = list(9,6,0,8)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_97 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 97)

perception_ratio_97_graph = grid.arrange(perception_ratio_97[[1]], perception_ratio_97[[2]],
                                         perception_ratio_97[[3]], perception_ratio_97[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (97% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 97.png", perception_ratio_97_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 98% threshold

padding_per_genes = list(9,6,0,8)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_98 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 98)

perception_ratio_98_graph = grid.arrange(perception_ratio_98[[1]], perception_ratio_98[[2]],
                                         perception_ratio_98[[3]], perception_ratio_98[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (98% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 98.png", perception_ratio_98_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


## Plotting and saving the result of the analysis for the 99% threshold

padding_per_genes = list(8,6,0,8)
names(padding_per_genes) = c("12S", "16S", "COI", "CytB")

perception_ratio_99 = plot_perception_analysis(perception_per_primers_final, colour_per_genes, 
                                               padding_per_genes, threshold = 99)

perception_ratio_99_graph = grid.arrange(perception_ratio_99[[1]], perception_ratio_99[[2]],
                                         perception_ratio_99[[3]], perception_ratio_99[[4]],
                                         ncol = 2, nrow = 2, widths = c(1,1.04), 
                                         bottom = textGrob("Number of sequences retrieved from the community", 
                                                           gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)),
                                         left = textGrob("Number of OTU perceived (99% threshold) / Number BIN community", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 10)))

ggsave("Graph/community_simulations/Perception analysis graph - 99.png", perception_ratio_99_graph, device = png,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, unit = "cm")


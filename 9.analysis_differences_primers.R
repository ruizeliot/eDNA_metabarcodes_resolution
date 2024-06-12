##### INITIALISATION #####

## Installing the functions designed for this study directly from Github 

#library(remotes) ; install_github("Eliot-RUIZ/eDNAevaluation", upgrade = "never")
#library(eDNAevaluation) -> NOT POSSIBLE NOW BECAUSE PACKAGE IS PRIVATE


## In case of error, even after updating R (see package "installr"), please install the function locally and load dependencies manually

miceadds::source.all("Functions")
library(tibble)
library(DECIPHER)
library(stringr)
library(ggplot2)
library(extrafont)
library(tidyverse)
library(robustlmm)
library(insight)
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

act_mitogenomes_bin_full_infos = tibble(merge(act_mitogenomes_bin_corrected, bin_database_corrected))
act_mitogenomes_bin_full_infos

mean_length_metabarcodes = tibble(read.csv("Data/primers_infos/Mean length data.csv"))
mean_length_metabarcodes %>% print(n = 22)

primers = sort(unique(filtered_metabarcodes_final$PRIMERS))
primers





##### MULTIPLE EVALUATION ON THE WHOLE DATASET #####

## Loading the result of intra-BIN and inter-BIN analyses

intra_bin_summary = readRDS("Output/intra_BIN/Intra-BIN analysis - Summary - All primers.rds")
intra_bin_summary[1:2]

inter_bin_summary = readRDS("Output/inter_BIN/Inter-BIN analysis - Summary - All primers.rds")
inter_bin_summary[1:2]


## Combining the results of all analyses in one single dataframe

for(i in 1:length(intra_bin_summary)) { 
  
  similarity_analysis_resolution_temporary = tibble(cbind(tibble(SIMILARITY = names(intra_bin_summary)[i]),
                                                          intra_bin_summary[[names(intra_bin_summary)[i]]], 
                                                          inter_bin_summary[[names(intra_bin_summary)[i]]][,-1]))
  
  similarity_analysis_resolution_temporary$CUMULATED_PERCENT_ERRORS_PER_ANALYSIS = 
    similarity_analysis_resolution_temporary$PERCENT_INTRA_ERRORS +
    similarity_analysis_resolution_temporary$PERCENT_INTER_ERRORS
  
  if(i == 1) similarity_analysis_resolution = similarity_analysis_resolution_temporary
  
  else similarity_analysis_resolution = rbind(similarity_analysis_resolution, similarity_analysis_resolution_temporary)
  
}


## Filtering and reordering the results

similarity_analysis_resolution = subset(similarity_analysis_resolution, select = 
                                          c("SIMILARITY", "PRIMERS", 
                                            "NUMBER_INTRA_ERRORS", "NUMBER_INTER_ERRORS",
                                            "PERCENT_INTRA_ERRORS", "PERCENT_INTER_ERRORS", 
                                            "CUMULATED_PERCENT_ERRORS_PER_ANALYSIS",
                                            "NUMBER_BIN_WITH_INTRA_ERRORS", "NUMBER_BIN_PAIRS_WITH_INTER_ERRORS",
                                            "PERCENT_BIN_WITH_INTRA_ERRORS", "PERCENT_BIN_PAIRS_WITH_INTER_ERRORS",
                                            "MEAN_NUMBER_INTRA_ERRORS_PER_BIN", "MEAN_NUMBER_INTER_ERRORS_PER_BIN_PAIRS",
                                            "MEAN_PERCENT_INTRA_ERRORS_PER_BIN", "MEAN_PERCENT_INTER_ERRORS_PER_BIN_PAIRS"))
similarity_analysis_resolution


## Saving results

write.csv(similarity_analysis_resolution, "Output/analysis_differences_primers/Multiple resolution analysis - All similarity thresholds.csv", row.names = F)
similarity_analysis_resolution = tibble(read.csv("Output/analysis_differences_primers/Multiple resolution analysis - All similarity thresholds.csv"))
similarity_analysis_resolution





##### MULTIPLE EVALUATION ON THE WHOLE DATABASE BUT WITH SUBSET OF ORDERS #####

## Preparing thresholds reference table

thresholds = 90:99
similarity_reference = tibble(data.frame(THRESHOLD = thresholds, INDEX = 1:length(thresholds)))
similarity_reference


## Choosing orders with sufficient number of BIN in each (> 100 BIN)

orders_to_study = act_mitogenomes_bin_corrected %>% group_by(ORDER, .drop = F) %>% 
  summarise(NUMBER_BIN_ORDER = n_distinct(BIN)) %>% filter(NUMBER_BIN_ORDER >= 100)
orders_to_study


## Running evaluations for all orders

start3 = Sys.time()
orders_analysis_resolution = multiple_resolution_analyses(metabarcodes_data = filtered_metabarcodes_final, 
                                                          infos_data = act_mitogenomes_bin_full_infos, 
                                                          elements_to_subset = orders_to_study$ORDER, 
                                                          similarity_threshold = 90:99, 
                                                          name_col_subset = "ORDER", 
                                                          name_col_taxa = "BIN", 
                                                          name_col_primer = "PRIMERS", 
                                                          name_col_accession = "ACCESSION", 
                                                          name_col_fragment = "FRAGMENT", 
                                                          # Change it with the path + name of your VSEARCH program
                                                          program_name = "vsearch-2.17.1-win-x86_64.exe")
end3 = Sys.time()
difftime(end3, start3) # Duration: 38mn
orders_analysis_resolution


## Saving results

write.csv(orders_analysis_resolution, "Output/analysis_differences_primers/Multiple resolution analysis - All major orders.csv", row.names = F)
orders_analysis_resolution = tibble(read.csv("Output/analysis_differences_primers/Multiple resolution analysis - All major orders.csv"))
orders_analysis_resolution




##### STATISTICAL ANALYSIS OF EFFECT OF SIMILARITY ON RESOLUTION AND IF THERE IS DIFFERENCE OF EFFECTS BETWEEN GENES #####

## Choosing model parameters

# Reference group: Taking COI metabarcodes as reference because they are in general the most resolutive and there is FishF1-FishR1, 
# the metabarcode used to define fishes' BIN.

stat_analysis_resolution = tibble(merge(mean_length_metabarcodes, similarity_analysis_resolution))

stat_analysis_resolution$GENE = factor(stat_analysis_resolution$GENE, levels = c("COI", "CytB", "12S", "16S"))

# Random effects (i.e. conditional variables of irrelevant importance here): 
# PRIMERS -> To take into account groups of repeated measures in similarity, and not study the variability of primers


## Running linear mixed-effects models

lmer_intra_bin_similarity = lmer(PERCENT_INTRA_ERRORS ~ GENE * SIMILARITY + 
                                   (1 | PRIMERS), data = stat_analysis_resolution)  # Intercept of PRIMERS can vary randomly
lmer_inter_bin_similarity = lmer(PERCENT_INTER_ERRORS ~ GENE * SIMILARITY + 
                                   (1 | PRIMERS), data = stat_analysis_resolution)  # Intercept of PRIMERS can vary randomly

shapiro.test(residuals(lmer_intra_bin_similarity)) # No normal distribution
shapiro.test(residuals(lmer_inter_bin_similarity)) # No normal distribution

plot(lmer_intra_bin_similarity, which = 1) # Heteroskedasticity + Outliers // Quite linear
plot(lmer_inter_bin_similarity, which = 1) # Heteroskedasticity + Outliers // Quite linear


## Running robust linear mixed-effects models

rlmer_intra_bin_similarity = rlmer(PERCENT_INTRA_ERRORS ~ GENE * SIMILARITY + 
                                     (1 | PRIMERS), data = stat_analysis_resolution)
rlmer_inter_bin_similarity = rlmer(PERCENT_INTER_ERRORS ~ GENE * SIMILARITY + 
                                     (1 | PRIMERS), data = stat_analysis_resolution)

summary(rlmer_intra_bin_similarity)
summary(rlmer_inter_bin_similarity)

plot.rlmerMod(rlmer_intra_bin_similarity, which = 2)[[1]] + stat_qq_line(aes(sample = resid(rlmer_intra_bin_similarity))) + theme_()
plot.rlmerMod(rlmer_inter_bin_similarity, which = 2)[[1]] + stat_qq_line(aes(sample = resid(rlmer_inter_bin_similarity))) + theme_()
# Clear absence of normality in the data -> more robustness is given to values in the middle of the plot

plot(rlmer_intra_bin_similarity, which = 3)[[1]] + theme_()
plot(rlmer_inter_bin_similarity, which = 3)[[1]] + theme_()
# Random effects follow a more normal distribution thus almost all weigths are equal to 1

plot(rlmer_intra_bin_similarity, which = 1)[[1]] + theme_()
plot(rlmer_inter_bin_similarity, which = 1)[[1]] + theme_()
# Spread of residuals increase towards rigth of plot = heteroskedasticity
# Many outliers with large residuals values are likely to cause this heteroskedasticity -> small weigths
# Appart from these few values, linearity seems to hold reasonably well


## Checking R squared of the model

R2_rlmer = function(rlmer_model){
  
  fixed_variance = get_variance_fixed(rlmer_model)
  
  random_variance = get_variance_random(rlmer_model)
  
  residual_variance = get_variance_residual(rlmer_model)
  
  marginal_R2 = fixed_variance / (fixed_variance + random_variance + residual_variance) # Nakagawa & Schielzeth (2013)
  
  conditional_R2 = (fixed_variance + random_variance) / (fixed_variance + random_variance + residual_variance) # Nakagawa & Schielzeth (2013)
  
  cat(paste0("Marginal R2 (fixed effects): ", round(marginal_R2, 4), " (", round(marginal_R2 * 100, 2), "% of variance explained)")) 
  
  cat("\n\n")
  
  cat(paste0("Conditional R2 (fixed + random effects): ", round(conditional_R2, 4), " (", round(conditional_R2 * 100, 2), "% of variance explained)"))
  
}

R2_rlmer(rlmer_intra_bin_similarity) # MR²: 74.77% vs CR²: 81.64%

R2_rlmer(rlmer_inter_bin_similarity) # MR²: 63.57% vs CR²: 90.05%


## Statistical analysis of fixed effects

stats.rlmerMod = function(object, level = 0.95) { # Function from Lockwood et al. (2021) supplementary material (Psychological Science)
  
  significance = function(p) {
    
    if(p >= 0.05 && p < 0.1)    significance = "."
    
    else if(p >= 0.01 && p < 0.05)    significance = "*"
    
    else if(p >= 0.001 && p < 0.01)    significance = "**"
    
    else if(p < 0.001)    significance = "***"
    
    else significance = "ns"
    
    significance
    
  }
  
  # Extract beta coefficients
  beta <- fixef(object)
  
  # Extract names of coefficients
  parm <- names(beta)
  
  # Extract standard errors for the coefficients
  se <- sqrt(diag(vcov(object)))
  
  # Set level of confidence interval
  conf.level <- qnorm((1 + level) / 2)
  
  # Calculate z value
  z = beta/se
  
  # Calculate CI and create table
  ctab <- cbind(round(beta, 5),
                round(beta - (conf.level * se), 5), 
                round(beta + (conf.level * se), 5),
                round(se, 5),
                round(z, 5),
                round(2*pnorm(-abs(z)), 5),
                sapply(2*pnorm(-abs(z)), significance))
  
  
  # label column names
  
  output = data.frame(apply(ctab[,1:6], 2, as.numeric), ctab[,7])
  
  rownames(output) = rownames(ctab)
  
  colnames(output) <- c('beta',
                        paste(100 * ((1 - level) / 2), '%'),
                        paste(100 * ((1 + level) / 2), '%'),
                        'SE',
                        'z',
                        'p',
                        ' ')
  
  # Output
  return(output)
  
}

stats.rlmerMod(rlmer_intra_bin_similarity)

stats.rlmerMod(rlmer_inter_bin_similarity)


# Effect of similarity on the reference (COI gene):

# Intra-BIN: Percent of wrongly discrimined cases per BIN for COI metabarcodes

# Intercept of COI is significantly different of 0 for COI, and homogeneous for 12S and 16S
# Intercept of CytB is significantly inferior than the intercept of COI

# Slopes are significantly different from 0 for all metabarcodes (COI different and all higher)
# Slopes are homogeneous among COI, 12S and 16S, but significantly higher than COI for CytB
# COI, 12S & 16S: increase by approximately 0.06 for every increase of 1% of similarity threshold
# CytB: increase by approximately 0.34 (0.06 + 0.28) for every increase of 1% of similarity threshold


# Inter-BIN: Percent of confounded sequences sharing different BIN for COI metabarcodes

# Intercept of COI significantly different of 0 and not significantly different than the intercept of CytB
# Intercept of 12S and 16S both significantly higher than the intercept of COI

# Slope of COI is significantly different from 0, and decreasesing by -0.09% for every increase of 1% of similarity threshold
# Slope of CytB is not significantly different than the slope for COI
# Slope of 12s and 16S are significantly different than the slope for COI
# 12S: decrease by approximately -0.24% (-0.09 - 0.15) for every increase of 1% of similarity threshold
# 16S: decrease by approximately -0.20% (-0.09 - 0.11) for every increase of 1% of similarity threshold


# The magnitude of the effect of changes of similarity is 0.5 time less important for intra-BIN than for inter-BIN resolution (for COI)
abs(-6.16811 * 0.06842) / abs(8.83808 * -0.08923)


## Plotting model predictions

similarity_rlmer_graph = tibble(rbind(cbind(stat_analysis_resolution, TYPE = "Intra-BIN", 
                                            PREDICTED = predict(rlmer_intra_bin_similarity)), 
                                      cbind(stat_analysis_resolution, TYPE = "Inter-BIN", 
                                            PREDICTED = predict(rlmer_inter_bin_similarity))))
similarity_rlmer_graph$TYPE = factor(similarity_rlmer_graph$TYPE, levels = c("Intra-BIN", "Inter-BIN"))

similarity_rlmer_plot = ggplot(similarity_rlmer_graph, aes(x = SIMILARITY, colour = GENE)) +
  facet_wrap(~ TYPE, labeller = labeller(TYPE = as_labeller(c("Inter-BIN" = "Inter-BIN~(R[~(c)]^{2}=={0.901})",
                                                              "Intra-BIN" = "Intra-BIN~(R[~(c)]^{2}=={0.816})"), label_parsed))) +
  geom_smooth(data = subset(similarity_rlmer_graph, TYPE == "Intra-BIN"), method = "lm", 
              aes(y = PREDICTED, fill = GENE), alpha = 0.1, lwd = 0.8) +
  geom_smooth(data = subset(similarity_rlmer_graph, TYPE == "Inter-BIN"), method = "lm", 
              aes(y = PREDICTED, fill = GENE), alpha = 0.1, lwd = 0.8) +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = c(0, 1, 2, 3), labels = function(x) paste0(x, "%")) + 
  scale_color_manual(values = c("#FF0000", "#00B6D8", "#00CF15", "#BC00B5")) +
  scale_fill_manual(values = c("#FF0000", "#00B6D8", "#00CF15", "#BC00B5")) +
  theme_() + labs(x = "Similarity", y = "Predicted misidentified cases in each analysis", fill = "REGION", colour = "REGION") +
  theme(strip.text = element_text(size = 16, family = "Segoe UI Semibold"), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16), axis.text = element_text(size = 13),
        legend.title = element_text(size = 14, family = "Segoe UI Semibold"),
        legend.text = element_text(size = 13))






##### STATISTICAL ANALYSIS OF EFFECT OF LENGTH ON RESOLUTION AND IF THERE IS DIFFERENCE OF EFFECTS BETWEEN GENES #####

## Choosing model parameters

# Random effects (i.e. conditional variables of irrelevant importance here): 
# SIMILARITY -> main factor impacting resolution as shown above -> randomized to study other smaller effects
# PRIMERS -> randomized to detach each primer from its length and gene and study their effect separately

# For each level of similarity, an evaluation was made for the same metabarcodes -> Crossed random effects -> (1 | SIMILARITY) + (1 | PRIMERS)


## Running linear mixed-effects models

lmer_intra_bin_length = lmer(PERCENT_INTRA_ERRORS ~ MEAN_LENGTH * GENE +         # Intercept of primers can vary randomly
                               (1 | SIMILARITY) + (1 | PRIMERS), data = stat_analysis_resolution) # Intercept of similarity can vary randomly
lmer_inter_bin_length = lmer(PERCENT_INTER_ERRORS ~ MEAN_LENGTH * GENE +         
                               (1 | SIMILARITY) + (1 | PRIMERS), data = stat_analysis_resolution) 

shapiro.test(residuals(lmer_intra_bin_length)) # No normal distribution
shapiro.test(residuals(lmer_inter_bin_length)) # No normal distribution

plot(lmer_intra_bin_length, which = 1) # Heteroskedasticity + Outliers // But linearity seems to be quite satisfied
plot(lmer_inter_bin_length, which = 1) # Heteroskedasticity + Outliers // But linearity seems to be quite satisfied 


## Running robust linear mixed-effects models

rlmer_intra_bin_length = rlmer(PERCENT_INTRA_ERRORS ~ MEAN_LENGTH * GENE +         
                                 (1 | SIMILARITY) + (1 | PRIMERS), data = stat_analysis_resolution)
rlmer_inter_bin_length = rlmer(PERCENT_INTER_ERRORS ~ MEAN_LENGTH * GENE +         
                                 (1 | SIMILARITY) + (1 | PRIMERS), data = stat_analysis_resolution)

summary(rlmer_intra_bin_length)
summary(rlmer_inter_bin_length)

plot.rlmerMod(rlmer_intra_bin_length, which = 2)[[1]] + stat_qq_line(aes(sample = resid(rlmer_intra_bin_length))) + theme_()
plot.rlmerMod(rlmer_inter_bin_length, which = 2)[[1]] + stat_qq_line(aes(sample = resid(rlmer_inter_bin_length))) + theme_()
# Clear absence of normality in the data -> more robustness is given to values in the middle of the plot

plot(rlmer_intra_bin_length, which = 3)[[1]] + theme_()
plot(rlmer_inter_bin_length, which = 3)[[1]] + theme_()
# Random effects follow a more normal distribution thus almost all weigths are equal to 1

plot(rlmer_intra_bin_length, which = 1)[[1]] + theme_()
plot(rlmer_inter_bin_length, which = 1)[[1]] + theme_()
# Spread of residuals increase towards rigth of plot = heteroskedasticity
# Many outliers with large residuals values are likely to cause this heteroskedasticity -> small weigths
# Appart from these few values, linearity seems to hold reasonably well


## Checking R squared of the model

R2_rlmer(rlmer_intra_bin_length) # MR²: 43.56% vs CR²: 82.56%

R2_rlmer(rlmer_inter_bin_length) # MR²: 23.07% vs CR²: 83.32%


## Statistical analysis of fixed effects for intra-BIN resolution

stats.rlmerMod(rlmer_intra_bin_length)

# Intercept: Wrongly discrimined sequences sharing same BIN are not significantly different than 0 for COI metabarcodes

# Main effect of length: No changes in number of wrongly discrimined sequences for COI metabarcodes depending on length

# Intercept of CytB: +1.38% [0.94 ; 1.83] of sequences wrongly compared to COI (and in general as the intercept is not significant)

# Not significant difference of slope and intercept compared to COI of 12S and 16S metabarcodes

# Even if length alone does not explain the intra-BIN resolution, longer CytB metabarcodes are 
# sligthly, but significantly, more resolutive than shorter ones (-0.002% [-0.001% ; -0.003%]) .


## Statistical analysis of fixed effects for inter-BIN resolution

stats.rlmerMod(rlmer_inter_bin_length)

# Intercept: Confounded sequences sharing different BIN are not significantly different than 0 for COI metabarcodes

# Main effect of length: No changes in number of percent of counfonded sequences for COI metabarcodes depending on length

# Genes do not explain any significant increase or decrease of the percentage for COI metabarcodes

# Length of metabarcodes in each genes do not explain any differences of percentage (with COI, also not similar to 0)


## Plotting of results

length_rlmer_graph = tibble(rbind(cbind(stat_analysis_resolution, TYPE = "Intra-BIN", PREDICTED = predict(rlmer_intra_bin_length)), 
                                  cbind(stat_analysis_resolution, TYPE = "Inter-BIN", PREDICTED = predict(rlmer_inter_bin_length))))
length_rlmer_graph$TYPE = factor(length_rlmer_graph$TYPE, levels = c("Intra-BIN", "Inter-BIN"))

length_rlmer_plot = ggplot(length_rlmer_graph, aes(x = MEAN_LENGTH, colour = GENE)) +
  facet_wrap(~ TYPE, labeller = labeller(TYPE = as_labeller(c("Inter-BIN" = "Inter-BIN~(R[~(c)]^{2}=={0.833})",
                                                              "Intra-BIN" = "Intra-BIN~(R[~(c)]^{2}=={0.826})"), label_parsed))) +
  geom_smooth(data = subset(length_rlmer_graph, TYPE == "Intra-BIN"), method = "lm", 
              aes(y = PREDICTED, fill = GENE), alpha = 0.1, lwd = 0.8) +
  geom_smooth(data = subset(length_rlmer_graph, TYPE == "Inter-BIN"), method = "lm",
              aes(y = PREDICTED, fill = GENE), alpha = 0.1, lwd = 0.8) +
  scale_y_continuous(breaks = c(0, 1, 2), labels = function(x) paste0(x, "%")) + 
  scale_color_manual(values = c("#FF0000", "#00B6D8", "#00CF15", "#BC00B5")) +
  scale_fill_manual(values = c("#FF0000", "#00B6D8", "#00CF15", "#BC00B5")) +
  theme_() + labs(x = "Mean length", y = "Predicted misidentified cases in each analysis", fill = "REGION", colour = "REGION") +
  theme(strip.text = element_text(size = 16, family = "Segoe UI Semibold"), axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16), axis.text = element_text(size = 13),
        legend.title = element_text(size = 14, family = "Segoe UI Semibold"),
        legend.text = element_text(size = 13))

combined_rlmer_plot = grid.arrange(similarity_rlmer_plot, length_rlmer_plot, ncol = 1, nrow = 2,
                                   left = textGrob("Predicted misidentified cases in each analysis", 
                                                   rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 16)))

ggsave("Graph/analysis_differences_primers/Mixed effect analysis - Similarity and length.png", combined_rlmer_plot,
       device = png, dpi = 600, width = 18.97063, height = 20.78666, units = "cm")






##### STATISTICAL ANALYSIS OF THE RESOLUTION ACCORDING TO THE DIFFERENT ORDERS #####

## Choosing model parameters

similarity_analysis_resolution$ORDER = "All"
similarity_analysis_resolution_rlmer = subset(similarity_analysis_resolution, select = 
                                                 c("ORDER", "SIMILARITY", "PRIMERS", "PERCENT_INTRA_ERRORS", "PERCENT_INTER_ERRORS"))
similarity_analysis_resolution_rlmer

orders_analysis_resolution_rlmer = subset(orders_analysis_resolution, select = 
           c("ORDER", "SIMILARITY", "PRIMERS", "PERCENT_INTRA_ERRORS", "PERCENT_INTER_ERRORS"))
orders_analysis_resolution_rlmer

stat_analysis_orders = rbind(similarity_analysis_resolution_rlmer, orders_analysis_resolution_rlmer)

stat_analysis_orders = tibble(merge(mean_length_metabarcodes, stat_analysis_orders))

levels_primers = unique(stat_analysis_orders[order(sapply(stat_analysis_orders$GENE, function(x) 
  which(x == c("COI", "CytB", "12S", "16S")))), ]$PRIMERS)
levels_primers 

stat_analysis_orders$PRIMERS = factor(stat_analysis_orders$PRIMERS, levels = levels_primers) # FishF1-FishR1 will serve as reference

c("All", orders_to_study$ORDER) ; stat_analysis_orders$ORDER = factor(stat_analysis_orders$ORDER, # "All" will serve as reference
                                                                      levels = c("All", orders_to_study$ORDER))

# Random effects (i.e. conditional variables of irrelevant importance here): 
# SIMILARITY -> main factor impacting resolution as shown above -> randomized to study other smaller effects


## Running linear mixed-effects models

lmer_intra_bin_orders = lmer(PERCENT_INTRA_ERRORS ~ ORDER + ORDER:PRIMERS + 
                               (1 | SIMILARITY), data = stat_analysis_orders) # Intercept of similarity can vary randomly
lmer_inter_bin_orders = lmer(PERCENT_INTER_ERRORS ~ ORDER + ORDER:PRIMERS +         
                               (1 | SIMILARITY), data = stat_analysis_orders) 

shapiro.test(residuals(lmer_intra_bin_orders)) # No normal distribution
shapiro.test(residuals(lmer_inter_bin_orders)) # No normal distribution

plot(lmer_intra_bin_orders, which = 1) # Heteroskedasticity + Outliers // But linearity seems to be quite satisfied
plot(lmer_inter_bin_orders, which = 1) # Heteroskedasticity + Outliers // But linearity seems to be quite satisfied 


## Running robust linear mixed-effects models

rlmer_intra_bin_orders = rlmer(PERCENT_INTRA_ERRORS ~ ORDER + ORDER:PRIMERS  +         
                                 (1 | SIMILARITY), data = stat_analysis_orders)
rlmer_inter_bin_orders = rlmer(PERCENT_INTER_ERRORS ~ ORDER + ORDER:PRIMERS +         
                                 (1 | SIMILARITY), data = stat_analysis_orders)

summary(rlmer_intra_bin_orders)
summary(rlmer_inter_bin_orders)

plot.rlmerMod(rlmer_intra_bin_orders, which = 2)[[1]] + stat_qq_line(aes(sample = resid(rlmer_intra_bin_orders))) + theme_()
plot.rlmerMod(rlmer_inter_bin_orders, which = 2)[[1]] + stat_qq_line(aes(sample = resid(rlmer_inter_bin_orders))) + theme_()
# Clear absence of normality in the data -> more robustness is given to values in the middle of the plot

plot(rlmer_intra_bin_orders, which = 3)[[1]] + theme_()
plot(rlmer_inter_bin_orders, which = 3)[[1]] + theme_()
# Random effects follow a more normal distribution thus almost all weigths are equal to 1

plot(rlmer_intra_bin_orders, which = 1)[[1]] + theme_()
plot(rlmer_inter_bin_orders, which = 1)[[1]] + theme_()
# Spread of residuals increase towards rigth of plot = heteroskedasticity
# Many outliers with large residuals values are likely to cause this heteroskedasticity -> small weigths
# Appart from these few values, linearity seems to hold reasonably well


## Checking R squared of the model

R2_rlmer(rlmer_intra_bin_orders) # MR²: 75.26% vs CR²: 86.04%

R2_rlmer(rlmer_inter_bin_orders) # MR²: 82.62% vs CR²: 94.15%


## Statistical analysis of fixed effects for intra-BIN resolution

intra_bin_orders_estimate = stats.rlmerMod(rlmer_intra_bin_orders)
intra_bin_orders_estimate
# FishF1/R1	is a good reference because no significant negative difference

analysis_interaction_rlmer = function(rlmer_pvalue, reference, orders){
  
  p_apa = function(p) {
    
    if(p >= 0.001 && p < 0.01)
      p = "< .01"
    
    else if(p < 0.001)
      p = "< .001"
    
    else 
      p = paste("=", substr(round(p, digits = 3), 2, 5))
    p
    
  }
  
  s = function(p) {
    
    if(p >= 0.05 && p < 0.1)    significance = " (.)"
    
    else if(p >= 0.01 && p < 0.05)    significance = " (*)"
    
    else if(p >= 0.001 && p < 0.01)    significance = " (**)"
    
    else if(p < 0.001)    significance = " (***)"
    
    else significance = " (ns)"
    
    significance
    
  }
  
  output_list = list()
  
  for(i in 1:length(orders)){
    
    order_df = rlmer_pvalue[grep(paste0(orders[i], ":"), rownames(rlmer_pvalue)), ]
    
    order_df = tibble(ORDER = orders[i], PRIMERS = word(rownames(order_df), 2, sep = fixed("PRIMERS")),
                      DIFFERENCE_REFERENCE_TEMPORARY = order_df[,1], DIFFERENCE_REFERENCE = paste0(round(order_df[,1], 2), " %"),
                      CONF_INT_DIFFERENCE = paste0("95% CI [", round(order_df[,2], 2), " %, ", round(order_df[,3], 2), " %]"),
                      PVALUE = paste0("p-value ", sapply(order_df[,6], p_apa), "              "),
                      PVALUE_TEMPORARY = paste0(sapply(order_df[,6], s), "     "),
                      MAX_UNSIGNIFICANT = ifelse(!any(order_df[,6] >= 0.05), NA, max(order_df[which(order_df[,6] >= 0.05), 1])))
    
    order_df = tibble(merge(reference, order_df))
    
    order_df = order_df[order(order_df$DIFFERENCE_REFERENCE_TEMPORARY), ]
    
    colnames(order_df)[10] = ""
    
    output_list[[i]] = subset(order_df, select = -DIFFERENCE_REFERENCE_TEMPORARY) 
    
  }
  
  names(output_list) = orders
  
  return(output_list)
  
}

intra_bin_orders_estimate[1:length(orders_to_study$ORDER) + 1, ]
# More inter-BIN errors for Cypriniformes (3%) than when considering the full dataset of this study for FishF1-R1

intra_bin_orders_interaction = analysis_interaction_rlmer(intra_bin_orders_estimate, mean_length_metabarcodes, orders_to_study$ORDER)

intra_bin_orders_interaction[[1]] %>% print(n = 21)

intra_bin_orders_interaction[[2]] %>% print(n = 21)

intra_bin_orders_interaction[[3]] %>% print(n = 21)

intra_bin_orders_interaction[[4]] %>% print(n = 21)

intra_bin_orders_interaction[[5]] %>% print(n = 21)

sapply(intra_bin_orders_interaction, function(x) round(mean(as.numeric(word(x$DIFFERENCE_REFERENCE, 1))), 2))
sapply(intra_bin_orders_interaction, function(x) round(max(as.numeric(word(x$DIFFERENCE_REFERENCE, 1))), 2))
# Decrease of intra-BIN resolution: In general very important for Siluriformes (up to 4%), around 1% for Cypriniformes and 
# Perciformes, around 0.5% for Gobiiformes and inexistant for Tetraodontiformes (except L2513)

rank_order_rlmer = function(data_list){
  
  for(i in 1:length(data_list)){
    
    primers_order = data_list[[i]]$PRIMERS
    
    names(primers_order) = 1:length(primers_order)
    
    primers_order = sort(primers_order)
    
    primers_order_df = tibble(primers_order, names(primers_order))
    
    colnames(primers_order_df) = c("PRIMERS", paste0("RANK_", toupper(names(data_list)[i])))
    
    if(i == 1) primers_order_all = primers_order_df
    
    else primers_order_all = tibble(merge(primers_order_all, primers_order_df))
    
  }
  
  primers_order_all$RANK_SUM = rowSums(apply(primers_order_all[,-1], 2, as.numeric))
  
  primers_order_all = primers_order_all[order(primers_order_all$RANK_SUM), ]
  
  return(primers_order_all)
  
}

rank_order_rlmer(intra_bin_orders_interaction) %>% print(n = 21)
rank_order_rlmer(intra_bin_orders_interaction[c(1,4)]) %>% print(n = 21)

# Including only Cypriniformes and Siluriformes for which differences are marked
# Minibar and Ac12S are the best markers for intra-BIN resolution (also the case for the remaining orders)
# L2513 is very good for both orders but the only one significantly different from FishF1/R1 for Tetraodontiformes
# Teleo1 for Gobiiformes and 16SFD for Perciformes are also the only significantly different from FishF1/R1


## Statistical analysis of fixed effects for inter-BIN resolution

inter_bin_orders_estimate = stats.rlmerMod(rlmer_inter_bin_orders)
inter_bin_orders_estimate
# FishF1/R1 is a good reference because no significant negative difference

inter_bin_orders_estimate[1:length(orders_to_study$ORDER) + 1, ]
# More inter-BIN errors for Gobiiformes (2.5%) and Tetraodontiformes (1.6%) than when considering the full dataset of this study for FishF1-R1

inter_bin_orders_interaction = analysis_interaction_rlmer(inter_bin_orders_estimate, mean_length_metabarcodes, orders_to_study$ORDER)

inter_bin_orders_interaction[[1]] %>% print(n = 21)

inter_bin_orders_interaction[[2]] %>% print(n = 21)

inter_bin_orders_interaction[[3]] %>% print(n = 21)

inter_bin_orders_interaction[[4]] %>% print(n = 21)

inter_bin_orders_interaction[[5]] %>% print(n = 21)

sapply(inter_bin_orders_interaction, function(x) round(mean(as.numeric(word(x$DIFFERENCE_REFERENCE, 1))), 2))
sapply(inter_bin_orders_interaction, function(x) round(max(as.numeric(word(x$DIFFERENCE_REFERENCE, 1))), 2))
# Generally more inter-BIN errors than intra-BIN errors and much higher maximum errors 
# Less differences between orders however, compared to intra-BIN errors


rank_order_rlmer(inter_bin_orders_interaction) %>% print(n = 21)

# The best metbarcodes for all orders are L14912 and L14841
# FishCB is also good, and is the best for Gobiiformes



## Preparing the tables for plotting the intra-BIN and inter-BIN resolution for each orders

graph_interaction_rlmer = function(rlmer_pvalue, reference, orders){
  
  output_list = list()
  
  for(i in 1:length(orders)){
    
    order_df = rlmer_pvalue[grep(paste0(orders[i], ":"), rownames(rlmer_pvalue)), ]
    
    output_list[[i]] = tibble(ORDER = orders[i], 
                              PRIMERS = word(rownames(order_df), 2, sep = fixed("PRIMERS")),
                              DIFFERENCE_REFERENCE = order_df[,1], 
                              CONF_INT_INF_DIFFERENCE = order_df[,2],
                              CONF_INT_SUPP_DIFFERENCE = order_df[,3],
                              PVALUE = ifelse(order_df[,6] >= 0.05, "unsignif", "signif"),
                              MAX_UNSIGNIFICANT = ifelse(!any(order_df[,6] >= 0.05), NA, 
                                                         max(order_df[which(order_df[,6] >= 0.05), 1])))
    
  }
  
  return(tibble(do.call(rbind, output_list)))
  
}

intra_rlmer_graph = graph_interaction_rlmer(intra_bin_orders_estimate, mean_length_metabarcodes, orders_to_study$ORDER)
intra_rlmer_graph

inter_rlmer_graph = graph_interaction_rlmer(inter_bin_orders_estimate, mean_length_metabarcodes, orders_to_study$ORDER)
inter_rlmer_graph

intra_rlmer_graph = setNames(intra_rlmer_graph, c(colnames(intra_rlmer_graph)[1:2], 
                                                  paste0("INTRA_", colnames(intra_rlmer_graph)[3:7])))

inter_rlmer_graph = setNames(inter_rlmer_graph, c(colnames(inter_rlmer_graph)[1:2], 
                                                  paste0("INTER_", colnames(inter_rlmer_graph)[3:7])))

intra_inter_rlmer_graph = tibble(merge(intra_rlmer_graph, inter_rlmer_graph))
intra_inter_rlmer_graph


## Plotting of results for intra-BIN and inter-BIN resolution

intra_inter_rlmer_graph$PRIMERS = as.character(intra_inter_rlmer_graph$PRIMERS)

intra_inter_rlmer_graph = subset(intra_inter_rlmer_graph, !(PRIMERS %in% c("L14735c2", "Fish2deg")))

intra_inter_rlmer_graph$PRIMERS = replace(intra_inter_rlmer_graph$PRIMERS, 
                                          which(intra_inter_rlmer_graph$PRIMERS %in% "Fish2b"), "Fish2b/deg")

intra_inter_rlmer_graph$PRIMERS = replace(intra_inter_rlmer_graph$PRIMERS, 
                                          which(intra_inter_rlmer_graph$PRIMERS %in% "L14735c"), "L14735c/c2")

intra_inter_rlmer_graph$PRIMERS = replace(intra_inter_rlmer_graph$PRIMERS, 
                                          which(intra_inter_rlmer_graph$PRIMERS %in% "FishF1-FishR1"), "FishF1-R1")

intra_inter_rlmer_graph$PRIMERS = factor(intra_inter_rlmer_graph$PRIMERS, levels = 
                                           names(sort(sapply(split(intra_inter_rlmer_graph, intra_inter_rlmer_graph$PRIMERS), function(x) 
                                             sum(c(x$INTRA_DIFFERENCE_REFERENCE, x$INTER_DIFFERENCE_REFERENCE))))))

orders_names_space = setNames(paste0(c("Cyprini", "Gobii", "Perci", "Siluri", "Tetraodonti"), "-\nformes (", 
                                     orders_to_study$NUMBER_BIN_ORDER, " BIN)"), orders_to_study$ORDER)

intra_inter_rlmer_graph$ORDER = factor(intra_inter_rlmer_graph$ORDER, levels = 
                                         orders_to_study[order(-orders_to_study$NUMBER_BIN_ORDER),]$ORDER)

intra_rlmer_plot = ggplot(intra_inter_rlmer_graph, aes(x = PRIMERS, y = INTRA_DIFFERENCE_REFERENCE, colour = INTRA_PVALUE)) +
  facet_wrap(~ ORDER, nrow = length(unique(intra_inter_rlmer_graph$ORDER)), scales = "free_y", 
             strip.position = "right", labeller = as_labeller(orders_names_space)) +
  geom_hline(aes(yintercept = 0), size = 1, color = "black") +
  geom_errorbar(aes(ymin = INTRA_CONF_INT_INF_DIFFERENCE, ymax = INTRA_CONF_INT_SUPP_DIFFERENCE), 
                width = 0.2, position = position_dodge(0.05)) + geom_point() + 
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(-1, (max(x) + 1) * 1.1)))), 
                     minor_breaks = seq(-10, 50, 1), labels = function(x) paste0(x, "%")) +
  # geom_rect(aes(xmin = 0, xmax = length(unique(PRIMERS)) + 1, ymin = -Inf, ymax = INTRA_MAX_UNSIGNIFICANT), fill = "green", alpha = 0.01) +
  # geom_rect(aes(xmin = 0, xmax = length(unique(PRIMERS)) + 1, ymin = INTRA_MAX_UNSIGNIFICANT, ymax = Inf), fill = "red", alpha = 0.01) +
  # geom_rect(aes(xmin = 0, xmax = length(unique(PRIMERS)) + 1, ymin = -Inf, ymax = MIN_UNSIGNIFICANT), fill = "red", alpha = 0.01) +
  scale_colour_manual(values = c("#e41a1c", "#4daf4a")) +
  theme_() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = unit(c(0.05,0.05,-0.3,0.05), "cm")), 
                   axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none",
                   strip.text = element_text(size = 8.5, family = "Segoe UI Semibold"),
                   panel.grid.major.y = element_line(size = 0.1, color = "darkgrey"),
                   panel.grid.minor.y = element_line(size = 0.1, color = "darkgrey"))

inter_rlmer_plot = ggplot(intra_inter_rlmer_graph, aes(x = PRIMERS, y = INTER_DIFFERENCE_REFERENCE, colour = INTER_PVALUE)) +
  facet_wrap(~ ORDER, nrow = length(unique(intra_inter_rlmer_graph$ORDER)), scales = "free_y", 
             strip.position = "right", labeller = as_labeller(orders_names_space)) +
  geom_hline(aes(yintercept = 0), size = 1, color = "black") +
  geom_errorbar(aes(ymin = INTER_CONF_INT_INF_DIFFERENCE, ymax = INTER_CONF_INT_SUPP_DIFFERENCE), 
                width = 0.2, position = position_dodge(0.05)) + geom_point() + 
  scale_y_continuous(breaks = function(x) unique(floor(pretty(seq(-1, (max(x) + 1) * 1.1)))), 
                     minor_breaks = seq(-10, 50, 1), labels = function(x) paste0(x, "%")) +
  # geom_rect(aes(xmin = 0, xmax = length(unique(PRIMERS)) + 1, ymin = -Inf, ymax = INTER_MAX_UNSIGNIFICANT), fill = "green", alpha = 0.01) +
  # geom_rect(aes(xmin = 0, xmax = length(unique(PRIMERS)) + 1, ymin = INTER_MAX_UNSIGNIFICANT, ymax = Inf), fill = "red", alpha = 0.01) +
  # geom_rect(aes(xmin = 0, xmax = length(unique(PRIMERS)) + 1, ymin = -Inf, ymax = MIN_UNSIGNIFICANT), fill = "red", alpha = 0.01) +
  scale_colour_manual(values = c("#e41a1c", "#4daf4a")) +
  theme_() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, margin = unit(c(0.05,0.05,-0.3,0.05), "cm")), 
                   axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none",
                   strip.text = element_text(size = 8.5, family = "Segoe UI Semibold"),
                   panel.grid.major.y = element_line(size = 0.1, color = "darkgrey"),
                   panel.grid.minor.y = element_line(size = 0.1, color = "darkgrey"))

combined_rlmer_order_plot = grid.arrange(intra_rlmer_plot, inter_rlmer_plot, ncol = 2, nrow = 1,
                                         left = textGrob("Predicted difference in proportion of errors compared to FishF1-R1", 
                                                         rot = 90, gp = gpar(fontfamily = "Segoe UI Semilight", fontsize = 11)))

ggsave("Graph/analysis_differences_primers/Mixed effect analysis - Orders.png", combined_rlmer_order_plot,
       dpi = 600, type = "cairo", width = 25.53229, height = 14.39333, units = "cm", device = png)


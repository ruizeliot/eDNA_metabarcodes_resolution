\name{multiple_resolution_analyses}
\title{
Intra-taxa and inter-taxa analyses per subset
}
\description{
This function automatically performs intra-taxa and inter-taxa comparisons (\code{vsearch_pairwise_similarity} & \code{intra_taxa_similarity}) and analyses (\code{intra_taxa_analyses} & \code{inter_taxa_analyses}) for every given subset of the dataset, and return the summary per primers for each similarity threshold.
}
\usage{
multiple_resolution_analyses(metabarcodes_data, 
                             infos_data, 
                             elements_to_subset, 
                             similarity_threshold,
                             name_col_subset, 
                             name_col_taxa, 
                             name_col_primer, 
                             name_col_accession, 
                             name_col_fragment, 
                             keep_combinations = NULL, 
                             program_name, 
                             verbose = F)
}
\arguments{
  \item{metabarcodes_data}{ A table containing all DNA fragments (columns \code{name_col_fragment}) with their accession numbers (columns \code{name_col_accession}), for each primer (columns \code{name_col_primer}). }
  \item{infos_data}{ A table in which there is the column \code{name_col_accession}, the column \code{name_col_subset} and the column \code{name_col_taxa}. }
  \item{elements_to_subset}{ If each subset can done based on one element only, a vector can be provided with the names of the element to subset. Otherwise, it must be a list containing the names for each subset in each elements of the list. If the function also needs to be run with the entire dataset, the name provided must be \code{"All"}. If a subset must correspond to all other elements that were present in the other subset, its name must be \code{"Others"}. }
  \item{similarity_threshold}{ Numeric vector containing one or more similarity thresholds (in percentage) for which to perform the analyses, for each subset. }
  \item{name_col_subset}{ Name of the column of \code{infos_data} containing the names of the subset to which each sequence belong in order to subset the dataset based on \code{elements_to_subset}. }
  \item{name_col_taxa}{ Name of the column of \code{infos_data} containing the names of the taxa to which each sequence belong. It must be the lowest taxonomic rank of the dataset (i.e. no subtaxa within these taxa). }
  \item{name_col_primer}{ Name of the column containing the name(s) of the primer(s) surrounding each DNA fragment in \code{metabarcodes_data}. }
  \item{name_col_accession}{ Name of the column containing the accession number of each DNA fragment in \code{metabarcodes_data}. A column with the same name must also be present in \code{infos_data} to associate informations to each comparisons. }
  \item{name_col_fragment}{ Name of the column containing the DNA fragment in \code{metabarcodes_data}. }
  \item{keep_combinations}{ A vector containing all comparisons to keep before merging informations (i.e. reduced number of comparisons = reduced merging time). Each elements must be the two accession numbers separated by \code{':'}, in both orders (e.g. \code{c("AP011938:KP993137", "KP993137:AP011938", "KT211434:AB238967", "AB238967:KT211434")}). }
  \item{program_name}{ Name of the VSEARCH program. Add .exe for Windows executables. If not located in current working directory, add absolute or relative path as a prefix of the program name. }
  \item{verbose}{ To print the steps carried out by the different functions executed. }
}
\value{
This function return all the summary tables from \code{intra_taxa_analyses} and \code{inter_taxa_analyses} per subset, and per similarity thresholds binded together in a single table. The columns \code{OVERALL_MEAN_SIMILARITY} and \code{OVERALL_SD_SIMILARITY} are removed from these table since the same statistic cannot be computed for inter-taxa analysis (i.e. usually only relevant comparisons are retrieved). On the contrary, the column \code{CUMULATED_PERCENT_ERRORS_PER_ANALYSIS} is added to represent the sum of the percentages of intra-errors and inter-errors.
}
\references{
ARTICLE SOON PUBLISHED

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahe, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: leliotruiz1@gmail.com.
}
\seealso{
\code{\link{vsearch_pairwise_similarity}}, \code{\link{intra_taxa_similarity}}, \code{\link{intra_taxa_analysis}}, \code{\link{inter_taxa_analysis}}
}
\examples{
## Creating random data

primers = c(rep("PRIMER_1", 36), rep("PRIMER_2", 36), rep("PRIMER_3", 36))

random_metabarcodes = sapply(1:36, function(x) 
  paste(sample(c("A","G","T","C"), 1000, replace = T), collapse = ""))

random_accession = sapply(1:36, function(x) # Format of NCBI accession numbers
  paste0(c(sample(LETTERS, 2, replace = T), 
           sample(1:9, 6, replace = T)), collapse = ""))

random_species = split(taxize::plantNames, stringr::word(taxize::plantNames, 1))

random_species = random_species[which(lengths(random_species) > 3)]

random_species = lapply(random_species, function(x) sample(x, 3, replace = T))

random_species = unname(unlist(random_species[sample(1:length(random_species), 3)]))

random_genus = stringr::word(random_species, 1)


## Assembling random data in tables

random_metabarcodes_data = tibble::tibble(PRIMERS = primers,
                                          ACCESSION = rep(random_accession, 3),
                                          FRAGMENT = rep(random_metabarcodes, 3))
random_metabarcodes_data

random_infos_data = tibble::tibble(ACCESSION = random_accession,
                                   SPECIES = rep(random_species, 4),
                                   GENUS = rep(random_genus, 4),
                                   ZONE = paste0("ZONE_", rep(1:4, each = 9)))
random_infos_data


## Running resolution analysis separately for different subsets

subset_list = list(DATASET = "All", SUBSET_1 = "ZONE_1", 
                   SUBSET_2 = c("ZONE_2", "ZONE_3"), SUBSET_3 = "Others")
subset_list

subsets_analyses = 
  multiple_resolution_analyses(metabarcodes_data = random_metabarcodes_data, 
                               infos_data = random_infos_data, 
                               elements_to_subset = subset_list, 
                               similarity_threshold = 40:45,
                               name_col_subset = "ZONE", 
                               name_col_taxa = "SPECIES", 
                               name_col_primer = "PRIMERS", 
                               name_col_accession = "ACCESSION",  
                               name_col_fragment = "FRAGMENT", 
                               # Change it with the path + name of your VSEARCH program
                               program_name = "vsearch-2.17.1-win-x86_64.exe")
subsets_analyses
}
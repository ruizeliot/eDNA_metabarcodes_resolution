\name{vsearch_pairwise_similarity}
\title{
Fast pairwise similarity using VSEARCH
}
\description{
Performs pairwise alignments of all sequences provided to compute the similarity for each comparison. Then, add all necessary informations associated to each sequences to filter irrelevant comparisons.
}
\usage{
vsearch_pairwise_similarity(metabarcodes_data, 
                            infos_data = NULL, 
                            min_similarity_threshold, 
                            max_similarity_threshold, 
                            number_rows_considered = c(0,0), 
                            name_col_primer, 
                            name_col_accession, 
                            name_col_fragment, 
                            keep_combinations = NULL, 
                            keep_col_infos_data = NULL, 
                            order_col_keep = NULL, 
                            program_name, 
                            verbose = T)
}
\arguments{
  \item{metabarcodes_data}{ A table containing all DNA fragments (columns \code{name_col_fragment}) with their accession numbers (columns \code{name_col_accession}), for each primer (columns \code{name_col_primer}). }
  \item{infos_data}{ If \code{keep_col_infos_data} is not \code{NULL}, a table in which there is the columns provided in \code{keep_col_infos_data} as well as a column \code{name_col_accession} to identify the corresponding comparisons. }
  \item{min_similarity_threshold}{ Filters all comparisons below the given similarity threshold. See Details. }
  \item{max_similarity_threshold}{ Filters all comparisons above the given similarity threshold. See Details. }
  \item{number_rows_considered}{ By default \code{c(0,0)}. See Details. }
  \item{name_col_primer}{ Name of the column containing the name(s) of the primer(s) surrounding each DNA fragment in \code{metabarcodes_data}. }
  \item{name_col_accession}{ Name of the column containing the accession number of each DNA fragment in \code{metabarcodes_data}. If \code{bin_data} are \code{infos_data} provided (i.e. not \code{NULL}), a column with the same name must also be present in these tables to associate BIN and/or taxonomic informations. }
  \item{name_col_fragment}{ Name of the column containing the DNA fragment in \code{metabarcodes_data}. }
  \item{keep_combinations}{ A vector containing all comparisons to keep before merging informations (i.e. reduced number of comparisons = reduced merging time). Each elements must be the two accession numbers separated by \code{':'}, in both orders (e.g. \code{c("AP011938:KP993137", "KP993137:AP011938", "KT211434:AB238967", "AB238967:KT211434")}). }
  \item{keep_col_infos_data}{ Names of the columns from \code{sequence_data} (if not \code{NULL}) to keep when associating informations to each comparisons. }
  \item{order_col_keep}{ Character vector with the order of the columns from \code{infos_data} in the table returned. As each row of the output corresponds to two sequences (comparisons), all columns provided in \code{keep_col_infos_data} must be given with the suffix \code{'1'} and \code{'2'}. }
  \item{program_name}{ Name of the VSEARCH program. Add .exe for Windows executables. If not located in current working directory, add absolute or relative path as a prefix of the program name. }
  \item{verbose}{ To print the steps carried out by VSEARCH for each primers (dynamic in R console and static in RStudio). }
}
\details{
For each sequences, VSEARCH first sorts all other sequences by number of shared k-mers in a table, and then performs pairwise alignments row-wise by descending order. The first element of the vector \code{number_rows_considered} controls the maximum number of rows (i.e. comparisons) passing the similarity acceptation criteria to consider before stopping the search. The second element controls the maximum number of comparisons not passing the acceptation criteria to consider before stopping the search. The default value in VSEARCH is \code{c(1,32)} (useful for assignation) but the default value in \code{vsearch_pairwise_similarity} is \code{c(0,0)} to perform all the possible pairwise alignments as it is necessary in most resolution analysis.

The arguments \code{min_similarity_threshold} and \code{max_similarity_threshold} allow to filter comparisons outside of the required range of similarity (i.e. similarity acceptation criteria). As this filtering is performed by VSEARCH outside of R, before writting the text file, it avoids dealing with very large tables containing many unnecessary comparisons (i.e. reading and merging of infos very long).
}
\value{
A \code{\link[tibble:tibble]tibble} per primer (list) containing in each row the similarity (\code{"SIMILARITY"}) between two sequences (\code{"ACCESSION1"} & \code{"ACCESSION1"}), as well as the columns provided in \code{order_col_keep}.

Caution: The comparisons are always duplicated in the output (e.g. \code{"AP011938:KP993137"} & \code{"KP993137:AP011938"}).
}
\references{
ARTICLE SOON PUBLISHED

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahe, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{intra_taxa_similarity}}, \code{\link{inter_taxa_analysis}}, \code{\link{multiple_resolution_analysis}}
}
\examples{
## Creating random data

primers = c(rep("PRIMER_1", 9), rep("PRIMER_2", 9), rep("PRIMER_3", 9))

random_metabarcodes = sapply(1:9, function(x) 
  paste(sample(c("A","G","T","C"), 1000, replace = T), collapse = ""))

random_accession = sapply(1:9, function(x) # Format of NCBI accession numbers
  paste0(c(sample(LETTERS, 2, replace = T), 
           sample(1:9, 6, replace = T)), collapse = ""))

random_species = sample(taxize::plantNames, 9)

random_genus = stringr::word(random_species, 1)


## Assembling random data in tables

random_metabarcodes_data = tibble::tibble(PRIMERS = primers,
                                          ACCESSION = rep(random_accession, 3),
                                          FRAGMENT = rep(random_metabarcodes, 3))
random_metabarcodes_data

random_infos_data = tibble::tibble(ACCESSION = random_accession,
                                   SPECIES = random_species,
                                   GENUS = random_genus)
random_infos_data


## Running VSEARCH

comparisons_not_filtered = 
  vsearch_pairwise_similarity(metabarcodes_data = random_metabarcodes_data, 
                              infos_data = random_infos_data,  
                              min_similarity_threshold = 0, 
                              max_similarity_threshold = 100, 
                              name_col_primer = "PRIMERS", 
                              name_col_accession = "ACCESSION", 
                              name_col_fragment = "FRAGMENT",
                              keep_col_infos_data = c("SPECIES", "GENUS"),
                              order_col_keep = c("ACCESSION1", "ACCESSION2", 
                                                   "SPECIES1", "SPECIES2",
                                                   "GENUS1", "GENUS2"),
                              # Change it with the path + name of your VSEARCH program
                              program_name = "vsearch-2.17.1-win-x86_64.exe")
comparisons_not_filtered

relevant_comparisons = c(paste0(random_accession[1:2], collapse = ":"),
                         paste0(random_accession[2:1], collapse = ":"),
                         paste0(random_accession[2:3], collapse = ":"),
                         paste0(random_accession[3:2], collapse = ":"))
relevant_comparisons

comparisons_not_filtered = 
  vsearch_pairwise_similarity(metabarcodes_data = random_metabarcodes_data, 
                              infos_data = random_infos_data,  
                              min_similarity_threshold = 0, 
                              max_similarity_threshold = 100, 
                              keep_combinations = relevant_comparisons,
                              name_col_primer = "PRIMERS", 
                              name_col_accession = "ACCESSION", 
                              name_col_fragment = "FRAGMENT",
                              keep_col_infos_data = c("SPECIES", "GENUS"),
                              order_col_keep = c("ACCESSION1", "ACCESSION2", 
                                                   "SPECIES1", "SPECIES2",
                                                   "GENUS1", "GENUS2"),
                              # Change it with the path + name of your VSEARCH program
                              program_name = "vsearch-2.17.1-win-x86_64.exe")
comparisons_not_filtered
}

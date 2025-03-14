\name{inter_taxa_analysis}
\title{
Inter-taxa resolution analysis
}
\description{
Takes the output of \code{\link{vsearch_pairwise_similarity}} and determines the number of comparisons of sequences within the same taxa that were wrongly discrimined using one or more similarity thresholds. 
}
\usage{
inter_taxa_analysis(inter_comparisons_list, 
                    infos_data, 
                    similarity_threshold, 
                    name_col_similarity, 
                    name_col_accession,
                    name_col_taxa,
                    verbose = T)
}
\arguments{
  \item{inter_comparisons_list}{ The list of comparisons (one table per primers) returned by \code{\link{vsearch_pairwise_similarity}} and filtered to retain only comparisons of sequences belonging to different taxa. Each tables must contain the column \code{name_col_similarity}, as well as two columns containing the taxa named like \code{name_col_taxa} with the suffixes \code{1} and \code{2} (e.g. \code{"TAXA1"} & \code{"TAXA2"}). Each element of the list must be named with the name of the primers surrounding the DNA fragment compared. }
  \item{infos_data}{ A table with the accession numbers of every sequences (\code{name_col_accession}) and the taxa to which they belong (\code{name_col_taxa}). }
  \item{similarity_threshold}{ Numeric vector containing one or more similarity thresholds (in percentage) for which to perform the analyses. }
  \item{name_col_similarity}{ Name of the column containing the pairwise similarity in percentage in each tables of \code{inter_comparisons_list}. }
  \item{name_col_accession}{ Name of the columns containing the accession number in \code{inter_comparisons_list} and \code{infos_data}. }
  \item{name_col_taxa}{ Name of the columns containing the names of the taxa in \code{inter_comparisons_list} and \code{infos_data}.
  \item{verbose}{ To print the different steps of the analysis and the current similarity threshold used. }}
}
\details{
The table provided in \code{infos_data} is used to compute the total number of comparisons of sequences belonging to two different taxa. Then, the function count the number of wrongly undiscrimined sequences belonging to different taxa at a certain similarity threshold, not taking into account the duplicated cases. This allows to compute the percentage of inter-taxa errors with the total number of comparisons per combination of different taxa established before. Finally, it computes summary statistics such as the overall percentage of inter-taxa errors, the percentage of combinations of taxa (i.e. pairs) having some inter-taxa errors and the mean percentage of inter-taxa errors per combination of different taxa.
}
\value{
  \item{SUMMARY}{ A table for each similarity threshold (list) with each row corresponding to the summary of the inter-taxa analysis for one primer. }
  \item{DETAILS_UNDISCRIMINED_PAIRS}{ A table for each similarity threshold (list) with each row giving the number of wrongly undiscrimined sequences belonging to the combination of two different taxa, and the similarity in those cases. }
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{vsearch_pairwise_similarity}}, \code{\link{intra_taxa_analysis}}
}
\examples{
## Creating random data

primers = c(rep("PRIMER_1", 18), rep("PRIMER_2", 18), rep("PRIMER_3", 18))

random_metabarcodes = sapply(1:18, function(x) 
  paste(sample(c("A","G","T","C"), 1000, replace = T), collapse = ""))

random_accession = sapply(1:18, function(x) # Format of NCBI accession numbers
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
                                   SPECIES = rep(random_species, 2),
                                   GENUS = rep(random_genus, 2))
random_infos_data


## Getting all comparisons from between two sequences from different species

all_comparisons_possible = 
  vsearch_pairwise_similarity(metabarcodes_data = random_metabarcodes_data, 
                              infos_data = random_infos_data,  
                              min_similarity_threshold = 0, 
                              max_similarity_threshold = 100, 
                              name_col_primer = "PRIMERS", 
                              name_col_accession = "ACCESSION", 
                              name_col_fragment = "FRAGMENT",
                              keep_col_infos_data = c("SPECIES", "GENUS"),
                              order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                                   "SPECIES1", "SPECIES2",
                                                   "GENUS1", "GENUS2"),
                              # Change it with the path + name of your VSEARCH program
                              program_name = "vsearch-2.17.1-win-x86_64.exe")
all_comparisons_possible

inter_species_comparisons = lapply(all_comparisons_possible, function(x) 
  x[which(x$SPECIES1 != x$SPECIES2),])
inter_species_comparisons


## Performing the inter-species analysis

inter_species_analysis = 
  inter_taxa_analysis(inter_comparisons_list = inter_species_comparisons, 
                      infos_data = random_infos_data, 
                      similarity_threshold = 40:45, 
                      name_col_similarity = "SIMILARITY", 
                      name_col_accession = "ACCESSION", 
                      name_col_taxa = "SPECIES")
inter_species_analysis
}

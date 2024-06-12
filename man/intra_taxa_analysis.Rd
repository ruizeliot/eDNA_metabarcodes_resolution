\name{intra_taxa_analysis}
\title{
Intra-taxa analysis
}
\description{
Takes the output of \code{\link{intra_taxa_similarity}} established for the lowest taxonomic rank of the dataset and determines the number of comparisons of sequences within the same taxa that were wrongly discrimined using one or more similarity thresholds. 
}
\usage{
intra_taxa_analysis(intra_comparisons_list, 
                    infos_data, 
                    similarity_threshold, 
                    name_col_similarity, 
                    name_col_accession, 
                    name_col_taxa, 
                    verbose = T)
}
\arguments{
  \item{intra_comparisons_list}{ The list of comparisons (one table per primers) returned by \code{\link{intra_taxa_similarity}} with the argument \code{name_col_lower_rank_taxa = NULL}. Each tables must contain the column \code{name_col_similarity}, as well as two columns containing the accession number named like \code{name_col_accession} with the suffixes \code{1} and \code{2} (e.g. \code{"ACCESSION1"} & \code{"ACCESSION2"}). Two similar columns but for \code{name_col_taxa} must also be present in each tables. Each element of the list must be named with the name of the primers surrounding the DNA fragment compared. }
  \item{infos_data}{ A table with the accession numbers of every sequences (\code{name_col_accession}) and the taxa to which they belong (\code{name_col_taxa}). }
  \item{similarity_threshold}{ Numeric vector containing one or more similarity thresholds (in percentage) for which to perform the analyses. }
  \item{name_col_similarity}{ Name of the column containing the pairwise similarity in percentage in each tables of \code{intra_comparisons_list}. }
  \item{name_col_accession}{ Name of the columns containing the accession number in \code{intra_comparisons_list} and \code{infos_data}. }
  \item{name_col_taxa}{ Name of the columns containing the names of the taxa in \code{intra_comparisons_list} and \code{infos_data}. }
  \item{verbose}{ To print the different steps of the analysis and the current similarity threshold used. }
}
\details{
The table provided in \item{infos_data} is used to get the total number of comparisons of sequences within each taxa. Then the function count the number of wrongly discrimined sequences at a certain similarity threshold, and compute the percentage of intra-taxa errors with the total number of comparisons established before. Finally, it computes summary statistics such as the overall percentage of intra-taxa errors, the percentage of taxa containing some intra-taxa errors, the mean percentage of intra-taxa errors per taxa or the overall mean pairwise similarity.

The fact that every comparisons are duplicated in the output of \code{\link{intra_taxa_similarity}} is taken into account to compute all statistics.
}
\value{
  \item{SUMMARY}{ A table for each similarity threshold (list) with each row corresponding to the summary of the intra-taxa analysis for one primer. }
  \item{DETAILS_SIMILARITY}{ A table with summary statistics of the similarity obtained for all comparisons of sequences within each taxa and for all primers. }
  \item{DETAILS_PER_TAXA}{ A table for each similarity threshold (list) with each row giving the number of wrongly discrimined sequences within a certain taxa and the similarity in those cases. }
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: leliotruiz1@gmail.com.
}
\seealso{
\code{\link{intra_taxa_similarity}}, \code{\link{inter_taxa_analysis}}
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


## Running VSEARCH for every species

intra_species_comparisons = 
  intra_taxa_similarity(metabarcodes_data = random_metabarcodes_data, 
                        infos_data = random_infos_data, 
                        name_col_taxa = "SPECIES", 
                        name_col_primer = "PRIMERS", 
                        name_col_accession = "ACCESSION", 
                        name_col_fragment = "FRAGMENT",
                        keep_col_infos_data = c("SPECIES", "GENUS"),
                        order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                             "SPECIES1", "SPECIES2",
                                             "GENUS1", "GENUS2"),
                        # Change it with the path + name of your VSEARCH program
                        program_name = "vsearch-2.17.1-win-x86_64.exe")
intra_species_comparisons


## Performing the intra-species analysis

intra_species_analysis = 
  intra_taxa_analysis(intra_comparisons_list = intra_species_comparisons, 
                      infos_data = random_infos_data, 
                      similarity_threshold = 40:45, 
                      name_col_similarity = "SIMILARITY", 
                      name_col_accession = "ACCESSION", 
                      name_col_taxa = "SPECIES")
intra_species_analysis
}
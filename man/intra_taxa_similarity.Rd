\name{intra_taxa_similarity}
\title{
Fast intra-taxa pairwise similarity using VSEARCH
}
\description{
Compute the similarity (pairwise alignments using VSEARCH) between all sequences belonging to different subtaxa within a taxa, for each taxa of a given dataset. If the taxa used is the lowest taxonomic rank of the dataset (i.e. \code{name_col_lower_rank_taxa = NULL}), all sequences within each taxa (e.g. BIN from BOLD, species, strains) are compared.
}
\usage{
intra_taxa_similarity(metabarcodes_data, 
                      infos_data, 
                      name_col_taxa, 
                      name_col_lower_rank_taxa = NULL, 
                      program_name, 
                      name_col_order, 
                      name_col_primer, 
                      name_col_accession, 
                      name_col_fragment, 
                      keep_col_infos_data = NULL, 
                      order_col_keeped = NULL,
                      verbose = T)
}
\arguments{
  \item{metabarcodes_data}{ A table containing all DNA fragments (columns \code{name_col_fragment}) with their accession numbers (columns \code{name_col_accession}), for each primer (columns \code{name_col_primer}). }
  \item{infos_data}{ A table in which there is the columns provided in \code{keep_col_infos_data} as well as a column \code{name_col_accession} to identify the corresponding comparisons. The columns \code{keep_col_infos_data} and \code{name_col_lower_rank_taxa} (if not \code{NULL}) must also be present in this table. }
  \item{name_col_taxa}{ Name of the column of \code{infos_data} containing the names of the taxa for which to perform the intra-taxa analysis. }
  \item{name_col_lower_rank_taxa}{ Name of the column of \code{infos_data} containing the names of the subtaxa to compare within each \code{name_col_taxa}. Set it as \code{NULL} if the \code{name_col_taxa} is already the lowest taxonomic rank of the dataset. }
  \item{name_col_primer}{ Name of the column containing the name(s) of the primer(s) surrounding each DNA fragment in \code{metabarcodes_data}. }
  \item{name_col_accession}{ Name of the column containing the accession number of each DNA fragment in \code{metabarcodes_data}. A column with the same name must also be present in \code{infos_data} to associate informations to each comparisons. }
  \item{name_col_fragment}{ Name of the column containing the DNA fragment in \code{metabarcodes_data}. }
  \item{keep_col_infos_data}{ Names of the columns from \code{infos_data} (if not \code{NULL}) to keep when associating informations to each comparisons. }
  \item{order_col_keeped}{ Character vector with the order of the columns from \code{infos_data} in the table returned. As each row of the output corresponds to two sequences (comparisons), all columns provided in \code{keep_col_infos_data} must be given with the suffix \code{'1'} and \code{'2'}. }
  \item{program_name}{ Name of the VSEARCH program. Add .exe for Windows executables. If not located in current working directory, add absolute or relative path as a prefix of the program name. }
  \item{verbose}{ To print the name of the taxa for which comparisons are currently computed, and the number of relevant comparisons that will be retrieved. }
}
\details{
The function first compute the total number of relevant comparisons to be done for each taxa to inform the user and remove subtaxa without comparisons (i.e. only one sequence). Then it runs VSEARCH using \code{vsearch_pairwise_similarity} independently for each taxa and for each primers, and removes irrelevant comparisons at each step (i.e. comparisons between sequences sharing the same subtaxa). No acceptation criteria based on similarity is set, and all possible pairwise alignements are performed to retrieve all relevant comparisons per taxa.
}
\value{
A \code{\link[tibble:tibble]tibble} per primer (list) containing in each row the similarity (\code{"SIMILARITY"}) between two sequences (\code{"ACCESSION1"} & \code{"ACCESSION1"}), as well as the columns provided in \code{order_col_keeped}.

Caution: The comparisons are always duplicated in the output (e.g. \code{"AP011938:KP993137"} and \code{"KP993137:AP011938"}).
}
\references{
ARTICLE SOON PUBLISHED

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahe, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: leliotruiz1@gmail.com.
}
\seealso{
\code{\link{vsearch_pairwise_similarity}}, \code{\link{intra_taxa_analysis}}, \code{\link{multiple_resolution_analyses}}
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


## Running VSEARCH

intra_species_comparisons = 
  intra_taxa_similarity(metabarcodes_data = random_metabarcodes_data, 
                        infos_data = random_infos_data, 
                        name_col_taxa = "SPECIES", 
                        name_col_lower_rank_taxa = NULL,
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

intra_genus_comparisons = 
  intra_taxa_similarity(metabarcodes_data = random_metabarcodes_data, 
                        infos_data = random_infos_data,  
                        name_col_taxa = "GENUS", 
                        name_col_lower_rank_taxa = "SPECIES",
                        name_col_primer = "PRIMERS", 
                        name_col_accession = "ACCESSION", 
                        name_col_fragment = "FRAGMENT",
                        keep_col_infos_data = c("GENUS", "SPECIES"),
                        order_col_keeped = c("ACCESSION1", "ACCESSION2", 
                                             "SPECIES1", "SPECIES2",
                                             "GENUS1", "GENUS2"),
                        # Change it with the path + name of your VSEARCH program
                        program_name = "vsearch-2.17.1-win-x86_64.exe")
intra_genus_comparisons
}

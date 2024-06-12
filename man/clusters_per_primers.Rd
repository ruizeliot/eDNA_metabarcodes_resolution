\name{clusters_per_primers}
\title{
Clustering per threshold and per primer
}
\description{
Takes the intra-taxa errors and inter-taxa errors obtained using \code{\link{intra_taxa_similarity}} and \code{\link{vsearch_pairwise_similarity}} to quickly construct a distance matrix per primer. A dendrogram per primer is then created based on these matrices to assign each sequences to a given cluster, using multiple levels of cutoff (i.e. similarity thresholds). The construction of the dendrograms and the clustering is performed by the function \code{\link[DECIPHER:IdClusters]{IdClusters}}.
}
\usage{
clusters_per_primers(infos_data, 
                     intra_comparisons_list, 
                     inter_comparisons_list, 
                     similarity_threshold, 
                     tree_method, 
                     name_col_similarity, 
                     name_col_accession, 
                     name_col_taxa, 
                     nb_processors = NULL, 
                     verbose = T)
}
\arguments{
  \item{infos_data}{ A table with the accession numbers of every sequences (\code{name_col_accession}) and the taxa to which they belong (\code{name_col_taxa}). }
  \item{intra_comparisons_list}{ The list of comparisons (one table per primers) returned by \code{\link{intra_taxa_similarity}} with the argument \code{name_col_lower_rank_taxa = NULL}. Each tables must contain the column \code{name_col_similarity}, as well as two columns containing the accession number named like \code{name_col_accession} with the suffixes \code{1} and \code{2} (e.g. \code{"ACCESSION1"} & \code{"ACCESSION2"}). Two similar columns but for \code{name_col_taxa} must also be present in each tables. Each element of the list must be named with the name of the primer(s) surrounding the DNA fragment compared. }
  \item{inter_comparisons_list}{ The list of comparisons (one table per primers) returned by \code{\link{vsearch_pairwise_similarity}} and filtered to retain only comparisons of sequences belonging to different taxa. Each tables must contain the column \code{name_col_similarity}, as well as two columns containing the taxa named like \code{name_col_taxa} with the suffixes \code{1} and \code{2} (e.g. \code{"TAXA1"} & \code{"TAXA2"}). Each element of the list must be named with the name of the primer(s) surrounding the DNA fragment compared. }
  \item{tree_method}{ Grouping method used to construct the trees. Currently, the methods \code{"complete"}, \code{"single"}, \code{"UPGMA"}, \code{"WPGMA"}, \code{"NJ"}, and \code{"inexact"} are supported. See the help of \code{\link[DECIPHER:IdClusters]{IdClusters}} to known more about each method. }
  \item{similarity_threshold}{ Numeric vector containing one or more similarity thresholds (in percentage) for which to perform clustering. }
  \item{name_col_similarity}{ Name of the column containing the pairwise similarity in percentage in each tables of \code{intra_comparisons_list} and \code{inter_comparisons_list}. }
  \item{name_col_accession}{ Name of the columns containing the accession number in \code{intra_comparisons_list}, \code{inter_comparisons_list} and \code{infos_data}. }
  \item{name_col_taxa}{ Name of the columns containing the names of the taxa in \code{intra_comparisons_list}, \code{inter_comparisons_list}, and \code{infos_data}. }
  \item{nb_processors}{ Number of processors to use when constructing trees and performing the clustering, or \code{NULL} to automatically detect and use all available processors. }
  \item{verbose}{ To print the different steps of the analysis and the progress bar of the clustering for each primer. }
}
\details{
As only comparisons with similarities falling within the range of thresholds provided are relevant to perform clustering using these cutoffs, constructing similarity matrices does not require knowing all the pairwise similarities between all sequences of the dataset. Therefore, the function first computes the theoric list of all possible comparisons, and then assign to each comparisons of sequences belonging to the same taxa with similarity not comprised within the range of thresholds provided a similarity value of 100 (i.e. no intra-taxa error because all sequences within the same taxa are all similar). The same is done for comparisons not corresponding to inter-taxa errors but the value of similarity is 0 (i.e. no inter-taxa error because all sequences belonging to different taxa are correctly discrimined since do not share any common patterns). All remaining comparisons are relevant for the clustering, and the accurate values of similarities obtained after pairwise alignment are added in the matrices. All of this allows to directly use the output \code{\link{intra_taxa_similarity}} and \code{\link{vsearch_pairwise_similarity}} established with the same similarity threshold(s). It also greatly diminish computation time since it prevent manipulating very large tables (i.e. only relevant inter-taxa comparisons).
}
\value{
A \code{\link[tibble:tibble]tibble} per primer (list) containing in each row the accession number of a column and its membership to the different clusters delineated at the different similarity thresholds (one column for each threshold).
}
\references{
ARTICLE SOON PUBLISHED

Wright, E. S. (2016). Using DECIPHER v2. 0 to analyze big biological sequence data in R. R J., 8(1), 352.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: leliotruiz1@gmail.com.
}
\seealso{
\code{\link{vsearch_pairwise_similarity}}, \code{\link{intra_taxa_similarity}}
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


## Running VSEARCH to get intra and inter-species comparisons 

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


## Getting the cluster for each similarity thresholds and each primers

clusters_per_primers = 
  clusters_per_primers(infos_data = random_infos_data, 
                       intra_comparisons_list = intra_species_comparisons, 
                       inter_comparisons_list = inter_species_comparisons, 
                       similarity_threshold = 40:45, 
                       tree_method = "NJ", 
                       name_col_similarity = "SIMILARITY",
                       name_col_accession = "ACCESSION", 
                       name_col_taxa = "SPECIES")
clusters_per_primers
}

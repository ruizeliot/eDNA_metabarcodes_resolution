\name{cluster_metabarcodes_decipher}
\title{
Clustering per threshold and per metabarcode using DECIPHER
}
\description{
Takes the intra-taxa errors and inter-taxa errors obtained using \code{\link{intra_taxa_similarity}} and \code{\link{vsearch_pairwise_similarity}} to quickly construct a pairwise distance matrix per metabarcode. The different matrix-based clustering methods implemented in the DECIPHER functions \code{\link[DECIPHER:IdClusters]{IdClusters}} (old package version) or \code{\link[DECIPHER:Treeline]{Treeline}} (new package version) can then be run one by one for each similarity threshold to avoid creating a new matrix for each method. The new DECIPHER function \code{\link[DECIPHER:Clusterize]{Clusterize}} for fast linear clustering that does not requires as input a distance-matrix can be run separately on a set of unaligned metabarcodes for various similarity thresholds as well.}
\usage{
cluster_metabarcodes_decipher(unaligned_seqs_per_primer_df = NULL, 
                              alignment_per_primer_df = NULL, 
                              infos_data = NULL, 
                              intra_comparisons_list = NULL, 
                              inter_comparisons_list = NULL, 
                              similarity_threshold, 
                              tree_methods, 
                              max_time_clustering = NULL,
                              name_col_accession, 
                              name_col_fragment = NULL, 
                              name_col_alignment = NULL, 
                              name_col_primer = NULL,
                              name_col_similarity = NULL, 
                              name_col_taxa = NULL, 
                              nb_threads = NULL, 
                              verbose = T)
}
\arguments{
  \item{unaligned_seqs_per_primer_df}{ A table containing unaligned fragments (\code{name_col_fragment}) for each metabarcode (\code{name_col_primer}). This table is mandatory if \code{"tree_methods"} is set to \code{"clusterize"}. }
  \item{alignment_per_primer_df}{ A table containing aligned fragments (\code{name_col_fragment}) for each metabarcode (\code{name_col_primer}). This table is mandatory if \code{"tree_methods"} is set to \code{"ML"}, \code{"ME"}, or \code{"MP"}. }
  \item{infos_data}{ A table with the accession numbers of every sequences (\code{name_col_accession}) and the taxa to which they belong (\code{name_col_taxa}). }
  \item{intra_comparisons_list}{ The list of comparisons (one table per primers) returned by \code{\link{intra_taxa_similarity}} with the argument \code{name_col_lower_rank_taxa = NULL}. Each tables must contain the column \code{name_col_similarity}, as well as two columns containing the accession number named like \code{name_col_accession} with the suffixes \code{1} and \code{2} (e.g. \code{"ACCESSION1"} & \code{"ACCESSION2"}). Two similar columns but for \code{name_col_taxa} must also be present in each tables. Each element of the list must be named with the name of the primer(s) surrounding the DNA fragment compared. }
  \item{inter_comparisons_list}{ The list of comparisons (one table per primers) returned by \code{\link{vsearch_pairwise_similarity}} and filtered to retain only comparisons of sequences belonging to different taxa. Each tables must contain the column \code{name_col_similarity}, as well as two columns containing the taxa named like \code{name_col_taxa} with the suffixes \code{1} and \code{2} (e.g. \code{"TAXA1"} & \code{"TAXA2"}). Each element of the list must be named with the name of the primer(s) surrounding the DNA fragment compared. }
  \item{similarity_threshold}{ Numeric vector containing one or more similarity thresholds (in percentage) for which to perform clustering. }
  \item{tree_methods}{ Grouping method used to construct the trees. Currently, the methods supported are \code{"complete"}, \code{"single"}, \code{"UPGMA"}, \code{"WPGMA"}, \code{"NJ"}, \code{"ML"}, \code{"ME"}, \code{"MP"}, and \code{"clusterize"}. See the help of \code{\link[DECIPHER:IdClusters]{IdClusters}}, \code{\link[DECIPHER:Treeline]{Treeline}} and \code{\link[DECIPHER:Clusterize]{Clusterize}} to known more about each method. }
  \item{max_time_clustering}{ Maximum number of hours the clustering algorithm is allowed to run for the methods \code{"ML"}, \code{"ME"}, or \code{"MP"}. }
  \item{name_col_accession}{ Name of the columns containing the accession number in \code{intra_comparisons_list}, \code{inter_comparisons_list} and \code{infos_data}. }
  \item{name_col_fragment}{ Name of the column containing unaligned metabarcodes in \code{unaligned_seqs_per_primer_df} for use with the method \code{"clusterize"}. }
  \item{name_col_alignment}{ Name of the column containing aligned metabarcodes in \code{alignment_per_primer_df} for use with the methods \code{"ML"}, \code{"ME"}, or \code{"MP"}. }
  \item{name_col_primer}{ Name of the column containing the names of metabarcodes in \code{unaligned_seqs_per_primer_df} or \code{alignment_per_primer_df}. }
  \item{name_col_similarity}{ Name of the column containing the pairwise similarity in percentage in each table of \code{intra_comparisons_list} and \code{inter_comparisons_list}. }
  \item{name_col_similarity}{ Name of the column containing the pairwise similarity in percentage in each table of \code{intra_comparisons_list} and \code{inter_comparisons_list}. }
  \item{name_col_taxa}{ Name of the columns containing the names of the taxa in \code{intra_comparisons_list}, \code{inter_comparisons_list}, and \code{infos_data}. }
  \item{nb_threads}{ Number of processors to use when constructing trees and performing the clustering, or \code{NULL} to automatically detect and use all available processors. }
  \item{verbose}{ To print the different steps of the analysis and the progress bar of the clustering for each metabarcode. }
}
\details{
As only comparisons with similarities falling within the range of thresholds provided are relevant to perform clustering using these cutoffs, constructing similarity matrices does not require knowing all the pairwise similarities between all sequences of the dataset. Therefore, the function first computes the theoric list of all possible comparisons, and then assign to each comparisons of sequences belonging to the same taxa with similarity not comprised within the range of thresholds provided a similarity value of 100 (i.e. no intra-taxa error because all sequences within the same taxa are all similar). The same is done for comparisons not corresponding to inter-taxa errors but the value of similarity is 0 (i.e. no inter-taxa error because all sequences belonging to different taxa are correctly discrimined since do not share any common patterns). All remaining comparisons are relevant for the clustering, and the accurate values of similarities obtained after pairwise alignment are added in the matrices. All of this allows to directly use the output \code{\link{intra_taxa_similarity}} and \code{\link{vsearch_pairwise_similarity}} established with the same similarity threshold(s). It also greatly diminish computation time since it prevent manipulating very large tables (i.e. only relevant inter-taxa comparisons).

The function \code{\link[DECIPHER:Clusterize]{Clusterize}} proposes an alternative by sorting unaligned metabarcodes by relatedness (k-mer similarity) to linearize the clustering problem and therefore provide clusters in a reasonable amount of time.
}
\value{
A \code{\link[tibble:tibble]tibble} per primer (list) containing in each row the accession number of a column and its membership to the different clusters delineated at the different similarity thresholds (one column for each threshold).
}
\references{
ARTICLE SOON PUBLISHED

Wright, E. S. (2016). Using DECIPHER v2. 0 to analyze big biological sequence data in R. R J., 8(1), 352.

Wright, E. S. (2024). Growing phylogenetic trees with TreeLine.

Wright, E. (2024). Accurately clustering biological sequences in linear time by relatedness sorting. Nature Communications, 15(1), 3047.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{vsearch_pairwise_similarity}}, \code{\link{intra_taxa_similarity}}, \code{\link{cluster_metabarcodes_mothur_swarm}}}
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
                        order_col_keep = c("ACCESSION1", "ACCESSION2", 
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
                              order_col_keep = c("ACCESSION1", "ACCESSION2", 
                                                   "SPECIES1", "SPECIES2",
                                                   "GENUS1", "GENUS2"),
                              # Change it with the path + name of your VSEARCH program
                              program_name = "vsearch-2.17.1-win-x86_64.exe")
all_comparisons_possible

inter_species_comparisons = lapply(all_comparisons_possible, function(x) 
  x[which(x$SPECIES1 != x$SPECIES2),])
inter_species_comparisons


## Getting the cluster for each similarity thresholds and each primers using hierarchical clustering methods

clusters_per_primers = 
  cluster_metabarcodes_decipher(infos_data = random_infos_data, 
                                intra_comparisons_list = intra_species_comparisons, 
                                inter_comparisons_list = inter_species_comparisons, 
                                similarity_threshold = 40:45, 
                                tree_method = c("complete", "single", "UPGMA", "WPGMA", "NJ"), 
                                name_col_similarity = "SIMILARITY",
                                name_col_accession = "ACCESSION", 
                                name_col_taxa = "SPECIES")
clusters_per_primers
}

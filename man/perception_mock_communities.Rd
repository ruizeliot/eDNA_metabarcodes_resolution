\name{perception_mock_communities}
\alias{perception_mock_communities}
\title{Diversity bias evolution per in-silico mock communities}
\description{
The \code{perception_mock_communities} function simulates in-silico sequencing outputs from a wide range of mock communities to evaluate how metabarcode taxonomic errors can influence OTU estimates compared to a known taxonomic baseline. The simulations aim to replicate realistic community compositions and sequencing scenarios by considering variability in sequence counts and taxonomic diversity (more or less redundancy).
}
\usage{
perception_mock_communities(clusters_taxa_list,
                            nb_random_community,
                            nb_seq_sampled,
                            nb_seq_per_taxa = NULL,
                            clusters_cols_prefix,
                            name_col_taxa, 
                            name_col_accession,
                            output_path)
}
\arguments{
\item{clusters_taxa_list}{A list of tibbles containing clustering data for each primer set, formatted similarly and with the same sequence name and taxa names. Each tibble must include an \code{ACCESSION} column, a \code{TAXA} column, and clustering results across similarity thresholds.}
\item{nb_random_community}{An integer specifying the number of random community simulations to perform for each interval of \code{nb_seq_sampled}.}
\item{nb_seq_sampled}{A numeric vector defining the total number of sequences sampled for each mock community.}
\item{nb_seq_per_taxa}{An optional numeric vector specifying the uniform distribution in which to sample randomly the diversity (see Details). If \code{NULL}, defaults to \code{c(seq(0.5, 1, by = 0.1), seq(1, 2, by = 0.2)[-1])}.}
\item{clusters_cols_prefix}{ Prefix common to all clustering columns that correspond to various similarity thresholds. }
\item{name_col_taxa}{ Name of the column containing taxonomic information for each sequence. }
\item{name_col_accession}{ Name of the columns containing the accession number for each sequence. }
\item{output_path}{Local path to the directory in which to save simulations results (if not existing, the directory will be created).}
}

\details{
Simulations of \code{nb_random_community} mock communities per intervals of \code{nb_seq_sampled} (total number of sequences in the community), with random diversity (randomly picked between 0.5x and 2x the current \code{nb_seq_sampled}). To add some extent of redundancy that is typically lacking in DNA reference databases in order to minimize potential risks of false-positive errors underestimation, only the diversity of a given mock community is again randomly sampled between 1 and the current value while the total number of sequences is maintained the same. To do so, the number of unique sequences within each of the remaining taxa to reach again the initial number of sequences. Then, the number of perceived OTU obtained at various similarity thresholds is compared to the known diversity for each mock community to assess if the diversity bias is higher (overestimation of biodiversity of false-positive errors) or lower than 100\% (underestimation of biodiversity due to a predominance of false-negative errors).
}

\value{
The function outputs CSV files for each number of sequences in the mock community, containing the following columns:
\itemize{
\item \code{SEQ}: Total number of sequences sampled.
\item \code{TAXA}: Number of unique taxa observed.
\item \code{SIMILARITY}: Similarity thresholds used in clustering.
\item Perceived OTU counts per primer set for each threshold.
}
After performing mock community simulation, the functions summarizes the number of taxa, OTU and diversity bias (mean + sd in each cases) per total number of sequence, per metabarcode and per similarity threshold in a single summary table.
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{cluster_metabarcodes_decipher}}, \code{\link{cluster_metabarcodes_mothur_swarm}}
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
                                tree_method = c("single"), 
                                name_col_similarity = "SIMILARITY",
                                name_col_accession = "ACCESSION", 
                                name_col_taxa = "SPECIES")
clusters_per_primers


## Associating the taxonomic information to the clusters

clusters_per_primers_species = lapply(1:length(clusters_per_primers), function(i)
  tibble(merge(random_infos_data[,1:2], clusters_per_primers[[i]], by = "ACCESSION")))
names(clusters_per_primers_species) = names(clusters_per_primers)
clusters_per_primers_species[1:2]


## Performing the mock community simulation to determine the average diversity bias per metabarcode

species_perception_single_link = perception_mock_communities(clusters_per_primers_species, nb_random_community = 10, 
                            nb_seq_sampled = 4:8, clusters_cols_prefix = "CLUSTERS_", 
                            name_col_taxa = "ACCESSION", name_col_accession = "SPECIES", output_path = "Output")
species_perception_single_link
}
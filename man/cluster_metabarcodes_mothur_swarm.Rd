\name{cluster_metabarcodes_mothur_swarm}
\title{
Clustering per threshold and per metabarcode using MOTHUR and SWARM
}
\description{
Clusters a set of aligned metabarcodes with all methods available in MOTHUR and/or using SWARM, for various similarity thresholds and metabarcodes.
}
\usage{
cluster_metabarcodes_mothur_swarm(alignment_per_primer_df,
                                  output_path, 
                                  dist_threshold,
                                  cutoff_dist_matrix,
                                  methods = "all", 
                                  name_col_accession, 
                                  name_col_alignment, 
                                  name_col_primers, 
                                  mothur_path, 
                                  vsearch_path, 
                                  swarm_path, 
                                  nb_threads = 1,
                                  verbose = T)
}
\arguments{
  \item{alignment_per_primer_df}{ A table with alignments (\code{name_col_alignment}) performed for each metabarcode (more precise than using a gene alignment for example) for each ID (\code{name_col_accession}). Unaligned metabarcodes can also be provided in \code{name_col_alignment} but it will only work with \code{methods = "swarm"} as SWARM perfoms its own alignment. }
    \item{output_path}{ The relative path from the current directory in which output files from MOTHUR can be saved. }
    \item{dist_threshold}{ One or more distance threshold (between 0 and 1) for which to perform the clustering, given as a numeric vector. }
    \item{cutoff_dist_matrix}{ Maximum distance after which pairwise distance does not have to be reported for a smaller file size. Usually, it is best to use a value slightly higher than the maximum \code{dist_threshold} even though larger values or even complete matrices (1) can avoid hierarchical methods to fail. }
    \item{methods}{ One or more method to use (see descriptions in the section \details). The default is to use all methods, which is called by "all". }
    \item{name_col_accession}{ Name of the column containing the accession number of each DNA fragment in \code{alignment_per_primer_df}. }
    \item{name_col_alignment}{ Name of the column containing aligned sequences for each metabarcode in \code{alignment_per_primer_df}. }
    \item{name_col_primers}{ Name of the column containing name of metabarcodes in \code{alignment_per_primer_df}. }
    \item{mothur_path}{ Name of the MOTHUR program. Add .exe for Windows executables. If not located in current working directory, add absolute or relative path as a prefix of the program name. }
    \item{vsearch_path}{ Name of the VSEARCH program. Add .exe for Windows executables. If not located in current working directory, add absolute or relative path as a prefix of the program name. }
    \item{swarm_path}{ Name of the SWARM program. Add .exe for Windows executables. If not located in current working directory, add absolute or relative path as a prefix of the program name. }
    \item{nb_threads}{ Number of threads to use for MOTHUR distance matrix computation and SWARM clustering. The default is 1. }
    \item{verbose}{ To print the steps carried out by VSEARCH for each primers (dynamic in R console and static in RStudio). }
}
\details{
The program MOTHUR implements various different methods to perform clustering: 
\itemize{
\item{\strong{Nearest neighbor clustering (\code{method = "nearest"}):}} Each of the sequences within an OTU are at most x percent (distance threshold) distant from the most similar sequence in the OTU.
\item{\strong{Furthest neighbor clustering (\code{method = "furthest"}):}} All of the sequences within an OTU are at most x percent (distance threshold) distant from all of the other sequences within the OTU.
\item{\strong{Average neighbor clustering (\code{method = "average"}):}} This method is a middle ground between the first two algorithms, but often cannot perform clustering under a given distance threshold dependent on aligned sequences: "If one cell has a distance of 0.03 and the cell it is getting merged with has a distance above 0.05 (chosen distance threshold) then the cutoff is reset to 0.03, because it’s not possible to merge at a higher level and keep all the data."
\item{\strong{OPTICLUST clustering (\code{method = "opti"}):}} OTUs are assembled using metrics to determine the quality of clustering (the default setting in MOTHUR).
\item{\strong{Abundance-based greedy clustering (\code{method = "agc"}):}} Abundance-based greedy clustering implemented in VSEARCH (option --sizeorder).
\item{\strong{Distance-based greedy clustering (\code{method = "dgc"}):}} Distance-based greedy clustering implemented in VSEARCH (option --sizeorder). It generally gives very similar outputs than the AGC method.
\item{\strong{SWARM clustering (\code{method = "swarm"}):}} Fast and robust method that recursively groups amplicons based on a number of number of
mismatches (substitutions, insertions or deletions) in the optimal pairwise global alignment performed by SWARM (aligned sequences are therefore unaligned, but the function still needs them for compatibility with other methods). The number of mismatches is computed from the length of the longest unaligned fragment multiplied by the distance threshold. The abundance of each fragment is set to 1 as the purpose of this function is to compare reference metabarcodes without taking into account their abundance in PCR.
}
}
\value{
A \code{\link[tibble:tibble]tibble} per metabarcode (list) containing in each row the accession number of a column and its membership to the different clusters delineated at the different similarity thresholds (one column for each threshold), and for each method selected (column METHOD).
}
\references{
ARTICLE SOON PUBLISHED

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahe, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584.

Schloss PD et al. 2009. Introducing mothur: Open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and Environmental Microbiology 75:7537–7541.

Westcott SL, Schloss PD. 2017. OptiClust, an Improved Method for Assigning Amplicon-Based Sequence Data to Operational Taxonomic Units. mSphere 2:e00073-17.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{cluster_metabarcodes_decipher}}
}
\examples{
## Creating random aligned metabarcodes

random_aligned_primer = tibble(PRIMERS = c(rep("METABARCODE1", 4), rep("METABARCODE2", 4)), 
                               ACCESSION = rep(c("AAAA0000", "AAAA0001", "AAAA0002", "AAAA0003"), 2),
                               ALIGNMENT_PER_PRIMER = c("TTTTGTGGGTACGGAGTCATTAAAAATCGTGACTCGATTACTCGATAAGTTCT------TATCGAGAACCTAGAATAGATATGCTCCGGTCGGTTTAGGAGAGCAAGTTTTTTGCCCAAAA",
                                                        "-TTTGTGGGTACTGGGTCATTACAAATCGTCTGTCGATTACTCGATAAGTTCATTAATTTATCGAGAACCAAGTATATATATGCTCCGGTCGGTTCAG--GCGCAAGTC---TGCCCACAA",
                                                        "-TTTGTGGGTACGGAGTCATTACAAATCGTGAGTCGATTACTCGATAAGTTCATTA---TATCGAGAACCTAGAATAGATGTTCTCCGGTGGTTTCAG--GAGCAATTT---TGCCCAAAA",
                                                        "-TTTGTGCGTACGGAGTCATTACTCAACGTGAGTCGATTACTCGATAAGTTCA----TTTATCGAGAACCTAGAATAGATATGCTGCGGTCGGCTCAG-TGTTCAAGAT---TGCCCAAA-",
                                                        "ATGCT-AGGCTAGGCTAAC---GTTGATCGTACGATGGAATGCTAGTGCTTGGC-TTAACTGA----TGGAATCGGTTGACCGTGACTGGAATCGCGTTAACGTAGC", 
                                                        "-TGCTTAGGCTAGGCTAAG-GAACCGA-TCGTATCTGGAATG-TAGTGCTTAGG-TTAACCGATCGATGGAATCGGT--ACCGTGACTGGATAGCGT-TAACGTAGC", 
                                                        "ATGCTTAGG-TAGGCTAAT--AACCGATCACGATCTGGAATGCTAGTGCTTGGT--TAACCGATCGATGG--TCGGTTGACCGTGACTGGAATCTAG-CGTTAGTAC", 
                                                        "ATGCTTAGGCTAGGCTACGTTAACCGATTACGATCTGGAATGCTAGTGCTTAGGCTTAACCGATCGATGGAATCGGTTGACCGTGACTGGAATCGCG-TAACGTAGC"))


## Perform clustering on the basis of the global alignment with SWARM and MOTHUR (involving VSEARCH for the "dgc" method here)

example_mothur_clusters = cluster_metabarcodes_mothur_swarm(alignment_per_primer_df = random_aligned_primer, 
                                                            dist_threshold = seq(0.01, 0.10, by = 0.01), 
                                                            # These are the only methods that will work for such small dataset
                                                            methods = c("swarm", "dgc"), 
                                                            cutoff_dist_matrix = 0.1,
                                                            name_col_accession = "ACCESSION", 
                                                            name_col_alignment = "ALIGNMENT_PER_PRIMER", 
                                                            name_col_primers = "PRIMERS",
                                                            # Saving of MOTHUR output in the directory "Output"
                                                            output_path = "Output",
                                                            # Change it with the path + name of your MOTHUR program
                                                            mothur_path = "mothur.exe", 
                                                            # Change it with the path + name of your VSEARCH program
                                                            vsearch_path = "vsearch-2.17.1-win-x86_64.exe",
                                                            # Change it with the path + name of your SWARM program
                                                            swarm_path = "swarm.exe")
example_mothur_clusters
}
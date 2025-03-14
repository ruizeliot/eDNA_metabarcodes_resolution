\name{compute_otu_errors_metabarcodes}
\title{
Clustering errors compared to known taxonomy
}
\description{
Compute false-negatives (intra) and false-positives (inter) errors in clusters obtained using various metabarcodes, methods and thresholds by comparison to a taxonomic reference baseline such as BOLD's BINs.
}
\usage{
compute_otu_errors_metabarcodes(clusters_taxa_list,
                                name_col_method,
                                name_col_taxa,
                                clusters_cols_prefix,
                                verbose = TRUE)
}
\arguments{
\item{clusters_taxa_list}{ A named list of data frames where each data frame corresponds to a metabarcode and contains information about clustering results for various thresholds (columns starting by \code{clusters_cols_prefix}) and methods (\code{name_col_method}), along with known taxonomic assignments (\code{name_col_taxa}). }
\item{name_col_method}{ Name of the column indicating the clustering method used. }
\item{name_col_taxa}{ Name of the column containing taxonomic information for each sequence. }
\item{clusters_cols_prefix}{ Prefix common to all clustering columns that correspond to various similarity thresholds. }
\item{verbose}{ Logical. If \code{TRUE}, the function prints progress information during execution. }
}
\details{
For each metabarcode, the function iterates through the specified clustering columns (corresponding to different similarity thresholds). It computes:
\itemize{
\item{\strong{False-positives errors:} The number and percentage of taxa assigned to multiple clusters within the same taxa.}
\item{\strong{False-negatives errors:} The number and percentage of clusters containing multiple taxa.}
}
The function aggregates results for each clustering method and similarity threshold, providing a comprehensive summary of clustering accuracy.
}
\value{
A list with two elements:
\itemize{
\item{\code{INTRA_ERRORS}: A list of \code{tibble} summarizing false-positives errors for each similarity threshold.}
\item{\code{INTER_ERRORS}: A list of \code{tibble} summarizing fals-negatives errors for each similarity threshold.}
}
Each \code{tibble} includes columns:
\itemize{
\item{\code{PRIMERS}: Name of the primer.}
\item{\code{METHOD}: Clustering method used.}
\item{\code{NUMBER_INTRA_ERRORS} / \code{NUMBER_INTER_ERRORS}: Number of observed errors.}
\item{\code{PERCENT_INTRA_ERRORS} / \code{PERCENT_INTER_ERRORS}: Percentage of observed errors.}
\item{\code{MEAN_NUMBER_INTRA_ERRORS_PER_BIN} & \code{MEAN_NUMBER_INTER_ERRORS_PER_CLUSTER}: Mean number of errors per taxa or cluster.}
\item{\code{MEAN_PERCENT_INTRA_ERRORS_PER_BIN} & \code{MEAN_PERCENT_INTER_ERRORS_PER_CLUSTER}: Mean percentage of errors per taxa or cluster.}
}
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{cluster_metabarcodes_decipher}}, \code{\link{cluster_metabarcodes_mothur_swarm}}
}
\examples{
library(tibble)

## Creating random data

accessions = c(rep("AB123456", 10), rep("CD789012", 10))
bins = c(rep("BOLD:AAA0001", 10), rep("BOLD:AAA0002", 10))
methods = c("NJ", "furthest", "nearest", "opti", "average", "WPGMA", "cluster", "complete", "dgc", "single")


## Creating a list of clusters indices per similarity thresholds, clustering methods and metabarcodes

clusters_taxa_list = list(
  `PRIMER_1` = tibble(ACCESSION = accessions, BIN = bins, PRIMERS = rep("PRIMER_1", 20),
                      CLUSTERS_99 = sample(1:500, 20, replace = TRUE), 
                      CLUSTERS_98 = sample(1:500, 20, replace = TRUE), METHOD = rep(methods, 2)),
  `PRIMER_2` = tibble(ACCESSION = accessions, BIN = bins, PRIMERS = rep("PRIMER_2", 20),
                      CLUSTERS_99 = sample(1:500, 20, replace = TRUE), 
                      CLUSTERS_98 = sample(1:500, 20, replace = TRUE), METHOD = rep(methods, 2)))
clusters_taxa_list


## Summary of intra-BIN (false-positives) and inter-BIN (false-negatives)

summary_intra_inter_errors = compute_otu_errors_metabarcodes(clusters_taxa_list, name_col_method = "METHOD", 
                                                             name_col_taxa = "BIN", clusters_cols_prefix = "CLUSTERS_")
summary_intra_inter_errors$INTRA_ERRORS
summary_intra_inter_errors$INTER_ERRORS
}
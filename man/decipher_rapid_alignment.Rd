\name{decipher_rapid_alignment}
\title{
Rapid alignment of DNA fragments per group using DECIPHER
}
\description{
This function automatise the construction of a chained guide tree (i.e. dendrogram) to fasten the alignment of long and numerous DNA sequences using the function \code{\link[DECIPHER:AlignSeqs]AlignSeqs} from the DECIPHER package (Wright, 2020). It also allows to repeat this operation for each group of metabarcodes (e.g. per primers or per genes), and save the alignments at each step.
}
\usage{
decipher_rapid_alignment(data, 
                         name_col_group, 
                         name_col_accession, 
                         name_col_fragment, 
                         progressive_saving = F, 
                         saving_name = NULL, 
                         nb_processors = NULL, 
                         verbose = T,
                         ...)
}
\arguments{
  \item{data}{ A table containing the columns \code{name_col_group}, \code{name_col_accession} and \code{name_col_fragment}. }
  \item{name_col_group}{ Name of the column in \code{data} containing the groups of fragments for which to perform the alignments separately. }
  \item{name_col_accession}{ Name of the column in \code{data} containing the accession numbers of each full sequence from NCBI. }
  \item{name_col_fragment}{ Name of the column in \code{data} containing the DNA fragments to align as strings (i.e. single character vectors). }
  \item{progressive_saving}{ To save the alignment of every sequences within the same group as a CSV file, before starting the alignment of the next group. It avoids restarting the alignment from scratch if it fails at one point for some reasons. }
  \item{nb_processors}{ Number of processors to use when aligning sequences, or \code{NULL} to automatically detect and use all available processors. }
  \item{verbose}{ To print the progress bar of the alignment for each group. }
  \item{...}{ Further arguments to be passed directly to \code{\link[DECIPHER:AlignSeqs]AlignSeqs}, excluding \code{myXStringSet}, \code{guideTree}, \code{iterations}, \code{refinements}, \code{guideTree}, \code{processors} and \code{verbose}. }
}
\value{
A \code{\link[tibble:tibble]tibble} containing in each row the accession numbers of the full DNA sequences along with the different aligned fragments (\code{name_col_fragment}) and the groups to which they belong (\code{name_col_group}).
}
\references{
ARTICLE SOON PUBLISHED

Wright, E. S. (2020). The art of multiple sequence alignment in R. The Biconductor, 29.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: leliotruiz1@gmail.com.
}
\seealso{
\code{\link{metabarcode_extraction}}
}
\examples{
## Creating semi-random sequences (conserved + random parts)

conserved_parts = lapply(1:3, function(x) sapply(1:10, function(y) 
  paste0(sample(c("AAA", "TTT", "GGG", "CCC"), 10, replace = T), collapse = "")))

random_parts = lapply(1:3, function(x) lapply(1:5, function(y) sapply(1:10, function(z) 
  paste0(sample(c("A","G","T","C"), sample(25:35, 1), replace = T), collapse = ""))))

semi_random_sequences = lapply(1:3, function(x) sapply(1:5, function(y) 
  paste0(c(rbind(conserved_parts[[x]], random_parts[[x]][[y]])), collapse = "")))

primers = c("PRIMER_1", "PRIMER_2", "PRIMER_3")

random_accession = sapply(1:5, function(x) # Format of NCBI accession numbers
  paste0(c(sample(LETTERS, 2, replace = T), 
           sample(1:9, 6, replace = T)), collapse = ""))


## Assembling random sequences in a table

semi_random_df = do.call(rbind, lapply(1:3, function(x) 
  tibble::tibble(PRIMERS = primers[x],
                 ACCESSION = random_accession,
                 FRAGMENT = semi_random_sequences[[x]])))
semi_random_df


## Aligning every fragments for each primer

alignments_per_primer = 
  decipher_rapid_alignment(data = semi_random_df, 
                           name_col_group = "PRIMERS", 
                           name_col_accession = "ACCESSION",
                           name_col_fragment = "FRAGMENT")
alignments_per_primer


## Visualising the alignments for each primers

DECIPHER::BrowseSeqs(Biostrings::DNAStringSet(
  subset(alignments_per_primer, PRIMERS == "PRIMER_1")$ALIGNED_SEQUENCES))

DECIPHER::BrowseSeqs(Biostrings::DNAStringSet(
  subset(alignments_per_primer, PRIMERS == "PRIMER_2")$ALIGNED_SEQUENCES))

DECIPHER::BrowseSeqs(Biostrings::DNAStringSet(
  subset(alignments_per_primer, PRIMERS == "PRIMER_3")$ALIGNED_SEQUENCES))
}
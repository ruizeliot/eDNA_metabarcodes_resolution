\name{view_metabarcodes}
\title{
Visual checking of automatically extracted metabarcodes
}
\description{
This function aims at easing the checking of the automatic metabarcodes extraction performed by \code{\link{metabarcode_extraction}}. For each primer and each different results of the primer search (i.e. none, one or both primers found), it shows the alignment between the primers, some automatically extracted metabarcodes, and possibly some manually extracted metabarcodes, on the original full DNA fragment. It also compares the length of extracted metabarcodes with the reference length(s).
}
\usage{
view_metabarcodes(data, 
                  manual_data = NULL, 
                  data_primer, 
                  fragment_data, 
                  nb_sequences = 1, 
                  show_all = F, 
                  size_margin = 15, 
                  name_col_fragment, 
                  name_col_accession, 
                  name_col_manual_primer_name = NULL, 
                  name_col_manual_sequences = NULL,
                  name_col_primer_name, 
                  name_col_orientation, 
                  name_col_forward_primer, 
                  name_col_reverse_primer, 
                  length1_col_name, 
                  length2_col_name = NULL, 
                  name_col_data_metabarcodes = "FRAGMENT"
                  name_col_data_primers = "PRIMERS", 
                  name_col_data_found = "PRIMER_FOUND", 
                  name_col_data_length = "LENGTH", 
                  name_col_data_length_with_primers = "LENGTH_WITH_PRIMERS")
}
\arguments{
  \item{data}{ A table, typically the output of \code{\link{metabarcode_extraction}}, containing the columns \code{name_col_accession}, \code{name_col_fragment}, \code{name_col_data_primers}, \code{name_col_col_data_found}, \code{name_col_data_length} and \code{name_col_data_length_with_primers}. The last four columns correspond by default to the names given by \code{\link{metabarcode_extraction}}. } 
  \item{manual_data}{ A table containing manually extracted metabarcodes in \code{name_col_manual_sequences} for each primers (\code{name_col_manual_primer_name}) and each manually extracted metabarcodes \code{name_col_accession}. If manual extraction of metabarcodes was not performed, let it set to \code{NULL} (default value). } 
  \item{data_primer}{ A table containing the names of each primer set (\code{name_col_primer_name}) and the sequences of the forward (\code{name_col_forward_primer}) and reverse primers (\code{name_col_reverse_primer}). It must also contain the column \code{name_col_orientation} to precise explicitely the order of the primers. A reference length for the fragment must also be given in \code{length1_col_name}, and a second can be provided in \code{length2_col_name}. }
  \item{fragment_data}{ A table containing the full unaligned DNA sequences in which metabarcodes were extracted (\code{name_col_fragment}) as well as their NCBI accession numbers (\code{name_col_accession}). } 
  \item{nb_sequences }{ Number of randomly selected accession numbers to show in the alignment. If \code{show_all = F}, it is not necessary to specify more than \code{4} if \code{manual_data = NULL}, or more than \code{2} if manually extracted are provided, because only the first and last 5 sequences (10 in total) of the alignment will be shown in the console. }
  \item{show_all}{ A boolean to indicate if the alignment should be viewed in an HTML browser (\code{\link[DECIPHER:BrowseSeqs]BrowseSeqs}) to check more than \code{2} or \code{4} sequences, and view the full length of the alignment. } 
  \item{size_margin}{ An integer indicating the number of nucleotides to print before the forward primer and after the reverse primer in the console. } 
  \item{name_col_fragment}{ Name of the column in \code{fragment_data} containing the full unaligned DNA sequences in which metabarcodes were extracted. } 
  \item{name_col_accession}{ Name of the column in \code{data}, \code{manual_data} (if not \code{NULL}) and \code{fragment_data} containing NCBI accession numbers.  } 
  \item{name_col_manual_primer_name}{ Name of the column in \code{manual_data} containing the names of the primer set used to manually extract each metabarcodes. It will not be used if \code{manual_data = NULL} (by default set to \code{NULL}). } 
  \item{name_col_manual_sequences}{ Name of the column in \code{manual_data} containing the manually extracted metabarcodes. It will not be used if \code{manual_data = NULL} (by default set to \code{NULL}). }
  \item{name_col_primer_name}{ Name of the column in \code{data_primer} containing the names of the primers for which to extract corresponding metabarcodes. }
  \item{name_col_orientation}{ Name of the column in \code{data_primer} explicitely specifying the order of each primer set. It can be \code{"2 x 1"} if the sequences in \code{name_col_forward_primer} correspond to reverse primers, or \code{"1 x 2"} otherwise. Indeed, primers given as forward primers in the original article turn out to be reverse primers sometimes (and inversely for reverse primers), and this argument allows to simply test which configuration is accurate. }
  \item{name_col_forward_primer}{ Name of the column in \code{data_primer} containing the sequences of all forward primers to use. }
  \item{name_col_reverse_primer}{ Name of the column in \code{data_primer} containing the sequences of all reverse primers to use. }
  \item{length1_col_name}{ Name of the column in \code{data_primer} containing a known reference length for the metabarcode alone (\code{remove_primers = T}) or for the metabarcode with primers (\code{remove_primers = F}). }
  \item{length2_col_name}{ A second optional argument (otherwise set to \code{NULL}) to specify a second column name of \code{data_primer} containing the known length of fragments from another reference. In this case, the mean between both lengths is taken, excluding \code{NA}. }
  \item{name_col_data_primers}{ Name of the column in \code{data} containing the names of the primers. Set by default with the column name of the table returned by \code{\link{metabarcode_extraction}}.   } 
  \item{name_col_data_found}{ Name of the column in \code{data} where the result of the primer search is stored. Set by default with the column name of the table returned by \code{\link{metabarcode_extraction}}. } 
  \item{name_col_data_length}{ Name of the column in \code{data} containing the length of the extracted fragment. Set by default with the column name of the table returned by \code{\link{metabarcode_extraction}}. } 
  \item{name_col_data_length_with_primers}{ Name of the column in \code{data} containing of the length of the extracted fragment with primers. Set by default with the column name of the table returned by \code{\link{metabarcode_extraction}}. }
  \item{name_col_data_metabarcodes}{ Name of the column in \code{data} containing of the extracted fragment. Set by default with the column name of the table returned by \code{\link{metabarcode_extraction}}. }
}
\value{
For each primer and each different results of the primer search (i.e. none, one or both primers found), it prints the alignment of all selected sequences as a \code{\link[Biostrings:DNAStringSet]DNAStringSet} in the console and shows a summary of fragment lengths compared with the known reference(s).
}
\references{
ARTICLE SOON PUBLISHED

Wright, E. S. (2016). Using DECIPHER v2.0 to analyze big biological sequence data in R. The R Journal, 8 (1), 352–359.
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{metabarcode_extraction}}
}
\examples{
## Creating a primer with 1 random mismatch far from the 3' end

mismatch_primer1 = lapply(1:10, function(x) sample(c("A", "T", "G", "C"), 1))

forward_primer1 = lapply(mismatch_primer1, function(x) 
  paste0(c(sample(c(x, rep("A", 9))), rep("A", 10)), collapse = ""))

reverse_primer1 = lapply(mismatch_primer1, function(x) 
  paste0(c(rep("T", 10), sample(c(x, rep("T", 9)))), collapse = ""))


## Creating a primer with 2 random mismatches and/or 1 indel far from the 3' end

mismatch_primer2 = lapply(1:10, function(x) sample(c("A", "T", "G", "C", ""), 1))

mismatch_primer2 = lapply(1:10, function(x) c(mismatch_primer1[[x]], mismatch_primer2[[x]]))

forward_primer2 = lapply(mismatch_primer2, function(x) 
  paste0(c(sample(c(x, rep("G", 8))), rep("G", 10)), collapse = ""))

reverse_primer2 = lapply(mismatch_primer2, function(x) 
  paste0(c(rep("C", 10), sample(c(x, rep("C", 8)))), collapse = ""))


## Creating metabarcodes with one random mismatch on the middle

metabarcodes_primer1 = lapply(1:10, function(x)
  paste0(c(rep("CGCG", 3), sample(c("A", "T"), 1), rep("CGCG", 3)), collapse = ""))


## Creating metabarcodes with many random mismatches and indels everywhere

metabarcodes_primer2 = sapply(1:10, function(x) 
  paste0(rep(c("AT", rbind(rep("AT", 7), sample(c("A", "T", "G", "C", ""), 7, replace = T)), 
               "AT"), 4), collapse = ""))


## Creating a gene with metabarcodes + primers placed in conserved areas

conserved_parts = sapply(1:7, function(y) 
  paste0(sample(c("AA", "TT", "GG", "CC"), 10, replace = T), collapse = ""))

random_parts = lapply(1:10, function(y) sapply(1:5, function(z) 
  paste0(rep(sample(c("A","G","T","C"), 20, replace = T), sample(1:3, 1)), collapse = "")))

conserved_parts_and_metabarcodes = lapply(1:10, function(y)
  c(conserved_parts[1], paste0(c(conserved_parts[2], "GGGGG",
                                 forward_primer1[[y]], 
                                 metabarcodes_primer1[[y]], 
                                 reverse_primer1[[y]],
                                 "GGGGG", conserved_parts[3]), collapse = ""),
    conserved_parts[4], paste0(c(conserved_parts[5], "AAAAA",
                                 forward_primer2[[y]], 
                                 metabarcodes_primer2[[y]], 
                                 reverse_primer2[[y]],
                                 "AAAAA", conserved_parts[6]), collapse = ""),
    conserved_parts[3]))

simulated_gene = sapply(1:10, function(y) 
  paste0(c(rbind(conserved_parts_and_metabarcodes[[y]], 
                 random_parts[[y]])), collapse = ""))


## Aligning all the sequences within this gene

random_accession = sapply(1:10, function(x) # Format of NCBI accession numbers
  paste0(c(sample(LETTERS, 2, replace = T), 
           sample(1:9, 6, replace = T)), collapse = ""))

gene_alignment_dna = DECIPHER::AlignSeqs(setNames(Biostrings::DNAStringSet(unlist(simulated_gene)), 
                                                   random_accession),
                                          iterations = 5, refinements = 3)
gene_alignment_dna

DECIPHER::BrowseSeqs(Biostrings::DNAStringSet(gene_alignment_dna))

gene_alignment = 
  tibble::tibble(ACCESSION = names(gene_alignment_dna),
                 ALIGNED_SEQUENCES = BiocGenerics::paste(gene_alignment_dna))
gene_alignment


## Extracting both primers from the gene

infos_primers = tibble::tibble(NAME = c("PRIMER1", "PRIMER2"), ORIENTATION = "1 x 2", 
                               KNOWN_LENGTH_METABARCODE = c(25,100),
                               FORWARD = c(paste0(rep("A", 20), collapse = ""), 
                                           paste0(rep("G", 20), collapse = "")), 
                               REVERSE = c(paste0(rep("A", 20), collapse = ""), 
                                           paste0(rep("G", 20), collapse = "")))
infos_primers

extracted_primers = 
  metabarcode_extraction(data_primer = infos_primers, 
                          data_aligned = gene_alignment, 
                          remove_primers = T, 
                          max_mismatch_whole_primer = 2, 
                          indels = T, 
                          name_col_alignement = "ALIGNED_SEQUENCES", 
                          name_col_accession = "ACCESSION", 
                          name_col_primer_name = "NAME", 
                          name_col_orientation = "ORIENTATION", 
                          name_col_forward_primer = "FORWARD", 
                          name_col_reverse_primer = "REVERSE",
                          length1_col_name = "KNOWN_LENGTH_METABARCODE", verbose = F)
extracted_primers


## Checking if the extraction was correct

correct_metabarcodes_data = 
  tibble::tibble(ACCESSION = rep(random_accession, 2),
                 PRIMERS = c(rep("PRIMER1", 10), rep("PRIMER2", 10)),
                 LENGTH = c(sapply(metabarcodes_primer1, nchar), 
                            sapply(metabarcodes_primer2, nchar)),
                 CORRECT_METABARCODES = c(unlist(metabarcodes_primer1), 
                                          unlist(metabarcodes_primer2)))
correct_metabarcodes_data

view_metabarcodes(data = extracted_primers, 
                  manual_data = correct_metabarcodes_data, 
                  data_primer = infos_primers, 
                  fragment_data = tibble::tibble(ACCESSION = random_accession,
                                                 UNALIGNED_GENE = unlist(simulated_gene)), 
                  nb_sequences = 2,
                  name_col_fragment = "UNALIGNED_GENE", 
                  name_col_accession = "ACCESSION", 
                  name_col_manual_primer_name = "PRIMERS", 
                  name_col_manual_sequences = "CORRECT_METABARCODES",
                  name_col_primer_name = "NAME", 
                  name_col_orientation = "ORIENTATION",
                  name_col_forward_primer = "FORWARD", 
                  name_col_reverse_primer = "REVERSE",
                  length1_col_name = "KNOWN_LENGTH_METABARCODE")

# This little example aims at showing the advantages and disadvantages of our method.
# Indeed, there is very few sequences and the region of PRIMER2 contains many mismatches
# and indels, with very little chances for PRIMER2 to match perfectly the primer region.
# In this context, automatic alignment often fails and manual correction are needed.
# It is why some metabarcodes between both PRIMER2 are sometimes a bit too short/long.

# With many sequences, manual correction would be too long, and the advantage of our 
# method is therefore to rely on well-aligned and not too much degenerated primers
# regions, that should represent the majority of sequences if primers are well designed,
# to still retrieve sequences poorly matching with primers and/or poorly aligned, by
# choosing the best positions for cutting the alignment.
}
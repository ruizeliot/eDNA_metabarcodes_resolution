\name{gene_extraction}
\title{
Extracts genes in a DNA sequence at given positions
}
\description{
This function extract genes in each DNA sequence (linear or circular) according to their orientation (i.e. 5' to 3' or 3' to 5') and their starting/ending positions, typically obtained from a GenBank file with \code{\link{gene_position}}.
}
\usage{
gene_extraction(position_data, 
                dna_data, 
                genes, 
                name_col_accession)
}
\arguments{
  \item{position_data}{ A table containing all starting and ending positions (numeric columns) for each genes to extract (i.e. no \code{NA}, \code{NaN} or \code{NULL}), as well as the accession numbers of each sequence. The columns containing the starting positions must be named \code{"STARTING_NAME_GENE"} (NAME_GENE provided in \code{genes}), those with ending positions \code{"ENDING_NAME_GENE"}, and those with the gene orientation \code{"ORIENTATION_NAME_GENE"}, which is what output \code{\link{gene_position}}. }
  \item{dna_data}{ A \code{\link[Biostrings:DNAStringSet]DNAStringSet} containing all DNA sequences (linear or circular) named with their accession numbers only. }
  \item{genes}{ Name of the genes to output in the columns of the final table with extracted genes. It must also be the second part the column names of \code{position_data} with positions and orientation. }
  \item{name_col_accession}{ Name of the column in \code{position_data} containing the accession numbers of each sequence present in \code{dna_data}. }
}
\value{
A \code{\link[tibble:tibble]tibble} containing in each row the accession numbers of the full DNA sequences along with the different genes (one column named \code{FRAGMENT_NAME_GENE} for each gene).
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{gene_position}}, \code{\link{decipher_rapid_alignment}}
}
\examples{
## Creating random data

random_sequences = sapply(1:2, function(x) 
  paste(sample(c("A","G","T","C"), sample(15500:16500, 1), replace = T), collapse = ""))

random_first_position_gene = lapply(1:2, function(x) 
  c(1, cumsum(sample(500:1500, 3))))

random_last_position_gene = lapply(1:2, function(x) 
  c(random_first_position_gene[[x]][-1] - 1, 
    sum(random_first_position_gene[[x]][4] + sample(500:1500, 1))))

random_accession = sapply(1:2, function(x) # Format of NCBI accession numbers
  paste0(c(sample(LETTERS, 2, replace = T), 
           sample(1:9, 6, replace = T)), collapse = ""))

random_species = sample(taxize::plantNames, 2)


## Assembling infos into a typical GenBank file from NCBI (as if read with fread)

typical_ncbi_genbank_file = data.table(do.call(rbind, lapply(1:2, function(x) 
  data.frame(V1 = c(paste0("LOCUS        ", random_accession[x], "      ", 
                           nchar(random_sequences[x]), " bp    DNA     circular"),
                    paste0("DEFINITION   ", random_species[x], 
                           " mitochondrion, complete genome."),
                    paste0("ACCESSION    ", random_accession[x]),
                    paste0("VERSION      ", random_accession[x], ".1"),
                    "... further infos in GENBANK ...",
                    "FEATURES             Location/Qualifiers",
                    paste0("source          1..", nchar(random_sequences[x])),
                    paste0('/organism="', random_species[x], '"'),
                    '/organelle="mitochondrion"',
                    '/mol_type="genomic DNA"',
                    c(rbind(paste0(c("rRNA", "gene", "mRNA", "CDS"),
                                   "            ", random_first_position_gene[[x]], 
                                   "..", random_last_position_gene[[x]]),
                            paste0('/gene="GENE', 1:4, '"'),
                            paste0('/product="detailed name of gene ', 1:4, '"'),
                            "... other infos about the gene ...")),
                    "... same for all others genes in the mitochondrial genome ...",
                    "ORIGIN",
                    "... full DNA sequence (nucleotides) ...",
                    "//",
                    " ")))))
typical_ncbi_genbank_file # Example for 4 genes in a complete mitochondrial genome


## Extracting the positions of each genes

name_genes = paste0("GENE", 1:4)
name_genes

full_name_genes = setNames(lapply(1:4, function(x) 
  paste0(c(paste0('GENE', x), paste0('gene ', x)), collapse = '"|')), name_genes)
full_name_genes

title_row_position = setNames(as.list(c("rRNA  ", "gene  ", "mRNA  ", "CDS  ")), name_genes)
title_row_position

extracted_gene_positions = 
  gene_position(position_gb = typical_ncbi_genbank_file, 
                genes = name_genes, 
                pattern_list = full_name_genes,
                type_list = title_row_position)
extracted_gene_positions


## Extracting each genes for the given positions

# As if it was the sequences from the FASTA file available on NCBI along with the GenBank file
dna_under_biostrings_format = Biostrings::DNAStringSet(random_sequences)
names(dna_under_biostrings_format) = random_accession
dna_under_biostrings_format 

extracted_genes = 
  gene_extraction(position_data = extracted_gene_positions, 
                  dna_data = dna_under_biostrings_format, 
                  genes = name_genes,
                  name_col_accession = "ACCESSION")
extracted_genes
}

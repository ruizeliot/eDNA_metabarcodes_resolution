\name{gene_position}
\title{
Extracts the position of genes in a GenBank file
}
\description{
This function takes a GenBank file (all informations about a sequence in NCBI) and extracts the starting and ending positions of certain genes if present. It also check if the orientation of the sequence (i.e. 5' to 3' or 3' to 5').
}
\usage{
gene_position(position_gb, 
              genes, 
              pattern_list, 
              type_list, 
              col_name = "V1")
}
\arguments{
  \item{position_gb}{ A GenBank file from NCBI, opened as a table in R. }
  \item{genes}{ Name of the genes to output in the final table with positions. The element of the lists given in \code{pattern_list} and \code{type_list} must also be named with these same names to correctly assign the searching informations to the correct genes. }
  \item{pattern_list}{ String with the full name(s) of the gene (often given in \code{/gene=} and/or \code{/product=}) in the GenBank file, in a list named with \code{genes}. Different names given to the same genes (gene names not standardised in NCBI) can be provided under this format \code{NAME1|NAME2|NAME3}. }
  \item{type_list}{ String with the name(s) on the left column on NCBI (not separated in a different column in the GenBank file) identifying the block of infos per genes, in a list named with \code{genes}. The position of the gene in the full sequence is always given in the first row of the block along with this name (same column in the GenBank file). Different names given to the same genes (not standardised) can be provided under this format \code{TYPE1|TYPE2|TYPE3}. }
  \item{col_name}{ Name of the unique column of the GenBank file from NCBI, opened as a table in R. It is \code{V1} by default because the function \code{\link[data.table:fread]fread}, which is useful for quickly opening large GenBank files, gives this name to the column by default. }
}
\details{
This function first searches the index of rows containing \code{type_list}, which is identify blocks of informations for multiple genes. Then, it searches the index of rows containing the full name(s) of the gene of interest \code{pattern_list} in the GenBank file. Finally, it selects the closest occurrence of \code{type_list} to an occurrence of \code{pattern_list} (see below of its works, using \code{tail}), extract the position of the gene, extract the accession number (row with the mention \code{"VERSION"}), and check the orientation (i.e. mention \code{"complement"} or not).

It is useful to add three spaces behind the names in \code{type_list}, since it select only occurrences in left column on NCBI (same column in the GenBank file). It is also useful to provide \code{pattern_list} with a \code{"} behind the name (e.g. \code{NAME1"|NAME2"|NAME3"}), since the gene names often appear in the block of infos within \code{/gene="NAME"} and/or \code{/product="NAME"}.
}
\value{
  \item{COMPLETE_POSITIONS}{A table containing only sequences with positions found for all required genes (not returned if empty). Along with the accession number of each sequences, the table contains three columns per genes, with the starting position of the genes, the ending position, and the orientation of the fragment. }
  \item{MISSING_POSITIONS}{Same table than \code{COMPLETE_POSITIONS} but for sequences missing at least one position for one gene (not returned if empty). }
}
\references{
ARTICLE SOON PUBLISHED

Annotated GenBank Record. (2021, January 12). NCBI. Retrieved March 28, 2022, from https://www.ncbi.nlm.nih.gov/genbank/samplerecord/
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{gene_extraction}}
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

title_row_position = setNames(as.list(c("rRNA   ", "gene   ", "mRNA   ", "CDS   ")), name_genes)
title_row_position

extracted_gene_positions = 
  gene_position(position_gb = typical_ncbi_genbank_file, 
                genes = name_genes, 
                pattern_list = full_name_genes,
                type_list = title_row_position)
extracted_gene_positions


## Demonstration of how the function works internally

# To split the dataframe after each occurence of the word "LOCUS" which is the first category
position_gb_list = split(typical_ncbi_genbank_file, 
                         findInterval(1:nrow(typical_ncbi_genbank_file), 
                                      which(grepl("LOCUS", typical_ncbi_genbank_file$V1)) + 1))[-1]
position_gb_list

# To extract the position of each row named "rRNA ", indicating that it is the 1st gene
index_gene = lapply(position_gb_list, function(x) grep(title_row_position[[1]], x$V1))
index_gene

# To get the position of each row containing the full names of the 1st gene
position = lapply(position_gb_list, function(x) grep(full_name_genes[[1]], x$V1))
position

# To get the last position in which the search_pattern has been found (tail 1 to get the last element)
tail(position[[1]], 1)

# To exclude all position in which "rRNA " occurs after the last position of the search_pattern
index_gene[[1]][index_gene[[1]] < tail(position[[1]], 1)]

# To get the closest position of "rRNA  to the search_pattern (last element of the above result)
tail(index_gene[[1]][index_gene[[1]] < tail(position[[1]], 1)], 1)

# It gives the position of the 1st gene in which is it easy to extract the start and end position
position_gb_list[[1]][tail(index_gene[[1]][index_gene[[1]] < tail(position[[1]], 1)], 1)]
}

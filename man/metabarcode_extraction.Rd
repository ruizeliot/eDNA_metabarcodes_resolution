\name{metabarcode_extraction}
\title{
Metabarcode extraction in aligned sequences for each primer
}
\description{
Extracts the metabarcodes between one or more primers in an alignment containing forward and reverse primers regions, using the same best starting and ending positions for every sequences per primers in the alignment. These combinations of best positions per primers are selected to give the as much fragments as possible with sizes closest to one or two references, and the best match of primers regions with forward and reverse primers.
}
\usage{
metabarcode_extraction(data_primer, 
                       data_aligned, 
                       remove_primers = T, 
                       max_mismatch_whole_primer, 
                       indels, 
                       name_col_alignment, 
                       name_col_accession, 
                       name_col_primer_name, 
                       name_col_orientation, 
                       name_col_forward_primer, 
                       name_col_reverse_primer,
                       length1_col_name, 
                       length2_col_name = NULL, 
                       verbose = T)
}
\arguments{
  \item{data_primer}{ A table containing the names of each primer set (\code{name_col_primer_name}) and the sequences of the forward (\code{name_col_forward_primer}) and reverse primers (\code{name_col_reverse_primer}). It must also contain the column \code{name_col_orientation} to precise explicitely the order of the primers. A reference length for the fragment must also be given in \code{length1_col_name}, and a second can be provided in \code{length2_col_name}. }
  \item{data_aligned}{ A table containing the aligned DNA sequences (\code{name_col_alignment}) and their NCBI accession numbers (\code{name_col_accession}). }
  \item{remove_primers}{ A boolean indicating if only the metabarcode must be returned with the forward and reverse primers surrounding it or not. }
  \item{max_mismatch_whole_primer}{ An integer giving the maximum number of mismatches allowed to hybridise each primer to each sequence with \code{\link[Biostrings:matchPattern]matchPattern}. Currently, only one integer is supported (i.e. same maximum number of mismatches allowed for all forward and reverse primers). }
  \item{indels}{ A boolean specifying if one insertion or deletion ("indels") is allowed to hybridise each primer to each sequence with \code{\link[Biostrings:matchPattern]matchPattern}. It is best to set it as \code{indels = T} if unsure about the presence of indels in at least one primer region, because only best local matches are reported in this case, avoiding treating numerous irrelevant positions. }
  \item{name_col_alignment}{ Name of the column in \code{data_aligned} containing each aligned sequences of fragments with all metabarcodes and primers regions. }
  \item{name_col_accession}{ Name of the column in \code{data_aligned} containing the NCBI accession numbers of each aligned sequences. }
  \item{name_col_primer_name}{ Name of the column in \code{data_primer} containing the names of the primers for which to extract corresponding metabarcodes. }
  \item{name_col_orientation}{ Name of the column in \code{data_primer} explicitely specifying the order of each primer set. It can be \code{"2 x 1"} if the sequences in \code{name_col_forward_primer} correspond to reverse primers, or \code{"1 x 2"} otherwise. Indeed, primers given as forward primers in the original article turn out to be reverse primers sometimes (and inversely for reverse primers), and this argument allows to simply test which configuration is accurate. }
  \item{name_col_forward_primer}{ Name of the column in \code{data_primer} containing the sequences of all forward primers to use. }
  \item{name_col_reverse_primer}{ Name of the column in \code{data_primer} containing the sequences of all reverse primers to use. }
  \item{length1_col_name}{ Name of the column in \code{data_primer} containing a known reference length for the metabarcode alone (\code{remove_primers = T}) or for the metabarcode with primers (\code{remove_primers = F}). }
  \item{length2_col_name}{ A second optional argument (otherwise set to \code{NULL}) to specify a second column name of \code{data_primer} containing the known length of fragments from another reference. In this case, the mean between both lengths is taken, excluding \code{NA}. }
  \item{verbose}{ A boolean indicating if the progression of the extraction (in percentage) should be shown for each primers or not. }
}
\details{
After storing the position of each nucleotides in the alignment, this function unalign each sequence to ease the search of primers with \code{\link[Biostrings:matchPattern]matchPattern}, allowing for a certain number of mismatches and/or one insertion or deletion for each forward and reverse primers. The starting and ending positions of the primers (\code{remove_primers = F}) or the metabarcode between them (\code{remove_primers = T}) is then retrieved. If one or both primers were found on multiple locations, each combinations of starting and ending positions are sorted by decreasing number of matched nucleotides by both primers. The length of the metabarcode region with or without primers (argument \code{remove_primers}) is then computed for each combinations of positions. The combination that gives the closest length to one known length (e.g. given by the authors that designed the primers) or the mean of two known lengths (e.g. two different sources), is then selected as the best combination for a given sequence, to avoid misplacing the primers. If multiple combinations give the same closest length to the reference (often occurs when there is a mismatch or one insertion/deletion), the combination that maximise the number of matched nucleotides by the primers is chosen. 

This operation is repeated for each sequences, to obtain the most probable starting and ending positions in the alignment for each of them. Then, the best starting and ending positions are chosen as the most frequent across all best positions for each sequences. These two positions are used to cut the alignment on the same positions for every sequences and therefore retrieves metabarcodes for all sequences, even if primers are poorly matching or poorly aligned. As, the function saves the results of the primer search for each sequence (i.e. one or both primers found or not found), sequences with too much degenerated primers regions can be easily identified. The function \code{\link{view_metabarcodes}} can then be used to check each of these cases, notably to detect if the metabarcodes are a bit too long or too short due to automatic alignment errors to further correct it manually, or remove wrong metabarcodes (minority of sequences if primers were well designed and the dataset contains enough sequences).
}
\value{
A \code{\link[tibble:tibble]tibble} containing in each row the metabarcodes extracted for each primer and for each aligned sequences. The column \code{PRIMER_FOUND} indicates if none, one or both primers were found with \code{\link[Biostrings:matchPattern]matchPattern}. The best starting and ending positions in the alignment are provided in the column \code{POSITION_IN_ALIGNMENT}. Finally, the length of the fragment with or without primers is also given (note that if \code{remove_primers = F}, both lengths will be the same because the metabarcode will be surrounded by primers).
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{view_metabarcodes}}
}
\examples{
## Simplified extraction process to show how the function works

# Redefining the internal function index_letter in the global environment
index_letter = function(splited_dna){
  
  letter = 0
  
  pos_list = list()
  
  for(r in 1:length(splited_dna)){
    
    if(splited_dna[r] == "-") pos_list[r] = 0
    
    else {
      
      if(letter == 0){
        
        letter = 1
        
        pos_list[r] = letter
        
      }
      
      else{
        
        letter = letter + 1
        
        pos_list[r] = letter
        
      }
      
    }
    
  }
  
  letter_positions = unlist(pos_list)
  
  names(letter_positions) = splited_dna
  
  return(letter_positions)
  
}

# Defining a nucleotide index
test_fragment_aligned = "-----GGGGGGGGGG-----AAATTTGGGCCCAAAAAAAAAA-----TTTTTTTTTTCCCGGGTTTAAA-----CCCCCCCCCC"
nucleotide_index = index_letter(strsplit(test_fragment_aligned, "")[[1]]) 
nucleotide_index

# Unaligning the sequences
test_fragment_unaligned = gsub("-", "", test_fragment_aligned) 
test_fragment_unaligned

#Searching starting position of metabarcode
starting_position = Biostrings::end(Biostrings::matchPattern(
  "AAATTTGGGCCC", test_fragment_unaligned, max.mismatch = 0, with.indels = F)) + 1 
starting_position

# Searching ending position of metabarcode
ending_position = Biostrings::start(Biostrings::matchPattern( 
  "CCCGGGTTTAAA", test_fragment_unaligned, max.mismatch = 0, with.indels = F)) - 1 
ending_position

# Finding to which position it corresponds in the alignment
which(nucleotide_index == starting_position) 
which(nucleotide_index == ending_position) 

# Extracting the aligned fragment along the same position for all sequence (only 1 here)
aligned_metabarcode = substr(test_fragment_aligned, which(nucleotide_index == starting_position), 
                             which(nucleotide_index == ending_position))
aligned_metabarcode 

# Unaligned the metabarcode
unaligned_metabarcode = gsub("-", "", aligned_metabarcode)
unaligned_metabarcode


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
                          name_col_alignment = "ALIGNED_SEQUENCES", 
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
# It is why some metabarcodes between both PRIMER2 are sometimes a bit too short/long
# in this random example.

# With many sequences, manual correction would be too long, and the advantage of our 
# method is therefore to rely on well-aligned and not too much degenerated primers
# regions, that should represent the majority of sequences if primers are well designed,
# to still retrieve sequences poorly matching with primers and/or poorly aligned, by
# choosing the best positions for cutting the alignment.
}
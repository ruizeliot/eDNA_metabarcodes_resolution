\name{accession_2_bin}
\alias{accession_2_bin}
\title{Retrieve BINs and taxonomy from BOLD from NCBI accession Numbers or BOLD IDs}
\description{
The \code{accession_2_bin} function searches for BOLD's BINs and associated taxonomy using NCBI accession numbers or BOLD IDs. It utilizes the \code{bold} R package for querying data.
}

\usage{
accession_2_bin(accession, 
                division_number = 10, 
                type = "accession", 
                final_saving = FALSE, 
                multiple_accession_list = FALSE)
}

\arguments{
\item{accession}{A vector of NCBI accession numbers or a named list of BOLD IDs.}
\item{division_number}{An integer indicating how many subsets to divide the data into for processing.}
\item{type}{A character string specifying the type of input, either \code{"accession"} or \code{"bold_id"}.}
\item{final_saving}{A logical value indicating whether to save the final output. Default is \code{FALSE}.}
\item{multiple_accession_list}{A logical value indicating if the input is a list of multiple accession sets. Default is \code{FALSE}.}
}

\value{
A \code{tibble} containing the BINs and associated taxonomy for each accession number or BOLD ID.
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\examples{
# Returning BOLD's BINs for each accession number as well as BOLD's taxonomy if attributed
fish_accessions = c("AB028664", "AB032554", "AB038381", "AB038556", "AB040676", 
                    "AB042861", "AB042952", "AB046473", "AB043025", "AB043068")
accession_2_bin(fish_accessions, division_number = 10, type = "accession", final_saving = F)

# Can also directly be used with a list of BOLD's IDs (for example returned by bold_identify for a sequence)
fish_bold_id_list = setNames(list(c("GBMNA17497-19", "GBMNA17487-19", "GBMNA17484-19"),
              c("NADBP039-16", "GBMIN133025-17", "GBMIN130151-17")), c("PP109429", "NN252689"))
accession_2_bin(accession = fish_bold_id_list, type = "bold_id", final_saving = F, multiple_accession_list = T)
}
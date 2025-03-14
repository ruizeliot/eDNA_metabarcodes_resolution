\name{id_to_classification}
\alias{id_to_classification}
\title{Retrieve taxonomic classification from TaxID}
\description{
The \code{id_to_classification} function performs consecutive searches of the classification for each TaxID in a local NCBI database. It utilizes the \code{classification} function from the \code{taxizedb} package to retrieve taxonomic classification data. The function extracts the required rank names, assigning \code{NA} where information is missing. Additionally, it retrieves a corrected species name and its rank for further filtering.
}

\usage{
id_to_classification(id, 
                     db = "ncbi")
}

\arguments{
\item{id}{A numeric vector of TaxIDs for which to retrieve taxonomic classification.}
\item{db}{A character string indicating the \code{taxizedb} local database to query. Default is \code{"ncbi"}.}
}

\details{
\code{id_to_classification} operates as follows:
\itemize{
\item Retrieves classification data for the provided TaxIDs using \code{taxizedb::classification}.
\item Extracts the names corresponding to each standard taxonomic rank (superkingdom, kingdom, phylum, class, order, family, genus).
\item Assigns \code{NA} if any taxonomic rank is missing.
\item Retrieves the species name and its rank using \code{taxizedb::taxid2name} and \code{taxizedb::taxid2rank}.
\item Optionally extracts the ID for the specified \code{require_id} rank.
\item Outputs the final result as a \code{tibble} with columns for each rank and the required ID if specified.
}
}

\value{
A \code{tibble} containing the following columns:
\itemize{
\item \code{ID}: The input TaxID.
\item \code{SPECIES}: The corresponding species name.
\item \code{RANK}: The taxonomic rank of the species.
\item \code{GENUS}, \code{FAMILY}, \code{ORDER}, \code{CLASS}, \code{PHYLUM}, \code{KINGDOM}, \code{SUPERKINGDOM}: Names corresponding to each taxonomic rank.
}
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{complete_taxonomy}}
}
\examples{
# Random species names
taxonomic_species = c("Panthera leo", "Macaca mulatta", "Orcinus orca", "Alocasia cuprea",
                      "Bembrops curvatura", "Psammoperca waigiensis")

# Adding NCBI classification from the TaxIDs
taxonomic_species_classified = id_to_classification(taxizedb::name2taxid(taxonomic_species))
taxonomic_species_classified
}
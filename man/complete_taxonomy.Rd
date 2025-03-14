\name{complete_taxonomy}
\alias{complete_taxonomy}
\title{Complete taxonomic information based on multiple taxonomic databases from TaxID}
\description{
Using NCBI TaxIDs the \code{complete_taxonomy} function completes missing taxonomic ranks by leveraging multiple taxonomic databases. It identifies missing ranks and uses the closest known inferior rank to find the corresponding upper taxon. The search sequence prioritizes GBIF, ITIS, Catalog of Life, World Flora Online (for non-fish sequences) since they can be installed locally with the function \code{\link[taxizedb:db_download]{db_download}}, and finally WoRMS due to its slower online query process.
}
\usage{
complete_taxonomy(data, 
                  BASE_TAXA = "automatic")
}
\arguments{
\item{data}{A table containing a column with TaxIDs named \code{ID} along with taxonomic information with the following required columns provided in this order: \code{ID, \code{SPECIES}, \code{RANK} (facultative column), \code{GENUS} (facultative column), \code{FAMILY}, \code{ORDER}, \code{CLASS}, \code{PHYLUM}, \code{KINGDOM}, \code{SUPERKINGDOM}.}}
\item{BASE_TAXA}{A character vector specifying the reference columns for missing data. Set to \code{"automatic"} for automatic detection based on standard taxonomic hierarchy.}
}
\details{
\code{complete_taxonomy} operates in the following manner:
\itemize{
\item Identifies columns with missing taxonomic information.
\item Searches for superior taxa starting with GBIF, followed by ITIS, Catalog of Life, World Flora Online, and finally WoRMS.
\item Uses the closest known inferior rank to guide the search for the missing higher taxon.
\item For WoRMS, it first retrieves the \code{AphiaID} and then queries the taxonomy associated with it.
\item If conflicts arise (e.g., multiple potential taxa), the function conserves all options and separates them with "or" to avoid splitting taxa into multiple synonyms.
}
}
\value{
A table formatted similarly than the input but with completed taxonomic information where possible. Columns with remaining \code{NA} values indicate unresolved taxonomy.
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\seealso{
\code{\link{id_to_classification}}
}
\examples{
# Random species names
taxonomic_species = c("Panthera leo", "Macaca mulatta", "Orcinus orca", "Alocasia cuprea",
                      "Bembrops curvatura", "Psammoperca waigiensis")

# Adding NCBI classification from the TaxIDs
taxonomic_species_classified = id_to_classification(taxizedb::name2taxid(taxonomic_species))
taxonomic_species_classified

# Adding NCBI classification from the TaxIDs
taxonomic_species_complete = complete_taxonomy(taxonomic_species_classified)
taxonomic_species_complete
}
\name{replace_values}
\alias{replace_values}
\title{Replace values in dataframe based on corrected data}
\description{
The \code{replace_values} function replaces values in a large dataframe using corrected values from an independent dataset. The function supports multiple variables and allows conditional replacement based on conflicts.
}

\usage{
replace_values(data_original, 
               data_model, 
               variables_original,
               variables_model = NULL, 
               id_original, 
               id_model = NULL,
               total_replacement)
}

\arguments{
\item{data_original}{A dataframe containing the original data.}
\item{data_model}{A dataframe containing the corrected values.}
\item{variables_original}{A character vector specifying the variable names in \code{data_original} to be replaced.}
\item{variables_model}{A character vector specifying the variable names in \code{data_model} to use for replacement. If \code{NULL}, defaults to \code{variables_original}.}
\item{id_original}{The column name in \code{data_original} representing unique identifiers.}
\item{id_model}{The column name in \code{data_model} representing unique identifiers. If \code{NULL}, defaults to \code{id_original}.}
\item{total_replacement}{A character vector specifying which columns require total replacement, or \code{"NA"} to keep original data if conflicts occur. Use \code{"all"} to replace all columns.}
}

\details{
\code{replace_values} operates as follows:
\itemize{
\item Matches rows in \code{data_original} and \code{data_model} based on the identifier columns.
\item For each matched row, replaces values in the specified columns according to the replacement rules.
\item If \code{total_replacement} is set to \code{"NA"}, original data is kept in case of conflict.
\item If set to \code{"all"}, all specified columns are replaced.
}
}

\value{
A \code{tibble} containing the original data with values replaced based on the corrected dataset. Columns with unresolved conflicts retain their original values if \code{total_replacement} is set to \code{"NA"}.
}
\references{
ARTICLE SOON PUBLISHED
}
\author{
Eliot RUIZ (MARBEC, Uni. Montpellier, France). Send me an email in case of bug: eliotruiz3@gmail.com.
}
\examples{
# Creating a random other column (any form of table is accepted anyway)
random_accession = sapply(1:6, function(x) # Format of NCBI accession numbers
  paste0(c(sample(LETTERS, 2, replace = T), 
           sample(1:9, 6, replace = T)), collapse = ""))

# Creating wrong species names (missing double letters)
species_wrong = c("Pantera leo", "Macaca mulata", "Orcinus orca", "Alocasia cuprea",
                      "Bembrops curvatura", "Psamoperca waigiensis")
species_id_df_wrong = tibble(TAXID = taxizedb::name2taxid(species_wrong), 
                             SPECIES = species_wrong, ACCESSION = random_accession)
species_id_df_wrong

# Creating the replacement table with right names
species_right = c("Panthera leo", "Macaca mulatta", "Psammoperca waigiensis")
species_id_df_replace = tibble(INITIAL_SPECIES = species_wrong[c(1,2,6)], SPECIES = species_right, 
       TAXID = taxizedb::name2taxid(species_right))
species_id_df_replace

# Replace values in the initial table without affecting other columns
replace_values(data_original = cbind(species_id_df_wrong, INITIAL_SPECIES = species_id_df_wrong$SPECIES), 
               data_model = species_id_df_replace, variables_original = c("TAXID", "SPECIES"), 
               id_original = "INITIAL_SPECIES", total_replacement = "all")
}
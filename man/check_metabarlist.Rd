% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_metabarlist.R
\name{check_metabarlist}
\alias{check_metabarlist}
\title{Checks if a list of tables is a \code{metabarlist}}
\usage{
check_metabarlist(metabarlist)
}
\arguments{
\item{metabarlist}{a \code{metabarlist} object}
}
\value{
\code{check_metabarlist} returns TRUE if \code{metabarlist} fulfills of conditions to be of class \code{metabarlist}` (see details). If not, the function throws a stop.
}
\description{
Tests if a list of tables contains correct information on MOTU abundances, MOTU characteristics, PCR characteristics, and sample characteristics to be classed as a \code{metabarlist} object.
}
\details{
The function \code{check_metabarlist} checks if \code{metabarlist} has all the properties required to be a well formed \code{metabarlist} object:

\itemize{
\item {\code{metabarlist} should be a list with four tables named
       `reads`, `motus`, `pcrs` and `samples`.}
\item {the `reads` table should be a numeric matrix.}
\item {the `motus`, `pcrs` and `samples` tables should be data frames.}
\item {the `reads` table row names should be exactly the same as row names
       in the `pcrs` table.}
\item {the `reads` table column names should be exactly the same as row names
       in the `motus` table.}
\item {the `pcrs` table should contain the mandatory columns required,
       i.e.  `sample_id`, `type` and `control_type`.}
\item {values in `type` should be properly defined:
       they can only be `sample` or `control`.}
\item {values in `control_type` should be properly defined:
       they can only be `sequencing`, `pcr`, `extraction`, `positive`, (or NA for samples).}
}

The function will stop if these basic criteria are not met.

\code{check_metabarlist} also issues warnings if tables lack columns which are non-mandatory for the \code{metabaR} package to run, but that are mandatory (or have a mandatory format) for specific functions (e.g., \code{\link{ggpcrplate}})
\itemize{
\item {the column `sequence` for the `motus` table.}
\item {the columns `tag_fwd`, `tag_rev`, `primer_fwd`, `primer_rev`,
       `plate_no`, `plate_col`, and `plate_row` for the `pcrs` table.}
}


In addition, the function issues a warning if any PCR or MOTU has a total count of 0 across the entire dataset.
}
\examples{

data(soil_euk)
check_metabarlist(soil_euk)

}
\seealso{
\code{\link{metabarlist_generator}}
}
\author{
Clément Lionnet & Frédéric Boyer & Lucie Zinger
}

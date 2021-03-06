% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcr_outlier.R
\name{pcr_outlier}
\alias{pcr_outlier}
\alias{filter_replicat}
\alias{coa_function}
\alias{bray_function}
\title{Identify non-conforming PCR replicates}
\usage{
pcr_outlier(
  metabarlist,
  FUN = bray_function,
  groups = metabarlist$pcrs$sample_id,
  graphics = FALSE
)

filter_replicat(sub_matrix, threshold)

coa_function(reads)

bray_function(reads)
}
\arguments{
\item{metabarlist}{a \code{metabarlist} object}

\item{FUN}{a function returning a distance matrix. The distance matrix should be a `dist` object which has the same dimensions as the number of rows in the `reads`table, i.e number of PCRs.}

\item{groups}{a vector containing the replicate identifier. The vector must have the same dimensions as the `PCRs` table from a \code{metabarlist} object. Default = metabarlist$pcrs$sample_id}

\item{graphics}{a boolean value to plot the distance densities for each iteration. Default = FALSE}

\item{sub_matrix}{a distance matrix for replicates comparisons}

\item{threshold}{a threshold below which a pcr is considered as outlier}

\item{reads}{a reads table from a \code{metabarlist} object}
}
\value{
a data frame with the replicate groups and a `replicating`column. If not possible, the function will terminate and return an error message.
}
\description{
Identify the non-conforming sample or control PCR replicates in the PCRs table from a \code{metabarlist} object.
Process numerous reiterations to compare distance densities within PCR replicates and between PCR replicates.
}
\details{
This function identifies non-conforming sample or control replicates.

The parameter `groups` defines groups of replicates. The vector should be arranged following the format of a `PCRs` table from the \code{metabarlist}.
Note: If the distance within replicates is higher than the distance between replicates, the function cannot return any result because all replicates are removed.
The parameter `FUN` defines the function used to compute the distance matrix. The function will return an object of class `dist` with the same length as the input table.
The default function use the `decostand` and `vegdist` functions from the `vegan` package to perform a correspondance analysis of the `reads` table from the \code{metabarlist}, and returns a distance matrix.
Default function detail:
bray_function <- function(reads) {
  distance_matrix <- vegdist(decostand(reads, method = 'total'), method='bray')
  return(distance_matrix)
}

When the `graphics` parameter is defined as True, a graphic is plotted with the density of distance within replicate and between replicate. The threshold is also plotted as a vertical line at the intersection of the two densities.

Note: In the circumstance where numerous sequencing projects have been pooled and analysed on the the same PCR plates, the function must be processed individually for each project, to avoid calculating distances which are meaningless.
}
\section{Functions}{
\itemize{
\item \code{pcr_outlier}: Identifying the non-replicating samples or controls in the table PCRs from a \code{metabarlist} object.

\item \code{filter_replicat}: recursive function to find the non replicating samples or controls

\item \code{coa_function}: distance function with ade4 package and coa analysis

\item \code{bray_function}: distance function with vegan package and Bray-Curtis distance
}}

\examples{
\dontrun{
data(soil_euk)

sample_subset <- subset_metabarlist(soil_euk, "pcrs",
                                    soil_euk$pcrs$type == "sample")
filter_replicat(sample_subset)
}
}
\author{
Frédéric Boyer & Clément Lionnet
}

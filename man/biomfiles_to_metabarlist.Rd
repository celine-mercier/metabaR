% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/biomfiles_to_metabarlist.R
\name{biomfiles_to_metabarlist}
\alias{biomfiles_to_metabarlist}
\title{Import BIOM and associated files to create a metabarlist object}
\usage{
biomfiles_to_metabarlist(
  file_biom,
  file_samples,
  file_pcrs = NULL,
  file_motus = NULL,
  ...
)
}
\arguments{
\item{file_biom}{path for the \code{BIOM} file. This is either a JSON formatted file
(biom file format version 1) or a HDF5 formatted file (\code{BIOM} file format
version 2 and 2.1), as described at \href{http://biom-format.org/}{http://biom-format.org/}.
This file should include at least MOTUs abundance data.
It may also store MOTUs and/or PCRs attributes data.
Mandatory fields for MOTUs and PCRs attributes data are described below.}

\item{file_samples}{path for the sample characteristics table.
The first column of this table should contain the sample names.}

\item{file_pcrs}{path for the PCRs characteristics table (e.g. tags, primers, plate wells, etc.),
if the \code{BIOM} file is missing these data. Mandatory fields:
(i) `sample_id`, i.e. the name of each sample. (ii) `type`,
i.e. the type of PCR; can be `sample` or `control`. (iii) `control_type`,
i.e. the type of control if applicable. Should be: 'NA' for samples,
`extraction` for extraction negative controls, `pcr` for pcr negative controls,
`sequencing` for sequencing negative controls (e.g. unused tag combinations),
and `positive` for positive controls. The first column of this table should
correspond to the names of the PCRs.}

\item{file_motus}{path for the MOTUs characteristics table (e.g. taxonomy, sequence, etc.),
if the \code{BIOM} file is missing these data. Rows of the table should
correspond to MOTUs, and the columns to their characteristics.
Mandatory fields: 'sequence', i.e. the most abundant sequence of the MOTU.
The first column of this table should contain MOTU names.}

\item{...}{other arguments to be pasted from \code{read.table}.}
}
\value{
a \code{metabarlist} object
}
\description{
Imports and formats \code{BIOM} and associated files to create a \code{metabarlist} object.
}
\details{
This function imports a \code{BIOM} file and associated files into \R to create a \code{metabarlist} object. Two files are required: a \code{BIOM} file, as well as a sample characteristics table. If the \code{BIOM} file does not contain PCRs and MOTUs attributes data, two other files containing these data are required. The files are imported in \R, included into a list of class \code{metabarlist} with the \code{\link{metabarlist_generator}} function, and congruencies between all tables are tested with the \code{\link{check_metabarlist}} function.
}
\examples{

\donttest{

dir <- tempdir()
url = "https://raw.githubusercontent.com/metabaRfactory/metabaR_external_data/master/"

litiere_euk_reads_hdf5_file = "litiere_euk_reads_hdf5.biom"
litiere_euk_reads_hdf5_url = paste(url, litiere_euk_reads_hdf5_file, sep="")
litiere_euk_reads_hdf5_path <- file.path(dir, litiere_euk_reads_hdf5_file)
download.file(litiere_euk_reads_hdf5_url, litiere_euk_reads_hdf5_path)

litiere_euk_motus_file = "litiere_euk_motus.txt"
litiere_euk_motus_url = paste(url, litiere_euk_motus_file, sep="")
litiere_euk_motus_path <- file.path(dir, litiere_euk_motus_file)
download.file(litiere_euk_motus_url, litiere_euk_motus_path)

litiere_euk_pcrs_file = "litiere_euk_pcrs.txt"
litiere_euk_pcrs_url = paste(url, litiere_euk_pcrs_file, sep="")
litiere_euk_pcrs_path = file.path(dir, litiere_euk_pcrs_file)
download.file(litiere_euk_pcrs_url, litiere_euk_pcrs_path)

litiere_euk_samples_file = "litiere_euk_samples.txt"
litiere_euk_samples_url = paste(url, litiere_euk_samples_file, sep="")
litiere_euk_samples_path = file.path(dir, litiere_euk_samples_file)
download.file(litiere_euk_samples_url, litiere_euk_samples_path)

soil_euk <- biomfiles_to_metabarlist(
 file_biom = litiere_euk_reads_hdf5_path,
 file_motus = litiere_euk_motus_path,
 file_pcrs = litiere_euk_pcrs_path,
 file_samples = litiere_euk_samples_path,
 sep = "\t")

}

}
\references{
\url{http://biom-format.org/}
}
\seealso{
\code{\link{check_metabarlist}}, \code{\link{metabarlist_generator}},
          \code{\link{obifiles_to_metabarlist}}, \code{\link{tabfiles_to_metabarlist}},
          or the \code{biomformat} package
}
\author{
Anne-Sophie Benoiston
}

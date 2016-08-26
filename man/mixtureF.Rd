\docType{data}
\name{mixtureF}
\alias{mixtureF}
\title{RnaQuin mixture F}
\format{Data frame for each RnaQuin sequin:
\itemize{
  \item Length: Length in nucleotide bases
  \item MXF: Input concentration in attomol/ul
}}
\usage{
data(mixtureF)
}
\value{
A data frame with columns defined in \code{Format}.
}
\description{
RnaQuin sequins are combined together across a range of concentration to
formulate a mixture for emulating expression level. The concentration of
each sequin within a mixture is specificed in a CSV mixture file.

This is a flat mixture (F) for RnaQuin, user can also download the file
directly from \code{s3.amazonaws.com/sequins/mixtures/MRN030_v001.csv}.
}
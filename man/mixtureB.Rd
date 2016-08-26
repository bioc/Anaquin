\docType{data}
\name{mixtureB}
\alias{mixtureB}
\title{RnaQuin mixture B}
\format{Data frame for each RnaQuin sequin:
\itemize{
  \item Length: Length in nucleotide bases
  \item MXB: Input concentration in attomol/ul
}}
\usage{
data(mixtureB)
}
\value{
A data frame with columns defined in \code{Format}.
}
\description{
RnaQuin sequins are combined together across a range of concentration to
formulate a mixture for emulating expression level. The concentration of
each sequin within a mixture is specificed in a CSV mixture file.

This is a staggered mixture (B) for RnaQuin, user can also download the
file directly from \code{s3.amazonaws.com/sequins/mixtures/MRN028_v001.csv}.
}
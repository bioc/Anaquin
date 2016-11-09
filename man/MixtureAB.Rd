\docType{data}
\name{MixtureAB}
\alias{MixtureAB}
\title{RnaQuin mixture AB}
\format{Data frame for each RnaQuin sequin:
\itemize{
  \item Length: Length in nucleotide bases
  \item Fold: Expected fold-change
}}
\usage{
data(MixtureAB)
}
\value{
A data frame with columns defined in \code{Format}.
}
\description{
RnaQuin sequins are combined together across a range of concentration to
formulate a mixture for emulating expression level. The concentration of
each sequin within a mixture is specificed in a CSV mixture file.

This dataset contains the fold-change between mixture A and B. Users can also
download the file from \code{http://www.sequin.xyz}.
}